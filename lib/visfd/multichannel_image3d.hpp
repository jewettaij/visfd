///   @file multichannel_image3d.hpp
///   @brief a class for storing images containing an arbitrary number of
///          channels.  (IE. multiple numbers per voxel, such as r,g,b values
///          or directional vector information.)
///   @author Andrew Jewett
///   @date 2018-2-26


#ifndef _MULTICHANNEL_IMAGE3D_HPP
#define _MULTICHANNEL_IMAGE3D_HPP



namespace visfd {


/// @class CompactMultiChannelImage3D
///
/// @brief   This class is useful IF you have a volumetric multi-channel image
///             (ie, multiple numbers are stored for every voxel)
///          AND you want to avoid allocating space for voxels you don't need.
///        The constructor will allocate the aaaafI[][][][] array which 
///        is public.  The constructor also accepts a "aaafMask" argument.
///        It will only allocate space for voxels if the corresponding entry
///        in the aaafMask[iz][iy][ix] array is non-zero.  If zero, then
///        aaaafI[iz][iy][ix]=nullptr.  If not, it points to an array of N numbers.
///        ("N" is the number of channels, an argument to the constructor.)
///        The resulting array, aaaafI, can be passed to any function
///        which expects multi-channel images (and accepts a aaafMask argument).
///        This includes all tensor-voting operations.
///         
/// @note  In practice, this often does not usually save that much space.
///        Keep in mind that, at a minimum, this data structure allocates a
///        pointer for every voxel in the image (8 bytes for a 64 bit system)
///        regardless of the contents of the aaafMask[][][] array.
///        However for the 6-channel images used in Tensor-Voting, the space
///        savings can be substantial.

template<typename Scalar>

class CompactMultiChannelImage3D
{

private:

  Scalar *afI;       // a 1-D array storing all the numbers 
  size_t n_good_voxels;
  int n_channels_per_voxel;
  int image_size[3];

public:

  Scalar ****aaaafI; // a 3D array of pointers into the afI array


  int
  nchannels() {
    return n_channels_per_voxel;
  }

  CompactMultiChannelImage3D(int set_n_channels_per_voxel) {
    n_channels_per_voxel = set_n_channels_per_voxel;
    n_good_voxels = 0;
    aaaafI = nullptr;
  }

  CompactMultiChannelImage3D(int set_n_channels_per_voxel,
                             int const set_image_size[3],
                             Scalar const *const *const *aaafMask = nullptr,
                             ostream *pReportProgress = nullptr  //!< print progress to the user?
                             )
  {
    n_channels_per_voxel = set_n_channels_per_voxel;
    n_good_voxels = 0;
    aaaafI = nullptr;
    Resize(image_size, aaafMask, pReportProgress);
  }

  void
  Resize(int const set_image_size[3],
         Scalar const *const *const *aaafMask = nullptr,
         ostream *pReportProgress = nullptr  //!< print progress to the user?
         )
  {
    if (aaaafI)
      Dealloc();
    Alloc(set_image_size, aaafMask, pReportProgress);
  }

private:

  void
  Alloc(int const set_image_size[3],
        Scalar const *const *const *aaafMask = nullptr,
        ostream *pReportProgress = nullptr  //!< print progress to the user?
        )
  {
    image_size[0] = set_image_size[0];
    image_size[1] = set_image_size[1];
    image_size[2] = set_image_size[2];

    if (pReportProgress)
      *pReportProgress
        << " -- Attempting to allocate space for a "
        << n_channels_per_voxel << "-channel image\n"
        << " -- (If this crashes your computer, find a computer with\n"
        << " --  more RAM and use \"ulimit\", OR use a smaller image.)\n";

    aaaafI = Alloc3D<Scalar*>(image_size);

    for (int iz = 0; iz < image_size[2]; iz++)
      for (int iy = 0; iy < image_size[1]; iy++)
        for (int ix = 0; ix < image_size[0]; ix++)
          aaaafI[iz][iy][ix] = nullptr;

    for (int iz = 0; iz < image_size[2]; iz++) {
      for (int iy = 0; iy < image_size[1]; iy++) {
        for (int ix = 0; ix < image_size[0]; ix++) {
          if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
            continue;
          n_good_voxels++;
        }
      }
    }
    afI = new Scalar[n_good_voxels * n_channels_per_voxel];
    int n = 0;
    for (int iz = 0; iz < image_size[2]; iz++) {
      for (int iy = 0; iy < image_size[1]; iy++) {
        for (int ix = 0; ix < image_size[0]; ix++) {
          if (aaafMask && (aaafMask[iz][iy][ix] == 0.0))
            continue;
          aaaafI[iz][iy][ix] =
            &(afI[n * n_channels_per_voxel]);
          n++;
        }
      }
    }

    if (pReportProgress)
      *pReportProgress
        << "        done\n" << endl;
  } //Alloc()

  void
  Dealloc()
  {
    delete [] afI;

    // Now delete the array of pointers (aaaafI, which pointed to afI)
    Dealloc3D(aaaafI);
    afI = nullptr;
    aaaafI = nullptr;
  }

public:

  // destructor, copy and move constructor, swap, and assignment operator

  ~CompactMultiChannelImage3D() {
    Dealloc();
  }

  CompactMultiChannelImage3D(const CompactMultiChannelImage3D<Scalar>& source) {
    Resize(source.image_size);
    // copy the numbers
    std::copy(source.afI, source.afI + n_good_voxels, afI);
    // copy the pointers
    for (int iz = 0; iz < image_size[2]; iz++) {
      for (int iy = 0; iy < image_size[1]; iy++) {
        for (int ix = 0; ix < image_size[0]; ix++) {
          // We want the pointers in aaaafI to point into the afI array
          // (instead of the source.afI array), so we copy the offsets.
          if (source.aaaafI[iz][iy][ix]) {
            Scalar* offset = source.aaaafI[iz][iy][ix] - source.afI;
            aaaafI[iz][iy][ix] = afI + offset;
          }
          else
            aaaafI[iz][iy][ix] = nullptr;
        }
      }
    }
  }

  void swap(CompactMultiChannelImage3D<Scalar> &other) {
    std::swap(n_good_voxels, other.n_good_voxels);
    std::swap(n_channels_per_voxel, other.n_channels_per_voxel);
    std::swap(image_size, other.image_size);
    std::swap(afI, other.afI);
    std::swap(aaaafI, other.aaaafI);
  }

  // Move constructor (C++11)
  CompactMultiChannelImage3D(CompactMultiChannelImage3D<Scalar>&& other) {
    this->swap(other);
  }

  // Using the "copy-swap" idiom for the assignment operator
  CompactMultiChannelImage3D<Scalar>&
    operator = (CompactMultiChannelImage3D<Scalar> source) {
    this->swap(source);
    return *this;
  }

}; //class CompactMultiChannelImage3D




} //namespace visfd



#endif //#ifndef _MULTICHANNEL_IMAGE3D_HPP
