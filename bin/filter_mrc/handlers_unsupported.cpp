#include <visfd.hpp>
using namespace visfd;
#include <feature.hpp>
#include <segmentation.hpp>
#include <draw.hpp>
#include "filter3d_variants.hpp"
#include "feature_variants.hpp"
#include "file_io.hpp"
#include "handlers_unsupported.hpp"

/// UNSUPPORTED CODE
/// THIS FILE IMPLEMENTS FEATURES WHICH ARE NO LONGER MAINTAINED AND TESTED.
/// THIS CODE WILL LIKELY BE DELETED IN THE FUTURE.




#ifndef DISABLE_DOGGXY

void
HandleDoggXY(const Settings &settings,
             MrcSimple &tomo_in,
             MrcSimple &tomo_out,
             MrcSimple &mask,
             float voxel_width[3])
{

  cerr << "filter_type = Difference-of-Generalized-Gaussians in the XY plane\n";
  // Generate a filter
  //
  //   h(x,y,z) = h_xy(x,y) * h_z(z)
  //
  //Take advantage of the fact that the filter we
  //are using (ie. the function we are convolving with the source image),
  //is the product of a function of X,Y, with a function of Z.
  //This makes it a "seperable" filter:  We can perform the filters in the Z
  //direction, followed by filtering the result in the X & Y directions.
  //This reduces the computation by a factor of O(filter.halfwidth[2]))
  //(A 1-D convolution followed by a 2-D convolution is much faster per 
  // voxel than a full 3-D convolution.)

  // First, generate the filter in the Z direction:

  int filter_truncate_halfwidthZ = -1;
  Filter1D<float, int> filterZ;
  if (settings.filter_truncate_ratio > 0.0)
    filter_truncate_halfwidthZ = floor(settings.width_a[2] *
                                       settings.filter_truncate_ratio);
  else
    filter_truncate_halfwidthZ = floor(settings.width_a[2] *
                                       sqrt(-2*log(settings.filter_truncate_threshold)));

  filterZ = GenFilterGauss1D(settings.width_a[2],
                             filter_truncate_halfwidthZ);

  float C;       // let the user know what C coefficient was used
  C = filterZ.afH[filter_truncate_halfwidthZ]; //(C=peak height located at the
                                              //   middle of the array)

  // Then generate the filter in the XY directions

  Filter2D<float, int> filterXY;
  float A, B;       // let the user know what A B coefficients were used

  filterXY = GenFilterDogg2D(settings.width_a,//"a" parameter in formula
                             settings.width_b,//"b" parameter in formula
                             settings.m_exp,  //"m" parameter in formula
                             settings.n_exp,  //"n" parameter in formula
                             settings.filter_truncate_ratio,
                             settings.filter_truncate_threshold,
                             &A,
                             &B);

  // Precompute the effect of the filter in the Z direction.

  // Create temporary 1-D arrays to perform the filter in the Z-direction:
  float *afIZorig = new float [tomo_in.header.nvoxels[2]];
  float *afIZnew  = new float [tomo_in.header.nvoxels[2]];
  float *afMask         = nullptr;
  if (mask.aaafI)
    afMask       = new float [tomo_in.header.nvoxels[2]];

  // Then apply the filter in the Z direction
  // (and store the filtered 3-D image in the original tomogram array)
  for (int ix = 0; ix < tomo_in.header.nvoxels[0]; ix++) {
    for (int iy = 0; iy < tomo_in.header.nvoxels[1]; iy++) {
      for (int iz = 0; iz < tomo_in.header.nvoxels[2]; iz++) {
        afIZorig[iz] = tomo_in.aaafI[iz][iy][ix];
        if (afMask)
          afMask[iz] = mask.aaafI[iz][iy][ix];
      }
      filterZ.Apply(tomo_in.header.nvoxels[2],
                    afIZorig,
                    afIZnew,
                    afMask,
                    settings.normalize_near_boundaries);

      //It would be wasteful to allocate a temporary array to store this
      //Instead store the result of the convolution in the original array:
      for (int iz = 0; iz < tomo_in.header.nvoxels[2]; iz++)
        tomo_in.aaafI[iz][iy][ix] = afIZnew[iz];
    } //for (int iy = 0; iy < tomo_in.header.nvoxels[1]; iy++) {
  } // for (int ix = 0; ix < tomo_in.header.nvoxels[0]; ix++) {
  delete [] afIZorig;
  delete [] afIZnew;
  if (afMask)
    delete [] afMask;

  // Now apply the filter in the X and Y directions:
  cerr << "  progress: processing plane#" << endl;
  for (int iz = 0; iz < tomo_in.header.nvoxels[2]; iz++) {
    cerr << "  " << iz+1 << " / " << tomo_in.header.nvoxels[2] << "\n";
    float **aafMaskXY = nullptr;
    if (mask.aaafI)
      aafMaskXY = mask.aaafI[iz];
    filterXY.Apply(tomo_in.header.nvoxels,
                   tomo_in.aaafI[iz],
                   tomo_out.aaafI[iz],
                   aafMaskXY,
                   false);
  }

  cerr << " Filter Used:\n"
    //" h(x,y,z) = (h_a(x,y) - h_b(x,y)) * C * exp(-0.5*(z/s)^2)\n"
    " h(x,y,z) = (h_a(x,y) - h_b(x,y)) * C * exp(-(z/s)^2)\n"
    " h_a(x,y) = A*exp(-((x/a_x)^2 + (y/a_y)^2)^(m/2))\n"
    " h_b(x,y) = B*exp(-((x/b_x)^2 + (y/b_y)^2)^(n/2))\n"
    "  ... where  A = " << A << "\n"
    "             B = " << B << "\n" 
    "             C = " << C << "\n" 
    "             m = " << settings.m_exp << "\n"
    "             n = " << settings.n_exp << "\n" 
    "   (a_x, a_y) = "
       << "(" << settings.width_a[0]
       << " " << settings.width_a[1] << ")\n"
    "   (b_x, b_y) = "
       << "(" << settings.width_b[0]
       << " " << settings.width_b[1] << ")\n"
    "            s = " << settings.width_a[2] << endl;
  cerr << " You can plot a slice of this function\n"
       << "     in the X direction using:\n"
    " draw_filter_1D.py -dogg " << A*C << " " << B*C
       << " " << settings.width_a[0] << " " << settings.width_b[0]
       << " " << settings.m_exp << " " << settings.n_exp << endl;
  if ((settings.width_a[1] != settings.width_a[0]) ||
      (settings.width_a[2] != settings.width_a[0])) {
    cerr << " and in the Y direction using:\n"
      " draw_filter_1D.py -dogg " << A*C << " " << B*C
         << " " << settings.width_a[1] << " " << settings.width_b[1]
         << " " << settings.m_exp << " " << settings.n_exp << endl;
    cerr << " and in the Z direction using:\n"
      " draw_filter_1D.py -gauss " << C   // * (A-B)   <--COMMENTING OUT (A-B)
         << " " << settings.width_a[2] << endl;
  }
} //HandleDoggXY()






#ifndef DISABLE_INTENSITY_PROFILES
void
HandleBlobRadialIntensity(const Settings &settings,
                          MrcSimple &tomo_in,
                          MrcSimple &tomo_out,
                          MrcSimple &mask,
                          float voxel_width[3])
{
  vector<array<float, 3> > sphere_centers; // coordinates for the center of each sphere (blob)
  vector<float> diameters;         // diameter of each sphere
  vector<float> scores;            // the "score" (contrast) of each blob

  assert(settings.in_coords_file_name != "");
  assert(settings.blob_profiles_file_name_base != "");

  ReadBlobCoordsFile(settings.in_coords_file_name,
                     &sphere_centers, // store blob coordinates here
                     &diameters,
                     &scores,
                     voxel_width[0],
                     settings.sphere_decals_diameter,
                     settings.sphere_decals_foreground,
                     settings.sphere_decals_scale);

  size_t Nblobs = sphere_centers.size();
  assert(Nblobs == diameters.size());
  assert(Nblobs == scores.size());

  vector<vector<float> > intensity_profiles(Nblobs); //!< store the intensity profiles here

  int image_size[3];
  for (int d = 0; d < 3; d++)
    image_size[d] = tomo_in.header.nvoxels[d];

  if ((sphere_centers.size() > 0) && (mask.aaafI != nullptr))
    DiscardMaskedBlobs(sphere_centers,
                       diameters,
                       scores,
                       mask.aaafI);

  cerr << "  creating intensity-vs-radius profiles for "
       << sphere_centers.size() << " blobs.\n" << endl;

  for (size_t i = 0; i < intensity_profiles.size(); i++) {

    int blob_effective_center[3];

    float radius_profile_width;
    #ifdef INTENSITY_PROFILE_RADIUS
    radius_profile_width = INTENSITY_PROFILE_RADIUS;
    #endif

    BlobIntensityProfile(image_size,
                         tomo_in.aaafI,
                         mask.aaafI,
                         sphere_centers[i],
                         diameters[i],
                         settings.blob_profiles_center_criteria,
                         intensity_profiles[i],
                         radius_profile_width,
                         blob_effective_center);

    stringstream intensity_vs_r_file_name_ss;
    intensity_vs_r_file_name_ss
      << settings.blob_profiles_file_name_base
      << "_" << i+1 << ".txt";
    cerr << "  creating \"" << intensity_vs_r_file_name_ss.str() << "\"" << endl;
    fstream f;
    f.open(intensity_vs_r_file_name_ss.str(), ios::out);
    if (! f)
      throw VisfdErr("Error: unable to open \""+
                     intensity_vs_r_file_name_ss.str()+"\" for writing.\n");
    for (int ir = 0; ir < intensity_profiles[i].size(); ir++) {
      f << ir*voxel_width[0] << " " << intensity_profiles[i][ir] << "\n";
    }
    f.close();




    // CRUFT ALERT!!  REMOVE THIS CRUFT IN THE #ifdef BELOW:
    // (This was something temporary for a side project we
    //  were working on in 2019-6.)
    #ifdef ENABLE_MEASURE_DISTANCE_TO_BOUNDARY

    int image_size[3];
    int sphere_centers_i[3];
    for (int d=0; d<3; d++) {
      image_size[d] = tomo_in.header.nvoxels[d];
      //convert the blob location from floats to ints:
      sphere_centers_i[d] = static_cast<int>(floor(sphere_centers[i][d]+0.5));
    }


    double ave_brightness = 0.0;
    double stddev_brightness = 0.0;
    {
      int Rsphere = ceil(diameters[i]/2); //radius of the sphere (surrounding the current blob)
      int ixs = sphere_centers_i[0];
      int iys = sphere_centers_i[1];
      int izs = sphere_centers_i[2];
      double sum = 0.0;
      size_t n_vox = 0;
      for (int jz = -Rsphere; jz <= Rsphere; jz++) {
        int izs_jz = izs + jz;
        for (int jy = -Rsphere; jy <= Rsphere; jy++) {
          int iys_jy = iys + jy;
          for (int jx = -Rsphere; jx <= Rsphere; jx++) {
            int ixs_jx = ixs + jx;
            if ((jx*jx + jy*jy + jz*jz) > Rsphere*Rsphere)
              continue;
            if (! ((0 <= ixs_jx) && (ixs_jx < image_size[0]) &&
                   (0 <= iys_jy) && (iys_jy < image_size[1]) &&
                   (0 <= izs_jz) && (izs_jz < image_size[2])))
              continue;
            if (mask.aaafI && (mask.aaafI[izs_jz][iys_jy][ixs_jx] == 0.0))
              continue;
            sum += tomo_in.aaafI[izs_jz][iys_jy][ixs_jx];
            n_vox++;
          }
        }
      }
      ave_brightness = sum / n_vox;
      double sum_sq = 0.0;
      for (int jz = -Rsphere; jz <= Rsphere; jz++) {
        int izs_jz = izs + jz;
        for (int jy = -Rsphere; jy <= Rsphere; jy++) {
          int iys_jy = iys + jy;
          for (int jx = -Rsphere; jx <= Rsphere; jx++) {
            int ixs_jx = ixs + jx;
            if ((jx*jx + jy*jy + jz*jz) > Rsphere*Rsphere)
              continue;
            if (! ((0 <= ixs_jx) && (ixs_jx < image_size[0]) &&
                   (0 <= iys_jy) && (iys_jy < image_size[1]) &&
                   (0 <= izs_jz) && (izs_jz < image_size[2])))
              continue;
            if (mask.aaafI && (mask.aaafI[izs_jz][iys_jy][ixs_jx] == 0.0))
              continue;
            sum_sq += SQR(tomo_in.aaafI[izs_jz][iys_jy][ixs_jx]-ave_brightness);
          }
        }
      }
      if (n_vox > 0) {
        double variance = sum_sq / n_vox;
        stddev_brightness = sqrt(variance);
      }
    } // calculate "stddev_brightness"


    // Calculate the distance of this blob (blob i) to
    // every type of object in the image.
    // Here, the mask image is assumed to store a segmentation of the image.
    // Every voxel in the mask image is an integer (a "label") indicating
    // to what (kind of) object that voxel belongs.
    // How many types of different objects are there in the image?
    size_t n_labels = 0;
    if (mask.aaafI) {
      mask.FindMinMaxMean(); //stores the maximum brightness voxel in "dmax"
      n_labels = static_cast<int>(floor(mask.header.dmax + 0.5)) + 1;
      if (n_labels > 1000000) //HACK to make sure we don't overflow memory
        n_labels = 1000000;
    }
    vector<float> min_dist_sq(n_labels,
                              std::numeric_limits<float>::infinity());
    // Now find the distance of blob i (sphere_centers_i) to each label type.
    // (I was too lazy to bother using an efficient algorithm here. 
    //  For most applications this is fast enough.)
    for (int iz = 0; iz < image_size[2]; iz++) {
      for (int iy = 0; iy < image_size[1]; iy++) {
        for (int ix = 0; ix < image_size[0]; ix++) {
          if (! mask.aaafI)
            continue;
          int label = static_cast<int>(floor(mask.aaafI[iz][iy][ix] + 0.5));
          if (label < n_labels) {
            float dist_sq = (SQR(ix-sphere_centers_i[0]) +
                             SQR(iy-sphere_centers_i[1]) +
                             SQR(iz-sphere_centers_i[2]));
            if ((min_dist_sq[label] < 0.0) || (dist_sq < min_dist_sq[label]))
              min_dist_sq[label] = dist_sq;
          }
        }
      }
    }

    // calculate the maximum slope of intensity-vs-r
    double max_slope = 0.0;
    for (int ip = 0; ip < intensity_profiles[i].size()-1; ip++) {
      double slope = intensity_profiles[i][ip] - intensity_profiles[i][ip+1];
      if ((ip==0) || (slope > max_slope))
        max_slope = slope;
    }

    // calculate the minimum and maximum values of intensity-vs-r
    double contrast_profile_min_max = 0.0;
    double profile_min = ave_brightness;
    double profile_max = ave_brightness;
    {
      if (intensity_profiles[i].size() > 0) {
        profile_max = intensity_profiles[i][0];
        for (int ip = 0; ip < intensity_profiles[i].size(); ip++) {
          if ((ip==0) || (intensity_profiles[i][ip] < profile_min))
            profile_min = intensity_profiles[i][ip];
        }
        contrast_profile_min_max = profile_max - profile_min;
      }
    }

    // now estimate the width of the peak:
    // I define the "width" as the distance from the peak's center where the
    // average the brightness at that distance from the center equals the
    // average brightness of all of the voxels within the entire blob (ie all
    // of the voxels within the sphere that circumscribes the entire blob).
    double peak_width = 0.0;
    {
      if (intensity_profiles[i].size() > 0) {
        // REMOVE THIS CRUFT:
        //int i_steepest_slope = 0;
        //double steepest_slope = 0.0;
        //for (int j = 1; j < intensity_profiles[i].size(); j++) {
        //  double slope = intensity_profiles[i][j]-intensity_profiles[i][j-1];
        //  if (slope < steepest_slope) {
        //    steepest_slope = slope;
        //    i_steepest_slope = j;
        //  }
        //}
        //peak_width = i_steepest_slope+0.5;

        // COMMENTING OUT:
        //The width of the peak is the radius where the brightness decays
        //to the average brightness surrounding the peak.
        //float peak_width_threshold = ave_brightness;

        //The width of the peak is the radius where the brightness decays below
        //a value equal to half way between the brightness of the peak
        //(profile_max) and the minimum of the profile curve (profile_min).

        float peak_width_threshold = 0.5*(profile_max + profile_min);
        int i_width = 0;

        for (int ip = 0; ip < intensity_profiles[i].size(); ip++) {
          if (intensity_profiles[i][ip] <= peak_width_threshold) {
            i_width = ip;
            break;
          }
        }
        if (i_width == 0)
          peak_width = 0.0;
        else {
          // interpolate between i_width-1 and i_width to find the
          // distance where the brightness-vs-i intersects "peak_width_threshold"
          double x = 0.0;
          // "x" is the fraction of the distance between i_width-1 and i_width
          // where the linear interpolation of the curve passes through the
          // "peak_width_threshold".  Solving this equation for x:
          //   (1-x) * intensity_profiles[i][i_width-1] +
          //     x   * intensity_profiles[i][i_width]
          //   = peak_width_threshold
          // solution -->
          x = ((peak_width_threshold -
                intensity_profiles[i][i_width-1])
               /
               (intensity_profiles[i][i_width] -
                intensity_profiles[i][i_width-1]));
          // peak_width is the "x" weighted average of i_width-1 and i_width
          peak_width = (i_width-1)*(1.0-x) + i_width*x;
        }
      }
    }
    cout << sphere_centers[i][0]<<" "<<sphere_centers[i][1]<<" "<<sphere_centers[i][2]
         << " " << diameters[i]
         << " " << scores[i]
         << " " << ave_brightness
         << " " << stddev_brightness
         << " " << intensity_profiles[i][0]  // <--brightest voxel's brightness
         << " " << peak_width
         << " " << max_slope
         << " " << contrast_profile_min_max
         << " " << blob_effective_center[0]<<" "<<blob_effective_center[1]<<" "<<blob_effective_center[2];

    // Now print the distance to various labelled objects in the mask image
    for (int il = 0; il < n_labels; il++) {
      cout << " " << ((min_dist_sq[il]!=std::numeric_limits<float>::infinity())
                      ?
                      (sqrt(min_dist_sq[il]) - diameters[i]/2)*voxel_width[0]
                      :
                      std::numeric_limits<float>::infinity());
      // note: sqrt(inf) should equal inf, but results are compiler dependent.
      //       The ? operator used above should insure that sqrt(inf) = inf
    }
    cout << "\n";

    #endif //#ifdef ENABLE_MEASURE_DISTANCE_TO_BOUNDARY <-- REMOVE THIS CRUFT!



  } //for (size_t i = 0; i < intensity_profiles.size(); i++)

} // HandleBlobIntensityProfiles()

#endif //#ifndef DISABLE_INTENSITY_PROFILES






#ifndef DISABLE_BOOTSTRAPPING

void
HandleBootstrapDogg(const Settings &settings,
                    MrcSimple &tomo_in,
                    MrcSimple &tomo_out,
                    MrcSimple &mask);
{
  cerr << "filter_type = Difference-of-Generalized-Gaussians with BOOTSTRAPPING\n"
       << "              (BOOTSTRAP_DOGG)\n";

  Filter3D<float, int> filter;
  float A, B;       // let the user know what A B coefficients were used

  filter = GenFilterDogg3D(settings.width_a,//"a" parameter in formula
                           settings.width_b,//"b" parameter in formula
                           settings.m_exp,  //"m" parameter in formula
                           settings.n_exp,  //"n" parameter in formula
                           settings.filter_truncate_ratio,
                           settings.filter_truncate_threshold,
                           &A,
                           &B,
                           &cerr);

  if (settings.bs_ntests > 0) {

    // IN OLDER VERSIONS OF THE CODE, I WOULD ONLY CONSIDER VOXELS WHOSE
    // FILTERED VALUE EXCEEDED SOME USER-SPECIFIED THRESHOLD.
    // NOW, WE USE THE "MASK" INSTEAD.  IN OTHER WORDS, THE USER IS
    // RESPONSIBLE FOR RUNNING A FILTER ON THE INPUT IMAGE (OR USING SOME
    // OTHER METHOD) TO SELECT OUT THE VOXELS THEY WANT US TO CONSIDER.
    // SAVE THESE VOXELS IN A DIFFERENT FILE, AND LOAD IT USING "-mask"
    //
    // (If the don't, this procedure will take a really, really long time.)
    //
    if (mask.aaafI == nullptr) {
      cerr << "WARNING: THIS FILTER IS VERY SLOW UNLESS YOU USE THE \"-mask\" ARGUMENT\n"
           << "         TO SPECIFY THE VOXELS YOU WANT TO CONSIDER.\n"
           << "         (BY DEFAULT, THIS PROGRAM USES ALL THE VOXELS).\n";
      for (int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
        for (int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
          for (int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
            vvGoodVoxels.push_back(ixiyiz);
            vCountFP.push_back(0);
          }
        }
      }
    }
    else { 
      for (int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
        for (int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
          for (int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
            if (mask.aaafI[ix][iy][iz] != 0.0) {
              vvGoodVoxels.push_back(ixiyiz);
              vCountFP.push_back(0);
            }
          }
        }
      }
    }
          
    //Initialize the random number generator we will need later:
    RANDOM_INIT(settings.bs_random_seed);

    for (int i_sim = 0; i_sim < settings.bs_ntests; i_sim++) {
      cerr <<
        "\n"
        " ---------------------------------------------------\n"
        "   False positives bootstrap simulation# " << isim+1
           << " / " << settings.bs_ntests << "\n"
           << endl;

      ScrambleImage(tomo_in.header.nvoxels,
                    tomo_in.aaafI,
                    settings.bs_scramble_radius,
                    aaafScrambled);

      for (i_vox = 0; i_vox < vvGoodVoxels.size(); i_vox++) {
        int ix = vGoodVoxels[i_vox][0];
        int iy = vGoodVoxels[i_vox][1];
        int iz = vGoodVoxels[i_vox][2];
        float g =
          filter.ApplyToVoxel(ix, iy, iz,
                              tomo_in.header.nvoxels,
                              //tomo_in.aaafI,
                              aaafScrambled,
                              mask.aaafI,
                              false);

        //if (g * settings.bs_threshold_sign >= //filter applied to scramble
        //    tomo_out.aaafI[iz][iy][ix]) //>= applied to orig image
        //  vCountFP[i_vox] += 1;

        if (g * settings.bs_threshold_sign >=
            settings.bs_threshold)
          vCountFP[i_vox] += 1;

      }
    } // for (int i_sim = 0; i_sim < settings.bs_ntests; i_sim++)

    // Now copy the measured probabilities into the output image:
    for (i_vox = 0; i_vox < vvGoodVoxels.size(); i_vox++) {
      int ix = vGoodVoxels[i_vox][0];
      int iy = vGoodVoxels[i_vox][1];
      int iz = vGoodVoxels[i_vox][2];
      // The tomogram should store the 1 - (probability of a false positive)
      // at this location.  In the code above, this probability
      // has been calculated by counting the number of times that a
      // (locally) scrambled image generated a signal as strong as the
      // result of applying the same filter to the unscrambled image there.
      tomo_out.aaafI[iz][iy][ix] =
        (1.0
         -
         static_cast<float>(vCountFP[i_vox]) / settings.fp_nsims);
    }

    ////Dealloc3D(tomo_in.header.nvoxels,
    ////          &aiCount,
    ////          &aaaiCount);

  } //if (settings.bs_ntests > 0)

} //HandleBootstrapDogg()


#endif //#ifndef DISABLE_BOOTSTRAPPING






#ifndef DISABLE_TEMPLATE_MATCHING


      // First, lets characterize the problem we are trying to solve
      // graphically and mathematically:
      //
      // 1-dimensional example:
      //
      //  I_i
      //            original image voxel intensities
      //  /|\                 __   
      //   |                 |  |   __
      //   |            __   |  |  |  |
      //   |           |  |  |  |__|  |
      //   |           |  |__|        |   __
      //   |           |              |  |  |
      //   |           |              |__|  |
      //   |      __   |                    |
      //   |     |  |__|                    |__
      //   |   __|                             |____
      //   |_____________________________________________\  i (voxel position)
      //        .. i0-2  i0-1  i0  i0+1  i0+2 ...        /
      //
      //
      // Notation:
      //  i = index (this integer corresponds to a voxel from the template)
      //      Don't worry about the fact that this is a 3-dimensional image.
      //      Pretend that images are 1-dimensional, and i is an integer.
      // I_ = the collection of intensities of voxels from the original image
      //    = (I_1, I_2, I_3, ... )  in vector notation
      //
      //
      // p_ = intensities of voxels from the image after translation by -i0
      //    = (  p_1,      p_2,      p_3,    ... )
      //    = ( I_(1-i0), I_(2-i0), I_(3-i0) ... )
      //
      // This will shift the image so that voxels that appear nearby i0 are
      // now located nearby 0.  (I.E, it will center the image at position i0)
      //
      //
      //  [Alternate notation:  p_ =  T(-i0) I_   
      //  where "T(-i0)" is the Translational Shift Operator acting on vector I_
      //   T(-i0) I_  =  (I_(1-i0), I_(2-i0), I_(3-i0), ... ) ]
      //
      //
      // Now, we want to compare this with the shape of another function, q_i
      // which stores the shape of template we are looking for in our image.
      // 
      //  q_i
      //
      //  /|\             template intensities  
      //   |                 __,--------.__              
      //   |              __/              \__
      //   |             /                    \
      //   |           _/                      \_
      //        ___,--'                          '--.____
      //                               
      //  /_________________________________________________\ i (voxel position)
      //  \            -4 -3 -2 -1  0  1  2  3  4           /
      //
      // q_ = the collection of intensities of voxels from the template
      //      which we are comparing with the image.
      //    = (..., q_-2, q_-1, q_0, q_1, q_2, q_3, ... )  in vector notation
      //
      // PROBLEM:
      //
      // Unfortunatley, the original image will be noisy, and the overall
      // brightness and contrast of the object will vary significantly from
      // image to image (and sometimes at different locations in the same image)
      // For example, we cannot say what the density (brightness) of, say,
      // a ribosome in a cell should be or even that it should be,
      // say 10% more dense than its surroundings.
      //
      // To see if the template fits the shape of the image at all,
      //   (...at this location in the image)
      // we slide and stretch the template function in the vertical
      // direction until it agrees with the original image as much as possible.
      // 
      // Mathematically, our goal is to solve for the scalars b,c which minimize
      // the sum-of-squares-difference between the intensities of the image and
      // the template, after scaling and intensity offsets have been applied
      //
      // Notation:
      // We can scale the intensity of voxels in the template
      // by multiplying the vector q_ by a scalar (which we call "c")
      // Likewise can offset intensity of voxels by adding a vector b*J
      // where "b" is a scalar
      //   and "J_" is a vector filled with 1s,  J_ = (1, 1, 1, ...)
      //
      // In this notation, we see that the goal is to choose b, c to minimize:
      //
      //   variance  =  |   p_   -   c*(q_ + b*J_) |^2
      //                                             
      //                   /|\      /|\     /|\
      //                    |        |       |
      //                    |        |       |
      //               translated  scaling intensity
      //                  image            offset
      //               (centered 
      //                 at i0)
      //
      // The "variance" is the mean squared difference between
      //   the intensity of the voxels in our source image,  and ...
      //   the intensity of our voxels in our template after shifting & scaling
      //
      // After we have solved for b and c, we will consider the template
      // to be a good match (at this location) if:
      //    1) the variance is small
      //    2) c is large
      //
      // Equivalently, we define the scoring function
      //
      //    c / sqrt(variance)
      //
      //    (A good match <--> lambda is small)
      //
      // First, we have to find an easy way to find the optimal values of b & c:
      //
      // --- Solution Method: ---
      //
      // First shift both the image and the template up and down, so that
      // have the same average value (considering voxels nearby)
      // Introduce :
      //   P_ = T(-i0) p_ - <p_> J  = ( p_1-<p_>, p_2-<p_>, p_3-<p_>, ... )
      //   Q_ =        q_ - <q_> J  = ( q_1-<q_>, q_2-<q_>, q_3-<q_>, ... )
      //
      // Now the goal is to find c to minimize
      //
      //    | P_ - c * Q_ |^2
      //
      // The optimal solution is to choose c so that p'_ - c*q'_
      // lies along the dotted line below
      //
      //              _.
      //        P_    /|
      //             / :
      //            /  :
      //           /   : P_ - c*Q_
      //          /    : (optimal c)
      //         /     : 
      //        /      :
      //       /       :
      //      /        :
      //     /         :
      //    /__________________\ Q_
      //                       /
      //
      //    __________\
      //              / c*Q_
      //                         
      //         where |c*Q_|^2 = <P_, P_> - <P_,Q_>^2/|Q_|^2
      //
      // SOLUTION:
      //
      //        c = <P_, Q_> / |Q_|^2
      //
      // variance = | P_  -  c Q_ |^2
      //          = < (P_ - c Q_), (P_ - c Q_) >
      //          = <P_, P_> - 2*c*<P_, Q_> + c^2*<Q_, Q_>
      //          = <P_, P_> - <P_,Q_>^2/|Q_|^2
      //
      // ...where we introduce the weighted inner product ("dot product")
      // between vectors P_ and Q_:
      //
      // <P_, Q_>  =  sum_i( w_i * P_i * Q_i )
      //
      // The norm (vector length) is also weighted by w_i:
      // | Q |^2 = <Q_, Q_> = sum_i( w_i * Q_i * Q_i )
      //
      //
      //
      // Note: Why do we use weighted inner products (w_i)?
      //
      // We only want to match the template to a region in the image locally,
      // ie., in the viscinity of i0 (see above)
      // To insure this, we weight all the terms in all the sumations by
      //       w_i  = a function of i which is large when i is near 0,
      //              and zero elsewhere
      //       w_i is usually a Gaussian whose width is at least as wide 
      //           as template we are fitting:
      //       w_i = A * exp(-(i/s)^2)
      //
      // ... and we only perform the sum in regions where w_i is non-zero



void
HandleTemplateGGauss(const Settings &settings,
                     MrcSimple &tomo_in,
                     MrcSimple &tomo_out,
                     MrcSimple &mask,
                     float voxel_width[3])
{

  // Filter weights, w_i:
  Filter3D<float, int>
    w = GenFilterGenGauss3D(settings.template_background_radius,
                            settings.n_exp,
                            //template_profile.halfwidth,
                            settings.filter_truncate_ratio,
                            settings.filter_truncate_threshold);

  // GenFilterGenGauss3D() creates normalized gaussians with integral 1.
  // That's not what we want for the weights, w_i:
  // The weights w_i should be 1 in the viscinity we care about, and 0 outside
  // So, we can undo this normalization by dividing all the w_i weights
  // by their maximum value max(w_i)    (at the central peak, at halfwidth).
  // This will mean the maximum w_i is 1, and the others decay to 0, as we want
  float wpeak = w.aaafH[0][0][0];
  w.MultiplyScalar(1.0 / wpeak);

  // Template function, Q:
  TemplateMatcher3D<float, int>
    Q = GenFilterGenGauss3D(settings.width_a, //"a" parameter in formula
                            settings.m_exp,   //"m" parameter in formula
                            w.halfwidth);

  // Q_ = q_ - <q_>
  float qave = Q.Average(w.aaafH);
  Q.AddScalar(-qave);

  // It's also convenient to scale the voxel intensities in Q_ so that the 
  // range of values, |Q_|^2, (=the standard deviation) is also unity.
  // This makes it meaningful to compare "c" with "sqrt(variance)" directly.

  // Now calculate Q_dot_Q
  float Q_dot_Q = Q.SumSqr(w.aaafH);
  Q.MultiplyScalar(1.0/sqrt(Q_dot_Q));
  Q_dot_Q = Q.SumSqr(w.aaafH); // = |Q_|^2  (same everywhere)
  assert(abs(Q_dot_Q - 1.0) < 0.001);  // (should equal 1 but may vary due to roundoff error of 32bit floats)


  cerr << " ------ Calculating the average of nearby voxels: ------\n";

  float template_background_sigma[3];
  for (int d = 0; d < 3; d++)
    template_background_sigma[d] = (settings.template_background_radius[d]
                                    /
                                    sqrt(3.0));

  // P = original image (after subtracting average nearby intensities):

  MrcSimple P = tomo_in; //this will take care of allocating the array

  // First, let's calculate the weighted average voxel intensity in the 
  // source image
  if (settings.n_exp == 2.0) {

    // then do it the fast way with seperable (ordinary) Gaussian filters
    ApplyGauss(tomo_in.header.nvoxels,
               tomo_in.aaafI,
               P.aaafI,    // <-- save result here
               mask.aaafI,
               template_background_sigma,
               settings.filter_truncate_ratio,
               settings.filter_truncate_threshold,
               settings.normalize_near_boundaries);
  }
  else {
    w.Apply(tomo_in.header.nvoxels,
            tomo_in.aaafI,
            P.aaafI,    // <-- save result here
            mask.aaafI,
            settings.normalize_near_boundaries,
            &cerr);
  }

  // Subtract the average value from the image intensity, and store in P:
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
        P.aaafI[iz][iy][ix] = tomo_in.aaafI[iz][iy][ix]-P.aaafI[iz][iy][ix];

  cerr << "\n"
    " ------ Calculating <P_, Q_> / |Q_|^2 ... -------\n" << endl;

  // Now calculate <P_, Q_>
  // Recall <P_, Q_> = sum_i (w_i * Q_i * P_i)
  // precompute w_i * Q_i and call it "Q_times_w_i"
  //        <P_, Q_> = sum_i (Q_times_w_i * P_i)
  // Then we can reuse the code we normally use to calculate sum_i(A_i B_i)
  Filter3D<float, int> Q_times_w = w;
  for (int jz=-w.halfwidth[2]; jz<=w.halfwidth[2]; jz++)
    for (int jy=-w.halfwidth[1]; jy<=w.halfwidth[1]; jy++)
      for (int jx=-w.halfwidth[0]; jx<=w.halfwidth[0]; jx++)
        Q_times_w.aaafH[jz][jy][jx] = 
          w.aaafH[jz][jy][jx] * Q.aaafH[jz][jy][jx];
               

  // Calculate P_dot_Q = <P_, Q_> = sum_i (Q_times_w_i * P_i)

  MrcSimple P_dot_Q = tomo_in;  //(note: this will allocate P_dot_Q.aaafI)

  Q_times_w.Apply(tomo_in.header.nvoxels,
                  P.aaafI,
                  P_dot_Q.aaafI, // <-- store result here
                  mask.aaafI,
                  false,
                  &cerr);


  //cerr << " ----- DEBUGGING: Calculating RMSE of fit directly -----\n"<< endl;

  //Q._ScanTemplateError(P.header.nvoxels,
  //                     P.aaafI,
  //                     tomo_out.aaafI,
  //                     P_dot_Q.aaafI,
  //                     w.aaafH,
  //                     2.0,
  //                     mask.aaafI,
  //                     &cerr);
  //for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
  //  for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
  //    for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
  //      tomo_out.aaafI[iz][iy][ix] = 
  //        P_dot_Q.aaafI[iz][iy][ix] / tomo_out.aaafI[iz][iy][ix];




  // Calculate <P_, P_>

  // Because memory is so precious, we must reuse arrays whenever we can.
  // So we will now use "P" (which used to denote voxel intensity)
  // to store the square of intensity instead
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
        P.aaafI[iz][iy][ix] *= P.aaafI[iz][iy][ix];

  cerr << "\n"
    " ------ Calculating <P_, P_> ... ------\n" << endl;

  //   Original code:  Commented out because it requires too much memory:
  //MrcSimple P_dot_P = tomo_in; //(this takes care of allocating the array)
  //   Memory efficient, ugly code:
  // We're done with tomo_in.  Re-use this array to store P_dot_P temporarilly,
  // but give it a new name, to make the code slightly less confusing to read:
  float ***P_dot_P_aaafI = tomo_in.aaafI;

  w.Apply(tomo_in.header.nvoxels,
          P.aaafI,
          //P_dot_P.aaafI,
          P_dot_P_aaafI,  // <-- store result here
          mask.aaafI,
          false,
          &cerr);

  // Now calculate "c", and "rmse" (sqrt(variance))
  // We can calculate this from the P_dot_Q and P_dot_P arrays we have now.
  // Save the result in tomo_out, and save to a file.

  // First, write out RMSE, root-mean-squared-error  (sqrt(variance))
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
        //     variance = <P_, P_> - <P_,Q_>^2 / |Q_|^2
        //              = <P_, P_> - <P_,Q_>^2 / <Q_,Q_>

        //float variance = (P_dot_P.aaafI[iz][iy][ix]
        float variance = (P_dot_P_aaafI[iz][iy][ix]
                          -
                          SQR(P_dot_Q.aaafI[iz][iy][ix]) / Q_dot_Q);

        //Optional:
        //Compensate for dividing by w.aaafH[][][] by "wpeak" earlier.
        //This enables us to interpret variance as RMSE, (a.k.a. root-
        //mean-squared-error.  The formula above only calculates the
        //"mean" if w_i are normalized, which they were before we divided
        //them all by wpeak.  (Also: This insures that this fast method is
        //equivalent to the "SlowDebug" method commented out above.)
        //(Note: It is important for other reasons that w_i (w.aaafH[][][])
        //       were not normalized until this point.)

        variance *= wpeak;

        if (variance < 0.0)
          variance = 0.0;

        tomo_out.aaafI[iz][iy][ix] = sqrt(variance);
      } //for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
    } //for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
  } //for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)

  tomo_out.FindMinMaxMean();
  string rmse_file_name;
  if ((EndsWith(settings.out_file_name, ".rec")) ||
      (EndsWith(settings.out_file_name, ".mrc")))
    rmse_file_name =
      (settings.out_file_name.substr(0,
                                     settings.out_file_name.length()-4)
       + string("_rmse.mrc"));
  else 
    rmse_file_name = settings.out_file_name + string("_rmse.mrc");
  tomo_out.Write(rmse_file_name);


  // Then write out the "c" parameters to a file
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
        //            c = <P_, Q_> / |Q_|^2
        //              = <P_, Q_> / <Q_,Q_>
        float c = P_dot_Q.aaafI[iz][iy][ix] / Q_dot_Q;
        tomo_out.aaafI[iz][iy][ix] = c;
      } //for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
    } //for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
  } //for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
  tomo_out.FindMinMaxMean();
  //No need to write this file out now. It happens later. Commenting out:
  //tomo_out.Write(settings.out_file_name + string("_c.rec"));
  

  // Optional: Create a (3D) image file containing the filter we are using
  //           multiplied by the weight function.  Q_times_w is the function
  //           that we are actually convolving with our original image.
  //           It's sometimes useful to visualize it, so I save it as a file
  MrcSimple Q_times_w_mrc(Q_times_w.array_size,
                          Q_times_w.aaafH);
  Q_times_w_mrc.Write("Q_times_w.mrc");
  // Do the same for the filter (Q) as well.
  MrcSimple Q_mrc(Q.array_size,
                  Q.aaafH);
  Q_mrc.Write("Q.mrc");


  //{//FOR DEBUGGING ONLY: GENERATE NOISY SIGNAL WITH KNOWN AMPLITUDE AND WIDTH
  //  int Q_rand_halfwidth[3];
  //  Q_rand_halfwidth[0] = floor(2 * w.halfwidth[0]);
  //  Q_rand_halfwidth[1] = floor(2 * w.halfwidth[1]);
  //  Q_rand_halfwidth[2] = floor(2 * w.halfwidth[2]);
  //  Filter3D<float, int>
  //    Q_rand = GenFilterGenGauss3D(settings.width_a,//"a" parameter in formula
  //                                 settings.m_exp,  //"m" parameter in formula
  //                                 Q_rand_halfwidth,
  //                                 static_cast<float*>(nullptr),
  //                                 &cerr);
  //  RANDOM_INIT();
  //  MrcSimple Q_rand_mrc(Q_rand.array_size,
  //                       Q_rand.aaafH);
  //  Q_rand_mrc.Rescale01();
  //  for (int jz=-Q_rand.halfwidth[2]; jz<=Q_rand.halfwidth[2]; jz++) {
  //    for (int jy=-Q_rand.halfwidth[1]; jy<=Q_rand.halfwidth[1]; jy++) {
  //      for (int jx=-Q_rand.halfwidth[0]; jx<=Q_rand.halfwidth[0]; jx++) {
  //        Q_rand_mrc.aaafI[jz][jy][jx] *= 0.5;
  //        Q_rand_mrc.aaafI[jz][jy][jx] += 0.2*RANDOM_GAUSSIAN();
  //      }
  //    }
  //  }
  //  Q_rand_mrc.FindMinMaxMean();
  //  Q_rand_mrc.Write("Q_rand.mrc");
  //}

} //HandleTemplateGGauss()






void
HandleTemplateGauss(const Settings &settings,
                    MrcSimple &tomo_in,
                    MrcSimple &tomo_out,
                    MrcSimple &mask,
                    float voxel_width[3])
{
  // Filter weights, w_i:
  Filter3D<float, int>
    w = GenFilterGenGauss3D(settings.template_background_radius,
                            static_cast<float>(2.0), //settings.template_background_exponent,
                            //template_profile.halfwidth,
                            settings.filter_truncate_ratio,
                            settings.filter_truncate_threshold);

  // GenFilterGenGauss() creates normalized gaussians with integral 1.
  // That's not what we want for the weights, w_i:
  // The weights w_i should be 1 in the viscinity we care about, and 0 outside
  // So, we can undo this normalization by dividing all the w_i weights
  // by their maximum value max(w_i)    (at the central peak, at halfwidth).
  // This will mean the maximum w_i is 1, and the others decay to 0, as we want
  float wpeak = w.aaafH[0][0][0];
  w.MultiplyScalar(1.0 / wpeak);

  // Template function, Q:
  Filter3D<float, int>
    q = GenFilterGenGauss3D(settings.width_a, //"a" parameter in formula
                            static_cast<float>(2.0), //settings.m_exp,
                            w.halfwidth);

  // Q_ = q_ - <q_>
  Filter3D<float, int> Q = q;
  float qave = q.Average(w.aaafH);
  Q.AddScalar(-qave);

  // It's also convenient to scale the voxel intensities in Q_ so that the 
  // range of values, |Q_|^2, (=the standard deviation) is also unity.
  // This makes it meaningful to compare "c" with "variance" directly.
  // (If c/sqrt(variance), then it means that the fit is

  // Now calculate Q_dot_Q
  float Q_dot_Q = Q.SumSqr(w.aaafH);
  Q.MultiplyScalar(1.0/sqrt(Q_dot_Q));

  // We might as well do the same thing for q and qave:
  q.MultiplyScalar(1.0/sqrt(Q_dot_Q));
  qave = q.Average(w.aaafH);
  Q_dot_Q = Q.SumSqr(w.aaafH); // = |Q_|^2  (same everywhere)
  assert(abs(Q_dot_Q - 1.0) < 0.001);  // (should equal 1 but may vary due to roundoff error of 32bit floats)



  cerr << "\n"
    " ------ Calculating the average of nearby voxels: ------\n";
  // P = original image (after subtracting average nearby intensities):

  MrcSimple P = tomo_in; //this will take care of allocating the array

  // First, let's calculate the weighted average voxel intensity in the 
  // source image
  //
  // THIS WORKS, BUT THERE IS A FASTER WAY.  COMMENTING OUT:
  //
  //w.Apply(tomo_in.header.nvoxels,
  //        tomo_in.aaafI,
  //        P.aaafI,   // <-- save result here
  //        mask.aaafI,
  //        true,
  //        &cerr);
  //
  // TRY THIS INSTEAD (works only when w is an ordinary Gaussian):

  float template_background_sigma[3];
  for (int d = 0; d < 3; d++)
    template_background_sigma[d] = (settings.template_background_radius[d]
                                    /
                                    sqrt(3.0));

  ApplyGauss(tomo_in.header.nvoxels,
             tomo_in.aaafI,
             P.aaafI,   // <-- save result here
             mask.aaafI,
             template_background_sigma, //width of Gaussian
             settings.filter_truncate_ratio,
             settings.filter_truncate_threshold,
             true,
             &cerr);

  // Subtract the average value from the image intensity, and store in P:
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
        P.aaafI[iz][iy][ix] = tomo_in.aaafI[iz][iy][ix]-P.aaafI[iz][iy][ix];



  cerr << "\n"
    " ------ Calculating <P_, Q_> / |Q_|^2 ... -------\n" << endl;




  //  REMOVE THIS CRUFT:
  //
  //// Now calculate <P_, Q_>
  //// Recall <P_, Q_> = sum_i (w_i * Q_i * P_i)
  //// precompute w_i * Q_i and call it "Q_times_w_i"
  ////        <P_, Q_> = sum_i (Q_times_w_i * P_i)
  //// Then we can reuse the code we normally use to calculate sum_i(A_i B_i)
  //Filter3D<float, int> Q_times_w = w;
  //for (int jz=-w.halfwidth[2]; jz<=w.halfwidth[2]; jz++)
  //  for (int jy=-w.halfwidth[1]; jy<=w.halfwidth[1]; jy++)
  //    for (int jx=-w.halfwidth[0]; jx<=w.halfwidth[0]; jx++)
  //      Q_times_w.aaafH[jz][jy][jx] =
  //        w.aaafH[jz][jy][jx] * Q.aaafH[jz][jy][jx];
               



  // Calculate P_dot_Q = <P_, Q_> = sum_i (Q_times_w_i * P_i)
  MrcSimple P_dot_Q = tomo_in; //(this takes care of allocating the array)



  // The original code works, but there's a faster way.  Commenting out:
  //Q_times_w.Apply(tomo_in.header.nvoxels,
  //                P.aaafI,
  //                P_dot_Q.aaafI,
  //                mask.aaafI,
  //                false,
  //                &cerr);
  // Instead, use "ApplyGauss()", a fast version which only works for
  // Gaussians.  First we have to calculate the width of that Gaussian...
  float radius_Q_times_w[3];
  float sigma_Q_times_w[3];
  for (int d=0; d < 3; d++) {
    // The product of two Gaussians (widths a, b) is a Gaussian
    // whose width is narrower than the width of either Gaussian
    // e^(-(r/a)^2) * e^(-(r/b)^2) = e^(-(r/R)^2),   R=sqrt(1/(1/a^2+1/b^2))
    radius_Q_times_w[d] =
      sqrt(1.0 / (1.0 / SQR(settings.width_a[d])
                  +
                  1.0 / SQR(settings.template_background_radius[d])));
    sigma_Q_times_w[d] = radius_Q_times_w[d] / sqrt(3.0);
  }
  // We can use the faster "ApplyGauss()" to perform the convolution
  ApplyGauss(P.header.nvoxels, // = tomo_in.header.nvoxels,
             P.aaafI,
             P_dot_Q.aaafI,        // <-- save result here
             mask.aaafI,
             sigma_Q_times_w,
             settings.filter_truncate_ratio,
             settings.filter_truncate_threshold,
             false, // Don't normalize. Gaussian will have peak height=1
             &cerr);
  // The function above convolves the image, P, with a Gaussian whose
  // central peak has a height of 1.  What we want to do instead is convolve
  // it with a Gaussian whose central peak has a height equal to the height
  // of the "Q" gaussian.  This number is in the middle of the Q.aaafH array
  float qpeak = q.aaafH[0][0][0];



  //  REMOVE THIS CRUFT:
  //float q_times_w_peak = Q_times_w.aaafH[0][0][0]
  ////should be the same, since w peaked at 1
  //assert((abs(qpeak-qave)-abs(q_times_w_peak)) < 0.001*abs(q_times_w_peak));


  // Now multiply all of the voxels in the colvoved image (P_dot_Q)
  // by "qpeak" to compensate for using the wrong normalization earlier:
  for(int iz=0; iz<P_dot_Q.header.nvoxels[2]; iz++)
    for(int iy=0; iy<P_dot_Q.header.nvoxels[1]; iy++)
      for(int ix=0; ix<P_dot_Q.header.nvoxels[0]; ix++)
        P_dot_Q.aaafI[iz][iy][ix] *= qpeak;
  // (The speed gains are worth the extra code complexity)

  // We forgot to take into account the fact that Q_ = q_ - <q_>
  // and Q_times_w is the product of Q with another Gaussian.
  // The result is:
  // Q_times_w is not a Gaussian, but the difference of two Gaussians,
  //    one with width "radius_Q_times_w" (which we computed already), and
  //    one with width "settings.template_background_radius"
  // We have convolved P with q_times_w (before subtracting <q_> = "qave")
  // and stored it in P_dot_Q.  Now convolve P with qave*w and subtract
  // it from the result we obtained for P_dot_Q.
  // Again, we can use the faster "ApplyGauss()" to perform the convolution:
  ApplyGauss(P.header.nvoxels, // = tomo_in.header.nvoxels,
             P.aaafI,
             tomo_in.aaafI,        // <-- save result here
             mask.aaafI,
             template_background_sigma,
             settings.filter_truncate_ratio,
             settings.filter_truncate_threshold,
             false,  // Don't normalize. Gaussian will have peak height=1
             &cerr);

  for(int iz=0; iz<P_dot_Q.header.nvoxels[2]; iz++)
    for(int iy=0; iy<P_dot_Q.header.nvoxels[1]; iy++)
      for(int ix=0; ix<P_dot_Q.header.nvoxels[0]; ix++)
        P_dot_Q.aaafI[iz][iy][ix] -= qave * tomo_in.aaafI[iz][iy][ix];



  // Calculate <P_, P_>
  // (This should be easier than calculating <P_, Q_>)

  // Because memory is so precious, we must reuse arrays whenever we can.
  // So we will now use "P" (which used to denote voxel intensity)
  // to store the square of intensity instead
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
        P.aaafI[iz][iy][ix] *= P.aaafI[iz][iy][ix];



  cerr << "\n"
    " ------ Calculating <P_, P_> ... ------\n" << endl;


  //   Original code:  Commented out because it requires too much memory:
  //MrcSimple P_dot_P = tomo_in; //(this takes care of allocating the array)
  //   Memory efficient, ugly code:
  // We're done with tomo_in.  Re-use this array to store P_dot_P temporarilly,
  // but give it a new name, to make the code slightly less confusing to read:
  float ***P_dot_P_aaafI = tomo_in.aaafI;

  // The code below works too, but there's a faster way.  COMMENTING OUT:
  //w.Apply(tomo_in.header.nvoxels,
  //        P.aaafI,
  //        //P_dot_P.aaafI,
  //        P_dot_P_aaafI,
  //        mask.aaafI,
  //        false,
  //        &cerr);
  // Instead, we can use "ApplyGauss()" to perform the convolution quickly

  ApplyGauss(P.header.nvoxels, // = tomo_in.header.nvoxels,
             P.aaafI,
             //P_dot_P.aaafI,
             P_dot_P_aaafI,        // <-- save result here
             mask.aaafI,
             template_background_sigma,  //width of "w" Gaussian
             settings.filter_truncate_ratio,
             settings.filter_truncate_threshold,
             false,  // don't normalize. The "w" weights have a peak of 1
             &cerr);

  // Now calculate "c", and "rmse" (sqrt(variance))
  // We can calculate this from the P_dot_Q and P_dot_P arrays we have now.
  // Save the result in tomo_out, and save to a file.

  // First, write out RMSE, root-mean-squared-error  (sqrt(variance))
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
        //     variance = <P_, P_> - <P_,Q_>^2 / |Q_|^2
        //              = <P_, P_> - <P_,Q_>^2 / <Q_,Q_>

        //float variance = (P_dot_P.aaafI[iz][iy][ix]
        float variance = (P_dot_P_aaafI[iz][iy][ix]
                          -
                          SQR(P_dot_Q.aaafI[iz][iy][ix]) / Q_dot_Q);

        //Optional:
        //Compensate for dividing by w.aaafH[][][] by "wpeak" earlier.
        //This enables us to interpret variance as RMSE, (a.k.a. root-
        //mean-squared-error.  The formula above only calculates the
        //"mean" if w_i are normalized, which they were before we divided
        //them all by wpeak.
        //(Note: It is important for other reasons that w_i (w.aaafH[][][])
        //       were not normalized until this point.)

        variance *= wpeak;

        if (variance < 0.0)
          variance = 0.0;

        tomo_out.aaafI[iz][iy][ix] = sqrt(variance);
      } //for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
    } //for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
  } //for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
  tomo_out.FindMinMaxMean();
  tomo_out.FindMinMaxMean();
  string rmse_file_name;
  if ((EndsWith(settings.out_file_name, ".rec")) ||
      (EndsWith(settings.out_file_name, ".mrc")))
    rmse_file_name =
      (settings.out_file_name.substr(0,
                                     settings.out_file_name.length()-4)
       + string("_rmse.mrc"));
  else 
    rmse_file_name = settings.out_file_name + string("_rmse.mrc");
  tomo_out.Write(rmse_file_name);

  // Then write out the "c" parameters to a file
  for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
    for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
      for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
        //            c = <P_, Q_> / |Q_|^2
        //              = <P_, Q_> / <Q_,Q_>
        float c = P_dot_Q.aaafI[iz][iy][ix] / Q_dot_Q;
        tomo_out.aaafI[iz][iy][ix] = c;
      } //for(int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
    } //for(int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
  } //for(int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
  tomo_out.FindMinMaxMean();
  //No need to write this file out now. It happens later. Commenting out:
  //tomo_out.Write(settings.out_file_name + string("_c.rec"));



  // Optional: Create a (3D) image file containing the filter we are using
  //           multiplied by the weight function.  Q_times_w is the function
  //           that we are actually convolving with our original image.
  //           It's sometimes useful to visualize it, so I save it as a file
  //MrcSimple Q_times_w_mrc(Q_times_w.array_size,
  //                        Q_times_w.aaafH);
  //Q_times_w_mrc.Write("Q_times_w.mrc");
  // Do the same for the filter (Q) as well.
  MrcSimple Q_mrc(Q.array_size,
                      Q.aaafH);
  Q_mrc.Write("Q.mrc");

} //HandleTemplateGauss()

#endif //#ifndef DISABLE_TEMPLATE_MATCHING




void
HandleDistanceToPoints(const Settings &settings,
                       MrcSimple &tomo_in,
                       MrcSimple &tomo_out,
                       MrcSimple &mask,
                       float voxel_width[3])
{
  fstream coords_file;
  coords_file.open(settings.in_coords_file_name, ios::in);
  if (! coords_file)
    throw VisfdErr("Error: unable to open \""+
                   settings.in_coords_file_name +"\" for reading.\n");
  vector<array<int,3> > crds;
  while (coords_file) {
    double x, y, z;
    coords_file >> x;
    coords_file >> y;
    coords_file >> z;
    double ix, iy, iz;
    ix = floor((x / voxel_width[0]) + 0.5);
    iy = floor((y / voxel_width[1]) + 0.5);
    iz = floor((z / voxel_width[2]) + 0.5);
    array<int, 3> ixiyiz;
    ixiyiz[0] = ix;
    ixiyiz[1] = iy;
    ixiyiz[2] = iz;
    crds.push_back(ixiyiz);
  }

  cerr << " ------ calculating distance to points in "
       << settings.in_coords_file_name << " ------\n"
       << endl;

  // At some point I was trying to be as general as possible and allowed
  // for the possibility that voxels need not be cubes (same width x,y,z)
  // Now, I realize that allowing for this possibility would slow down some
  // calculations considerably, so I just assume cube-shaped voxels:
  float voxel_width_ = voxel_width[0];
  assert((voxel_width[0] == voxel_width[1]) &&
         (voxel_width[1] == voxel_width[2]));

  for (int iz=0; iz<tomo_out.header.nvoxels[2]; iz++) {
    cerr << "processing Z-plane " << iz+1 << " / "
         << tomo_out.header.nvoxels[2] << "\n";
    for (int iy=0; iy<tomo_out.header.nvoxels[1]; iy++) {
      for (int ix=0; ix<tomo_out.header.nvoxels[0]; ix++) {
        int rminsq_int = SQR(tomo_out.header.nvoxels[0] +
                             tomo_out.header.nvoxels[1] +
                             tomo_out.header.nvoxels[2]);

        if (mask.aaafI && (mask.aaafI[iz][iy][ix] == 0.0))
          continue;

        // Loop over all of the points and find the nearest one.
        // Calculate the distance to that point, and store that in the tomogram.
        // (Note: This is a terribly inefficient way to do this.
        //        A paint-bucket-like-tool using expanding spheres would
        //        be much faster and visit each voxel only once on average.)
        for (vector<array<int,3> >::const_iterator pxyz=crds.begin();
             pxyz != crds.end();
             pxyz++) {
          int rsqi=SQR(ix-(*pxyz)[0])+SQR(iy-(*pxyz)[1])+SQR(iz-(*pxyz)[2]);
          if (rsqi < rminsq_int)
            rminsq_int = rsqi;
        }
        float rmin = sqrt(rminsq_int * SQR(voxel_width_));
        tomo_out.aaafI[iz][iy][ix] = rmin;
      }
    }
  }
} //HandleDistanceToPoints()




#endif  //#ifndef DISABLE_DOGGXY
