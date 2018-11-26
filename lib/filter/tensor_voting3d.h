

template<class RealNum, class Integer>


class TV3Dslow::public Filter3D {

private:

  RealNum ***aaaf3Displacement[3];
  RealNum *af3Displacement[3];

public:

  TV3Dslow():Filter3D() {
    af3Displacement = NULL;
    aaaf3Displacement = NULL;
  }

  void SetSigma(RealNum sigma, RealNum filter_cutoff_ratio=2.5)
  {
    RealNum sigmas[3] = {sigma, sigma, sigma};
    // Precompute the Gaussian as a function of r, store in ::aaafH
    // ("GenFilterGenGauss3D()" is defined in filter3d.h)
    *(static_cast<Filter3D*>(this)) =
      GenFilterGenGauss3D(sigmas, 2, filter_cutoff_ratio);
    
    AllocNormalized();
    PrecalcNormalized(aaaf3Displacement, halfwidth);
  }

  ~TV3Dslow() {
    DeallocNormalized();
  }

  RealNum
  TVApplyStickToVoxel(Integer ix,
                      Integer iy,
                      Integer iz,
                      Integer const size_source[3],
                      RealNum const *const *const *aaafSource, //!< saliency (score) of each voxel (usually based on Hessian eigenvalues)
                  //CompactMultiChannelImage3D<RealNum> *pvSource, //!< vector associated with each voxel
                      RealNum const *const *const *aaa3fSource[3], //!< vector associated with each voxel
                      RealNum ***aaa33fDest[3][3], //!< votes will be collected here
                      RealNum const *const *const *aaafMaskSource,  //!< ignore voxels in source where mask==0
                      bool detect_curves_not_surfaces = false,
                      RealNum *pDenominator = NULL) const

  {
    RealNum g = 0.0;
    RealNum denominator = 0.0;

    for (Integer jz=-halfwidth[2]; jz<=halfwidth[2]; jz++) {
      Integer iz_jz = iz-jz;
      if ((iz_jz < 0) || (size_source[2] <= iz_jz))
        continue;

      for (Integer jy=-halfwidth[1]; jy<=halfwidth[1]; jy++) {
        Integer iy_jy = iy-jy;
        if ((iy_jy < 0) || (size_source[1] <= iy_jy))
          continue;

        for (Integer jx=-halfwidth[0]; jx<=halfwidth[0]; jx++) {
          Integer ix_jx = ix-jx;
          if ((ix_jx < 0) || (size_source[0] <= ix_jx))
            continue;

          RealNum filter_val = aaafH[jz][jy][jx];

          if (aaafMaskSource) {
            filter_val *= aaafMaskSource[iz_jz][iy_jy][ix_jx];
            if (filter_val == 0.0)
              continue;
          }
          //Note: The "filter_val" also is needed to calculate
          //      the denominator used in normalization.
          //      It is unusual to use a mask unless you intend
          //      to normalize the result later, but I don't enforce this

          RealNum radial_contribution =
            filter_val * aaafSource[iz_jz][iy_jy][ix_jx];

          RealNum const *r = aaaf3Displacement[iz_jz][iy_jy][ix_jx];
          RealNum const *n = aaa3fSource[jz][jy][jx];

          //
          //   .
          //   :\
          //   : \
          // ->:  \<- theta
          //   :   \
          //   :    \ 
          //   :     \
          //   :      \
          //   :       \  
          //   ^        \
          //   |         \ 
          // n |          \
          //   |        .-' r
          //   |,,..--''
          //  0
          //
          // "n" = the direction of the stick tensor.  If the object being 
          //       detected is a surface, n is perpendicular to that surface.
          //       If the object is a curve, then n points tangent to the curve.
          // "r" = the position of the vote receiver relative to the voter
          //
          // "theta" = the angle of r relative to the plane perpendicular to n.
          //           (NOT the angle of r relative to n.  This is a confusing
          //            convention, but this is how it is normally defined.)
          
          RealNum sintheta = DotProduct3(r, n); //(sin() not cos(), see diagram)
          RealNum sinx2 = sintheta * 2.0;
          RealNum sin2 = costheta * costheta;
          RealNum cos2 = 1.0 - sin2;
          RealNum angle_depencece2 = cos2;
          if (detect_curves_not_surfaces) {
            // If we are detecting 1-dimensional curves (polymers, etc...)
            // instead of 2-dimensional surfaces (membranes, ...)
            // then the stick direction is assumed to be along the 
            // of the curve, (as opposed to perpendicular to the 2D surface).
            // As such
            angle_dependence2 = sin2;
          }
          RealNum decay_function = radial_contribution;
          switch(n) {
          case 2:
            decay_function *= angle_dependence2;
            break;
          case 4:
            decay_function *= angle_dependence2 * angle_dependence2;
            break;
          case:
            //RealNum angle_dependence = sqrt(angle_dependence2);
            decay_function *= pow(angle_dependence2, 0.5*n);
            break;
          }

          RealNum n_rotated[3];
          for (Integer d=0; d<3; d++) {
            if (detect_curves_not_surfaces)
              n_rotated[d] = n[d] - sinx2*r[d];
            else
              n_rotated[d] = sinx2*r[d] - n[d]
          }

          RealNum tensor_vote[3][3];
          for (Integer di=0; di<3; di++) {
            for (Integer dj=0; dj<3; dj++) {
              if (di <= dj) {
                tensor_vote[di][dj] = (decay_function *
                                       n_rotated[d1] * n_rotated[d2]);
                aaa33fDest[iz][iy][ix][di][dj] += tensor_vote[di][dj];
              }
            }
          }

          if (pDenominator)
            denominator += filter_val;
        }
      }
    }

    if (pDenominator)
      *pDenominator = denominator;

    return g;
  } // TVApplyToVoxel()





  void
  _TVDenseStickSlow(Integer const image_size[3], //!< source image size
                    RealNum const *const *const *aaafSource, //!< saliency (score) of each voxel (usually based on Hessian eigenvalues)
                  //CompactMultiChannelImage3D<RealNum> *pvSource, //!< vector associated with each voxel
                    RealNum const *const *const *aaa3fSource[3], //!< vector associated with each voxel
                    RealNum ***aaa33fDest[3][3], //!< votes will be collected here
                    RealNum const *const *const *aaafMaskSource,  //!< ignore voxels in source where mask==0
                    RealNum const *const *const *aaafMaskDest,  //!< don't cast votes wherever mask==0
                    RealNum ***aaafDenominator = NULL,
                    RealNum sigma,  //!< Gaussian width of influence
                    Integer exponent, //!< angle dependence
                    RealNum truncate_ratio=2.5,  //!< how many sigma before truncating?
                    bool detect_curves_not_surfaces = false,
                    ostream *pReportProgress = NULL  //!< print progress to the user?
                    )
  {
    assert(aaafSource);
    assert(aaa3fSource);
    assert(aaa33fDest);

    Integer truncate_halfwidth = floor(sigma * truncate_ratio);

    //assert(pv);
    //assert(pV->nchannels() == 3);
    //pV->Resize(image_size, aaafMaskSource, pReportProgress);

    if (pReportProgress)
      *pReportProgress << "  progress: processing plane#" << endl;

    // The mask should be 1 everywhere we want to consider, and 0 elsewhere.
    // Multiplying the density in the tomogram by the mask removes some of 
    // the voxels from consideration later on when we do the filtering.
    // (Later, we will adjust the weight of the average we compute when we
    //  apply the filter in order to account for the voxels we deleted now.)

    for (Integer iz=0; iz<size_source[2]; iz++) {

      if (pReportProgress)
        *pReportProgress << "  " << iz+1 << " / " << size_source[2] << "\n";

      #pragma omp parallel for collapse(2)
      for (Integer iy=0; iy<size_source[1]; iy++) {

        for (Integer ix=0; ix<size_source[0]; ix++) {

          // Calculate the effect of the filter on
          // the voxel located at position ix,iy,iz

          if ((aaafMaskDest) && (aaafMaskDest[iz][iy][ix] == 0.0)) {
            for (Integer di=0; di<3; di++)
              for (Integer dj=0; dj<3; dj++)
                aaafDest[iz][iy][ix][di][dj] = 0.0;
            continue;
          }

          TVApplyStickToVoxel(ix, iy, iz,
                              size_source,
                              aaafSource,
                              aaaf3Source,
                              aaafMaskSource,
                              (aaafDenominator
                               ? &(aaafDenominator[iz][iy][ix])
                               : NULL));
        }
      }
    }
  } //_TVDenseStickSlow()

private:

  void DeallocNormalized() {
    if (aaaf3Displacement) {
      //shift pointers back to normal
      aaaf3Displacement -= halfwidth[2];
      for (Integer iz = 0; iz < array_size[2]; iz++) {
        aaaf3Displacement[iz] -= halfwidth[1];
        for (Integer iy = 0; iy < array_size[1]; iy++) {
          aaaf3Displacement[iz][iy] -= halfwidth[0];
        }
      }
      Dealloc3D(array_size,
                &af3Displacement,
                &aaaf3Displacement);
    }
  }
  void AllocNormalized() {
    if (aaaf3Displacement)
      Dealloc3D(array_size,
                &af3Displacement,
                &aaaf3Displacement);
    Alloc3D(array_size,
            &af3Displacement,
            &aaaf3Displacement);
    //shift pointers to enable indexing from i = -halfwidth .. +halfwidth
    aaaf3Displacement += halfwidth[2];
    for (Integer iz = 0; iz < array_size[2]; iz++) {
      aaaf3Displacement[iz] += halfwidth[1];
      for (Integer iy = 0; iy < array_size[1]; iy++) {
        aaaf3Displacement[iz][iy] += halfwidth[0];
      }
    }
  }

  void PrecalcNormalized(RealNum ***aaaf3Displacement[3]) {
    // pre-compute the normalized radius unit vector
    for (Integer iz = -halfwidth[2]; iz <= halfwidth[2]; iz++) {
      for (Integer iy = -halfwidth[1]; iy <= halfwidth[1]; iy++) {
        for (Integer ix = -halfwidth[0]; ix <= halfwidth[0]; ix++) {
          RealNum length = sqrt(ix*ix + iy*iy + iz*iz);
          if (length == 0)
            length = 1.0;
          aaaf3Displacement[iz][iy][ix][0] = ix / length;
          aaaf3Displacement[iz][iy][ix][1] = iy / length;
          aaaf3Displacement[iz][iy][ix][2] = iz / length;
        }
      }
    }
  }

}; // class TV3Dslow
