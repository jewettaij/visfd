#include <cassert>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <random>
#include <ctime>              // required for "struct timezone" and
                              // "gettimeofday()" used to set the randomseed
using namespace std;

#ifndef DISABLE_OPENMP
#include <omp.h>       // (OpenMP-specific)
#endif

//#include <boost/math/special_functions/gamma.hpp>

#include <err_report.hpp>
#include <mrc_simple.hpp>
#include <filter3d.hpp>
#include "settings.hpp"


// (Note: For gcc version 4.8.3, you must compile using: g++ -std=c++11)


int main(int argc, char **argv) {
  try {

    Settings settings; // parse the command-line argument list from the shell
    settings.ParseArgs(argc, argv);


    #ifndef DISABLE_OPENMP
    #pragma omp parallel
    {
      int rank, nthr;
      rank = omp_get_thread_num();
      //cerr << "rank=" << rank << endl;
      if (rank == 0) {
        nthr = omp_get_num_threads();
        cerr << "   (Using " << nthr << " threads (cpu cores).  You can change this using the \"-np n\"\n"
             << "    argument, or by setting the OMP_NUM_THREADS environment variable.)" << endl;
      }
    }
    #else
    cerr << " (Serial version)" << endl;
    #endif //#ifndef DISABLE_OPENMP


    // Read the input tomogram
    MrcSimple tomo_in;

    if (settings.in_file_name != "") {
      cerr << "Reading tomogram \""<<settings.in_file_name<<"\"" << endl;
      tomo_in.Read(settings.in_file_name, false);
      // (Note: You can also use "tomo_in.Read(cin);" or "cin >> tomo;")
      tomo_in.PrintStats(cerr);      //Optional (display the tomogram size & format)
    }
    else {
      assert((settings.image_size[0] > 0) &&
             (settings.image_size[1] > 0) &&
             (settings.image_size[2] > 0));
      tomo_in.Resize(settings.image_size);
    }

    // ---- mask ----

    // Optional: if there is a "mask", read that too
    MrcSimple mask;
    if (settings.mask_file_name != "") {
      cerr << "Reading mask \""<<settings.mask_file_name<<"\"" << endl;
      mask.Read(settings.mask_file_name, false);
      if ((mask.header.nvoxels[0] != tomo_in.header.nvoxels[0]) ||
          (mask.header.nvoxels[1] != tomo_in.header.nvoxels[1]) ||
          (mask.header.nvoxels[2] != tomo_in.header.nvoxels[2]))
        throw InputErr("Error: The size of the mask image does not match the size of the input image.\n");
    }

    // ---- Voxel size? ----


    float voxel_width[3] = {1.0, 1.0, 1.0};
    if (settings.voxel_width > 0.0) {
      // Did the user manually specify the width of each voxel?
      voxel_width[0] = settings.voxel_width;
      voxel_width[1] = settings.voxel_width;
      voxel_width[2] = settings.voxel_width;
    }
    else {
      // Otherwise, infer it from the header of the MRC file
      voxel_width[0] = tomo_in.header.cellA[0]/tomo_in.header.nvoxels[0];
      voxel_width[1] = tomo_in.header.cellA[1]/tomo_in.header.nvoxels[1];
      voxel_width[2] = tomo_in.header.cellA[2]/tomo_in.header.nvoxels[2];
      if (settings.voxel_width_divide_by_10) {
        voxel_width[0] *= 0.1;
        voxel_width[1] *= 0.1;
        voxel_width[2] *= 0.1;
      }
      cerr << "voxel width in physical units = ("
           << voxel_width[0] << ", "
           << voxel_width[1] << ", "
           << voxel_width[2] << ")\n";
    }

    if ((voxel_width[0] <= 0.0) ||
        (voxel_width[1] <= 0.0) ||
        (voxel_width[2] <= 0.0))
      throw InputErr("Error in tomogram header: Invalid voxel width(s).\n"
                     "Use the -w argument to specify the voxel width.");

    if ((voxel_width[0] != voxel_width[1]) ||
        (voxel_width[1] != voxel_width[2]))
      throw InputErr("Error in tomogram header: Unequal voxel widths in the x, y and directions.\n"
                     "Use the -w argument to specify the voxel width.");


    if (settings.in_coords_file_name != "") {
      for (int iz = 0; iz < tomo_in.header.nvoxels[2]; iz++)
        for (int iy = 0; iy < tomo_in.header.nvoxels[1]; iy++)
          for (int ix = 0; ix < tomo_in.header.nvoxels[0]; ix++)
            tomo_in.aaafI[iz][iy][ix] = 0.0;
      fstream coords_file;
      coords_file.open(settings.in_coords_file_name.c_str(), ios::in);
      if (! coords_file)
        throw InputErr("Error: unable to open \""+
                       settings.in_coords_file_name +"\" for reading.\n");
      while (coords_file) {
        float x, y, z;
        coords_file >> x;
        coords_file >> y;
        coords_file >> z;
        int ix, iy, iz;
        ix = static_cast<int>(x / voxel_width[0]);
        iy = static_cast<int>(y / voxel_width[1]);
        iz = static_cast<int>(z / voxel_width[2]);
        if (((0 <= ix) && (ix < tomo_in.header.nvoxels[0])) &&
            ((0 <= iy) && (iy < tomo_in.header.nvoxels[1])) &&
            ((0 <= iz) && (iz < tomo_in.header.nvoxels[2])))
          tomo_in.aaafI[iz][iy][ix] = 1.0;
      }
    } //if (settings.in_coords_file_name != "") {



    double vol_total = settings.compartment_volume;
    if (settings.compartment_volume < 0) {
      if (mask.aaafI) {
        for (int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
          for (int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
            for (int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
              vol_total += mask.aaafI[iz][iy][ix];
      }
      else {
        vol_total = (tomo_in.header.nvoxels[0] *
                     tomo_in.header.nvoxels[1] *
                     tomo_in.header.nvoxels[2]);
      }
      vol_total *= (voxel_width[0]*voxel_width[1]*voxel_width[2]);
    }

    if (settings.num_particles < 0) {
      // If the user did not explicitly specify the number of particles
      // (ie "ribosomes") present in the original image, we can attempt
      // to infer it from the image itself.  To do this, I assume that
      // the image is filled with voxels containing 0s everywhere except
      // at the center of where each particle is located, where the voxel is 1.
      // Consequently, the number of non-zero voxels (not excluded by the mask)
      // should equal the sum of all of the voxel entries.
      // (This should also equal the number of non-zero voxels.  However
      //  sometimes these images are themselves slightly blurred to make it
      //  easier to see each particle.  This blurring should not effect the sum)
      double sum = 0.0;
      for (int iz=0; iz<tomo_in.header.nvoxels[2]; iz++) {
        for (int iy=0; iy<tomo_in.header.nvoxels[1]; iy++) {
          for (int ix=0; ix<tomo_in.header.nvoxels[0]; ix++) {
            if (mask.aaafI)
              sum += tomo_in.aaafI[iz][iy][ix] * mask.aaafI[iz][iy][ix];
            else
              sum += tomo_in.aaafI[iz][iy][ix];
          }
        }
      }
      settings.num_particles = sum;
    }


    // DEBUG: The next if-then statement is for debugging only.
    if (settings.randomize_input_image) {
      size_t nparticles = floor(settings.num_particles);
      size_t nvoxels;
      if (mask.aaafI) {
        for (int iz=0; iz<tomo_in.header.nvoxels[2]; iz++)
          for (int iy=0; iy<tomo_in.header.nvoxels[1]; iy++)
            for (int ix=0; ix<tomo_in.header.nvoxels[0]; ix++)
              if (mask.aaafI[iz][iy][ix] != 0.0)
                nvoxels += 1;
      }
      else {
        nvoxels = (tomo_in.header.nvoxels[0] *
                   tomo_in.header.nvoxels[1] *
                   tomo_in.header.nvoxels[2]);
      }
      vector<bool> random_bit_list(nvoxels, false);
      for (size_t i=0; i < nparticles; i++)
        random_bit_list[i] = true;

      long random_seed = settings.random_seed;
      if (random_seed <= 0) {
        random_seed = time(NULL);
        cerr << "(random_seed = " << random_seed << ")" << endl;
      }
      shuffle(random_bit_list.begin(), random_bit_list.end(),
              default_random_engine(random_seed));
      size_t i=0;
      for (int iz = 0; iz < tomo_in.header.nvoxels[2]; iz++) {
        for (int iy = 0; iy < tomo_in.header.nvoxels[1]; iy++) {
          for (int ix = 0; ix < tomo_in.header.nvoxels[0]; ix++) {
            tomo_in.aaafI[iz][iy][ix] = 0.0;
            if ((! mask.aaafI) || mask.aaafI[iz][iy][ix] != 0.0) {
              tomo_in.aaafI[iz][iy][ix] = random_bit_list[i];
              i++;
            }
          }
        }
      }
    } //if (settings.randomize_input_image)



    // First we must allocate a new array to store the filtered image.
    MrcSimple tomo_out = tomo_in;


    for (int i_sig = 0; i_sig < settings.vfSigma.size(); ++i_sig) {

      // ---- main calculation: -----

      //A Gaussian filter is a (weighted) average over a volume of nearby voxels
      //Below we estimate the effective volume of those nearby voxels, by taking
      //advantage of the fact that Gaussian is normalized (it's integral = 1).
      //Normalization effectively reduces the height of the filter by the number
      //of voxels which effectively belong to the filter.
      // Example: a cube filter
      // If we were averaging the brightness of the 27 voxels
      // belonging to a 3*3*3 cube surrounding the voxel of interest
      // then the effective volume of the filter is is 3*3*3=27, and
      // weight of each voxel's contribution to the average is 1/27.  Similarly,
      // for a more general filter, we can interpret the peak height of the
      // filter as 1/volume. Below, we calculate the height of the 3D Gaussian's
      // peak after normalization.  Normally the 3D Gaussian is:
      //   (2*pi*sigma^2)^(-3/2) * exp(-0.5*(x^2+y^2+z^2)/sigma^2)
      // and the peak height is (2*pi*sigma^2)^(-3/2),
      // However this formula does not work for small sigma << 1.
      // For small sigma, you have to normalize by computing a discrete sum
      // (over a finite range).  We do this below (in a non-intuitive way).

      float sigma = settings.vfSigma[i_sig];

      sigma /= voxel_width[0]; //scale by voxel_width (if specified)


      // How wide is the filter window (over what size window do we integrate)?
      if (settings.filter_truncate_ratio <= 0) {
        assert(settings.filter_truncate_threshold > 0.0);
        settings.filter_truncate_ratio = sqrt(-2*log(settings.filter_truncate_threshold));
      }
      int filter_truncate_halfwidth = floor(sigma *
                                            settings.filter_truncate_ratio);

      float gauss_peak_height_3D;

      {
        // GenFilterGauss1D generates a 1D discrete normalized Gaussian (located
        // in the afH[]) member. afH[0] is the height at the central peak. 
        // However, for now I just use them to estimate the height of the 
        // Gaussian peak.  (See below)
        Filter1D<float, int> filter1D =
          GenFilterGauss1D(sigma, filter_truncate_halfwidth);

        // A 3D Gaussian blur can be performed by expoiting the fact that
        // multidimensional Gaussians are "seperable" filters:  You can blur
        // in the X direction, then in the Y direcion, then in the Z direction.
        // The result is the same that you would get from a 3D filter Gaussian
        // which is the product of Gaussians in the X,Y,Z directions.
        // So the peak height is the product of the 3 1D Gaussian peak heights.
        gauss_peak_height_3D =
          (filter1D.afH[0] *
           filter1D.afH[0] *
           filter1D.afH[0]);
      }
    
      // The Gaussian peak's height is 1/volume.

      float volume_gaussian_bin = 1.0 / gauss_peak_height_3D;
      volume_gaussian_bin *= (voxel_width[0]*voxel_width[1]*voxel_width[2]);

      float num_bins = vol_total / volume_gaussian_bin;

      // Now blur the original image (if the user did not do it already)
      if ((! settings.precomputed_gaussian_blur) ||
          (settings.vfSigma.size() > 1)) {

        ApplyGauss(tomo_in.header.nvoxels,
                   tomo_in.aaafI,
                   tomo_out.aaafI, //<-store resulting image here
                   mask.aaafI,
                   sigma,
                   filter_truncate_halfwidth,
                   true,
                   &cerr);

        // Densities everywhere are assumed to be in physical units
        // (ie. of 1/Angstroms^3),  NOT   1/voxels^3
        // The Gaussian blur computes densities in 1/voxels^3
        // To compensate for this, divide the densities by voxel_width^3
        MultiplyScalarArr(static_cast<float>(
                          1.0/(voxel_width[0]*voxel_width[1]*voxel_width[2])),
                          tomo_out.header.nvoxels,
                          tomo_out.aaafI);
      }

      // After bluring the image, find the lowest density:
    
      float extreme_density;
      int afXextreme[3] = {-1, -1, -1};
      float global_minima;
      float global_maxima;

      // Find the voxels with the minima and maxima intensities.
      // Discard global minima or maxima which lie on the boundary
      // The safe, careful way to do this is to only consider
      // voxels which are surrounded by other voxels and are
      // local minima or maxima.
      for (int iz=0; iz < tomo_out.header.nvoxels[2]; iz++) {
        for (int iy=0; iy < tomo_out.header.nvoxels[1]; iy++) {
          for (int ix=0; ix < tomo_out.header.nvoxels[0]; ix++) {

            bool is_local_minima = true;
            bool is_local_maxima = true;

            if (! settings.extrema_on_boundary) {
              if ((ix == 0) || (ix ==  tomo_out.header.nvoxels[0]-1) ||
                  (ix == 0) || (ix ==  tomo_out.header.nvoxels[0]-1) ||
                  (ix == 0) || (ix ==  tomo_out.header.nvoxels[0]-1)) {
                is_local_minima = false;
                is_local_maxima = false;
              }
            }

            float center_val = tomo_out.aaafI[iz][iy][ix];
            for (int jz=-1; jz<=1; jz++) {
              for (int jy=-1; jy<=1; jy++) {
                for (int jx=-1; jx<=1; jx++) {
                  if (! settings.extrema_on_boundary) {
                    if (mask.aaafI && (mask.aaafI[iz+jz][iy+jy][ix+jx] == 0.0))
                    {
                      is_local_minima = false;
                      is_local_maxima = false;
                    }
                  }
                  if (tomo_out.aaafI[iz+jz][iy+jy][ix+jx] <= center_val)
                    is_local_minima = false;
                  if (tomo_out.aaafI[iz+jz][iy+jy][ix+jx] >= center_val)
                    is_local_maxima = false;
                }
              }
            }
            if ((center_val < global_minima) || (afXextreme[0] == -1)) {
              global_minima = center_val;
              if (settings.use_min_density) {
                extreme_density = global_minima;
                afXextreme[0] = ix;
                afXextreme[1] = iy;
                afXextreme[2] = iz;
              }
            }
            if ((center_val > global_maxima) || (afXextreme[0] == -1)) {
              global_maxima = center_val;
              if (! settings.use_min_density) {
                extreme_density = global_maxima;
                afXextreme[0] = ix;
                afXextreme[1] = iy;
                afXextreme[2] = iz;
              }
            }
          } //for (int ix=1; ix < tomo_out.header.nvoxels[0]-1; ix++)
        } //for (int iy=1; iy < tomo_out.header.nvoxels[1]-1; iy++)
      } //for (int iz=1; iz < tomo_out.header.nvoxels[2]-1; iz++)


      // Did we fail to find any local minima or maxima densities?
      if (afXextreme[0] == -1) {
        string extrema_type = "minina";
        if (! settings.use_min_density)
          extrema_type = "maxima";
        stringstream msg_ss;
        msg_ss << "Error: There are no local density " << extrema_type << "\n"
               << "       in the image (at scale sigma = " << sigma << ")\n"
               << "       Aborting...\n";
        throw InputErr(msg_ss.str());
      }

      float ave_density = settings.num_particles / vol_total;


      // Suppose you divide the volume of the region (ie. cell, compartment...)
      // into equal size "bins" of equal size ("volume_gaussian_bin").
      // A fixed number of particles lie within this volume.
      // ...What is the probability that the number of particles in a bin
      //    does not exceed "extreme_density"?
      //  This is given by the (cumulative) Poisson distribution:
      //     prob = lambda^k * exp(-lambda) / k!
      //          = lambda^k * exp(-lambda) / Gamma(k+1)  <-continuous version
      // "k" is the number of particles in this Gaussian-shaped "bin"
      long double k = extreme_density * volume_gaussian_bin;
      // "lambda" is the expected number of particles in this "bin"
      long double lambda = ave_density * volume_gaussian_bin;


      long double prob_cdf_obvious_way = 0.0;
      if (settings.use_min_density) {
        // Calculate the culmulative probability of seeing this many
        // particles in the bin or less.  First try it the obvious way:

        for (long i=0; i <= floor(k); i++)
          prob_cdf_obvious_way += pow(lambda, i) * exp(-lambda) / tgamma(i+1.0);
      }
      else {
      // Calculate the culmulative probability of seeing this many
      // particles in the bin or more.  First try it the obvious way:
        for (long i=0; i < floor(k); i++)
          prob_cdf_obvious_way += pow(lambda, i) * exp(-lambda) / tgamma(i+1.0);
        prob_cdf_obvious_way = 1.0 - prob_cdf_obvious_way;
      }

      cerr << "####################################" << endl;
      cerr << "## DEBUG MESSAGES (PLEASE IGNORE) ##" << endl;
      cerr << "## location (in voxels) = ("
           << afXextreme[0] << ","
           << afXextreme[1] << ","
           << afXextreme[2] << ")" << endl;
      cerr << "## num_in_bin = " << k << endl;
      cerr << "## num_expected_in_bin = " << lambda << endl;
      cerr << "## prob_cdf_obvious_way = " << prob_cdf_obvious_way << endl;
      cerr << "####################################" << endl;

      long double prob_cdf;

      // However, in our case the number of particles in this "bin" is not (k)
      // necessarily an integer.  In that case use the continuous version of
      // (the cumulative distribution) of the Poisson distributio.
      // COMMENTING OUT: This requires the BOOST libraries
      // COMMENTING OUT: if (use_min_density)
      // COMMENTING OUT:   prob_cdf = gamma_q(k+1,lambda);
      // COMMENTING OUT: else
      // COMMENTING OUT:   prob_cdf = 1.0 - gamma_q(k,lambda);
      //    (gamma_q() is defined in boost)

      // The cumulative distribution of the Poisson distribution
      // is also given by the upper incomplete Gamma function
      // COMMENTING OUT: USE THE BOOST LIBRARY
      //double prob_cdf = gamma_q(floor(k)+1, lambda);
      // Unfortunately the  BOOST library can not be easily distributed with
      // this code.  (Even if I throw away most ofthe BOOST code using the
      // "bcp" utility, the remaining BOOST code exceeds the size of this
      // the code for this project by a factor of 50.)
      // Hence we will use the integer approximation ("obvious_way"):

      prob_cdf = prob_cdf_obvious_way;

      // ...Now what is the probability that the number of particles
      //    in ANY OF THE BINS does not exceed "extreme_density"?

      long double prob_total = 1.0 - pow((1.0 - prob_cdf), num_bins);

      double effective_bin_size = (pow(volume_gaussian_bin, 1.0/3) *
                                   voxel_width[0]);
      cout << prob_total
           << " " << extreme_density
           << " " << afXextreme[0]
           << " " << afXextreme[1]
           << " " << afXextreme[2]
           << " " << effective_bin_size
           << endl;

      // Discussion:
      //      What do we expect this number to be?
      // What is the expected value of "prob_total" for a randomly distributed
      // particles in the same number of bins.
      //
      // Answer: 0.5
      //
      // The proof is not specific to the kind of distribution we are using.
      // (a Poisson distribution)
      //
      // Proof:
      //
      // Let "c" be the cumulative probability distrubution for the entire
      // system (in this "prob_total", but it could also be "prob_cdf")
      //
      // Let p(x) be the probability density of measuring x.
      //
      // Let C(x) be the cumulative distribution of measuring something <= x.
      // (in this case the probability that none of the bins have
      //  a density exceeding x)  This means that:
      //
      // C(X) = \int_{-\infty}^X  dx  p(x)
      //
      // Let C^{-1}(c)  denote the inverse of C(x)
      // ==> C^{-1}(c) = x   if   C(x) = c
      //
      // What is the average value of c?
      //
      // <c> = \int_{-infty}^{\infty} * c(x) * p(x) * dx
      //
      //     = \int_0^1         c *  p(x(c)) * dx(c)/dc * dc
      //
      //     = \int_0^1 c * p( C^{-1}(c) ) * (d/dc) C^{-1}(c)  * dc
      //
      //     = \int_0^1 c * p( C^{-1}(c) ) * ( 1 / (d/dx) C( C^{-1}(c) ) ) * dc
      //
      //     = \int_0^1 c * p( C^{-1}(c) ) * ( 1 / p( C^{-1}(c) ) )  * dc
      //
      //     = \int_0^1 c * dc
      //
      //     = 0.5
      //
      // However this number will be lower if you repeat this multiple times
      // (using different bin-sizes ("sigma" values), for example)
      // because by repeating, you give it more opportunities to find
      // a more extreme value.  (Here we scan over a range of "sigma" values.)

    } // for (int i_sig = 0; i_sig < settings.vfSigma.size(); i_sig++) {...


    if ((settings.out_file_name != "") && (settings.vfSigma.size() == 1))
    {
      // Write the image file containing the filtered version for that sigma
      cerr << "writing a tomogram containing the most recently calculated density cloud"
           << endl;

      tomo_out.Write(settings.out_file_name);
      // (You can also use "file_stream << tomo_out;")
    }

  } //try {
  catch (InputErr& e) {
    cerr << "\n" << e.what() << endl;
    exit(1);
  }

} // main()

