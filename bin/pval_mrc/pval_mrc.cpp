#include <cassert>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

#include <boost/math/special_functions/gamma.hpp>

#include <err_report.h>
#include <mrc_simple.h>
#include <filter3d.h>
#include "settings.h"

// (Note: For gcc version 4.8.3, you must compile using: g++ -std=c++11)


int main(int argc, char **argv) {
  try {

    Settings settings; // parse the command-line argument list from the shell
    settings.ParseArgs(argc, argv);

    // Read the input tomogram
    cerr << "Reading tomogram \""<<settings.in_file_name<<"\"" << endl;
    MrcSimple tomo_in;
    tomo_in.Read(settings.in_file_name, false);
    // (Note: You can also use "tomo_in.Read(cin);" or "cin >> tomo;")
    tomo_in.PrintStats(cerr);      //Optional (display the tomogram size & format)

    // ---- mask ----

    // Optional: if there is a "mask", read that too
    MrcSimple mask;
    if (settings.mask_file_name != "") {
      cerr << "Reading mask \""<<settings.mask_file_name<<"\"" << endl;
      mask.Read(settings.mask_file_name, false);
      if ((mask.header.nvoxels[0] != tomo_in.header.nvoxels[0]) ||
          (mask.header.nvoxels[1] != tomo_in.header.nvoxels[1]) ||
          (mask.header.nvoxels[2] != tomo_in.header.nvoxels[2]))
        throw InputErr("Error: The size of the mask does not match the size of the tomogram.\n");
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



    
    float vol_total = settings.compartment_volume;
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
        float gauss_peak_height_3D =
          (filter1D.afH[0] *
           filter1D.afH[0] *
           filter1D.afH[0]);
      }
    
      // The Gaussian peak's height is 1/volume.

      float volume_gaussian_bin = 1.0 / gauss_peak_height_3D;

      float num_bins = vol_total / volume_gaussian_bin;

      // Now blur the original image (if the user did not do it already)
      if ((! settings.precomputed_gaussian_blur) ||
          (settings.vfSigma.size() > 1)) {

        ApplyGauss3D(tomo_in.header.nvoxels,
                     tomo_in.aaafI,
                     tomo_out.aaafI,
                     mask.aaafI,
                     sigma,
                     filter_truncate_halfwidth,
                     true,
                     &cerr);
      }

      // After bluring the image, find the lowest density:
    
      float extreme_density;
      if (settings.use_min_density)
        extreme_density = _MinArr(tomo_out.header.nvoxels,
                                  tomo_out.aaafI,
                                  mask.aaafI);
      else
        extreme_density = _MaxArr(tomo_out.header.nvoxels,
                                  tomo_out.aaafI,
                                  mask.aaafI);
      

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
      float k = extreme_density * volume_gaussian_bin;
      // "lambda" is the expected number of particles in this "bin"
      float lambda = ave_density * volume_gaussian_bin;


      double prob_cdf_obvious_way = 0.0;
      if (settings.use_min_density)
        // Calculate the culmulative probability of seeing this many
        // particles in the bin or less.  First try it the obvious way:
        for (long i=0; i <= floor(k); i++)
          prob_cdf_obvious_way += pow(lambda, i) * exp(-lambda) / tgamma(i+1.0);
      else {
      // Calculate the culmulative probability of seeing this many
      // particles in the bin or more.  First try it the obvious way:
        for (long i=0; i < floor(k); i++)
          prob_cdf_obvious_way += pow(lambda, i) * exp(-lambda) / tgamma(i+1.0);
        prob_cdf_obvious_way = 1.0 - prob_cdf_obvious_way;
      }

      cerr << "####################################" << endl;
      cerr << "## DEBUG MESSAGES (PLEASE IGNORE) ##" << endl;
      cerr << "## num_in_bin = " << k << endl;
      cerr << "## prob_cdf_obvious_way = " << prob_cdf_obvious_way << endl;
      cerr << "####################################" << endl;

      double prob_cdf;

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

      double prob_total = 1.0 - pow((1.0 - prob_cdf), num_bins);

      cout << pow(volume_gaussian_bin, 1.0/3) * voxel_width[0]
           << " " << prob_total << endl;

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

    } // for (int i_sig = 0; i_sig < vfSigma.size(); i_sig++) {...

  } //try {
  catch (InputErr& e) {
    cerr << "\n" << e.what() << endl;
    exit(1);
  }

} // main()

