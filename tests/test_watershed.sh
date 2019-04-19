#!/usr/bin/env bash

SIGMA=120

test_watershed() {
    cd tests/

    # 1D examples:
    # find the number of maxima in the image exceeding a certain threshold
    ../bin/filter_mrc/filter_mrc -w 1 -in test_1d_example.rec -find-maxima test_1d_example_maxima1.txt -maxima-threshold 1200
    N_MAXIMA_STEP1=`wc test_1d_example_maxima1.txt | awk '{print $1}'`
    # generate an image of solid non-overlapping white spheres at each maxima
    ../bin/filter_mrc/filter_mrc -w 1 -in test_1d_example.rec -out test_1d_example_spheres.rec -draw-spheres test_1d_example_maxima1.txt -diameters 3 -foreground 1 -background 0 -spheres-shell-ratio 1
    # find maxima in this image (which has many voxels with the same brightness)
    ../bin/filter_mrc/filter_mrc -w 1 -in test_1d_example_spheres.rec -find-maxima test_1d_example_maxima2.txt -maxima-threshold 0.5
    N_MAXIMA_STEP2=`wc test_1d_example_maxima2.txt | awk '{print $1}'`
    assertTrue "Failure: -find-maxima fails on images having multiple voxels with the same brightness." "[ $N_MAXIMA_STEP1 -eq $N_MAXIMA_STEP2 ]"


    # 3D examples:
    IN_FNAME_BASE="test_blob_detect"
    OUT_FNAME_BASE="test_image_watershed"
    OUT_FNAME_REC=${OUT_FNAME_BASE}.rec
    ../bin/filter_mrc/filter_mrc -w 19.2 -mask ${IN_FNAME_BASE}_mask.rec -i ${IN_FNAME_BASE}.rec -o ${IN_FNAME_BASE}_gauss_${SIGMA}.rec -gauss ${SIGMA}
    assertTrue "Failure during Gaussian blur:  File \"test_blob_detect_gauss_${SIGMA}.rec\" not created" "[ -s test_blob_detect_gauss_${SIGMA}.rec ]"

    # 1D examples:

    # Test -find-minima
    ../bin/filter_mrc/filter_mrc -w 19.2 -mask ${IN_FNAME_BASE}_mask.rec -i ${IN_FNAME_BASE}_gauss_${SIGMA}.rec -find-minima ${IN_FNAME_BASE}_gauss_${SIGMA}_minima.txt -o ${IN_FNAME_BASE}_gauss_${SIGMA}_minima.rec
    assertTrue "Failure: -find-minima failed.  File \"${IN_FNAME_BASE}_gauss_${SIGMA}_minima.txt\" not created" "[ -s ${IN_FNAME_BASE}_gauss_${SIGMA}_minima.txt ]"
    assertTrue "Failure: -find-minima failed to make an image.  File \"${IN_FNAME_BASE}_gauss_${SIGMA}_minima.rec\" not created" "[ -s ${IN_FNAME_BASE}_gauss_${SIGMA}_minima.txt ]"
    # Count the number of dark minima in the image
    N_MINIMA=`wc ${IN_FNAME_BASE}_gauss_${SIGMA}_minima.txt | awk '{print $1}'`
    N_MINIMA_IMAGE=`../bin/print_mrc_stats/print_mrc_stats ${IN_FNAME_BASE}_gauss_${SIGMA}_minima.rec | grep "minimum brightness" | awk '{print -$3}'`
    assertTrue "Failure: The number of minima reported by -find-minima is not consistent with the image that was created." "[ $N_MINIMA -eq $N_MINIMA_IMAGE ]"

    # Test watershed segmentation
    ../bin/filter_mrc/filter_mrc -w 19.2 -mask ${IN_FNAME_BASE}_mask.rec -i ${IN_FNAME_BASE}_gauss_${SIGMA}.rec -out ${OUT_FNAME_REC} -watershed minima >& test_log_e.txt
    assertTrue "Failure during watershed segmentation:  File \"${OUT_FNAME_REC}\" not created" "[ -s ${OUT_FNAME_REC} ]"
    N_BASINS=`grep 'Number of basins found: ' < test_log_e.txt | awk '{print $5}'`
    assertTrue "Failure: No basins found." "[ $N_BASINS -gt 0 ]"
    N_BASINS_IMAGE=`../bin/print_mrc_stats/print_mrc_stats ${OUT_FNAME_REC} | grep "maximum brightness" | awk '{print $3}'`
    assertTrue "Failure: The number of watershed basins in the image is not consistent with the message printed to the terminal." "[ $N_BASINS -eq $N_BASINS_IMAGE ]"
    assertTrue "Failure: The number of watershed basins does not equal the number of minima." "[ $N_BASINS -eq $N_MINIMA ]"


    # Now invert the image brightnesses (photo negative, dark <--> bright)

    ../bin/filter_mrc/filter_mrc -w 19.2 -mask ${IN_FNAME_BASE}_mask.rec -i ${IN_FNAME_BASE}_gauss_${SIGMA}.rec -out ${IN_FNAME_BASE}_gauss_${SIGMA}_inv.rec -invert

    # Test -find-maxima
    ../bin/filter_mrc/filter_mrc -w 19.2 -mask ${IN_FNAME_BASE}_mask.rec -i ${IN_FNAME_BASE}_gauss_${SIGMA}_inv.rec -find-maxima ${IN_FNAME_BASE}_gauss_${SIGMA}_inv_maxima.txt -o ${IN_FNAME_BASE}_gauss_${SIGMA}_inv_maxima.rec
    assertTrue "Failure: -find-maxima failed.  File \"${IN_FNAME_BASE}_gauss_${SIGMA}_maxima.txt\" not created" "[ -s ${IN_FNAME_BASE}_gauss_${SIGMA}_inv_maxima.txt ]"
    assertTrue "Failure: -find-maxima failed to make an image.  File \"${IN_FNAME_BASE}_gauss_${SIGMA}_inv_maxima.rec\" not created" "[ -s ${IN_FNAME_BASE}_gauss_${SIGMA}_inv_maxima.txt ]"

    N_MAXIMA=`wc ${IN_FNAME_BASE}_gauss_${SIGMA}_inv_maxima.txt | awk '{print $1}'`
    N_MAXIMA_IMAGE=`../bin/print_mrc_stats/print_mrc_stats ${IN_FNAME_BASE}_gauss_${SIGMA}_inv_maxima.rec | grep "maximum brightness" | awk '{print $3}'`
    assertTrue "Failure: The number of minima does not equal the number of maxima of the image negative." "[ $N_MINIMA -eq $N_MAXIMA ]"
    assertTrue "Failure: The number of maxima reported by -find-maxima is not consistent with the image that was created." "[ $N_MAXIMA -eq $N_MAXIMA_IMAGE ]"
    assertTrue "Failure: The number of maxima reported by -find-maxima is not consistent with the image that was created." "[ $N_MAXIMA -eq $N_BASINS ]"

    # Test watershed segmentation on the image negative:
    ../bin/filter_mrc/filter_mrc -w 19.2 -mask ${IN_FNAME_BASE}_mask.rec -i ${IN_FNAME_BASE}_gauss_${SIGMA}_inv.rec -out ${OUT_FNAME_REC} -watershed maxima >& test_log_e.txt
    assertTrue "Failure during watershed segmentation:  File \"${OUT_FNAME_REC}\" not created" "[ -s ${OUT_FNAME_REC} ]"
    # Count the number of bright maxima in the inverted (negative) image
    N_BASINS_INV=`grep 'Number of basins found: ' < test_log_e.txt | awk '{print $5}'`
    assertTrue "Failure: Watershed is inconsistent when applied to images with inverted brightness." "[ $N_BASINS -eq $N_BASINS_INV ]"

    # Test "-connect" to see if it behaves like "-watershed maxima"
    
    ../bin/filter_mrc/filter_mrc -w 19.2 -mask ${IN_FNAME_BASE}_mask.rec -i ${IN_FNAME_BASE}_gauss_${SIGMA}_inv.rec -out ${OUT_FNAME_REC} -connect 36.75 >& test_log_e.txt
    N_BASINS_INV=`grep 'Number of clusters found: ' < test_log_e.txt | awk '{print $5}'`
    assertTrue "Failure: Either \"-connect\" argument or the \"-invert\" argument is behaving in an new and unnexpected way, or the extrema finding algorithm is failing" "[ $N_BASINS_INV -eq 2 ]"

    # Delete temporary files:
    rm -rf "${OUT_FNAME_REC}" test_blob_detect_gauss_* test_log_e.txt test_1d_example_*

  cd ../
}

. shunit2/shunit2

