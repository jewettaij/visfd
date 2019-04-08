#!/usr/bin/env bash

SEP=1
THRESH=-90

test_blob_detection() {
  cd tests/
    # remove low frequencies from the image (not necessary for this tiny image)
    ../bin/filter_mrc/filter_mrc -w 19.6 -mask test_blob_detect_mask.rec -i test_blob_detect.rec -o test_blob_detect_dog_0_500.rec -dog 0 500

    assertTrue "-dog argument failed.  File test_blob_detect_dog_0_500.rec file was not created" "[ -s test_blob_detect_dog_0_500.rec ]"

    # use clipping to remove extremely bright or dark voxels
    ../bin/filter_mrc/filter_mrc -w 19.6 -mask test_blob_detect_mask.rec -i test_blob_detect_dog_0_500.rec -o test_blob_detect_dog_0_500_cl_-1.3_1.3.rec -cl -1.3 1.3

    assertTrue "-cl argument failed.  File test_blob_detect_dog_0_500_cl_-1.3_1.3.rec was not created" "[ -s test_blob_detect_dog_0_500_cl_-1.3_1.3.rec ]"

    #../bin/filter_mrc/filter_mrc -w 19.6 -mask test_blob_detect_mask.rec -in test_blob_detect_dog_0_500.rec -blob minima test_blobs.txt 160.0 280.0 1.01

    ../bin/filter_mrc/filter_mrc -w 19.6 -mask test_blob_detect_mask.rec -in test_blob_detect.rec -blob minima test_blobs.txt 160.0 280.0 1.01

    assertTrue "blob detection failed.  File test_blobs.txt was not created" "[ -s test_blobs.txt ]"

    ../bin/filter_mrc/filter_mrc -w 19.6 -mask test_blob_detect_mask.rec -in test_blob_detect.rec -discard-blobs test_blobs.txt test_blobs_sep_${SEP}_thresh_${THRESH}.txt -blob-separation ${SEP} -minima-threshold ${THRESH}

    assertTrue "non-max suppression failed.  File test_blobs_sep_${SEP}_thresh_${THRESH}.txt was not created" "[ -s test_blobs_sep_${SEP}_thresh_${THRESH}.txt ]"

    # Create an image with each blob represented by a single voxel of
    # brightness 1 surrounded by other voxels of brightness 0
    ../bin/filter_mrc/filter_mrc -w 19.6 -mask test_blob_detect_mask.rec -in test_blob_detect.rec -out test_blob_detect_results.rec -spheres test_blobs_sep_${SEP}_thresh_${THRESH}.txt -sphere-shell-ratio 0.1 -sphere-background 0 -sphere-foreground 1 -sphere-radii 0
    # Note
    # To superimpose the spheres onto the original image, get rid of these args:
    # -sphere-background 0 -sphere-foreground 1 -sphere-radii 0
    # ...and add this argument:
    # -sphere-background-scale 0.2
    
    assertTrue "visualization failed.  File test_blob_detect_results.rec was not created" "[ -s test_blob_detect_results.rec ]"

    NBLOBS_IN_LIST=`wc test_blobs_sep_${SEP}_thresh_${THRESH}.txt | awk '{print $1}'`
    NBLOBS_IN_IMAGE=`../bin/sum_voxels/sum_voxels -mask test_blob_detect_mask.rec test_blob_detect_results.rec`
    
    assertTrue "visualization failed: number of non-zero voxels != number of blobs" "[ $NBLOBS_IN_LIST -eq $NBLOBS_IN_IMAGE ]"

    rm -rf test_blob_detect_dog_0_500.rec test_blob_detect_dog_0_500_cl_-1.3_1.3.rec test_blobs*.txt test_blob_detect_results.rec

  cd ../
}

. shunit2/shunit2
