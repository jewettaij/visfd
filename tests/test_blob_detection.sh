#!/usr/bin/env bash

!!!!!!!!!! AS OF 2019-2-27, this is not ready yet !!!!!!!!!!

test_blob_detection() {
  cd tests/
    # remove low frequencies from the image (not necessary for this tiny image)
    ../bin/filter_mrc/filter_mrc -w 19.6 -i test_blob_detection.rec -o test_blob_detection_dog_0_500.rec -dog 0 500

    assertTrue "-dog argument failed.  File test_blob_detection_dog_0_500.rec file was not created" "[ -s test_blob_detection_dog_0_500.rec ]"

    # use clipping to remove extremely bright or dark voxels
    ../bin/filter_mrc/filter_mrc -w 19.6 -i test_blob_detection_dog_0_500.rec -o test_blob_detection_dog_0_500_cl_-1.3_1.3.rec -cl -1.3 1.3

    assertTrue "-cl argument failed.  File test_blob_detection_dog_0_500_cl_-1.3_1.3.rec was not created" "[ -s test_blob_detection_dog_0_500_cl_-1.3_1.3.rec ]"

    ../bin/filter_mrc/filter_mrc -w 19.6 -in test_blob_detection_dog_0_500.rec -blob test_blobs 160.0 280.0 1.01

    assertTrue "blob detection failed.  File test_blobs.minima.txt was not created" "[ -s test_blobs.minima.txt ]"

    bin/filter_mrc/filter_mrc -w 19.6 -in test_blob_detection_dog_0_500.rec -discard-blobs test_blobs.minima.txt test_blobs_minima_sep_0.9.txt -blob-separation 0.9 -minima-threshold -90

    
    THRESH=-55
    
   rm -rf test_image_fluct.rec
  cd ../
}

. shunit2/shunit2
