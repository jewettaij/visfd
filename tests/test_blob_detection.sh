#!/usr/bin/env bash


test_blob_detection() {
  cd tests/
    # remove low frequencies from the image (not necessary for this tiny image)
    ../bin/filter_mrc/filter_mrc -w 19.6 -i test_blob_detection.rec -o test_blob_detection_dog_0_500.rec -dog 0 500

    assertTrue "-dog argument failed.  File test_blob_detection_dog_0_500.rec file was not created" "[ -s test_blob_detection_dog_0_500.rec ]"

    # use clipping to remove extremely bright or dark voxels
    ../bin/filter_mrc/filter_mrc -w 19.6 -i test_blob_detection_dog_0_500.rec -o test_blob_detection_dog_0_500_cl_-1.3_1.3.rec -cl -1.3 1.3

    assertTrue "-cl argument failed.  File test_blob_detection_dog_0_500_cl_-1.3_1.3.rec was not created" "[ -s test_blob_detection_dog_0_500_cl_-1.3_1.3.rec ]"

    ../bin/filter_mrc/filter_mrc -w 19.6 -in test_blob_detection_dog_0_500.rec -blob minima test_blobs.txt 160.0 280.0 1.01

    assertTrue "blob detection failed.  File test_blobs.minima.txt was not created" "[ -s test_blobs.minima.txt ]"

    ../bin/filter_mrc/filter_mrc -w 19.6 -in test_blob_detection.rec -discard-blobs test_blobs.txt test_blobs_minima_sep_0.9.txt -blob-separation 0.9 -minima-threshold -80

    assertTrue "non-max suppression failed.  File test_blobs_minima_sep_0.9.txt was not created" "[ -s test_blobs_minima_sep_0.9.txt ]"

    ../bin/filter_mrc/filter_mrc -w 19.6 -in test_blob_detection.rec -out test_blob_detection_results.rec -spheres test_blobs_minima_sep_0.9.txt -sphere-shell-ratio 1 -sphere-background-scale 0.1 -sphere-radii 20.0 
    
    assertTrue "visualization failed.  File test_blob_detection_results.rec was not created" "[ -s test_blob_detection_results.rec ]"
    
    rm -rf test_blob_detection_dog_0_500.rec test_blob_detection_dog_0_500_cl_-1.3_1.3.rec test_blobs*.txt test_blob_detection_results.rec

  cd ../
}

. shunit2/shunit2
