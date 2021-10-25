#!/usr/bin/env bash

SEP=1.1
THRESH=-90   #(Note: This threshold depends on the "-truncate" and/or
              #       "-truncate-threshold" arguments.)

test_blob_detection() {
  cd tests/
    # remove low frequencies from the image (not necessary for this tiny image)
    ../bin/filter_mrc/filter_mrc -w 19.6 -mask test_blob_detect_mask.rec -in test_blob_detect.rec -o test_blob_detect_dog_0_500.rec -dog 0 500

    assertTrue "-dog argument failed.  File test_blob_detect_dog_0_500.rec file was not created" "[ -s test_blob_detect_dog_0_500.rec ]"

    # use clipping to remove extremely bright or dark voxels
    ../bin/filter_mrc/filter_mrc -w 19.6 -mask test_blob_detect_mask.rec -in test_blob_detect_dog_0_500.rec -o test_blob_detect_dog_0_500_cl_-1.3_1.3.rec -cl -1.3 1.3

    assertTrue "-cl argument failed.  File test_blob_detect_dog_0_500_cl_-1.3_1.3.rec was not created" "[ -s test_blob_detect_dog_0_500_cl_-1.3_1.3.rec ]"

    #../bin/filter_mrc/filter_mrc -w 19.6 -mask test_blob_detect_mask.rec -in test_blob_detect_dog_0_500.rec -blob minima test_blobs.txt 160.0 280.0 1.02

    ../bin/filter_mrc/filter_mrc -w 19.6 -mask test_blob_detect_mask.rec -in test_blob_detect.rec -blob minima test_blobs.txt 160.0 280.0 1.01

    assertTrue "blob detection failed.  File test_blobs.txt was not created" "[ -s test_blobs.txt ]"

    ../bin/filter_mrc/filter_mrc -w 19.6 -mask test_blob_detect_mask.rec -in test_blob_detect.rec -discard-blobs test_blobs.txt test_blobs_sep_${SEP}_thresh_${THRESH}.txt -blob-separation ${SEP} -minima-threshold ${THRESH}

    assertTrue "non-max suppression failed.  File test_blobs_sep_${SEP}_thresh_${THRESH}.txt was not created" "[ -s test_blobs_sep_${SEP}_thresh_${THRESH}.txt ]"
    NBLOBS_IN_LIST=`wc test_blobs_sep_${SEP}_thresh_${THRESH}.txt | awk '{print $1}'`
    assertTrue "Blob detection failed. The number of non-redundant blobs should be 2" "[ $NBLOBS_IN_LIST -eq 2 ]"


    # Create an image with each blob represented by a single voxel of
    # brightness 1 surrounded by other voxels of brightness 0

    ../bin/filter_mrc/filter_mrc -w 19.6 -mask test_blob_detect_mask.rec -in test_blob_detect_dog_0_500_cl_-1.3_1.3.rec -out test_blob_detect_results.rec -draw-spheres test_blobs_sep_${SEP}_thresh_${THRESH}.txt -background 0 -foreground 1 -sphere-radii 0

    # Note:
    # To superimpose the spheres onto the original image, get rid of these args:
    # -background 0 -foreground 1 -sphere-radii 0
    # ...and add this argument:
    # -background-auto -background-scale 0.2
    
    assertTrue "visualization failed.  File test_blob_detect_results.rec was not created" "[ -s test_blob_detect_results.rec ]"

    NBLOBS_IN_IMAGE=`../bin/sum_voxels/sum_voxels -mask test_blob_detect_mask.rec test_blob_detect_results.rec`
    
    assertTrue "visualization failed: number of non-zero voxels != number of blobs" "[ $NBLOBS_IN_LIST -eq $NBLOBS_IN_IMAGE ]"

    # Test supervised learning of blob threshold parameter (single image)
    ../bin/filter_mrc/filter_mrc -w 19.6 -mask test_blob_detect_mask.rec -in test_blob_detect.rec -discard-blobs test_blobs.txt test_blobs_sep_${SEP}_SUPERVISED.txt -blob-separation ${SEP} -auto-thresh score -supervised test_supervised_pos.txt test_supervised_neg.txt >& test_log_e.txt

    assertTrue "\"-supervised\" non-max suppression failed: File test_blobs_sep_${SEP}_SUPERVISED.txt was not created" "[ -s test_blobs_sep_${SEP}_SUPERVISED.txt ]"

    NBLOBS_SUPERVISED_SINGLE=`wc test_blobs_sep_${SEP}_SUPERVISED.txt | awk '{print $1}'`

    assertTrue "\"-supervised\" non-max suppression failed: Number of remaining blobs == 0" "[ $NBLOBS_SUPERVISED_SINGLE -gt 0 ]"

    THRESH_SUPERVISED_SINGLE=`grep 'threshold upper bound:' < test_log_e.txt | awk '{print $4}'`

    assertTrue "-supervised failed: threshold not found" "[ -n $THRESH_SUPERVISED_SINGLE ]"

    assertTrue "-supervised failed: threshold is non-sensical" "[ $THRESH_SUPERVISED_SINGLE != inf ]"

    assertTrue "-supervised failed: threshold is non-sensical" "[ $THRESH_SUPERVISED_SINGLE != -inf ]"

    # Test supervised learning of blob threshold parameter (multiple images)

    # set up the input files
    ../bin/filter_mrc/filter_mrc -w 19.6 -mask test_blob_detect_mask.rec -in test_blob_detect.rec -discard-blobs test_blobs.txt test_blobs_sep_${SEP}.txt -blob-separation ${SEP}

    echo "test_supervised_pos.txt test_supervised_neg.txt test_blobs_sep_${SEP}.txt" > test_supervised_multi.txt
    echo "test_supervised_pos.txt test_supervised_neg.txt test_blobs_sep_${SEP}.txt" >> test_supervised_multi.txt

    # run the learning algorithm to choose the threshold
    ../bin/filter_mrc/filter_mrc -w 19.6 -in test_blob_detect.rec -auto-thresh score -supervised-multi test_supervised_multi.txt >& test_log_e.txt

    THRESH_SUPERVISED_MULTI=`grep 'threshold upper bound:' < test_log_e.txt | awk '{print $4}'`

    assertTrue "-supervised-multi failed: threshold not found" "[ -n $THRESH_SUPERVISED_MULTI ]"

    # Since we used the same input files we did before for a single image
    # (duplicated twice), we should expect to get the same threshold

    assertTrue "-supervised-multi failed: threshold not consistent" "[ $THRESH_SUPERVISED_MULTI == $THRESH_SUPERVISED_SINGLE ]"

    # cleanup

    rm -rf test_blob_detect_dog_0_500.rec test_blob_detect_dog_0_500_cl_-1.3_1.3.rec test_blobs*.txt test_blob_detect_results.rec test_supervised_multi.txt test_log_e.txt deleteme*

    cd ../
}

. shunit2/shunit2
