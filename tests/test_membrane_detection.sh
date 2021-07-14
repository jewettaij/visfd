#!/usr/bin/env bash

test_membrane_detection() {
    cd tests/
    OUT_FNAME_BASE="test_image_membrane_surface_55_tv_4_4_clust_1e+09_45deg"
    OUT_FNAME_REC=${OUT_FNAME_BASE}.rec
    OUT_FNAME_NORMALS=${OUT_FNAME_BASE}.ply
    ../bin/filter_mrc/filter_mrc -w 19.2 -in test_image_membrane.rec -out ${OUT_FNAME_REC} -surface minima 55 -tv 4 -tv-angle-exponent 4 -bin 2 -save-progress test_image_membrane
    ../bin/filter_mrc/filter_mrc -w 19.2 -in test_image_membrane.rec -out ${OUT_FNAME_REC} -surface minima 55 -tv 4 -tv-angle-exponent 4 -bin 2 -load-progress test_image_membrane -connect 1e+09 -connect-angle 30 -surface-normals-file ${OUT_FNAME_NORMALS} -select-cluster 1 >& test_log_e.txt
    assertTrue "Failure during membrane detection, tensor-voting, or voxel clustering:  File \"${OUT_FNAME_REC}\" not created" "[ -s ${OUT_FNAME_REC} ]"
    N_CLUSTERS=`grep 'Number of clusters found:' < test_log_e.txt | awk '{print $5}'`
    assertTrue "Failure: No membrane clusters detected." "[ $N_CLUSTERS -gt 0 ]"
    # Sum the number of voxels in the largest surface.
    # These are voxels in the image we just created whose brightness equals 1.
    SUM_VOXELS=`../bin/sum_voxels/sum_voxels -thresh4 0.98 0.99 1.01 1.02 "${OUT_FNAME_REC}"`
    assertTrue "Failure during membrane detection. File \"${OUT_FNAME_REC}\" has a weird checksum." "[ $SUM_VOXELS -gt 50 ]"
    rm -rf "${OUT_FNAME_REC}" "${OUT_FNAME_NORMALS}" test_log_e.txt *tensor_?.rec
  cd ../
}

. shunit2/shunit2
