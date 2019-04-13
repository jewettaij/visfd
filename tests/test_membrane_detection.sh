#!/usr/bin/env bash

test_membrane_detection() {
    cd tests/
    OUT_FNAME_BASE="test_image_membrane_surface_55_tv_4_4_clust_1e+09_0.707_0.707_0.707_0.707"
    OUT_FNAME_REC=${OUT_FNAME_BASE}.rec
    OUT_FNAME_NORMALS=${OUT_FNAME_BASE}.ply
    ../bin/filter_mrc/filter_mrc -w 19.2 -i test_image_membrane.rec -out ${OUT_FNAME_REC} -surface minima 55 -surface-tv 4 -surface-tv-angle-exponent 4 -connect-vector-saliency 0.707 -connect-vector-neighbor 0.707 -connect-tensor-saliency 0.707 -connect-tensor-neighbor 0.707 -connect 1e+09 -surface-normals-file ${OUT_FNAME_NORMALS} -select-cluster 1 >& test_log_e.txt
    assertTrue "Failure during membrane detection, tensor-voting, or voxel clustering:  File \"${OUT_FNAME_REC}\" not created" "[ -s ${OUT_FNAME_REC} ]"
    N_CLUSTERS=`grep 'Number of clusters found:' < test_log_e.txt | awk '{print $6}'`
    assertTrue "Failure: No membrane clusters detected." "[ $N_CLUSTERS -gt 0 ]"
    SUM_VOXELS=`../bin/sum_voxels/sum_voxels "${OUT_FNAME_REC}"`
    assertTrue "Failure during membrane detection. File \"${OUT_FNAME_REC}\" has a weird checksum." "[ $SUM_VOXELS -gt 6000 ]"
    rm -rf "${OUT_FNAME_REC}" "${OUT_FNAME_NORMALS}" test_log_e.txt
  cd ../
}

. shunit2/shunit2

