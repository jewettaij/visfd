#!/usr/bin/env bash

test_membrane_detection() {
  cd tests/
    OUT_FNAME="test_image_membrane_planar_55_tv_4_4_clust_1e+09_0.707_0.707_0.707_0.707.rec"
    assertTrue "Failure during membrane detection, tensor-voting, or voxel clustering:  File \"$OUT_FNAME\" not created" "[ -s $OUT_FNAME ]"
    SUM_VOXELS=`bin/sum_voxels/sum_voxels $OUT_FNAME`
    assertTrue "Failure during membrane detection, tensor-voting, or voxel clustering." "[ $SUM_VOXELS -eq 7687 ]"
    rm -rf OUT_FNAME
  cd ../
}

. shunit2/shunit2

