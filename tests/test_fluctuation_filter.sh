#!/usr/bin/env bash

test_fluctuation_filter() {
  cd tests/
    ../bin/filter_mrc/filter_mrc -i test_image_membrane.rec -mask-rect 1 14 2 14 2 14 -out test_image_fluct.rec -fluct 60

    assertTrue "test_image_fluct.rec file not created" "[ -s test_image_fluct.rec ]"
    rm -rf test_image_fluct.rec
  cd ../
}

. shunit2/shunit2
