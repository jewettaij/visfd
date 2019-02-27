#!/usr/bin/env bash

test_fluctuation_filter() {
  cd tests/
    assertTrue "test_image_fluct.rec file not created" "[ -s test_image_fluct.rec ]"
    rm -rf test_image_fluct.rec
  cd ../
}

. shunit2/shunit2
