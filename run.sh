#!/bin/bash

if test -n "$9"; then
  nMin="$1"
  nMax="$2"
  dataPointNum="$3"
  J_ratios_path="$4"
  Ts_path="$5"
  isBeta="$6"
  start="$7"
  end="$8"
  flags="$9"
  saveTo_path="${10}"
else
  nMin="6"
  nMax="12"
  dataPointNum="200"
  J_ratios_path="/home/mmaschke/BA_Code/Data/args/J_ratios.txt"
  Ts_path="/home/mmaschke/BA_Code/Data/args/Ts.txt"
  isBeta="0"
  start="0"
  end="3"
  flags="1000000"
  saveTo_path="/home/mmaschke/BA_Code/Data"
fi
./cmake-build-release/SpinChainED $nMin $nMax $dataPointNum $J_ratios_path $Ts_path $isBeta $start $end $flags $saveTo_path
#python3 ./src/plotting.py $nMin $nMax $dataPointNum $J_ratios_path $Ts_path $flags $saveTo_path
