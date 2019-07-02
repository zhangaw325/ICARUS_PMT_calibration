#!/usr/bin/env bash

# loop over all waveform directories
echo "Starting analysis"
while read p; do
  # pass directory name to analyze script
  echo "Analyzing $p"
  echo $p | ./analyze_cl.sh
  # echo $p
done < ${1}.WaveformDirectories.txt # %1 refers to first argument of command

./fix_led_names.sh

echo "Deleting filelist files"
rm filelist*PMT*.txt

echo "Done"
