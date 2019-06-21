#!/usr/bin/env bash

# loop over all waveform directories
echo "Starting analysis"
while read p; do
  # pass directory name to analyze script
  echo "Analyzing $p"
  echo $p | ./analyze_cl.sh
  # echo $p
done < WaveformDirectories.txt

./fix_led_names.sh
