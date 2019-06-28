#!/usr/bin/env bash

echo "Beginning creation of waveform analysis files"

i=0;
while read line; do
	$line > ${i}.WaveformDirectories.txt
	echo "$line > ${i}.WaveformDirectories.txt"
	((i++))
done < WaveformDirectories.txt

echo "Done"