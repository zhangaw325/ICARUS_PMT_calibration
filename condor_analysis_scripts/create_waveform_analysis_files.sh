#!/usr/bin/env bash

echo "Beginning creation of waveform analysis files"

i=0;
while read line; do
	mkdir ${i}_files
	echo $line > ${i}_files/${i}.WaveformDirectories.txt
	echo "$line > ${i}_files/${i}.WaveformDirectories.txt"
	cp analyze_cl.sh ${i}_files/
	cp ProcessData_cl_threshold.py ${i}_files/
	cp fix_led_names.sh ${i}_files/
	((i++))
done < WaveformDirectories.txt

echo "Done"
