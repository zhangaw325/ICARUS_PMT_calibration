#!/usr/bin/env bash

read -p "Enter PMT chimney you wish to compile into spreadsheet form: " chimney
echo $chimney

# find LED off files
offFiles=(${chimney}_PMT_*_LedOff_output.txt)
numFiles=${#offfiles[@]}
echo "$numFiles LED off files found"

if [ "$numFiles" -ne "9" ]
do
	echo "Improper number of LED off files found."
	exit 1
fi

i=0
out_index=0
PMT_nums=()
out_nums=()

for file in "${offFiles[@]}"; do
	echo $file

	# record PMT number ordering
	chimney=$(echo $file | cut -f1 -d_)
	pmt0=$(echo $file | cut -f3 -d_)
	pmt1=$(echo $file | cut -f4 -d_)
	pmt2=$(echo $file | cut -f5 -d_)
	pmt3=$(echo $file | cut -f6 -d_)
	voltRaw=$(echo $file | cut -f7 -d_)

	voltage=$(echo $voltRaw | cut -f1 -dV)

	let "imod=$i%3"
	if [ "$imod" -eq "0" ]
		let "$out_index=$i"
		PMT_nums[$i]=pmt0
		let "j=$i+1"
		PMT_nums[$j]=pmt1
		let "j=$i+2"
		PMT_nums[$j]=pmt2
		let "j=$i+3"
		PMT_nums[$j]=pmt3
	fi

	#------------------------------
	# store dark pulses to output
	#------------------------------

	# read lines from file
	line0=$(sed '3q;d' $file)
	line1=$(sed '4q;d' $file)
	line2=$(sed '5q;d' $file)
	line3=$(sed '6q;d' $file)

	pulse0=$(echo $line0 | cut -f2 -d$':')
	pulse1=$(echo $line1 | cut -f2 -d$':')
	pulse2=$(echo $line2 | cut -f2 -d$':')
	pulse3=$(echo $line3 | cut -f2 -d$':')

	# output in proper order
	out_nums[$out_index]=$pulse0
	((out_index+=3))
	out_nums[$out_index]=$pulse1
	((out_index+=3))
	out_nums[$out_index]=$pulse2
	((out_index+=3))
	out_nums[$out_index]=$pulse3
	((out_index+=3))

	# increment index
	((i++))
done

# function from StackExchange
function join_by { local IFS="$1"; shift; echo "$*"; }

# print output and finish up
echo "PMT numbers are ${PMT_nums[*]}"

pmt_string=$(join_by _ "${PMT_nums[@]}")
filename="${chimney}_darkpulses_${pmt_string}"

printf "%s\n" "${out_nums[@]}" > ${filename}.txt
