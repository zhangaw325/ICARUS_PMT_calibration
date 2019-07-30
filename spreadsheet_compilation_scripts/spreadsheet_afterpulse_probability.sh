#!/usr/bin/env bash

read -p "Enter PMT chimney you wish to compile into spreadsheet form: " chimney
echo $chimney

# find LED on files
onFiles=(${chimney}_PMT_*_LedOn_output.txt)
numFiles=${#onFiles[@]}
echo "$numFiles LED on files found"

if [ "$numFiles" -ne "9" ]
then
	echo "Improper number of LED on files found."
	exit 1
fi

i=0
out_index=0
PMT_nums=()
out_nums=()

for ind in {0..29}; do
	out_nums[$ind]=0
done

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
	then	
		let "out_index=$i*4"
		let "j=4*$i/3"
		PMT_nums[$j]=$pmt0
		let "j=$j+1"
		PMT_nums[$j]=$pmt1
		let "j=$j+1"
		PMT_nums[$j]=$pmt2
		let "j=$j+1"
		PMT_nums[$j]=$pmt3
	fi

	if [ "$imod" -eq "1" ]
	then
		let "out_index=($i*4)-3"
	fi

	if [ "$imod" -eq "2" ]
	then
		let "out_index=($i*4)-6"
	fi

	#------------------------------
	# store dark pulses to output
	#------------------------------

	# read lines from file
	line0=$(sed '7q;d' $file)
	line1=$(sed '8q;d' $file)
	line2=$(sed '9q;d' $file)
	line3=$(sed '10q;d' $file)

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
echo "Pulse values are ${out_nums[*]}"

pmt_string=$(join_by _ "${PMT_nums[@]}")
filename="${chimney}_afterpulseprobs_${pmt_string}"

printf "%s\n" "${out_nums[@]}" > ${filename}.txt
