#!/usr/bin/env bash

read -p "Enter PMT chimney you wish to compile into spreadsheet form: " chimney
echo $chimney

# find LED off files
offFiles=(${chimney}_PMT_*_afterpulse.txt)
numFiles=${#offFiles[@]}
echo "$numFiles LED off files found"

if [ "$numFiles" -ne "3" ]
then
	echo "Improper number of LED off files found."
	exit 1
fi

i=0
out_index=0
PMT_nums=()

probs_out=()
locs_out=()

probs=()
first_peak=()
first_peak_error=()
first_chi=()
first_ndf=()
second_peak=()
second_peak_error=()
second_chi=()
second_ndf=()

for ind in {0..34}; do
	probs[$ind]="--"
	first_peak[$ind]="--"
	first_peak_error[$ind]="--"
	first_chi[$ind]="--"
	first_ndf[$ind]="--"
	second_peak[$ind]="--"
	second_peak_error[$ind]="--"
	second_chi[$ind]="--"
	second_ndf[$ind]="--"
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

	if [ "$i" -eq "0" ]
	then
		PMT_nums[$i]=$pmt0
		let "j=$i+1"
		PMT_nums[$j]=$pmt1
		let "j=$i+2"
		PMT_nums[$j]=$pmt2
		let "j=$i+3"
		PMT_nums[$j]=$pmt3
	fi
	if [ "$i" -eq "1" ]
	then
		let "j=$i*4"
		PMT_nums[$j]=$pmt0
		let "j=($i*4)+1"
		PMT_nums[$j]=$pmt1
		let "j=($i*4)+2"
		PMT_nums[$j]=$pmt2
		let "j=($i*4)+3"
		PMT_nums[$j]=$pmt3
	fi
	if [ "$i" -eq "2" ]
	then
		let "j=$i*8"
		PMT_nums[$j]=$pmt0
		let "j=($i*8)+1"
		PMT_nums[$j]=$pmt1
		let "j=($i*8)+2"
		PMT_nums[$j]=$pmt2
		let "j=($i*8)+3"
		PMT_nums[$j]=$pmt3
	fi

	((out_index = $i*12 + 2))

	#------------------------------
	# store afterpulses to output
	#------------------------------

	# only highest-voltage PMTs were analyzed

	# read lines from file
	prob_lines=()
	loc_lines=()
	for ind in 0 1 2 3; do
		prob_line_num=0
		((prob_line_num = 2*${ind}+1))
		prob_lines[$ind]=$(sed "${prob_line_num}q;d" $file)

		loc_line_num=0
		((loc_line_num = 2*${ind}+2))
		loc_lines[$ind]=$(sed "${loc_line_num}q;d" $file)
	done

	# read probability lines
	for ind in 0 1 2 3; do
		((place_index = $out_index + $ind*3))
		probs[$place_index]=$(echo ${prob_lines[$ind]} | cut -f5 -d' ')
	done
	# echo $probs

	# read peak lines
	for ind in 0 1 2 3; do
		((place_index = $out_index + $ind*3))
		first_peak[$place_index]=$(echo ${loc_lines[$ind]} | cut -f5 -d' ')
		first_peak_error[$place_index]=$(echo ${loc_lines[$ind]} | cut -f6 -d' ')
		first_chi[$place_index]=$(echo ${loc_lines[$ind]} | cut -f7 -d' ')
		first_ndf[$place_index]=$(echo ${loc_lines[$ind]} | cut -f8 -d' ')
		second_peak[$place_index]=$(echo ${loc_lines[$ind]} | cut -f9 -d' ')
		second_peak_error[$place_index]=$(echo ${loc_lines[$ind]} | cut -f10 -d' ')
		second_chi[$place_index]=$(echo ${loc_lines[$ind]} | cut -f11 -d' ')
		second_ndf[$place_index]=$(echo ${loc_lines[$ind]} | cut -f12 -d' ')
	done

	# increment index
	((i++))
done

# function from StackExchange
function join_by { local IFS="$1"; shift; echo "$*"; }

# print output and finish up
echo "PMT numbers are ${PMT_nums[*]}"

pmt_string=$(join_by _ "${PMT_nums[@]}")
filename="${chimney}_afterpulse_${pmt_string}.txt"

for ind in {0..35}; do
        line_to_append=$(printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" "${probs[$ind]}" "${first_peak[$ind]}" "${first_peak_error[$ind]}"\
        					"${first_chi[$ind]}" "${first_ndf[$ind]}" "${second_peak[$ind]}" "${second_peak_error[$ind]}"\
        					"${second_chi[$ind]}" "${second_ndf[$ind]}")

        if [ "$ind" -eq "0" ]
        then
                echo -e $line_to_append > $filename
        else
                echo -e $line_to_append >> $filename
        fi
done
