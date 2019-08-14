#!/usr/bin/env bash

echo "Starting analysis."

echo "$ROOTSYS"
echo "$PATH"

# Loop over result.root files in sets of 3
while read l1; do

	#save info into variables
	chimney=$(echo $l1 | cut -f1 -d_)
	pmt0=$(echo $l1 | cut -f3 -d_)
	pmt1=$(echo $l1 | cut -f4 -d_)
	pmt2=$(echo $l1 | cut -f5 -d_)
	pmt3=$(echo $l1 | cut -f6 -d_)
	volt1Raw=$(echo $l1 | cut -f7 -d_)
	ledRaw=$(echo $l1 | cut -f8 -d_)

	volt1=$(echo $volt1Raw | cut -f1 -dV)
	led=${ledRaw:3}

	read l2

	volt2Raw=$(echo $l2 | cut -f7 -d_)
	volt2=$(echo $volt2Raw | cut -f1 -dV)

	read l3

	volt3Raw=$(echo $l3 | cut -f7 -d_)
	volt3=$(echo $volt3Raw | cut -f1 -dV)

	echo "Fitting Chimney: $chimney; PMTs: $pmt0,$pmt1,$pmt2,$pmt3; Voltages: $volt1,$volt2,$volt3"

	ROOTSYS=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.16.00/x86_64-ubuntu18-gcc73-opt PATH=$ROOTSYS/bin:$PATH root -l -b <<EOF
        .x FitChargeDistributions.C("$chimney",$pmt0,$pmt1,$pmt2,$pmt3,$volt1,$volt2,$volt3,true);
        .q;
EOF
done < HistogramsToFit.txt

echo "Analysis finished."
