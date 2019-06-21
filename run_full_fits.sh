#!/usr/bin/env bash

echo "Starting analysis."

# Loop over result.root files in sets of 3
while read l1; do
	read l2
	read l3
	echo "line1 is $l1"
done < HistogramsToFit.txt

echo "Analysis finished."