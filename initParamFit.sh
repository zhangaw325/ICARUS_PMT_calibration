# Runs the initial parameter fit and creates a compiled .pdf file
# Assumes that there is only one .csv file in the given directory,
#  which holds the initial parameter information

#!/usr/bin/env bash

echo "Starting fit analysis..."

echo " Warning: All pdf, txt, and output root files from previous fits will be deleted. Proceed? (Y/N)"
read input
if [[ $input == "Y" || $input == "y" ]]; then
    #remove files to preserve versions
    rm *.pdf *fit.txt *gain.root
    #file name formatting
    csvName=`ls *.csv | head -1`
    pdfName="${csvName/csv/pdf}"   
    ROOTSYS=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.16.00/x86_64-ubuntu18-gcc73-opt PATH=$ROOTSYS/bin:$PATH root -l -b <<EOF
        .x FitChargeDistributions_InitParam.C("$csvName");
        .q;
EOF
    pdfunite *.pdf $pdfName
    echo "Analysis finished."
    echo "Output to: $pdfName"
    evince $pdfName &
else
    echo "Analysis cancelled."
fi
