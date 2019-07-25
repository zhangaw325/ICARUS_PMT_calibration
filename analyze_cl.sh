#!/usr/bin/env bash

#save current directory (or source directory of bash script within the script itself)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
#echo $DIR

# dirDataName="USER INPUT"
#prompt user for name of directory with .wfm files and store
read -p "Enter directory name you wish to anaylze: " dirDataName

#change directory to aiwu database
#ls all of the full-path names of the .wfm files
#save it to a .txt file in the currentDir location
cd /home/azhang/ICARUS/PMT/Data201905/$dirDataName
find ~+ -type f -name "*.wfm" > $DIR/filelist_$dirDataName.txt
#echo $DIR/filelist_$dirDataName.txt 

cd $DIR
#echo $DIR

#run python script
LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.16.00/x86_64-ubuntu18-gcc73-opt/lib/ PYTHONPATH=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.16.00/x86_64-ubuntu18-gcc73-opt/lib python ProcessData_cl_threshold.py $dirDataName > ${dirDataName}_output.txt

#print message if successful
echo "$dirDataName.root created"
