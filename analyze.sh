#!/usr/bin/env bash

#save current directory (or source directory of bash script within the script itself)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

#prompt user for name of directories with .wfm files and store
read -p "Enter directory name you wish to anaylze (ex: A10,1,2,3,4,On): " info

#put the list of directories into an array
shopt -s dotglob
shopt -s nullglob

#save information into variables
chimney=$(echo $info | cut -f1 -d,)
pmt1=$(echo $info | cut -f2 -d,)
pmt2=$(echo $info | cut -f3 -d,)
pmt3=$(echo $info | cut -f4 -d,)
pmt4=$(echo $info | cut -f5 -d,)
led=$(echo $info | cut -f6 -d,)

#find relevant data files based on input

if [ "$led" == "Off" ]
then
    dirArray=(/home/azhang/ICARUS/PMT/Data201905/${chimney}_PMT_${pmt1}_${pmt2}_${pmt3}_${pmt4}_*Off)
    led=false
else
    dirArray=(/home/azhang/ICARUS/PMT/Data201905/${chimney}_PMT_${pmt1}_${pmt2}_${pmt3}_${pmt4}_*On)
    dirArray+=(/home/azhang/ICARUS/PMT/Data201905/${chimney}_PMT_${pmt1}_${pmt2}_${pmt3}_${pmt4}_*ON) 
    led=true
fi
echo -e "PROCESSING: Chimney:$chimney\tPMT1:$pmt1\tPMT2:$pmt2\tPMT3:$pmt3\tPMT4:$pmt4\tLED:$led"

i=1

#loop through each voltage file and process data
for fullDir in "${dirArray[@]}"; do
    #change directory to aiwu database 
    cd $fullDir
    dir="$(basename $fullDir)"
    #ls all of the full-path names of the .wfm files
    #save it to a .txt file in the currentDir location
    find ~+ -type f -name "*.wfm" > $DIR/filelist_$dir.txt

    #move back to original directory
    cd $DIR    
    #run python script
    python ProcessData_cl.py $dir > $DIR/$dir.txt
    
    #fix to change the python output files from  LedON to LedOn
    if [ "$(echo $dir | cut -f8 -d_)" == "LedON" ]
    then
	mv filelist_"$dir".txt filelist_"$(echo $dir | cut -f1-7 -d_)"_LedOn.txt
	mv "$dir"_result.root "$(echo $dir | cut -f1-7 -d_)"_LedOn_result.root
	mv "$dir".txt "$(echo $dir | cut -f1-7 -d_)"_LedOn.txt
    fi
    
    #get voltage info
    volt[i]=$(echo $dir | cut -f7 -d_)
    volt[i]=${volt[i]%?}
    echo -e "\tRaw data processed for ${volt[i]}V"
    ((i++))
done

if [ "$led" == "true" ]
then
root -l -b <<EOF
     .x FitChargeDistributions.C("$chimney",$pmt1,$pmt2,$pmt3,$pmt4,1440,1470,1500,true);
     .q;
EOF
fi
