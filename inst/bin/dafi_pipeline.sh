#!/usr/bin/env bash

cd ~/work
cwd=$(pwd)
numOfCores=$(grep -c ^processor /proc/cpuinfo)
Label=${PWD##*/}
 run_FCSTrans2TXT.R FCS
 mkdir Preprocessed
 mv *.txt Preprocessed
 mkdir Gated
 cd Gated
 for file in $cwd/Preprocessed/*
 do
 filepath=$(basename $file )
 extension="${filepath##*.}"
 filename="${filepath%.*}"
 echo Processing File:  $filename
 mkdir $filename
 cd $filename
 dafi_gating $file $cwd/config/inclusion.config $cwd/config/exclusion.config 500 500 $numOfCores 2
 cd ..
 done
 parsePercentage.awk */DU*percentage.txt > Batch_percentages.txt
 parsePercentage.awk */DU*events.txt > Batch_events.txt
 cd ..
