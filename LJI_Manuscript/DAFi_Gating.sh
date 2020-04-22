#!/usr/bin/env bash

Project="$PWD"
Label=${PWD##*/}
CustomPlotConfigFile=$Project"/config/customplot.config"
HeaderList=$Project"/config/header.lst"
HeaderReplaceFile=$Project"/config/header.txt"

export OMP_NUM_THREADS=8
mkdir $Project/Gated_TXT
mkdir $Project/Plots_PNG

cd $Project/Gated_TXT
for f in $Project/Preprocessed_TXT/*.txt
do
	filepath=$(basename $f )
	extension="${filepath##*.}"
	filename="${filepath%.*}"
	echo Gating File:  $filename
	mkdir $filename
	cd $filename
	$Project/bin/dafi_gating $f $Project/config/$Label".config" $Project/config/$Label"_exclude.config" 200 500 8 2
	cd ..
done
$Project/bin/parsePercentage.awk */DU*percentage.txt > Batch_percentages.txt 
$Project/bin/parsePercentage.awk */DU*events.txt > Batch_events.txt 
cd ..
if [ -f "$CustomPlotConfigFile" ]
then
	while read p; do
		PlotLabel=$(echo $p | cut -f1 -d ';')
		CustomConfig=$(echo $p | cut -f2 -d ';')
		echo Processing Custom Plots for $CustomConfig
		mkdir $Project"/Plots_PNG/"$PlotLabel
		time $Project/bin/DAFi-plotting.py -f 'Gated_TXT/*/fl*.txt' -s 0.5 -p 8 --gridshow --titleshow --sort $CustomConfig
		mv *.png $Project"/Plots_PNG/"$PlotLabel
	done < $CustomPlotConfigFile
fi
echo Processing Gating Plots
time $Project/bin/DAFi-plotting.py -f 'Gated_TXT/*/fl*.txt' --config $Project/config/$Label".config" -s 0.5 -p 24 --gridshow --titleshow --sizeofaxislabels 14 --titleshow --sizeofpoplabels 10 --showmultigates --showparent --sort
mkdir $Project/Plots_PNG/Gates
mv *.png $Project/Plots_PNG/Gates
