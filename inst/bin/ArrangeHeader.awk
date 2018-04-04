#!/usr/bin/awk -f

BEGIN{FS="\t";OFS="\t"} 

FNR == 1 { ++fIndex } 
fIndex == 1{split($0,header,",")} 
fIndex == 2{Swap[$1]=$2;next} 
(fIndex > 2 && FNR ==1){
	for (i=1;i<=NF;i++){
		if($(i) in Swap){
			$(i)=Swap[$(i)]};
			loc[$(i)]=i
		}
	} 
(fIndex > 2){
	printf "%s", $(loc[header[1]]);
	for (i=2;i<=length(header);i++){
		if(header[i] in loc){
			printf "\t%s", $(loc[header[i]])
		}
	};
	printf "\n"
}
