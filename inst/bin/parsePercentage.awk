#!/usr/bin/awk -f
# LJI_Control_Donor_conf?_pop6/Tasks/*/DUnSup_pop_percentage.txt

BEGIN{
	IFS="\t";
	header="";
}

{
	split(FILENAME, path, "/");
	currentName=path[1];
	if(!files[currentName]++){header=header"\t"currentName};
	if(matrix[$1]==""){
		matrix[$1]=$2;			
	}else{
		matrix[$1]=matrix[$1]"\t"$2}
}

END{
	
	print header;
	n=length(matrix)-1;
	for (i=1; i<=n; i++){
		print i"\t"matrix[i];
	}
}
