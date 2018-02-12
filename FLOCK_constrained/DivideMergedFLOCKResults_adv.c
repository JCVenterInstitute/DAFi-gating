/*********************************************************************************************************************************
This program needs the following inputs:
a) the FLOCK results from a file merged from the filtered data of a set of files
b) The set of filtered files
c) the ORIGINAL set of files before filtering
d) An integer specifying which cell population to output for visualization (e.g., by DAG)

This program outputs:
a) a percentage table (number of events in each cell population divided by the number of events in the ORIGINAL file)
with rows as the cell populations (number of cell populations is the number of rows) 
and the columns the file_IDs (number of files is the number of columns)
b) a folder contains the cell population for each of the filtered file (in csv format for DAG visualization of each population)
c) a folder contains the flock_results.txt for each of the original file (make the data discarded the white background (cluster_ID=1) for Python visualization of the overview of populations)

The use of the program is:
DivideMergedFLOCKResults_adv flock_results.txt folder_filtered_data folder_original_data Population_ID

In the output, a TXT table named folder_percentage.txt will be created.
A folder folder_name_population_ID will be created. Inside the folder you will get
File1Name_filt.csv
File2Name_filt.csv

A folder folder_name_flock_results will be created. Inside the folder you will get
File1_flock_results.txt
File2_flock_results.txt


It is important to note that the filtered file can be empty. In such a case, all percentage values will be zero, the population for that
specific file will also be an empty file.

The utility of including the original file is only to get the size of the original file for calculating the percentage

The way to map the original file and the filtered file is to rely on the filename, the filtered file must be named as originalname_filt.txt
while the original file is originalname.txt

The way to identify the filtered data in the merged flock_results.txt is relying on two lines of the same values.
That is, when the filtered file has less than 2 lines, all the output values will be zero.

There needs a step to read through the merged flock_results.txt to get the number of clusters so that the memory can be allocated for
the output percentage table.

Date: April 22, 2016 
Author: Yu "Max" Qian, Ph.D.
Contact: mqian@jcvi.org or qianyu_cs@yahoo.com

***********************************************************************************************************************************/


#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
//#include <unistd.h> //this is for Unix/Linux
#include <assert.h>
#include <direct.h>

#define LINE_LEN 2048
#define WORD_LEN 512
#define DEBUG 1
#define MAX_VALUE 1000000000
#define WHITE_COLOR 1

void getfileinfo(FILE *f_src, int *file_Len, int *num_dm, char *name_string) //this function can read in the situation of the FLOCK result file, the file needs to have a header
{
  char src[LINE_LEN];
  char current_name[WORD_LEN];
  char prv;

  int num_rows=0;
  int num_columns=0;
  int ch='\n';
  int prev='\n';
  int time_pos=0;
  int i=0;
  int j=0;

  src[0]='\0';
  fgets(src, LINE_LEN, f_src);
 
 
  name_string[0]='\0';
  current_name[0]='\0';
  prv='\n';

  while ((src[i]==' ') || (src[i]=='\t')) //skip space and tab characters before the header
    i++;

  while ((src[i]!='\0') && (src[i]!='\n') && (src[i]!='\r')) //repeat until the end of the line
	{

	  current_name[j]=src[i];
		
      if (((src[i]==',') && (prv!=',')) || ((src[i]=='\t') && (prv!='\t')))//a complete word
        {
          current_name[j]='\0';
			
          
          num_columns++; 
          
          strcat(name_string,current_name); 

		  if (src[i]==',')
			strcat(name_string,",");
		  else
			strcat(name_string,"\t");
            

          current_name[0]='\0';
          j=0;			
        }		
		
     if ((src[i]!=',') && (src[i]!='\t'))
        j++;
		
      prv=src[i];
      i++;
	} //end of while
	
  //name_string[j]='\0';
  if ((prv!=',') && (prv!='\t'))//the last one hasn't been retrieved
    {
      current_name[j]='\0';
      
      num_columns++;
      strcat(name_string,current_name);

	 //strcat(name_string,",");
      
     
    }
  
  //start computing # of rows

  while ((ch = fgetc(f_src))!= EOF )
    {
      if (ch == '\n')
        {
          ++num_rows;
        }
      prev = ch;
    }
  if (prev!='\n') //if the last row doesn't have a new line character, we still need to count it in
    ++num_rows;
	
  *file_Len=num_rows;
  *num_dm=num_columns; 

  
}


/************************************* Read the source file into uncomp_data **************************************/
void readsource(FILE *f_src, int file_Len, int num_dm, int **uncomp_data)
{
  
  int index=0;

  int i=0;
  int j=0;
  int t=0;

 
  char src[LINE_LEN];
  char xc[WORD_LEN];

  src[0]='\0';
  fgets(src,LINE_LEN, f_src); //skip the first line about parameter names

  while (!feof(f_src) && (index<file_Len)) //index = 0, 1, ..., file_Len-1
    {
      src[0]='\0';	    
      fgets(src,LINE_LEN,f_src);
      i=0;
     
						    
      for (t=0;t<num_dm;t++) 
      {
          xc[0]='\0';
          j=0;
          while ((src[i]!='\0') && (src[i]!='\n') && (src[i]!=' ') && (src[i]!=',') && (src[i]!='\t'))
          {
              xc[j]=src[i];
              i++;
              j++;
          }
		
          xc[j]='\0';	    
          i++;

          uncomp_data[index][t]=atoi(xc);
       }	
      
      	
      index++;     	
      //fprintf(fout_ID,"%s",src);
    } //end of while
	
  /*if (DEBUG)
    {
      printf("the last line of the source data is:\n");
      for (j=0;j<num_dm;j++)
        printf("%d ",uncomp_data[index-1][j]);
      printf("\n");
    }
	*/
}


void main (int argc, char **argv)
{
  
  char tmpbuf[WORD_LEN];
  char directory_name[WORD_LEN];
  char orig_directory[WORD_LEN];
  char f_name[WORD_LEN];
  char pop_id_char[WORD_LEN];
  char file_orig_name[WORD_LEN];
  char csvfoldername[WORD_LEN];
  char orig_resultfoldername[WORD_LEN];
  char result_file_name[WORD_LEN];
  char percentage_name[WORD_LEN];
  char src[LINE_LEN];
  char para_name_string[LINE_LEN];
  char temp_para_name_string[LINE_LEN];
  char note_file_list[LINE_LEN];
  char command_line[LINE_LEN];

  FILE *f_flock_result; //input file of merged flock result
  FILE *f_file_list; //file list in the folder
  FILE *f_src; 
  FILE *f_orig_src;

  FILE *f_csvresult; //the output file for each population in a filtered data file
  FILE *f_pcnt; //the final percentage value output file
  FILE *f_results; //the flock result file output for each original file
  
  int specified_population_id=0; //the population id that needs to be output

  int num_files=0; //number of files in the filtered/original folder
  int file_index=0;
  int temp=0;
  int data_kept=0;
 
  int file_Len=0; //number of events in the flock_results.txt
  int temp_file_Len=0;  //temp number of events in the filtered file
  int temp_file_size=0;  //temp number of events in the original file

  int num_dm=0; //number of dimensions in the flock_results.txt
  int temp_num_dm=0; //temp number of dimensions in the original/filtered data, which should be equal to num_dm-2
 
  int num_clusters_identified=0; //total number of clusters in the flock_results.txt

  int i=0;
  int j=0;
  int t=0;
  int p=0;
  int file_ID=0;
  int total=0;

  int found=0;
  int fileISempty=0;
  
  char **file_name;

  int **result_data;
  int **temp_data;
  int **orig_temp_data;
  int **new_temp_data;


  double **percentage_value; //the percentage table with num_clusters_identified * num_files

  int *cluster_id;
  int *size_file; //number of events in the original file, size_file[num_files]
  int *size_population; //for a filtered data file, the number of events of each population/cluster, cluster_size[num_clusters_identified]
  

  _strtime_s(tmpbuf,sizeof(char)*WORD_LEN);
  printf( "Starting time:\t\t\t\t%s\n", tmpbuf );
  _strdate_s(tmpbuf,sizeof(char)*WORD_LEN);
  printf( "Starting date:\t\t\t\t%s\n", tmpbuf );

  if (argc!=5)
  {
	  fprintf(stderr,"Incorrect number of input parameters!\n"); 
      fprintf(stderr,"usage:\n");
      fprintf(stderr,"DivideMergedFLOCKResults_adv FLOCK_Result_File FilteredFileFolderName OriginalFileFolderName Population_ID\n");
	  fprintf(stderr,"The number of files and the file names under both directories must be exactly the same (filtered file should be origname_filt.txt). The first folder is the filtered data; the other set is the original files.\n");
	  exit(0);
  }


  f_flock_result=fopen(argv[1],"r");

  directory_name[0]='\0'; //filtered_data;
  strcpy(directory_name,argv[2]);

  orig_directory[0]='\0'; //original data
  strcpy(orig_directory,argv[3]);

  pop_id_char[0]='\0';
  strcpy(pop_id_char,argv[4]);
  specified_population_id=atoi(argv[4]);

 
  getfileinfo(f_flock_result, &file_Len, &num_dm, para_name_string); //file_Len is the total number of events in the merged file flock_results.txt

  if (DEBUG)
  {
      printf("header of FLOCK result file is %s\n",para_name_string);
	  printf("FLOCK result file size is %d; number of dimensions is %d\n", file_Len, num_dm);
  }

  rewind(f_flock_result);

  result_data = (int **)malloc(sizeof(int*)*file_Len);
  memset(result_data,0,sizeof(int*)*file_Len);
	
  for (i=0;i<file_Len;i++)
  {
		result_data[i]=(int *)malloc(sizeof(int)*num_dm);
		memset(result_data[i],0,sizeof(int)*num_dm);
  }

  readsource(f_flock_result, file_Len, num_dm, result_data);  //note that the populatin_ID starts from 1


  fclose(f_flock_result);

  //now get the number of clusters from the result_data, which is a combination of the data file with two additional columns: event_id and population_id
  
  num_clusters_identified=0;

  for (i=0;i<file_Len;i++)
  {
	  j=result_data[i][num_dm-1]; //the last column, which is the population_id, starting from 1
	  if (j>num_clusters_identified)
		  num_clusters_identified=j;
  }

  //now read the file names in the folder

	note_file_list[0]='\0';
  	strcat(note_file_list,orig_directory);
	strcat(note_file_list,"_file_list");
	
	command_line[0]='\0';
	strcpy(command_line,"dir /b "); //for Unix: strcpy(command_line,"ls ");
	strcat(command_line,orig_directory); //read the original file folder
	strcat(command_line," > ");
	strcat(command_line,note_file_list);
	
	if (DEBUG)
		printf("executing command: %s\n",command_line);
	system(command_line);
		
	f_file_list=fopen(note_file_list,"r");

	while (!feof(f_file_list)) //since the file is generated by ls command, this way can be used as we know there is no new line at the end of the file
	{
		src[0]='\0';
		fgets(src,LINE_LEN,f_file_list);
		num_files++;
	}

	num_files--; //remove the newline at the end of the file

	if (DEBUG)
		printf("Number of Files in Original File Folder %s is %d\n", orig_directory, num_files);
	//Begin to read the file list and save into file_name
	
	percentage_value = (double **)malloc(sizeof(double*)*num_clusters_identified); //this table will be the most important output, some of the values will be zero
	memset(percentage_value,0,sizeof(double*)*num_clusters_identified);

	for (i=0;i<num_clusters_identified;i++)
	{
		percentage_value[i]=(double*)malloc(sizeof(double)*num_files);
		memset(percentage_value[i],0,sizeof(double)*num_files);
	}

	size_population=(int *)malloc(sizeof(int)*num_clusters_identified); //number of events of each population in each file
	memset(size_population,0,sizeof(int)*num_clusters_identified);

	size_file=(int*)malloc(sizeof(int)*num_files);
	memset(size_file,0,sizeof(int)*num_files);

	file_name = (char **)malloc(sizeof(char*)*num_files);
	memset(file_name,0,sizeof(char*)*num_files);

	for (i=0;i<num_files;i++)
	{
		file_name[i]=(char *)malloc(sizeof(char)*WORD_LEN);
		memset(file_name[i],0,sizeof(char)*WORD_LEN);
	}
	
	file_index=0;
	rewind(f_file_list); //reset the file pointer

	while (file_index<num_files)
	{
		file_name[file_index][0]='\0';
		fgets(file_name[file_index],LINE_LEN,f_file_list);

		i=0;

		while ((file_name[file_index][i]!='\r') && (file_name[file_index][i]!='\n'))
			i++;
	
		file_name[file_index][i]='\0';
		
		if ((DEBUG) && ((file_index==0) || (file_index==(num_files-1))))
		{
			printf("file_name[%d] is |%s|\n",file_index,file_name[file_index]);
		}
		file_index++;
	}
		
	fclose(f_file_list);

	csvfoldername[0]='\0';
	strcpy(csvfoldername,directory_name);
	strcat(csvfoldername,"_Pop");
	strcat(csvfoldername,pop_id_char);

	_mkdir(csvfoldername);

	orig_resultfoldername[0]='\0';
	strcpy(orig_resultfoldername,directory_name);
	strcat(orig_resultfoldername,"_flockresult");

	_mkdir(orig_resultfoldername);

	_chdir(orig_directory); //enter the original data folder

	for (file_ID=0;file_ID<num_files;file_ID++)
	{
		f_src=fopen(file_name[file_ID],"r");

		getfileinfo(f_src, &temp_file_size, &temp_num_dm, temp_para_name_string);  //the original data is [size_file[file_ID]][temp_num_dm]
		size_file[file_ID]=temp_file_size;

		fclose(f_src);

		if (size_file[file_ID]<=0)
		{
			printf("There is an original data file that is empty, exit...\n");
			exit(0);
		}

	}

	
	//now we got the total number of events in each original file, in size_file[file_ID]

	//now we begin to process each filtered file, the first thing is to change the file_name

	_chdir(".."); 

	if (DEBUG)
		printf("begin to read file\n");

	for (file_ID=0;file_ID<num_files;file_ID++)
	{

		_chdir(directory_name); //enter the filtered data folder

		fileISempty=0;
		t=0;

		for (i=0;i<num_clusters_identified;i++)
			size_population[i]=0;
		
		f_name[0]='\0';
		while (file_name[file_ID][t]!='.') //the original file name without extension
		{
			f_name[t]=file_name[file_ID][t];
			t++;
		}
		f_name[t]='\0';

		file_orig_name[0]='\0';
		strcpy(file_orig_name,f_name);

		strcat(f_name,"_filt.txt");

		temp_para_name_string[0]='\0';
		temp_file_Len=0;
		temp_num_dm=0;

		f_src=fopen(f_name,"r");

		getfileinfo(f_src,&temp_file_Len,&temp_num_dm,temp_para_name_string); //note temp_file_Len changes for each filtered file
		
		if (DEBUG)
			printf("Number of Events in Filtered File %s is %d with number of dimensions %d\n",f_name, temp_file_Len,temp_num_dm);

		total=0;

		if (temp_file_Len<2)  //the filtered file is empty, all output should be zero
		{
			fileISempty=1;
			fclose(f_src);
		}
		else
		{
			fileISempty=0;

			if (temp_file_Len>file_Len)
			{
				printf("There is a file with more events than that of the result file, which is not ok, exiting...\n");
				exit(0);
			}

			if (temp_num_dm>num_dm)
			{
				printf("There is a file with more dimensions than that of the result file, which is not ok, exiting...\n");
				exit(0);
			}

			rewind(f_src);

			temp_data = (int **)malloc(sizeof(int*)*temp_file_Len);
			memset(temp_data,0,sizeof(int*)*temp_file_Len);
	
			for (i=0;i<temp_file_Len;i++)
			{
				temp_data[i]=(int *)malloc(sizeof(int)*temp_num_dm);
				memset(temp_data[i],0,sizeof(int)*temp_num_dm);
			}

			readsource(f_src,temp_file_Len,temp_num_dm,temp_data); //temp_data[temp_file_Len,temp_num_dm] stores the filted data

			fclose(f_src);

			cluster_id=(int *)malloc(sizeof(int)*temp_file_Len); //population_ID of each row of the filtered data
			memset(cluster_id,0,sizeof(int)*temp_file_Len);

			//read the first and the second line to identify the individual file in the merged file for getting cluster_id, which is the only thing we can rely on.

			t=0;
		
			for (i=0;i<file_Len;i++)  //go through the flock_results.txt to find the two rows that map with the filtered file
			{
				found=1;

				for (j=0;j<temp_num_dm;j++)
				{
					if (result_data[i][j]!=temp_data[0][j]) //compare with the first row of the filtered file
					{
						found=0;
						break;					
					}
				}
			
			
				if (1==found) //found the first row that maps
				{
					if (i>=(file_Len-1)) //if it is already the last row, that means failed to find the next one
					{
						found=0;
						break;			
					}
					else
					{
						for (j=0;j<temp_num_dm;j++)
						{
							if (result_data[i+1][j]!=temp_data[1][j]) //compare with the second row of the filtered file, if this also passes, it is the location
							{
								found=0;
								break;							
							}
						}
					}
				}
				else //not found
					continue;

				if (0==found)
					continue;
				else
				{
					t=i;
					break;
				}
			} //end i=0 to file_Len-1 for searching in the flock_results.txt

			if (0==found)
			{
				printf("Cannot find the filtered file in the merged result, which is not ok, exiting...\n");
				exit(0);
			}
			else
				printf("The %d line in the flock_results.txt maps the filtered file %d\n",t, file_ID+1);

		
			for (i=0;i<temp_file_Len;i++)  //now we can get the population_ID of each row in the filtered file
			{
				cluster_id[i]=result_data[t][num_dm-1]-1; //the last value is the cluster ID, but needs to minus 1 to get the ID

				size_population[cluster_id[i]]++;

				total++;

				t++;
			}

			

		} //end of else if temp_file_Len>=2

		if (DEBUG)
		{
			printf("File %s has %d events\n",f_name,total);
			printf("The original File has %d events\n",size_file[file_ID]);
		}

		for (i=0;i<num_clusters_identified;i++)
		{
			if (1==fileISempty)
				percentage_value[i][file_ID]=0.0;
			else
				percentage_value[i][file_ID]=((double)size_population[i]*100.0)/(double)size_file[file_ID];
		}
					
		//now we need to generate each population file (.csv) in the filtered data
		_chdir("..");

		_chdir(csvfoldername);

		f_csvresult=fopen(f_name,"w");

		fprintf(f_csvresult,"%s\n",temp_para_name_string);

		if (0==fileISempty) //when the filtered data file is not empty
		{
			for (i=0;i<temp_file_Len;i++)
			{
				if ((cluster_id[i]+1)==specified_population_id)  //the specified id starts from 1 while the cluster_id starts from 0
				{
					for (j=0;j<temp_num_dm;j++)
						if (j==(temp_num_dm-1))
							fprintf(f_csvresult,"%d\n",temp_data[i][j]);
						else
							fprintf(f_csvresult,"%d\t",temp_data[i][j]);
				}
				else
					continue;
			}
		}
		fclose(f_csvresult);

		//now we need to move the original folder to generate the flock_results.txt for each of the original file
		_chdir("..");
		_chdir(orig_directory);

		temp_file_size=size_file[file_ID];

		orig_temp_data = (int **)malloc(sizeof(int*)*temp_file_size);
		memset(orig_temp_data,0,sizeof(int*)*temp_file_size);
	
		for (i=0;i<temp_file_size;i++)
		{
			orig_temp_data[i]=(int *)malloc(sizeof(int)*temp_num_dm);
			memset(orig_temp_data[i],0,sizeof(int)*temp_num_dm);
		}

		new_temp_data = (int **)malloc(sizeof(int*)*temp_file_size);
		memset(new_temp_data,0,sizeof(int*)*temp_file_size);
	
		for (i=0;i<temp_file_size;i++)
		{
			new_temp_data[i]=(int *)malloc(sizeof(int)*num_dm); //new_temp_data has two more dimensions than orig_temp_data
			memset(new_temp_data[i],0,sizeof(int)*num_dm);
		}
		
		f_orig_src=fopen(file_name[file_ID],"r"); //this requires the filename between the original file folder and the new file folder is the SAME!
		readsource(f_orig_src,temp_file_size,temp_num_dm,orig_temp_data);
		fclose(f_orig_src);

		//enter the result folder to write the flock_result file
		_chdir("..");
		_chdir(orig_resultfoldername);
	
		result_file_name[0]='\0';
		strcat(result_file_name,file_orig_name);
		strcat(result_file_name,"_results.txt");
		f_results=fopen(result_file_name,"w"); //remember, for the python visualization code, color=1 means white background, so the cluster_ID needs to be changed; also the order of events needs to sorted based on color_ID

		fprintf(f_results,"%s\tEvent\tPopulation\n",temp_para_name_string);

		t=0;
		temp=0;

		for (i=0;i<temp_file_size;i++) //for each row of the original file
		{
			data_kept=1;

			if (1==fileISempty)
				data_kept=0;
			else
			{
				for (p=0;p<temp_file_Len;p++) //for each row of this filtered file
				{
					data_kept=1;
					for (j=0;j<temp_num_dm;j++)
					{
						if (orig_temp_data[i][j]!=temp_data[p][j])
						{
							data_kept=0;
							break;
						}
					}
					if (data_kept==1) //found 1 match
						break;
				} //end for p = each row of the filtered file
			}
			
			if (data_kept==0) //this row doesn't have a cluster_id
			{
				
				for (j=0;j<temp_num_dm;j++)
					new_temp_data[temp][j]=orig_temp_data[i][j];

				new_temp_data[temp][num_dm-2]=temp;
				new_temp_data[temp][num_dm-1]=WHITE_COLOR;	
				temp++;							
			}
			else  //this row has a cluster_id[p], and therefore needs to be put at the end of the output file so that they won't be covered by the background white color when being visualized
			{
				t++;
				for (j=0;j<temp_num_dm;j++)
					new_temp_data[temp_file_size-t][j]=orig_temp_data[i][j];

				new_temp_data[temp_file_size-t][num_dm-2]=temp_file_size-t+1;

				if ((cluster_id[p]+1)!=WHITE_COLOR)
					new_temp_data[temp_file_size-t][num_dm-1]=cluster_id[p]+1; //start from 1
				else
					new_temp_data[temp_file_size-t][num_dm-1]=num_clusters_identified+1; //the original "WHITE_COLOR" is changed to the largest color ID plus 1
				
			}			
		} //end for i = each row of the original file

		
		for (i=0;i<temp_file_size;i++)
		 {
			for (j=0;j<num_dm;j++)
			{
				if (j==num_dm-1)
					fprintf(f_results,"%d\n",new_temp_data[i][j]);
				else
					fprintf(f_results,"%d\t",new_temp_data[i][j]);
			}	
		 }
		fclose(f_results);



		//now we move back to the file folder
		
		_chdir("..");

		if (temp_file_Len>=2) 
		{
			for (i=0;i<temp_file_Len;i++)
				free(temp_data[i]);
			
			free(temp_data);
			free(cluster_id);
		}
		
		for (i=0;i<temp_file_size;i++)
		{
			free(orig_temp_data[i]);
			free(new_temp_data[i]);
		}
		
		free(orig_temp_data);
		free(new_temp_data);
				
		
	} //end of file_ID=0 to num_files-1 for identifying each filtered data file in the flock_results.txt

	percentage_name[0]='\0';
	strcpy(percentage_name,directory_name);
	strcat(percentage_name,"_All_Percentage.txt");
	f_pcnt=fopen(percentage_name,"w");
	fprintf(f_pcnt,"Population_ID\t");

	for (i=0;i<num_files;i++)
	{
		if (i==num_files-1)
			fprintf(f_pcnt,"%s\n",file_name[i]);
		else
			fprintf(f_pcnt,"%s\t",file_name[i]);
	}

	for (i=0;i<num_clusters_identified;i++)
	{
		fprintf(f_pcnt,"Pop_%d\t",i+1);
		for (j=0;j<num_files;j++)
		{
			if (j==num_files-1)
				fprintf(f_pcnt,"%.4f\n",percentage_value[i][j]);
			else
				fprintf(f_pcnt,"%.4f\t",percentage_value[i][j]);								
		}
	}
	fclose(f_pcnt);

	for (i=0;i<file_Len;i++)
		free(result_data[i]);
	free(result_data);

	free(size_file);
	free(size_population);

	for (i=0;i<num_files;i++)
		free(file_name[i]);
	free(file_name);

	for (i=0;i<num_clusters_identified;i++)
		free(percentage_value[i]);
	free(percentage_value);

	 _strtime_s(tmpbuf,sizeof(char)*WORD_LEN);
	printf( "Ending time:\t\t\t\t%s\n", tmpbuf );
	_strdate_s(tmpbuf,sizeof(char)*WORD_LEN);
	printf( "Ending date:\t\t\t\t%s\n", tmpbuf );

}