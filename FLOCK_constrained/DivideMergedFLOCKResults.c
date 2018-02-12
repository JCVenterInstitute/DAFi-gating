/*********************************************************************************************************************************
This program serves as one step in the cross-sample mapping pipeline. It maps FLOCK results from merged aligned data to the
original individual files.

This program reads in a FLOCK result file and two folders. One folder (after alignment) contains the aligned individual files
(Tab-delimited TXT files), and the other contains the original individual files (before alignment).
The program outputs in the before alignment folder a set of sub-folders, each contains FLOCK results for one original 
individual file. These result files can be copied to the visualization interface for generating dot plots and summary 
statistics.

The file names and number of files in both folders need to be exactly the same.

Another input to this program is the number of clusters in the FLOCK results from merged aligned data.

The algorithm procedure is:

Read the FLOCK result file into memory, and the header. The result file is for a concatenation of all aligned individual files.

Enter the folder containing the individual aligned files, reads the number of files
Initiate a loop for each file
Open the file, reads the first and the second rows
Extract the population_ID from the FLOCK result file by searching for the first and the second rows.

Enter the folder containing the individual unaligned files, use the population_ID to generate other result files (MFI, profile,
percentage etc.), output files of each file into a different folder.

Date: May 15, 2015 (written in one day!!)
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
	
  if (DEBUG)
    {
      printf("the last line of the source data is:\n");
      for (j=0;j<num_dm;j++)
        printf("%d ",uncomp_data[index-1][j]);
      printf("\n");
    }
}

void show(int **Matrix, int *cluster_id, int file_Len, int k, int num_dm, char *name_string)
{
	int situ1=0;
	int situ2=0;

	int i=0;
	int id=0;
	int j=0;
	int info_id=0;
	int nearest_id=0;
	int insert=0;
	int temp=0;
	int m=0;
	int n=0;
	int t=0;
	
	int *size_c;
	


	int **size_mybound_1;
	int **size_mybound_2;
	int **size_mybound_3;
	int **size_mybound_0;

	float interval=0.0;

	float *big;
	float *small;


	float **center;
	float **mybound;
	
	int **prof; //prof[i][j]=1 means population i is + at parameter j
	
	FILE *fpcnt_id; //proportion id
	//FILE *fcent_id; //center_id, i.e., centers of clusters within the original data
	FILE *fprof_id; //profile_id

	big=(float *)malloc(sizeof(float)*num_dm);
	memset(big,0,sizeof(float)*num_dm);

	small=(float *)malloc(sizeof(float)*num_dm);
	memset(small,0,sizeof(float)*num_dm);

	for (i=0;i<num_dm;i++)
	{
		big[i]=0.0;
		small[i]=(float)MAX_VALUE;
	}
	
	
	
	size_c=(int *)malloc(sizeof(int)*k);
	memset(size_c,0,sizeof(int)*k);

	center=(float**)malloc(sizeof(float*)*k);
	memset(center,0,sizeof(float*)*k);
	for (i=0;i<k;i++)
	{
		center[i]=(float*)malloc(sizeof(float)*num_dm);
		memset(center[i],0,sizeof(float)*num_dm);
	}

	

	mybound=(float**)malloc(sizeof(float*)*num_dm);
	memset(mybound,0,sizeof(float*)*num_dm);
	for (i=0;i<num_dm;i++) //there are 3 mybounds for 4 categories
	{
		mybound[i]=(float*)malloc(sizeof(float)*3);
		memset(mybound[i],0,sizeof(float)*3);
	}

	

	prof=(int **)malloc(sizeof(int*)*k);
	memset(prof,0,sizeof(int*)*k);
	for (i=0;i<k;i++)
	{
		prof[i]=(int *)malloc(sizeof(int)*num_dm);
		memset(prof[i],0,sizeof(int)*num_dm);
	}

	

	for (i=0;i<file_Len;i++)
	{
		id=cluster_id[i]-1; //this cluster_ID is not from zero but from 1, because it is read from the FLOCK result TXT file. But here we need it to start from 0.
		for (j=0;j<num_dm;j++)
		{
			center[id][j]=center[id][j]+(float)Matrix[i][j];
			if (big[j]<(float)Matrix[i][j])
				big[j]=(float)Matrix[i][j];
			if (small[j]>(float)Matrix[i][j])
				small[j]=(float)Matrix[i][j];
		}
		
		size_c[id]++;		
	}

	

	for (i=0;i<k;i++)
		for (j=0;j<num_dm;j++)
		{			
			if (size_c[i]!=0)
				center[i][j]=(center[i][j]/(float)(size_c[i]));
			else
				center[i][j]=0;	
		}

	

	for (j=0;j<num_dm;j++)
	{
		interval=((big[j]-small[j])/(float)4);
		//printf("interval[%d] is %f\n",j,interval);
		for (i=0;i<3;i++)
			mybound[j][i]=small[j]+((float)(i+1)*interval);
	}
	

	size_mybound_0=(int **)malloc(sizeof(int*)*k);
	memset(size_mybound_0,0,sizeof(int*)*k);
	
	for (i=0;i<k;i++)
	{
		size_mybound_0[i]=(int*)malloc(sizeof(int)*num_dm);
		memset(size_mybound_0[i],0,sizeof(int)*num_dm);		
	}

	size_mybound_1=(int **)malloc(sizeof(int*)*k);
	memset(size_mybound_1,0,sizeof(int*)*k);
	
	for (i=0;i<k;i++)
	{
		size_mybound_1[i]=(int*)malloc(sizeof(int)*num_dm);
		memset(size_mybound_1[i],0,sizeof(int)*num_dm);		
	}

	size_mybound_2=(int **)malloc(sizeof(int*)*k);
	memset(size_mybound_2,0,sizeof(int*)*k);
	
	for (i=0;i<k;i++)
	{
		size_mybound_2[i]=(int*)malloc(sizeof(int)*num_dm);
		memset(size_mybound_2[i],0,sizeof(int)*num_dm);	
	}

	size_mybound_3=(int **)malloc(sizeof(int*)*k);
	memset(size_mybound_3,0,sizeof(int*)*k);

	for (i=0;i<k;i++)
	{
		size_mybound_3[i]=(int*)malloc(sizeof(int)*num_dm);
		memset(size_mybound_3[i],0,sizeof(int)*num_dm);
	}
	
	for (i=0;i<file_Len;i++)
		for (j=0;j<num_dm;j++)
			{
				if ((float)Matrix[i][j]<mybound[j][0])// && ((Matrix[i][j]-small[j])>0)) //the smallest values excluded
					size_mybound_0[cluster_id[i]-1][j]++;
				else
				{
					if ((float)Matrix[i][j]<mybound[j][1])
						size_mybound_1[cluster_id[i]-1][j]++;
					else
					{
						if ((float)Matrix[i][j]<mybound[j][2])
							size_mybound_2[cluster_id[i]-1][j]++;
						else
							//if (Matrix[i][j]!=big[j]) //the biggest values excluded
								size_mybound_3[cluster_id[i]-1][j]++;
					}

				}
			}

	fprof_id=fopen("profile.txt","w");
	fprintf(fprof_id,"Population_ID\t");
	fprintf(fprof_id,"%s\n",name_string);
	
	for (i=0;i<k;i++)
	{
		fprintf(fprof_id,"%d\t",i+1); //i changed to i+1 to start from 1 instead of 0: April 16, 2009
		for (j=0;j<num_dm;j++)
		{
			
			if (size_mybound_0[i][j]>size_mybound_1[i][j])
				situ1=0;
			else
				situ1=1;
			if (size_mybound_2[i][j]>size_mybound_3[i][j])
				situ2=2;
			else
				situ2=3;

			if ((situ1==0) && (situ2==2))
			{
				if (size_mybound_0[i][j]>size_mybound_2[i][j])
					prof[i][j]=0;
				else
					prof[i][j]=2;
			}
			if ((situ1==0) && (situ2==3))
			{
				if (size_mybound_0[i][j]>size_mybound_3[i][j])
					prof[i][j]=0;
				else
					prof[i][j]=3;
			}
			if ((situ1==1) && (situ2==2))
			{
				if (size_mybound_1[i][j]>size_mybound_2[i][j])
					prof[i][j]=1;
				else
					prof[i][j]=2;
			}
			if ((situ1==1) && (situ2==3))
			{
				if (size_mybound_1[i][j]>size_mybound_3[i][j])
					prof[i][j]=1;
				else
					prof[i][j]=3;
			}
			
			//begin to output profile
			if (j==num_dm-1)
			{
				if (prof[i][j]==0)
					fprintf(fprof_id,"1\n");
				if (prof[i][j]==1)
					fprintf(fprof_id,"2\n");
				if (prof[i][j]==2)
					fprintf(fprof_id,"3\n");
				if (prof[i][j]==3)
					fprintf(fprof_id,"4\n");
			}
			else
			{
				if (prof[i][j]==0)
					fprintf(fprof_id,"1\t");
				if (prof[i][j]==1)
					fprintf(fprof_id,"2\t");
				if (prof[i][j]==2)
					fprintf(fprof_id,"3\t");
				if (prof[i][j]==3)
					fprintf(fprof_id,"4\t");
			}
		}
	}
	fclose(fprof_id);

	///////////////////////////////////////////////////////////
	

	fpcnt_id=fopen("percentage.txt","w");
	fprintf(fpcnt_id,"Population_ID\tPercentage\n");

	for (t=0;t<k;t++)
	{
		fprintf(fpcnt_id,"%d\t%.4f\n",t+1,(float)size_c[t]*100.0/(float)file_Len);	//t changed to t+1 to start from 1 instead of 0: April 16, 2009									
	}
	fclose(fpcnt_id);

	free(big);
	free(small);
	free(size_c);

	for (i=0;i<k;i++)
	{
		free(center[i]);
		free(prof[i]);
		free(size_mybound_0[i]);
		free(size_mybound_1[i]);
		free(size_mybound_2[i]);
		free(size_mybound_3[i]);
	}
	free(center);
	free(prof);
	free(size_mybound_0);
	free(size_mybound_1);
	free(size_mybound_2);
	free(size_mybound_3);

	for (i=0;i<num_dm;i++)
		free(mybound[i]);
	free(mybound);
	
}

void ID2Center_all(int **data_in, int file_Len, int num_dm, int num_clust, int *cluster_ID, float **population_center)
{
  int i=0;
  int j=0;
  int ID=0;
  int *size_c;



  size_c=(int *)malloc(sizeof(int)*num_clust);
  memset(size_c,0,sizeof(int)*num_clust);

  for (i=0;i<num_clust;i++)
    for (j=0;j<num_dm;j++)
      population_center[i][j]=0;

	
  for (i=0;i<file_Len;i++)
    {
         ID=cluster_ID[i]-1;
			
         if (ID==-1)
         {
            //printf("ID==-1! in ID2Center_all\n");
            //exit(0);
			fprintf(stderr,"Incorrect file format or input parameters (resulting in incorrect population IDs)\n"); //modified on July 23, 2010
			exit(0);
         }

         for (j=0;j<num_dm;j++)
           population_center[ID][j]=population_center[ID][j]+(float)data_in[i][j];
			
         size_c[ID]++;        
    }
	
 
  for (i=0;i<num_clust;i++)
    {
      for (j=0;j<num_dm;j++)
        if (size_c[i]!=0)
          population_center[i][j]=(population_center[i][j]/(float)(size_c[i]));
        else
		{
          population_center[i][j]=0;
		  //printf("size_c[%d]=0 in ID2center_all\n",i);
		}
    }


  free(size_c);

}

void main (int argc, char **argv)
{
  
  char tmpbuf[WORD_LEN];
  char directory_name[WORD_LEN];
  char orig_directory[WORD_LEN];
  char f_name[WORD_LEN];

  FILE *f_flock_result;
  FILE *f_file_list;
  FILE *f_src;
  FILE *f_orig_src;

  FILE *f_out; //coordinates
  FILE *f_cid; //population-ID of events
  FILE *f_ctr; //centroids of populations
  FILE *f_results; //coordinates file event and population column
  FILE *f_mfi; //added April 16, 2009 for mean fluorescence intensity
  FILE *f_parameters; //number of bins and density calculated by
                      //the algorithm. Used to update the database
  FILE *f_properties; //Properties file used by Image generation software

  char note_file_list[LINE_LEN];
  char command_line[LINE_LEN];

  char src[LINE_LEN];
  char para_name_string[LINE_LEN];
  char temp_para_name_string[LINE_LEN];

  int num_total_clusters=0;
  int num_files=0;
  int file_index=0;
  int num_dm=0;
  int file_Len=0;
  int temp_num_dm=0;
  int temp_file_Len=0;
  int i=0;
  int j=0;
  int t=0;
  int file_ID=0;
  int found=0;
  
  char **file_name;
  int **result_data;
  int **temp_data;
  int **orig_temp_data;
  float **population_center;

  int *cluster_id;

  int min = 999999;
  int max = 0;

  _strtime_s(tmpbuf,sizeof(char)*WORD_LEN);
  printf( "Starting time:\t\t\t\t%s\n", tmpbuf );
  _strdate_s(tmpbuf,sizeof(char)*WORD_LEN);
  printf( "Starting date:\t\t\t\t%s\n", tmpbuf );

  if (argc!=5)
  {
	  fprintf(stderr,"Incorrect number of input parameters!\n"); 
      fprintf(stderr,"usage:\n");
      fprintf(stderr,"DivideMergedFLOCKResults FLOCK_Result_File FolderNameContainingAlignedIndividualFiles Number_Of_Clusters_In_FLOCK_Result FolderNameContainingOriginalIndividualFiles\n");
	  fprintf(stderr,"The number of files and the file names under both directories must be exactly the same. One set is the original, while the other set is after alignment.\n");
	  exit(0);
  }


  f_flock_result=fopen(argv[1],"r");

  directory_name[0]='\0';
  strcpy(directory_name,argv[2]);

  num_total_clusters=atoi(argv[3]);

  orig_directory[0]='\0';
  strcpy(orig_directory,argv[4]);

  getfileinfo(f_flock_result, &file_Len, &num_dm, para_name_string);

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

  readsource(f_flock_result, file_Len, num_dm, result_data);


  fclose(f_flock_result);

  //now reads the file names in the folder

	note_file_list[0]='\0';
  	strcat(note_file_list,directory_name);
	strcat(note_file_list,"_file_list");
	
	command_line[0]='\0';
	strcpy(command_line,"dir /b "); //for Unix: strcpy(command_line,"ls ");
	strcat(command_line,directory_name);
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
		printf("Number of Files in Folder %s is %d\n", directory_name, num_files);
	//Begin to read the file list and save into file_name
	
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

	//now we begin to process each file
	

	if (DEBUG)
		printf("begin to read file\n");

	for (file_ID=0;file_ID<num_files;file_ID++)
	{

		_chdir(directory_name);

		f_src=fopen(file_name[file_ID],"r");

		t=0;
		f_name[0]='\0';
		while (file_name[file_ID][t]!='.')
		{
			f_name[t]=file_name[file_ID][t];
			t++;
		}
		f_name[t]='\0';

		temp_para_name_string[0]='\0';
		temp_file_Len=0;
		temp_num_dm=0;

		getfileinfo(f_src,&temp_file_Len,&temp_num_dm,temp_para_name_string);
		
		printf("temp_file_Len is %d; temp_num_dm is %d\n",temp_file_Len,temp_num_dm);

		if (temp_file_Len<2)
		{
			printf("There is a file with no or only 1 event, which is not ok, exiting...\n");
			exit(0);
		}

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

		readsource(f_src,temp_file_Len,temp_num_dm,temp_data);

		fclose(f_src);

		cluster_id=(int *)malloc(sizeof(int)*temp_file_Len);
		memset(cluster_id,0,sizeof(int)*temp_file_Len);

		//read the first and the second line to identify the individual file in the merged file for getting cluster_id, which is the only thing we can rely on.

		t=0;
		
		for (i=0;i<file_Len;i++)
		{
			found=1;

			for (j=0;j<temp_num_dm;j++)
			{
				if (result_data[i][j]!=temp_data[0][j]) //compare with the first row of the file
				{
					found=0;
					break;
					
				}
			}
			
			
			if (1==found)
			{
				if (i>=file_Len-1)
				{
					found=0;
					break;
					
				}
				else
				{
					for (j=0;j<temp_num_dm;j++)
					{
						if (result_data[i+1][j]!=temp_data[1][j]) //compare with the second row of the file, if this also passes, it is the location
						{
							found=0;
							break;
							
						}
					}
				}
			}
			else
				continue;

			if (0==found)
				continue;
			else
			{
				t=i;
				break;
			}
		}

		if (0==found)
		{
			printf("Cannot find the file in the result, which is not ok, exiting...\n");
			exit(0);
		}
		else
			printf("t=%d\n",t);

		
		for (i=0;i<temp_file_Len;i++)
		{
			cluster_id[i]=result_data[t][num_dm-1]; //the last value is the cluster ID
			t++;
		}
				
		//now we need to move the original folder to get the data
		_chdir("..");
		_chdir(orig_directory);

		orig_temp_data = (int **)malloc(sizeof(int*)*temp_file_Len);
		memset(orig_temp_data,0,sizeof(int*)*temp_file_Len);
	
		for (i=0;i<temp_file_Len;i++)
		{
			orig_temp_data[i]=(int *)malloc(sizeof(int)*temp_num_dm);
			memset(orig_temp_data[i],0,sizeof(int)*temp_num_dm);
		}
		
		f_orig_src=fopen(file_name[file_ID],"r");
		readsource(f_orig_src,temp_file_Len,temp_num_dm,orig_temp_data);
		fclose(f_orig_src);

		//use f_name to create a new folder

		_mkdir(f_name);

		_chdir(f_name);
		//printf("before show\n");
		//write results in 
		show(orig_temp_data,cluster_id,temp_file_Len,num_total_clusters,temp_num_dm,temp_para_name_string);

		//printf("completing show\n");

		population_center=(float **)malloc(sizeof(float*)*num_total_clusters);
		memset(population_center,0,sizeof(float*)*num_total_clusters);
		for (i=0;i<num_total_clusters;i++)
		{
			population_center[i]=(float*)malloc(sizeof(float)*temp_num_dm);
			memset(population_center[i],0,sizeof(float)*temp_num_dm);
		 }

		ID2Center_all(orig_temp_data,temp_file_Len,temp_num_dm,num_total_clusters,cluster_id,population_center);

		//printf("completing ID2Center\n");
  

		f_cid=fopen("population_id.txt","w");		
		f_out=fopen("coordinates.txt","w");
		f_results=fopen("flock_results.txt","w");

		for (i=0;i<temp_file_Len;i++)
			fprintf(f_cid,"%d\n",cluster_id[i]); //this cluster_id is from the result file, no need to be added by 1

		fprintf(f_out,"%s\n",temp_para_name_string);
		fprintf(f_results,"%s\tEvent\tPopulation\n",temp_para_name_string);
	
		for (i=0;i<temp_file_Len;i++)
		{
			for (j=0;j<temp_num_dm;j++)
			{
				if (orig_temp_data[i][j] < min) {
					min = orig_temp_data[i][j];
				}
				if (orig_temp_data[i][j] > max) {
					max = orig_temp_data[i][j];
				}
				if (j==temp_num_dm-1)
				{
					fprintf(f_out,"%d\n",orig_temp_data[i][j]);
					fprintf(f_results,"%d\t",orig_temp_data[i][j]);
				}
				else
				{
					fprintf(f_out,"%d\t",orig_temp_data[i][j]);
					fprintf(f_results,"%d\t",orig_temp_data[i][j]);
				}
			}
			fprintf(f_results,"%d\t",i + 1);
			fprintf(f_results,"%d\n",cluster_id[i]); //this cluster_id is from the result file, no need to be added by 1
		}

		fclose(f_cid);
		fclose(f_out);
		fclose(f_results);

		f_parameters=fopen("parameters.txt","w");
		fprintf(f_parameters,"Number_of_Bins\t0\n"); //cross-sample comparison
		fprintf(f_parameters,"Density\t0\n"); //crosss-sample comparison
		fprintf(f_parameters,"Min\t%d\n",min);
		fprintf(f_parameters,"Max\t%d\n",max);
		fclose(f_parameters);

		f_properties=fopen("fcs.properties","w");
		fprintf(f_properties,"Bins=0\n"); //cross-sample comparison
		fprintf(f_properties,"Density=0\n"); //cross-sample comparison
		fprintf(f_properties,"Min=%d\n",min);
		fprintf(f_properties,"Max=%d\n",max);
		fprintf(f_properties,"Populations=%d\n",num_total_clusters);
		fprintf(f_properties,"Events=%d\n",file_Len);
		fprintf(f_properties,"Markers=%d\n",num_dm);
		fclose(f_properties);

		f_ctr=fopen("population_center.txt","w");
		f_mfi=fopen("MFI.txt","w");

		for (i=0;i<num_total_clusters;i++) {
	
			fprintf(f_ctr,"%d\t",i+1);  //i starts from 1 instead of 0
			fprintf(f_mfi,"%d\t",i+1);

			for (j=0;j<temp_num_dm;j++) {
				if (j==temp_num_dm-1)
				{
					fprintf(f_ctr,"%.0f\n",population_center[i][j]);
					fprintf(f_mfi,"%.0f\n",population_center[i][j]);
				}
				else
				{
					fprintf(f_ctr,"%.0f\t",population_center[i][j]);
					fprintf(f_mfi,"%.0f\t",population_center[i][j]);
				}
			}
		}

  	
		fclose(f_mfi);
		fclose(f_ctr);


		//now we move back to the file folder
		
		_chdir("..");
		_chdir("..");

		for (i=0;i<temp_file_Len;i++)
		{
			free(temp_data[i]);
			free(orig_temp_data[i]);
		}
		free(temp_data);
		free(orig_temp_data);

		free(cluster_id);
		
		for (i=0;i<num_total_clusters;i++)
		{
			free(population_center[i]);
		}
		free(population_center);

	}

	for (i=0;i<file_Len;i++)
		free(result_data[i]);
	free(result_data);

	for (i=0;i<num_files;i++)
		free(file_name[i]);
	free(file_name);

	 _strtime_s(tmpbuf,sizeof(char)*WORD_LEN);
	printf( "Ending time:\t\t\t\t%s\n", tmpbuf );
	_strdate_s(tmpbuf,sizeof(char)*WORD_LEN);
	printf( "Ending date:\t\t\t\t%s\n", tmpbuf );

}