#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <direct.h>
#include <string.h>
#include <sys/stat.h>
#include <ctype.h>

#define LINE_LEN 4096
#define DEBUG 1

int getfileinfo(char *file_name, int *file_Len, int *num_dm, char *name_string)
{
  char src[LINE_LEN];
  char current_name[64];
  char prv;

  FILE *f_src;
  //int time_ID=-1;
  int name_string_len=0;

  int num_rows=0;
  int num_columns=0;
  int ch='\n';
  int prev='\n';
  //int time_pos=0;
  int i=0;
  int j=0;

  f_src=fopen(file_name,"r"); 
  src[0]='\0';
  fgets(src, LINE_LEN, f_src);

  name_string[0]='\0';
  current_name[0]='\0';
  prv='\n';

  // skip space and tab characters
  while ((src[i]==' ') || (src[i]=='\t'))
    i++;

  // repeat until the end of line is reached
  while ((src[i]!='\0') && (src[i]!='\n') && (src[i]!='\r'))
    {
      current_name[j]=src[i];
		
      if ((src[i]=='\t') && (prv!='\t')) //a complete word
        {
          current_name[j]='\0';
			
         
       //   if (0!=_stricmp(current_name,"Time"))
        //    {
              num_columns++; //num_columns does not inlcude the column of Time
         //     time_pos++;
              strcat(name_string,current_name); 
              strcat(name_string,"\t");
         //   }
         // else 
         //   {
         //     time_ID=time_pos;
         //   }
         

          //num_columns++;
          //strcat(name_string,current_name); 
         // strcat(name_string,"\t");

          current_name[0]='\0';
          j=0;			
        }		
		
      if ((src[i]=='\t') && (prv=='\t')) //a duplicate tab or space
        {
          current_name[0]='\0';
          j=0;
        }
		
      if (src[i]!='\t')
        j++;
		
      prv=src[i];
      i++;
    }
	
  if (prv!='\t') //the last column hasn't been retrieved
    {
      current_name[j]='\0';
     
   //   if (0!=_stricmp(current_name,"Time"))
  //      {
          num_columns++;
          strcat(name_string,current_name);
   //       time_pos++;
   //     }
   //   else
   //     {
   //       time_ID=time_pos;
	//	  name_string_len=strlen(name_string);
	//	  name_string[name_string_len-1]='\0';
     //   }
      
      //num_columns++;
      //strcat(name_string,current_name);
    }
  if (DEBUG==1)
    {
      //printf("time_ID is %d\n",*time_ID);
      printf("name_string is %s\n",name_string);
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
  if (prev!='\n')
    ++num_rows;

   
  *file_Len=num_rows;
  *num_dm=num_columns; 

  printf("original file size is %d; number of dimensions is %d\n", *file_Len, *num_dm);
  fclose(f_src);

  return 0;
}

void readsource(char *file_name, int file_Len, int num_dm, int **uncomp_data)
{
  //int time_pass=0; //to mark whether the time_ID has been passed
  int index=0;

  int i=0;
  int j=0;
  int t=0;

  char src[LINE_LEN];
  char xc[LINE_LEN/10];

  FILE *f_src;
  
  f_src=fopen(file_name,"r");

  src[0]='\0';
  fgets(src,LINE_LEN, f_src); //skip the first line about parameter names

  while (!feof(f_src) && (index<file_Len)) //index = 0, 1, ..., file_Len-1
    {
      src[0]='\0';	    
      fgets(src,LINE_LEN,f_src);
      i=0;
      //time_pass=0;
						
      
      for (t=0;t<num_dm;t++) 
      {
          xc[0]='\0';
          j=0;
          while ((src[i]!='\0') && (src[i]!='\n') && (src[i]!='\r') && (src[i]!=' ') && (src[i]!='\t'))
          {
             xc[j]=src[i];
             i++;
             j++;
          }
		
          xc[j]='\0';	    
          i++;
		  

		 // if (time_pass!=time_ID)
		  	uncomp_data[index][t]=atof(xc);
		//  else
		//	  t--;

		//  time_pass++;
			
       }	
      
       index++;     	
      //fprintf(fout_ID,"%s",src);
    } //end of while
	
  if (DEBUG == 1)
    {
      printf("the last line of the source data is:\n");
      for (j=0;j<num_dm;j++)
        printf("%d ",uncomp_data[index-1][j]);
      printf("\n");
    }
   fclose(f_src);
}

void readmapping(char *correctsequence,int real_num_dm, int *mapping_sequence)
{
	int i=0;
	int j=0;
	int temp=0;
	
	char position[LINE_LEN];
	
	
	for (temp=0;temp<real_num_dm;temp++)
	{
		j=0;
		position[j]='\0';

		while ((correctsequence[i]!='\r') && (correctsequence[i]!='\n') && (correctsequence[i]!='\0') && (correctsequence[i]!='\t'))
		{
			position[j]=correctsequence[i];
			i++;
			j++;
		}
		position[j]='\0';
		mapping_sequence[temp]=atoi(position);

		if (DEBUG)
			printf("mapping_sequence[%d]=%d\n",temp,mapping_sequence[temp]);

		i++;
		
	}

}

void main (int argc, char **argv)
{
	
	FILE *f_file_list;
	FILE *f_stats;
	FILE *f_out;
	FILE *f_header;

	char **file_name;
	
	//float **input_data;
	int **input_data; 

	int *mapping_sequence;

	char directory_name[LINE_LEN];
	char output_directory[LINE_LEN];
	char command_line[LINE_LEN];
	char src[LINE_LEN];
	char note_file_list[LINE_LEN];
	char stats_name[LINE_LEN];
	char name_string[LINE_LEN];
	char clean_file_name[LINE_LEN];
	char correctheader[LINE_LEN];
	char correctsequence[LINE_LEN];
	char header_file_name[LINE_LEN];
	char sequence_file_name[LINE_LEN];
	
	int file_index=0;
	int num_files=0;	
	int num_events=0;
	//int time_id=-1;
	int header_provided=0;
	int real_num_dm=1;

	int file_Len=0;
	int num_dm=0;

	int i=0;
	int j=0;
	int t=0;



	if ((argc!=2) && (argc!=5))
	{
		fprintf(stderr,"Incorrect number of input parameters!\n");
		printf("usage: ReplaceHeaderRemoveTime TXT_File_Directory [FileContainingHeader] [FileContainingMarkerSequence] [RealNumDimension]\n");
		exit(0);
	}

	directory_name[0]='\0';
	strcpy(directory_name,argv[1]);

	if (argc==5)
	{
		header_file_name[0]='\0';
		sequence_file_name[0]='\0';
		strcpy(header_file_name,argv[2]);
		strcpy(sequence_file_name,argv[3]);
		real_num_dm=atoi(argv[4]);
		

		if (DEBUG)
			printf("real_num_dm=%d\n",real_num_dm);

		header_provided=1;

		f_header=fopen(header_file_name,"r");
		correctheader[0]='\0';
		fgets(correctheader,LINE_LEN,f_header);
		fclose(f_header);

		f_header=fopen(sequence_file_name,"r");
		correctsequence[0]='\0';
		fgets(correctsequence,LINE_LEN,f_header);
		fclose(f_header);
	}

	command_line[0]='\0';
	strcpy(command_line,"ls ");
	strcat(command_line,directory_name);
	strcat(command_line," > ");
	strcat(command_line,directory_name);
	strcat(command_line,"_file_list");
	if (DEBUG)
		printf("executing command: %s\n",command_line);
	system(command_line);

	note_file_list[0]='\0';
	strcat(note_file_list,directory_name);
	strcat(note_file_list,"_file_list");
	
	f_file_list=fopen(note_file_list,"r");

	
	while (!feof(f_file_list)) 
	{
		src[0]='\0';
		fgets(src,LINE_LEN,f_file_list);
		num_files++;
	}

	num_files--; //remove the newline at the end of the file

	if (DEBUG)
		printf("num_files is %d\n",num_files);
	//Begin to read the file list and save into file_name
	
	file_name = (char **)malloc(sizeof(char*)*num_files);
	memset(file_name,0,sizeof(char*)*num_files);
	
	for (i=0;i<num_files;i++)
	{
		file_name[i]=(char *)malloc(sizeof(char)*LINE_LEN);
		memset(file_name[i],0,sizeof(char)*LINE_LEN);
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

	i=0;
		
	if (DEBUG)
		printf("begin to read file\n");
	//Begin to read files


	stats_name[0]='\0';
	strcat(stats_name,directory_name);
	strcat(stats_name,"_fcs_stats.txt");

	output_directory[0]='\0';
	strcat(output_directory,directory_name);
	strcat(output_directory,"_out");

	f_stats=fopen(stats_name,"w");
	fprintf(f_stats,"FileName\tNumEvents\tNumDimensions\n");

	mkdir(output_directory, S_IRWXU | S_IRWXG | S_IRWXO);

	mapping_sequence= (int *)malloc(sizeof(int)*real_num_dm);
	memset(mapping_sequence,0,sizeof(int)*real_num_dm);

	if (header_provided)
		readmapping(correctsequence,real_num_dm,mapping_sequence);


	for (i=0;i<num_files;i++)
	{
		chdir(directory_name);
		
		getfileinfo(file_name[i], &file_Len, &num_dm, name_string);

		fprintf(f_stats,"%s\t%d\t%d\n", file_name[i], file_Len, num_dm);

		input_data = (int **)malloc(sizeof(int*)*file_Len);
		memset(input_data,0,sizeof(int*)*file_Len);
		for (j=0;j<file_Len;j++)
		{
			input_data[j]=(int *)malloc(sizeof(int)*num_dm);
			 memset(input_data[j],0,sizeof(int)*num_dm);
		}

		
		readsource(file_name[i], file_Len, num_dm, input_data);

		

		chdir("..");
		chdir(output_directory);

		clean_file_name[0]='\0';
		strcpy(clean_file_name,file_name[i]);
		
		f_out=fopen(clean_file_name,"w");
		
		if (header_provided)
		{
			fprintf(f_out,correctheader);
			//fputc('\n',f_out);

			for (j=0;j<file_Len;j++)
			{
				for (t=0;t<real_num_dm;t++)
				{
					if (t==real_num_dm-1)
						fprintf(f_out,"%d\n",input_data[j][mapping_sequence[t]]);
					else
						fprintf(f_out,"%d\t",input_data[j][mapping_sequence[t]]);
				}
			}
		}
		else
		{
			fprintf(f_out,name_string);
			//fputc('\n',f_out);

			for (j=0;j<file_Len;j++)
			{
				for (t=0;t<num_dm;t++)
				{
					if (t==num_dm-1)
						fprintf(f_out,"%d\n",input_data[j][t]);
					else
						fprintf(f_out,"%d\t",input_data[j][t]);
				}
			}
		}

		
		

		fclose(f_out);
		
		for (j=0;j<file_Len;j++)
			free(input_data[j]);
		free(input_data);

		

		chdir("..");		
		
	}

	
	fclose(f_stats);
	free(mapping_sequence);
	
	for (t=0;t<num_files;t++)
		free(file_name[t]);
	free(file_name);

	
}