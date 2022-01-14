#include <time.h>
#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
//#include <direct.h>
#include <string.h>
#include <sys/stat.h>
#include <ctype.h>

#define LINE_LEN 4096
#define DEBUG 1

int getfileinfo(char *file_name, int *file_Len)
{
  char src[LINE_LEN];

  FILE *f_src;

  int num_rows=0;
  int ch='\n';
  int prev='\n';
  //int time_pos=0;

  f_src=fopen(file_name,"r");
  src[0]='\0';
  fgets(src, LINE_LEN, f_src);


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

   fclose(f_src);

  return 0;
}

int getkeywordinfo(char *file_name, int num_rows, int num_keywords, FILE *f_out)
{

	int i=0;
	int j=0;
	int t=1;
	int found=0;
	int found_n=0;
	int found_s=0;

	char src[LINE_LEN];

	char keyword_name[LINE_LEN];
	char desired_name_N[LINE_LEN];
	char desired_name_S[LINE_LEN];

	FILE *f_src;

	f_src=fopen(file_name,"r");

	for (t=1;t<=num_keywords;t++)
	{
		desired_name_N[0]='\0';
		desired_name_S[0]='\0';

		sprintf(desired_name_N,"$P%dN",t);
		sprintf(desired_name_S,"$P%dS",t);

		//printf("desired_name_N is |%s|",desired_name_N);
		rewind(f_src);

		found_n=0;
		found_s=0;

		for (i=0;i<num_rows;i++)
		{

			src[0]='\0';
			fgets(src,LINE_LEN,f_src);

			found=0;

			keyword_name[0]='\0';

			j=0;

			for (j=0;j<6;j++) //check the first 6 characters
			{
				if ((src[j]!='=') && (src[j]!='\0') && (src[j]!='\r') && (src[j]!='\n'))
					keyword_name[j]=src[j];
				else
					break;
			}
			keyword_name[j]='\0';


			if (src[j]=='=') //found a keyword
			{
				found=1; //found a keyword
				//printf("keyword_name is |%s|",keyword_name);
			}

			j++; //move to the next character in src

			if (found==1)
			{
				if (strcmp(keyword_name,desired_name_N)==0)
				{
					found_n=1;

					while ((src[j]!='\r') && (src[j]!='\n') && (src[j]!='\0'))
					{
						fputc(src[j],f_out);
						j++;
					}

				}

				if (strcmp(keyword_name,desired_name_S)==0)
				{
					found_s=1;
					fputc('_',f_out);

					while ((src[j]!='\r') && (src[j]!='\n') && (src[j]!='\0'))
					{
						fputc(src[j],f_out);
						j++;
					}

					if (t>=num_keywords)
						fputc('\n',f_out);
					else
						fputc('\t',f_out);
				}
			} //end for found==1

			if ((found_n==1) && (found_s==1)) //if both are found, no need to continue searching
				break;
		} //end for i=0 to num_rows

		if (found_n==0) //there is no $PnN, there must be no $PnS either, print NA
		{
			if (t>=num_keywords)
				fprintf(f_out,"NA\n");
			else
				fprintf(f_out,"NA\t");
		}
		if ((found_n==1) && (found_s==0)) //there is $PnN but no $PnS, add the separator after $PNN
		{
			if (t>=num_keywords)
				fputc('\n',f_out);
			else
				fputc('\t',f_out);
		}

	} //end for t=1 to num_keywords

	fclose(f_src);

	return 0;
}

void main (int argc, char **argv)
{

	FILE *f_file_list;

	FILE *f_out;

	char **file_name;


	char directory_name[LINE_LEN];
	char output_name[LINE_LEN];

	char command_line[LINE_LEN];
	char src[LINE_LEN];
	char note_file_list[LINE_LEN];


	int file_index=0;
	int num_files=0;

	int num_keywords=1; //number of keywords to read

	int i=0;
	int j=0;
	int t=0;

	int num_rows;

	if (argc!=3)
	{
		fprintf(stderr,"Incorrect number of input parameters!\n");
		printf("usage: ReadKeyWords KeyWord_File_Directory MaxNumberOfKeyWords\n");
		exit(0);
	}

	directory_name[0]='\0';
	strcpy(directory_name,argv[1]);

	num_keywords=atoi(argv[2]);

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


	output_name[0]='\0';
	strcat(output_name,directory_name);
	strcat(output_name,"_keywords.txt");

	f_out=fopen(output_name,"w");

	chdir(directory_name);

	for (i=0;i<num_files;i++)
	{
		getfileinfo(file_name[i], &num_rows);

		//printf("Number of Rows is %d\n",num_rows);

		fprintf(f_out,"%s\t",file_name[i]);

		getkeywordinfo(file_name[i], num_rows, num_keywords, f_out);
	}
	fclose(f_out);

	for (t=0;t<num_files;t++)
		free(file_name[t]);
	free(file_name);

}
