/************************************************************************************************************************
	
	FLOCK data clustering with user specified constraints about boundaries of positive and negative
	
	Author: Yu "Max" Qian, Ph.D., mqian@jcvi.org or qianyu_cs@yahoo.com
	
	Copyright: Author and J. Craig Venter Institute
	
	Development time: August 25, 2015
	
	Usage: FLOCK_constrained data_file num_bins density_threshold max_num_clusters profile_spec_file pop_proportion_threshold pop_range_distance_threshold
	
	data_file is a tab delimited file with compensated and transformed values.
	num_bins: 6-30
	density_threshold: 3-50
	max_num_clusters: upper bound of the number of clusters

	Setting num_bins or density_threshold 0, the program will automatically decide the thresholds

	Recommend use: fixing the num_bins and upper bound of number of clusters, but let the program to
	select the density threshold

	profile_spec_file: a two column tab delimited file: dimension ID in the data file, separating_coordinate.
	
	For example, if the data file header is CD3, CD4, CD25, smaller than 50 is CD3-; larger than 50 is CD3+
	1	50
	2	75
	3	80

	The algorithm is:
	Step 1: Run FLOCK to identify all clusters in full dimensional space
	Step 2: For each FLOCK population 
	For each dimension
	Calcluate its proportion outside each dimension separatining_coordinate
	Calculate the value range of the outside part on each dimension
	Check whether the proportion is smaller than the proportion_threshold AND the range of the dimension is smaller than the range_distance_threshold
	if yes, check the next cell population
	otherwise, generate a new centroid (the centroid of the outside events) on this dimension; check the next dimension
	Step 3: Run K-means based on the new set of centroids; Output the final result
	
	A new profile file with the negative/positive for each population on each dimension will be generated.
***********************************************************************************************************************/
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
//#include <direct.h>
#include <unistd.h>
#include <assert.h>


#define DEBUG 0
#define LINE_LEN 2048
#define FILE_NAME_LEN 128
#define PARA_NAME_LEN 64
#define MAX_VALUE 1000000000
#define MIN_GRID 6
#define MAX_GRID 30
#define NUM_CLUSTER_FOLD 50 //number of dense grids cannot be more than 20 times of the maximum number of clusters specified by the user
#define NORM_METHOD 3 //2 if z-score; 0 if no normalization; 1 if min-max; 3 if 0-4095 rescaling 
#define KMEANS_TERM 100
#define UPPER_NUM_REGION 1000000  //upper bound of number of hyper regions
#define DEFAULT_DENSITY_T 3
#define MAX_D_T 50
#define DEFAULT_NUM_POP 30
#define LOWER_NUM_CLUSTER 2
#define UPPER_NUM_CLUSTER 999 //upper bound of number of data clusters specified by the user
//#define DAG_BIN 200
#define UPPER_VALUE 4095
//#define PROP_FILTER_T 0.3

int find_connected(int **G, int num_dense_grids, int ndim, int *grid_clusterID);

void get_spec_num(FILE *f_spec, int *num_rows)
{
	int number_rows=0;
	char line[LINE_LEN];
	
	line[0]='\0';

	while (fgets(line,LINE_LEN,f_spec)!=NULL) 
	{
		number_rows++;
		line[0]='\0';
	}

	//printf("Number of Dimensions to be profiled is %d\n",number_rows);
	*num_rows=number_rows;
}

void get_spec_info(FILE *f_spec, int *filtered_d, int *filtered_low, int num_rows)
{
	int temp=0;
	char line[LINE_LEN];
	
	line[0]='\0';

	while ((fgets(line,LINE_LEN,f_spec)!=NULL) && (temp<num_rows))
	{
		sscanf(line,"%d\t%d\n",&filtered_d[temp],&filtered_low[temp]);
		line[0]='\0';
		temp++;
	}

	//if (temp==num_rows)
	//	printf("The last row of the spec file is %d\t%d\n", filtered_d[temp-1],filtered_low[temp-1]);

}


/************* Read basic info of the source file ****************************/
void getfileinfo(FILE *f_src, int *file_Len, int *num_dm, char *name_string, int *time_ID)
{
  char src[LINE_LEN];
  char current_name[64];
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

  if ((src[0]=='F') && (src[1]=='C') && (src[2]=='S'))
	{
		fprintf(stderr,"the correct input format is a tab-delimited txt file, instead of FCS file.\n");
		abort();
	}

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
			
          if (0!=strcmp(current_name,"Time"))
            {
              num_columns++; //num_columns does not inlcude the column of Time
              time_pos++;
              strcat(name_string,current_name); 
              strcat(name_string,"\t");
            }
          else
            {
              *time_ID=time_pos;
            }
          
          
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
	
  if (prv!='\t') //the last one hasn't been retrieved
    {
      current_name[j]='\0';
      
      if (0!=strcmp(current_name,"Time"))
        {
          num_columns++;
          strcat(name_string,current_name);
          time_pos++;
        }
      else
        {
          *time_ID=time_pos;
        }
      
      
    }
  if (DEBUG==1)
    {
      printf("time_ID is %d\n",*time_ID);
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

  //added on July 23, 2010
  if (num_rows<50)
  {
    fprintf(stderr,"Number of events in the input file is too few and should not be processed!\n"); //modified on July 23, 2010
	abort();
  }
  
  *file_Len=num_rows;
  *num_dm=num_columns; 

  printf("original file size is %d; number of dimensions is %d\n", *file_Len, *num_dm);
}



/************************************* Read the source file into uncomp_data **************************************/
void readsource(FILE *f_src, int file_Len, int num_dm, double **uncomp_data, int time_ID)
{
  int time_pass=0; //to mark whether the time_ID has been passed
  int index=0;

  int i=0;
  int j=0;
  int t=0;

  char src[LINE_LEN];
  char xc[LINE_LEN/10];

  src[0]='\0';
  fgets(src,LINE_LEN, f_src); //skip the first line about parameter names

  while (!feof(f_src) && (index<file_Len)) //index = 0, 1, ..., file_Len-1
    {
      src[0]='\0';	    
      fgets(src,LINE_LEN,f_src);
      i=0;
      time_pass=0;
						
      if (time_ID==-1)
        {
          for (t=0;t<num_dm;t++) //there is no time_ID
            {
              xc[0]='\0';
              j=0;
              while ((src[i]!='\0') && (src[i]!='\n') && (src[i]!=' ') && (src[i]!='\t'))
                {
                  xc[j]=src[i];
                  i++;
                  j++;
                }
		
              xc[j]='\0';	    
              i++;

              uncomp_data[index][t]=atof(xc);
            }	
        }
      else
        {
          for (t=0;t<=num_dm;t++) //the time column needs to be skipped, so there are num_dm+1 columns
            {
              xc[0]='\0';
              j=0;
              while ((src[i]!='\0') && (src[i]!='\n') && (src[i]!=' ') && (src[i]!='\t'))
                {
                  xc[j]=src[i];
                  i++;
                  j++;
                }
		
              xc[j]='\0';	    
              i++;

              if (t==time_ID)
                {
                  time_pass=1;
                  continue;
                }
				
              if (time_pass)
                uncomp_data[index][t-1]=atof(xc);
              else
                uncomp_data[index][t]=atof(xc);
            }
        }        	
      index++;     	
      //fprintf(fout_ID,"%s",src);
    } //end of while
	
  if (DEBUG == 1)
    {
      printf("the last line of the source data is:\n");
      for (j=0;j<num_dm;j++)
        printf("%f ",uncomp_data[index-1][j]);
      printf("\n");
    }
}


/**************************************** Normalization ******************************************/
void tran(double **orig_data, int file_Len, int num_dm, int norm_used, double **matrix_to_cluster, double *largest_value, double *smallest_value)
{
  int i=0;
  int j=0;

  //double biggest=-MAX_VALUE;
  //double smallest=MAX_VALUE;

  double *aver; //average of each column
  double *std; //standard deviation of each column

  aver=(double*)malloc(sizeof(double)*file_Len);
  memset(aver,0,sizeof(double)*file_Len);

  std=(double*)malloc(sizeof(double)*file_Len);
  memset(std,0,sizeof(double)*file_Len);	
		
  for (j=0;j<num_dm;j++)
  {
       largest_value[j]=-MAX_VALUE;
	   smallest_value[j]=MAX_VALUE;
       for (i=0;i<file_Len;i++)
       {
            if (orig_data[i][j]>largest_value[j])
                largest_value[j]=orig_data[i][j];
            if (orig_data[i][j]<smallest_value[j])
                smallest_value[j]=orig_data[i][j];
       }
  }

  if (norm_used==2) //z-score normalization
    {
      for (j=0;j<num_dm;j++)
        {
          aver[j]=0;
          for (i=0;i<file_Len;i++)
            aver[j]=aver[j]+orig_data[i][j];
          aver[j]=aver[j]/(double)file_Len;

          std[j]=0;
          for (i=0;i<file_Len;i++)
            std[j]=std[j]+(orig_data[i][j]-aver[j])*(orig_data[i][j]-aver[j]);
          std[j]=sqrt(std[j]/(double)file_Len);
			
          for (i=0;i<file_Len;i++)
            matrix_to_cluster[i][j]=(orig_data[i][j]-aver[j])/std[j];  //z-score normalization
        }
    }

  if (norm_used==1) //0-1 min-max normalization
    {
      for (j=0;j<num_dm;j++)
        {
            for (i=0;i<file_Len;i++)
            {
              if (largest_value[j]==smallest_value[j])
                matrix_to_cluster[i][j]=0.5;
              else
                matrix_to_cluster[i][j]=(orig_data[i][j]-smallest_value[j])/(largest_value[j]-smallest_value[j]);
            }
        }
    }

  if (norm_used==3) //0-1 0-4095 normalization
    {
      for (j=0;j<num_dm;j++)
            for (i=0;i<file_Len;i++)
                matrix_to_cluster[i][j]=orig_data[i][j]/(double)UPPER_VALUE;
    }

  if (norm_used==0) //no normalization
    {
      for (i=0;i<file_Len;i++)
        for (j=0;j<num_dm;j++)
          matrix_to_cluster[i][j]=orig_data[i][j];
    }

  free(aver);
  free(std);
}



/********************************************** RadixSort *******************************************/
/* Perform a radix sort using each dimension from the original data as a radix.
 * Outputs:
 * sorted_seq   -- a permutation vector mapping the ordered list onto the original data.
 *                  (sorted_seq[i] -> index in the original data of the ith element of the ordered list)
 * grid_ID      -- mapping between the original data and the "grids" (see below) found as a byproduct
 *                  of the sorting procedure.
 * num_nonempty -- the number of grids that occur in the data (= the number of distinct values assumed
 *                  by grid_ID)
 */

void radixsort_flock(int **position,int file_Len,int num_dm,int num_bin,int *sorted_seq,int *num_nonempty,int *grid_ID)
{
  int i=0;
  int length=0; 
  int start=0;
  int prev_ID=0;
  int curr_ID=0;
	
  int j=0;
  int t=0;
  int p=0;
  int loc=0;
  int temp=0;
  int equal=0;
	
  int *count; //count[i]=j means there are j numbers having value i at the processing digit
  int *index; //index[i]=j means the starting position of grid i is j
  int *cp; //current position
  int *mark; //mark[i]=0 means it is not an ending point of a part, 1 means it is (a "part" is a group of items with identical bins for all dimensions)
  int *seq; //temporary sequence

  count=(int*)malloc(sizeof(int)*num_bin);
  memset(count,0,sizeof(int)*num_bin);

  cp=(int*)malloc(sizeof(int)*num_bin);
  memset(cp,0,sizeof(int)*num_bin);

  index=(int*)malloc(sizeof(int)*num_bin); // initialized below

  seq=(int*)malloc(sizeof(int)*file_Len);
  memset(seq,0,sizeof(int)*file_Len);
	
  mark=(int*)malloc(sizeof(int)*file_Len);
  memset(mark,0,sizeof(int)*file_Len);
	
  for (i=0;i<file_Len;i++)
    {
      sorted_seq[i]=i;
      mark[i]=0;
      seq[i]=0;
    }
  for (i=0;i<num_bin;i++)
    {
      index[i]=0;
      cp[i]=0;
      count[i]=0;
    }

  for (j=0;j<num_dm;j++)
    {
      if (j==0) //compute the initial values of mark
        {
          for (i=0;i<file_Len;i++)
            count[position[i][j]]++; // initialize the count to the number of items in each bin of the 0th dimension

          index[0] = 0;
          for (i=0;i<num_bin-1;i++)
            {
              index[i+1]=index[i]+count[i];  //index[k]=x means k segment starts at x (in the ordered list)
              if ((index[i+1]>0) && (index[i+1]<=file_Len))
                {
                  mark[index[i+1]-1]=1; // Mark the end of the segment in the ordered list
                }
              else
                {
                  printf("out of myboundary for mark at index[i+1]-1.\n");
                }
            }
          mark[file_Len-1]=1;
			
          for (i=0;i<file_Len;i++)
            {
              /* Build a permutation vector for the partially ordered data.  Store the PV in sorted_seq */
              loc=position[i][j];
              temp=index[loc]+cp[loc]; //cp[i]=j means the offset from the starting position of grid i is j 
              sorted_seq[temp]=i;  //sorted_seq[i]=temp is also another way to sort
              cp[loc]++;
            }
        }
      else
        {
          //reset count, index, loc, temp, cp, start, and length
          length=0;
          loc=0;
          temp=0;
          start=0;
          for (p=0;p<num_bin;p++)
            {
              cp[p]=0;
              count[p]=0;
              index[p]=0;
            }

          for (i=0;i<file_Len;i++)
            {
              int iperm = sorted_seq[i]; // iperm allows us to traverse the data in sorted order.
              if (mark[i]!=1)
                {
                  /* Count the number of items in each bin of
                     dimension j, BUT we are going to reset at the end
                     of each "part".  Thus, the total effect is to do
                     a sort by bin on the jth dimension for each group
                     of data that has been identical for the
                     dimensions processed up to this point.  This is
                     the standard radix sort procedure, but doing it
                     this way saves us having to allocate buckets to
                     hold the data in each group of "identical-so-far"
                     elements. */
                  count[position[iperm][j]]++;  //count[position[i][j]]++;
                  length++;                     // This is the total length of the part, irrespective of the value of the jth component
                                                // (actually, less one, since we don't increment for the final element below)
                }
              if (mark[i]==1)
                {
                  //length++;
                  count[position[iperm][j]]++;//count[position[i][j]]++;  //the current point must be counted in
                  start=i-length; //this part starts from start to i: [start,i]
                  /* Now we sort on the jth radix, just like we did for the 0th above, but we restrict it to just the current part.
                     This would be a lot more clear if we broke this bit of code out into a separate function and processed recursively,
                     plus we could multi-thread over the parts.  (Hmmm...)
                  */
                  index[0] = start; // Let index give the offset within the whole ordered list.
                  for (t=0;t<num_bin-1;t++)
                    {
                      index[t+1]=index[t]+count[t];
						
                      if ((index[t+1]<=file_Len) && (index[t+1]>0))
                        {
                          mark[index[t+1]-1]=1; // update the part boundaries to include the differences in the current radix.
                        }
						
                    }
                  mark[i]=1;

                  /* Update the permutation vector for the current part (i.e., from start to i).  By the time we finish the loop over i
                     the PV will be completely updated for the partial ordering up to the current radix. */
                  for (t=start;t<=i;t++)
                    {
                      loc=position[sorted_seq[t]][j];//loc=position[t][j];
                      temp=index[loc]+cp[loc];
                      if ((temp<file_Len) && (temp>=0)) 
                        {
                          // seq is a temporary because we have to maintain the old PV until we have finished this step.
                          seq[temp]=sorted_seq[t];  //sorted_seq[i]=temp is also another way to sort
                          cp[loc]++;
                        }
                      else
                        {
                          printf("out of myboundary for seq at temp.\n");
                        }
                    }

                  for (t=start;t<=i;t++)
                    {
                      // copy the temporary back into sorted_seq.  sorted_seq is now updated for radix j up through
                      // entry i in the ordered list.
                      if ((t>=0) && (t<file_Len))
                        sorted_seq[t]=seq[t];
                      else
                        printf("out of myboundary for seq and sorted_seq at t.\n");
                    }
                  //reset count, index, seq, length, and cp
                  length=0;
                  loc=0;
                  temp=0;
                  for (p=0;p<num_bin;p++)
                    {
                      cp[p]=0;
                      count[p]=0;
                      index[p]=0;
                    }
                }
            }//end for i
        }//end else
    }//end for j

  /* sorted_seq[] now contains the ordered list for all radices.  mark[] gives the boundaries between groups of elements that are
     identical over all radices (= dimensions in the original data) (although it appears we aren't going to make use of this fact) */
  
  for (i=0;i<file_Len;i++)
    grid_ID[i]=0; //in case the initial value hasn't been assigned
  *num_nonempty=1; //starting from 1!	

  /* assign the "grid" identifiers for all of the data.  A grid will be what we were calling a "part" above.  We will number them
     serially and tag the *unordered* data with the grid IDs.  We will also count the number of populated grids (in general there will
     be many possible combinations of bin values that simply never occur) */
  
  for (i=1;i<file_Len;i++)
    {
      equal=1;
      prev_ID=sorted_seq[i-1];
      curr_ID=sorted_seq[i];
      for (j=0;j<num_dm;j++)
        {
          if (position[prev_ID][j]!=position[curr_ID][j])
            {	
              equal=0;  //not equal
              break;
            }
        }
		
      if (equal)
        {
          grid_ID[curr_ID]=grid_ID[prev_ID];
        }
      else
        {
          *num_nonempty=*num_nonempty+1;
          grid_ID[curr_ID]=grid_ID[prev_ID]+1;
        }
      //all_grid_vol[grid_ID[curr_ID]]++;
    }

  //free memory
  free(count);
  free(index);	
  free(cp);	
  free(seq);
  free(mark); 
  
}

/********************************************** Compute Position of Events ************************************************/
void compute_position(double **data_in, int file_Len, int num_dm, int num_bin, int **position, double *interval)
{
  /* What we are really doing here is binning the data, with the bins
     spanning the range of the data and number of bins = num_bin */
  int i=0;
  int j=0;

  double *small; //small[j] is the smallest value within dimension j
  double *big; //big[j] is the biggest value within dimension j
		
  small=(double*)malloc(sizeof(double)*num_dm);
  memset(small,0,sizeof(double)*num_dm);

  big=(double*)malloc(sizeof(double)*num_dm);
  memset(big,0,sizeof(double)*num_dm);
	
	
  for (j=0;j<num_dm;j++)
    {
      big[j]=-MAX_VALUE;
      small[j]=MAX_VALUE;
      for (i=0;i<file_Len;i++)
        {
          if (data_in[i][j]>big[j])
            big[j]=data_in[i][j];

          if (data_in[i][j]<small[j])
            small[j]=data_in[i][j];
        }
		
      interval[j]=(big[j]-small[j])/(double)num_bin;	//interval is computed using the biggest value and smallest value instead of the channel limit
      /* XXX: I'm pretty sure the denominator of the fraction above should be num_bin-1. */
	  /* I don't think so: num_bin is the number of bins */
    }
    
  for (j=0;j<num_dm;j++)
  {
     for (i=0;i<file_Len;i++)
     {	
        if (data_in[i][j]>=big[j])
           position[i][j]=num_bin-1;
        else
        {
           position[i][j]=(int)((data_in[i][j]-small[j])/interval[j]); //position[i][j]=t means point i is at the t grid of dimensional j
           if ((position[i][j]>=num_bin) || (position[i][j]<0))
           {
               //printf("position mis-computed in density analysis!\n");
               //exit(0);
			   fprintf(stderr,"Incorrect input file format or input parameters (number of bins overflows)!\n"); //modified on July 23, 2010
				abort();
			   
           }
        }
     }
  }


  free(small);
  free(big);
}

/********************************************** select_bin to select the number of bins **********************************/
//num_bin=select_bin(normalized_data, file_Len, num_dm, MIN_GRID, MAX_GRID, position, sorted_seq, all_grid_ID, &num_nonempty);
/* Determine the number of bins to use in each dimension.  Additionally sort the data elements according to the binned
 * values, and partition the data into "grids" with identical (binned) values.  We try progressively more bins until we
 * maximize a merit function, then return the results obtained using the optimal number of bins. 
 *
 * Outputs:
 * position     -- binned data values
 * sorted_seq   -- permutation vector mapping the ordered list to the original data
 * all_grid_ID  -- grid to which each data element was assigned.
 * num_nonempty -- number of distinct values assumed by all_grid_ID
 * interval     -- bin width for each data dimension
 * return value -- the number of bins selected.
 */

int select_bin(double **normalized_data, int file_Len, int num_dm, int min_grid, int max_grid, int **position, int *sorted_seq, 
                int *all_grid_ID, int *num_nonempty, double *interval, int user_num_bin)
{
 
  int num_bin=0;
  int select_num_bin=0;
  int m=0;
  int n=0;
	
  int i=0;
  int bin_scope=0;
  int temp_num_nonempty=0;

  int *temp_grid_ID;
  int *temp_sorted_seq;
  int **temp_position;

  //sorted_seq[i]=j means the event j ranks i

  double temp_index=0;
  double *bin_index;	
  double *temp_interval;

 
	
  temp_grid_ID=(int *)malloc(sizeof(int)*file_Len);
  memset(temp_grid_ID,0,sizeof(int)*file_Len);
	
  temp_sorted_seq=(int *)malloc(sizeof(int)*file_Len);
  memset(temp_sorted_seq,0,sizeof(int)*file_Len);

  temp_position=(int **)malloc(sizeof(int*)*file_Len);
  memset(temp_position,0,sizeof(int*)*file_Len);
  for (m=0;m<file_Len;m++)
    {
      temp_position[m]=(int*)malloc(sizeof(int)*num_dm);
      memset(temp_position[m],0,sizeof(int)*num_dm);
    }

  temp_interval=(double*)malloc(sizeof(double)*num_dm);
  memset(temp_interval,0,sizeof(double)*num_dm);

  bin_scope=max_grid-min_grid+1;
  bin_index=(double *)malloc(sizeof(double)*bin_scope);
  memset(bin_index,0,sizeof(double)*bin_scope);

  i=0;

  for (num_bin=min_grid;num_bin<=max_grid;num_bin++)
    {
      /* compute_position bins the data into num_bin bins.  Each
         dimension is binned independently.

         Outputs:
         temp_position[i][j] -- bin for the jth component of data element i.
         temp_interval[j]    -- bin-width for the jth component
      */
      compute_position(normalized_data, file_Len, num_dm, num_bin, temp_position, temp_interval);
      radixsort_flock(temp_position,file_Len,num_dm,num_bin,temp_sorted_seq,&temp_num_nonempty,temp_grid_ID);

      /* our figure of merit is the number of non-empty grids divided by number of bins per dimension.
         We declare victory when we have found a local maximum */
      bin_index[i]=((double)temp_num_nonempty)/((double)num_bin);
	  if ((double)(temp_num_nonempty)>=(double)(file_Len)*0.95)
		  break;
      if ((bin_index[i]<temp_index) && (user_num_bin==0))
         break;
	  if ((user_num_bin==num_bin-1) && (user_num_bin!=0))
		 break;
      
      /* Since we have accepted this trial bin, copy all the temporary results into
         the output buffers */
      memcpy(all_grid_ID,temp_grid_ID,sizeof(int)*file_Len);
      memcpy(sorted_seq,temp_sorted_seq,sizeof(int)*file_Len);
      memcpy(interval,temp_interval,sizeof(double)*num_dm);
		
      for (m=0;m<file_Len;m++)
        for (n=0;n<num_dm;n++)
          position[m][n]=temp_position[m][n];

      temp_index=bin_index[i];
      select_num_bin=num_bin;
      num_nonempty[0]=temp_num_nonempty;
      i++;
    }

   if ((select_num_bin<min_grid) || (select_num_bin>max_grid))
  {
    fprintf(stderr,"Number of events collected is too few in terms of number of markers used. The file should not be processed!\n"); //modified on Nov 4, 2010
	exit(0);
  }
	
  if (temp_index==0)
  {
	 fprintf(stderr,"Too many dimensions with too few events in the input file, or a too large number of bins used.\n"); //modified on July 23, 2010
	 abort();
  }
 

  free(temp_grid_ID);
  free(temp_sorted_seq);
  free(bin_index);
  free(temp_interval);
	
  for (m=0;m<file_Len;m++)
    free(temp_position[m]);
  free(temp_position);

  return select_num_bin; 
}

/************************************* Select dense grids **********************************/
// compute num_dense_grids, num_dense_events, dense_grid_reverse, and all_grid_vol
// den_cutoff=select_dense(file_Len, all_grid_ID, num_nonempty, &num_dense_grids, &num_dense_events, dense_grid_reverse);
/*
 * Prune away grids that are insufficiently "dense" (i.e., contain too few data items)
 *
 * Outputs:
 * num_dense_grids    -- number of dense grids
 * num_dense_events   -- total number of data items in all dense grids
 * dense_grid_reverse -- mapping from list of all grids to list of dense grids.
 * return value       -- density cutoff for separating dense from non-dense grids.
 */

int select_dense(int file_Len, int *all_grid_ID, int num_nonempty, int *num_dense_grids, int *num_dense_events, int *dense_grid_reverse, int den_t_event, int max_num_clusters)
	//max_num_clusters will be used to limit the number of dense hyper grids when selecting the density threshold
{
  

  int i=0;
  int vol_ID=0;
  int biggest_size=0; //biggest grid_size, to define grid_size_index
  //int biggest_index=0;
  //int actual_threshold=0; //the actual threshold on grid_size, e.g., 1 implies 1 event in the grid
  //int num_remain=0; //number of remaining grids with different density thresholds
  int temp_num_dense_grids=0;
  int temp_num_dense_events=0;
	
  int *grid_size_index;
  int *all_grid_vol;
  //int *grid_density_index;

  //double den_average=0;
 // double avr_index=0;
 

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute all_grid_vol
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  all_grid_vol=(int *)malloc(sizeof(int)*num_nonempty);
  memset(all_grid_vol,0,sizeof(int)*num_nonempty);

  /* Grid "volume" is just the number of data contained in the grid. */
  for (i=0;i<file_Len;i++)
    {
      vol_ID=all_grid_ID[i]; //vol_ID=all_grid_ID[sorted_seq[i]];
      all_grid_vol[vol_ID]++;  //all_grid_vol[i]=j means grid i has j points
    }

 
	
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute grid_size_index (histogram of grid sizes)
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  for (i=0;i<num_nonempty;i++)
    {
      if (biggest_size<all_grid_vol[i])
        {
          biggest_size=all_grid_vol[i];
        }
    }

  //added on July 23, 2010
  if (biggest_size<DEFAULT_DENSITY_T) //there is no hyper grid with 3 or more events
  {
	 fprintf(stderr,"Too many dimensions with too few events in the input file, or a too large number of bins used.\n"); //modified on July 23, 2010
	 abort();
  }

  grid_size_index=(int*)malloc(sizeof(int)*biggest_size);
  memset(grid_size_index,0,sizeof(int)*biggest_size);
	
  for (i=0;i<num_nonempty;i++)
    {
      grid_size_index[all_grid_vol[i]-1]++; //grid_size_index[0]=5 means there are 5 grids having size 1
    }



  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute den_cutoff
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  //grid_density_index=(int *)malloc(sizeof(int)*(biggest_size-2));//from event 2 to biggest_size-1, i.e., from 0 to biggest_size-3
  //memset(grid_density_index,0,sizeof(int)*(biggest_size-2));

  if (den_t_event==0)  //user doesn't define the density threshold
  {
	  for (i=2;i<biggest_size;i++) //the grid with 1 or 2 event will be skipped, i.e., grid_density_index[0] and [1] won't be defined
	  {
		  if ((grid_size_index[i]<(NUM_CLUSTER_FOLD*max_num_clusters)) && (grid_size_index[i]>0)) 
		  {
			den_t_event=i+1;
			break;
		  }		  
	  }
  }

  if (den_t_event==0) //if something is wrong, e.g., all hyper-grids only have 1 event
	  den_t_event=DEFAULT_DENSITY_T;

  for (i=0;i<num_nonempty;i++)
	  if (all_grid_vol[i]>=den_t_event) //the biggest number of events is allowed; otherwise, should be all_grid_vol[i]>den_t_event
		temp_num_dense_grids++;

  if (temp_num_dense_grids<DEFAULT_DENSITY_T)
  {
	  //modified on July 23, 2010
	  //printf("a too high density threshold is set! Please decrease your density threshold.\n");
	  //exit(0);
	  fprintf(stderr,"a too high density threshold is set! Please decrease your density threshold.\n"); //modified on July 23, 2010
	  abort();
  }

  if (temp_num_dense_grids>UPPER_NUM_REGION)
  {
	  //modified on July 23, 2010
	  //printf("a too low density threshold is set! Please increase your density threshold.\n");
	  //exit(0);
	  fprintf(stderr,"a too low density threshold is set! Please increase your density threshold.\n"); //modified on July 23, 2010
	  abort();
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Compute dense_grid_reverse
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  temp_num_dense_grids=0;
  temp_num_dense_events=0;
  for (i=0;i<num_nonempty;i++)
    {
      dense_grid_reverse[i]=-1;
		
      if (all_grid_vol[i]>=den_t_event) 
        {			
          dense_grid_reverse[i]=temp_num_dense_grids;  //dense_grid_reverse provides a mapping from all nonempty grids to dense grids.
          temp_num_dense_grids++;
          temp_num_dense_events+=all_grid_vol[i];						
        }
    }

  num_dense_grids[0]=temp_num_dense_grids;
  num_dense_events[0]=temp_num_dense_events;	

  free(grid_size_index);
  free(all_grid_vol);

  return den_t_event;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute densegridID_To_gridevent and eventID_To_denseventID
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	grid_To_event(file_Len, dense_grid_reverse, all_grid_ID, eventID_To_denseventID, densegridID_To_gridevent);
/*
 * Filter the data so that only the data belonging to dense grids is left
 *
 * Output:
 * eventID_To_denseeventID   -- mapping from original event ID to ID in list containing only events contained in dense grids.
 * densegridID_To_gridevent  -- mapping from dense grids to prototype members of the grids.
 *
 */

void grid_To_event(int file_Len, int *dense_grid_reverse, int *all_grid_ID, int *eventID_To_denseventID, int *densegridID_To_gridevent)
{
  int i=0;
  int dense_grid_ID=0;
  int dense_event_ID=0;

  for (i=0;i<file_Len;i++)
    {
      dense_grid_ID=dense_grid_reverse[all_grid_ID[i]];
      eventID_To_denseventID[i]=-1;
      if (dense_grid_ID!=-1) //for point (i) belonging to dense grids 
        {			
          eventID_To_denseventID[i]=dense_event_ID;
          dense_event_ID++;
		
          if (densegridID_To_gridevent[dense_grid_ID]==-1) //for point i that hasn't been selected
            densegridID_To_gridevent[dense_grid_ID]=i; //densegridID_To_gridevent maps dense_grid_ID to its representative point			
        }		
    }
	
	
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute dense_grid_seq
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	generate_grid_seq(file_Len, num_dm, sorted_seq, num_dense_grids, densegridID_To_gridevent, position, dense_grid_rank, dense_grid_seq);
/* Construct a table of binned data values for each dense grid.
 *
 * Output:
 *
 * dense_grid_seq  -- table of binned data values for each dense grid (recall that all members of a grid have identical binned data values)
 */

void generate_grid_seq(int num_dm, int num_dense_grids, int *densegridID_To_gridevent, int **position, int **dense_grid_seq)
{  

  int i=0;
  int j=0;
  int ReventID=0; //representative event ID of the dense grid

  for (i=0;i<num_dense_grids;i++)
    {
      ReventID = densegridID_To_gridevent[i];

      for (j=0;j<num_dm;j++)
        dense_grid_seq[i][j]=position[ReventID][j];
    }
}

//compare two vectors
int compare_value(int num_dm, int *search_value, int *seq_value)
{
  int i=0;

  for (i=0;i<num_dm;i++)
    {
      if (search_value[i]<seq_value[i])
        return 1;
      if (search_value[i]>seq_value[i])
        return -1;
      if (search_value[i]==seq_value[i])
        continue;
    }
  return 0;
}

//binary search the dense_grid_seq to return the dense grid ID if it exists
int binary_search(int num_dense_grids, int num_dm, int *search_value, int **dense_grid_seq)
{

  int low=0;
  int high=0;
  int mid=0;
  int comp_result=0;
  int match=0;
  //int found=0;
	
  low = 0;
  high = num_dense_grids-1;
    
  while (low <= high) 
    {
      mid = (int)((low + high)/2);
	
      comp_result=compare_value(num_dm, search_value,dense_grid_seq[mid]);
	
		
      switch(comp_result)
        {
        case 1:
          high=mid-1;
          break;
        case -1:
          low=mid+1;
          break;
        case 0:
          match=1;
          break;
        }
      if (match==1)
        break;
    }
	


  if (match==1)
    {
      return mid;
    }
  else
    return -1;   
}


/********************************************** Computing Centers Using IDs **********************************************/

void ID2Center(double **data_in, int file_Len, int num_dm, int *eventID_To_denseventID, int num_clust, int *cluster_ID, double **population_center)
{
  int i=0;
  int j=0;
  int ID=0;
  int eventID=0;
  int *size_c;



  size_c=(int *)malloc(sizeof(int)*num_clust);
  memset(size_c,0,sizeof(int)*num_clust);

  for (i=0;i<num_clust;i++)
    for (j=0;j<num_dm;j++)
      population_center[i][j]=0;

  for (i=0;i<file_Len;i++)
    {
      eventID=eventID_To_denseventID[i];

      if (eventID!=-1) //only events in dense grids count
        {
          ID=cluster_ID[eventID];
			
          if (ID==-1)
            {
              //modified on July 23, 2010
			  //printf("ID==-1! in ID2Center\n");
              //exit(0);
			  fprintf(stderr,"Incorrect file format or input parameters (no dense regions found!)\n"); //modified on July 23, 2010
			  abort();
            }

          for (j=0;j<num_dm;j++)
            population_center[ID][j]=population_center[ID][j]+data_in[i][j];
			
          size_c[ID]++;				
        }
    }
	

  for (i=0;i<num_clust;i++)
    {
      for (j=0;j<num_dm;j++)
        if (size_c[i]!=0)
          population_center[i][j]=(population_center[i][j]/(double)(size_c[i]));
        else
          population_center[i][j]=0;
    }

  free(size_c);

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Compute Population Center with all events
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ID2Center_all(double **data_in, int file_Len, int num_dm, int num_clust, int *cluster_ID, double **population_center)
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
         ID=cluster_ID[i];
			
         if (ID==-1)
         {
            //commented on July 23, 2010
			//printf("ID==-1! in ID2Center_all\n");
            //exit(0);
			fprintf(stderr,"Incorrect file format or input parameters (resulting in incorrect population IDs)\n"); //modified on July 23, 2010
			abort();
         }

         for (j=0;j<num_dm;j++)
           population_center[ID][j]=population_center[ID][j]+data_in[i][j];
			
         size_c[ID]++;        
    }
	

  for (i=0;i<num_clust;i++)
    {
      for (j=0;j<num_dm;j++)
        if (size_c[i]!=0)
          population_center[i][j]=(population_center[i][j]/(double)(size_c[i]));
        else
          population_center[i][j]=0;
    }

  free(size_c);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Merge neighboring grids to clusters
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int merge_grids(double **normalized_data, double *interval, int file_Len, int num_dm, int num_bin, int **position, int num_dense_grids, 
                 int *dense_grid_reverse, int **dense_grid_seq, int *eventID_To_denseventID, int *densegridID_To_gridevent, int *all_grid_ID,
                 int *cluster_ID, int *grid_ID, int *grid_clusterID)
{
  

  int i=0;
  int j=0;
  int t=0;
  int p=0;
  int num_clust=0;
  int ReventID=0;
  int denseID=0;
  int neighbor_ID=0;
  //int temp_grid=0;

  int *grid_value;
  int *search_value;
	
  int **graph_of_grid; //the graph constructed for dense grids: each dense grid is a graph node

  double real_dist=0;
  double **norm_grid_center;
	
  norm_grid_center=(double **)malloc(sizeof(double*)*num_dense_grids);
  memset(norm_grid_center,0,sizeof(double*)*num_dense_grids);
	
  for (i=0;i<num_dense_grids;i++)
    {
      norm_grid_center[i]=(double *)malloc(sizeof(double)*num_dm);
      memset(norm_grid_center[i],0,sizeof(double)*num_dm);
    }

  for (i=0;i<file_Len;i++)
    {
      denseID=eventID_To_denseventID[i];
      if (denseID!=-1) //only dense events can enter
        {
          grid_ID[denseID]=dense_grid_reverse[all_grid_ID[i]];
			
          if (grid_ID[denseID]==-1)
            {
              fprintf(stderr,"Incorrect input file format or input parameters (no dense region found)!\n"); //modified on July 23, 2010
              abort();
            }			
        }
    }

	
  /* Find centroid (in the normalized data) for each dense grid */
  ID2Center(normalized_data,file_Len,num_dm,eventID_To_denseventID,num_dense_grids,grid_ID,norm_grid_center);	
 
  //printf("pass the grid ID2 center\n"); //commmented on July 23, 2010

	
  graph_of_grid=(int **)malloc(sizeof(int*)*num_dense_grids);
  memset(graph_of_grid,0,sizeof(int*)*num_dense_grids);
  for (i=0;i<num_dense_grids;i++)
    {
      graph_of_grid[i]=(int *)malloc(sizeof(int)*num_dm);
      memset(graph_of_grid[i],0,sizeof(int)*num_dm);		
		

      for (j=0;j<num_dm;j++)
        graph_of_grid[i][j]=-1;
    }	
	
  grid_value=(int *)malloc(sizeof(int)*num_dm);
  memset(grid_value,0,sizeof(int)*num_dm);

  search_value=(int *)malloc(sizeof(int)*num_dm);
  memset(search_value,0,sizeof(int)*num_dm);

  
  for (i=0;i<num_dense_grids;i++)
    {
      ReventID=densegridID_To_gridevent[i];

      for (j=0;j<num_dm;j++)
        {
          grid_value[j]=position[ReventID][j];
          
        }
     

      /* For each dimension, find the neighbor, if any, that is equal in all other dimensions and 1 greater in
         the chosen dimension.  If the neighbor's centroid is not too far away, add it to this grid's neighbor
         list. */
      for (t=0;t<num_dm;t++)
        {
          for (p=0;p<num_dm;p++)
            search_value[p]=grid_value[p];

          if (grid_value[t]==num_bin-1)
            continue;

          search_value[t]=grid_value[t]+1; //we only consider the neighbor at the bigger side

          neighbor_ID=binary_search(num_dense_grids, num_dm, search_value, dense_grid_seq);
			
          if (neighbor_ID!=-1) 
            {
              real_dist=norm_grid_center[i][t]-norm_grid_center[neighbor_ID][t];
	
              if (real_dist<0)
                real_dist=-real_dist;
				
              if (real_dist<2*interval[t])
                graph_of_grid[i][t]=neighbor_ID;			
            }
        }
      grid_clusterID[i]=i; //initialize grid_clusterID
    } 
	
	
  //graph constructed
  //DFS-based search begins

  /* Use a depth-first search to construct a list of connected subgraphs (= "clusters").
     Note our graph as we have constructed it is a DAG, so we can use that to our advantage
     in our search. */
  //  num_clust=dfs(graph_of_grid,num_dense_grids,num_dm,grid_clusterID);
  num_clust=find_connected(graph_of_grid, num_dense_grids, num_dm, grid_clusterID);

	
  //computes grid_ID and cluster_ID
  for (i=0;i<file_Len;i++)
    {
      denseID=eventID_To_denseventID[i];
      if (denseID!=-1) //only dense events can enter
	  {
        cluster_ID[denseID]=grid_clusterID[grid_ID[denseID]];
		//if (cluster_ID[denseID]==-1)
		//	printf("catch you!\n");
	  }
    }
	
  free(search_value);
  free(grid_value);

  for (i=0;i<num_dense_grids;i++)
    {
      free(graph_of_grid[i]);
      free(norm_grid_center[i]);
    }
  free(graph_of_grid);
  free(norm_grid_center);

  return num_clust;
}

/********************************************* Merge Clusters to Populations *******************************************/
// This is the function future work can be on because it is about how to cluster the about 500 points accurately

int merge_clusters(int num_clust, int num_dm, double *interval, double **cluster_center, int *cluster_populationID, int max_num_pop)
{
  int num_population=0;
  int temp_num_population=0;

  int i=0;
  int j=0;
  int t=0;
  int merge=0;
  int smid=0;
  int bgid=0;
  double merge_dist=1.1; //initial value of merge_dist*interval should be slightly larger than the bin width

  int *map_ID;

  double diff=0;  

  map_ID=(int*)malloc(sizeof(int)*num_clust);
  memset(map_ID,0,sizeof(int)*num_clust);

  while ((num_population>max_num_pop) || (num_population<=1))
  {
  
	  if (num_population<=1)
	  	  merge_dist=merge_dist-0.1;

	  if (num_population>max_num_pop)
          merge_dist=merge_dist+0.1;
	  
	    
	 for (i=0;i<num_clust;i++)
		cluster_populationID[i]=i;

    for (i=0;i<num_clust-1;i++)
    {
      for (j=i+1;j<num_clust;j++)
        {
          merge=1;
			
          for (t=0;t<num_dm;t++)
            {
              diff=cluster_center[i][t]-cluster_center[j][t];
				
              if (diff<0)
                diff=-diff;
              if (diff>(merge_dist*interval[t]))
                merge=0;
            }

          if ((merge) && (cluster_populationID[i]!=cluster_populationID[j]))
            {
              if (cluster_populationID[i]<cluster_populationID[j])  //they could not be equal
                {
                  smid = cluster_populationID[i];
                  bgid = cluster_populationID[j];
                }
              else
                {
                  smid = cluster_populationID[j];
                  bgid = cluster_populationID[i];
                }
              for (t=0;t<num_clust;t++)
                {
                  if (cluster_populationID[t]==bgid)
                    cluster_populationID[t]=smid;
                }
            }
        }
    }

  

  for (i=0;i<num_clust;i++)
    map_ID[i]=-1;

  for (i=0;i<num_clust;i++)
    map_ID[cluster_populationID[i]]=1;

  num_population=0;
  for (i=0;i<num_clust;i++)
    {
      if (map_ID[i]==1)
        {
          map_ID[i]=num_population;
          num_population++;
        }
    }

  if ((temp_num_population>max_num_pop) && (num_population==1))
	  break;
  else
	  temp_num_population=num_population;

  if (num_clust<=1)
	break;
  }

  for (i=0;i<num_clust;i++)
    cluster_populationID[i]=map_ID[cluster_populationID[i]];
	
  free(map_ID);

  return num_population; 
}

///////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double kmeans(double **Matrix, int k, double kmean_term, int file_Len, int num_dm, int *shortest_id, double **center)
{
	
	int i=0;
	int j=0;
	int t=0;
	int random=0;
	int random1=0;
	int random2=0;
	int times=0;
	int dist_used=0; //0 is Euclidean and 1 is Pearson
	int random_init=0; //0: not use random seeds; 
	int real_Len=0;
	int skipped=0;
	
 	int *num;  //num[i]=t means the ith cluster has t points
	
	double vvv=1.0; // the biggest variation;
	double distance=0.0;
	double xv=0.0;
	double variation=0.0;

	double mean_dx=0;
	double mean_dy=0;
	double sum_var=0;
	double dx=0;
	double dy=0;
	double sd_x=0;
	double sd_y=0;	
	double diff=0;
	double distortion=0;

	double shortest_distance;

	double *temp_center;	
			
	double **sum;	
 
	temp_center = (double *)malloc(sizeof(double)*num_dm);
	memset(temp_center,0,sizeof(double)*num_dm);

	if (random_init)
	{
		for (i=0;i<k;i++)
		{	
			random1=rand()*rand();
			random2=abs((random1%5)+1);
			for (t=0;t<random2;t++)
				random2=random2*rand()+rand();
	
			random=abs(random2%file_Len);
			for (j=0;j<num_dm;j++)
				center[i][j]=Matrix[random][j];				
		}
	}


	num = (int *)malloc(sizeof(int)*k);
	memset(num,0,sizeof(int)*k);

	sum = (double **)malloc(sizeof(double*)*k);
	memset(sum,0,sizeof(double*)*k);
	for (i=0;i<k;i++)
	{
		sum[i] = (double *)malloc(sizeof(double)*num_dm);
		memset(sum[i],0,sizeof(double)*num_dm);
	}


	times=0;
	real_Len=0;

	while (((vvv>kmean_term) && (kmean_term<1)) || ((times<kmean_term) && (kmean_term>=1)))
	{
		for (i=0;i<k;i++)
		{
			num[i]=0;
			for (j=0;j<num_dm;j++)
				sum[i][j]=0.0;  
		}
		
		for (i=0;i<file_Len;i++)  //for each data point i, we compute the distance between Matrix[i] and center[j]
		{		
			skipped = 0;
			shortest_distance=MAX_VALUE;
			for (j=0;j<k;j++)  //for each center j
			{
				distance=0.0;
								
				if (dist_used==0)  //Euclidean distance
				{
					for (t=0;t<num_dm;t++) //for each dimension here num_dm is always 1 as we consider individual dimensions
					{
					
						diff=center[j][t]-Matrix[i][t];
					
						diff=diff*diff;  
						
						distance = distance+diff; //here we have a weight for each dimension
					}
				}
				else  //pearson correlation
				{
					mean_dx=0.0;
					mean_dy=0.0;
					sum_var=0.0;
					dx=0.0;
					dy=0.0;
					sd_x=0.0;
					sd_y=0.0;
					for (t=0;t<num_dm;t++)
					{
						mean_dx+=center[j][t];
						mean_dy+=Matrix[i][t];
					}
					mean_dx=mean_dx/(double)num_dm;
					mean_dy=mean_dy/(double)num_dm;
			
					for (t=0;t<num_dm;t++)
					{
						dx=center[j][t]-mean_dx;
						dy=Matrix[i][t]-mean_dy;
						sum_var+=dx*dy;
						sd_x+=dx*dx;
						sd_y+=dy*dy;
					}
					if (sqrt(sd_x*sd_y)==0)
						distance = 1.0;
					else
						distance = 1.0 - (sum_var/(sqrt(sd_x*sd_y))); // distance ranges from 0 to 2;
					//printf("distance=%f\n",distance);			
				}	//pearson correlation ends 				

				if ((distance<shortest_distance) && (skipped == 0))
				{
					shortest_distance=distance;						
					shortest_id[i]=j;  
				}			
			}//end for j
				real_Len++;
				num[shortest_id[i]]=num[shortest_id[i]]+1;
				for (t=0;t<num_dm;t++)
					sum[shortest_id[i]][t]=sum[shortest_id[i]][t]+Matrix[i][t];		
		}//end for i
	/* recompute the centers */
	//compute_mean(group);		
		vvv=0.0;
		for (j=0;j<k;j++)
		{
			memcpy(temp_center,center[j],sizeof(double)*num_dm);
			variation=0.0;
			if (num[j]!=0)
			{
				for (t=0;t<num_dm;t++)
				{
					center[j][t]=sum[j][t]/(double)num[j];
					xv=(temp_center[t]-center[j][t]);
					variation=variation+xv*xv;
				}
			}
			
			if (variation>vvv)
				vvv=variation;  //vvv is the biggest variation among the k clusters;			
		}
	//compute_variation;
		times++;
	} //end for while
	

	free(num);

	for (i=0;i<k;i++)
		free(sum[i]);
	free(sum);
	free(temp_center);


	return distortion;
		
}

//////////////////////////////////////////////////////
/*************************** Show *****************************/
void show(double **Matrix, int *cluster_id, int file_Len, int k, int num_dm, char *name_string, double *big, double *small, int *constrained_d, int *constrained_line, int num_rows)
{
	int situ1=0;
	int situ2=0;

	int i=0;
	int id=0;
	int j=0;
	int t=0;
	
	int *size_c;
	


	int **size_mybound_1;
	int **size_mybound_2;
	int **size_mybound_3;
	int **size_mybound_0;

	double interval=0.0;

	//double *big;
	//double *small;


	//double **center;
	double **mybound;
	
	int **prof; //prof[i][j]=1 means population i is in the first quantile at parameter j
	
	int **size_f;
	double **percent_f;

	FILE *fpcnt_id; //proportion id
	FILE *fprof_id; //profile_id

	FILE *fnegpos_id; //negative or positive

	//big=(double *)malloc(sizeof(double)*num_dm);
	//memset(big,0,sizeof(double)*num_dm);

	//small=(double *)malloc(sizeof(double)*num_dm);
	//memset(small,0,sizeof(double)*num_dm);

	//for (i=0;i<num_dm;i++)
	//{
	//	big[i]=0.0;
	//	small[i]=(double)MAX_VALUE;
	//}
	
	
	size_c=(int *)malloc(sizeof(int)*k);
	memset(size_c,0,sizeof(int)*k);

	/*center=(double**)malloc(sizeof(double*)*k);
	memset(center,0,sizeof(double*)*k);
	for (i=0;i<k;i++)
	{
		center[i]=(double*)malloc(sizeof(double)*num_dm);
		memset(center[i],0,sizeof(double)*num_dm);
	}*/

	mybound=(double**)malloc(sizeof(double*)*num_dm);
	memset(mybound,0,sizeof(double*)*num_dm);
	for (i=0;i<num_dm;i++) //there are 3 mybounds for 4 categories
	{
		mybound[i]=(double*)malloc(sizeof(double)*3);
		memset(mybound[i],0,sizeof(double)*3);
	}

	prof=(int **)malloc(sizeof(int*)*k);
	memset(prof,0,sizeof(int*)*k);
	for (i=0;i<k;i++)
	{
		prof[i]=(int *)malloc(sizeof(int)*num_dm);
		memset(prof[i],0,sizeof(int)*num_dm);
	}

	size_f=(int **)malloc(sizeof(int*)*k);
	memset(size_f,0,sizeof(int*)*k);
	for (i=0;i<k;i++)
	{
		  size_f[i]=(int *)malloc(sizeof(int)*num_rows);
		  memset(size_f[i],0,sizeof(int)*num_rows);
	}

	for (i=0;i<k;i++)
		for (j=0;j<num_rows;j++)
			size_f[i][j]=0;

	for (i=0;i<file_Len;i++)
		  for (j=0;j<num_rows;j++)
			  if (Matrix[i][constrained_d[j]-1]>(double)constrained_line[j]) //on dimension j the event i is positive
				  size_f[cluster_id[i]][j]++;

	for (i=0;i<file_Len;i++)
	{
		id=cluster_id[i];
		//for (j=0;j<num_dm;j++)
		//{
		//	center[id][j]=center[id][j]+Matrix[i][j];
		//}
		
		size_c[id]++;		
	}

	percent_f=(double **)malloc(sizeof(double*)*k);
	  memset(percent_f,0,sizeof(double*)*k);
	  for (i=0;i<k;i++)
	  {
		  percent_f[i]=(double*)malloc(sizeof(double)*num_rows);
		  memset(percent_f[i],0,sizeof(double)*num_rows);
	  }

	
	 for (i=0;i<k;i++)
		for (j=0;j<num_rows;j++)
		{
			if (size_c[i]!=0)
				percent_f[i][j]=(double)size_f[i][j]/(double)size_c[i];
			else
				percent_f[i][j]=1.0;		
		}

	fnegpos_id=fopen("Final_Profile_NegPos.txt","w");

	fprintf(fnegpos_id,"FLOCK_Population_ID_Profile\t%s\n",name_string);

	for (i=0;i<k;i++)
	  {
	  	  fprintf(fnegpos_id,"%d\t",i+1);

		  for (j=0;j<num_rows;j++) //this requires num_rows==num_dm, i.e., 
			  //the user needs to specify the same number of constraints in the same order as in the file header; i.e., the user should specify in the profile constraint file the columns in an order from 0 to num_dm-1 
		  {
			if (j==num_rows-1)
			{
				if (percent_f[i][j]>0.5)
					fprintf(fnegpos_id,"+\n");
				else
					fprintf(fnegpos_id,"-\n");
			}
			else
			{
				if (percent_f[i][j]>0.5)
					fprintf(fnegpos_id,"+\t");
				else
					fprintf(fnegpos_id,"-\t");
			}
		  }
	  }

	fclose(fnegpos_id);

	for (i=0;i<k;i++)
	{
		free(size_f[i]);
		free(percent_f[i]);	
	}

	free(size_f);
	free(percent_f);


	/*for (i=0;i<k;i++)
		for (j=0;j<num_dm;j++)
		{			
			if (size_c[i]!=0)
				center[i][j]=(center[i][j]/(double)(size_c[i]));
			else
				center[i][j]=0;	
		}*/

	for (j=0;j<num_dm;j++)
	{
		interval=((big[j]-small[j])/4.0);
		//printf("interval[%d] is %f\n",j,interval);
		for (i=0;i<3;i++)
			mybound[j][i]=small[j]+((double)(i+1)*interval);
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
				if (Matrix[i][j]<mybound[j][0])// && ((Matrix[i][j]-small[j])>0)) //the smallest values excluded
					size_mybound_0[cluster_id[i]][j]++;
				else
				{
					if (Matrix[i][j]<mybound[j][1])
						size_mybound_1[cluster_id[i]][j]++;
					else
					{
						if (Matrix[i][j]<mybound[j][2])
							size_mybound_2[cluster_id[i]][j]++;
						else
							//if (Matrix[i][j]!=big[j]) //the biggest values excluded
								size_mybound_3[cluster_id[i]][j]++;
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
		fprintf(fpcnt_id,"%d\t%.4f\n",t+1,(double)size_c[t]*100.0/(double)file_Len);	//t changed to t+1 to start from 1 instead of 0: April 16, 2009									
	}
	fclose(fpcnt_id);

	//free(big);
	//free(small);
	free(size_c);

	for (i=0;i<k;i++)
	{
		//free(center[i]);
		free(prof[i]);
		free(size_mybound_0[i]);
		free(size_mybound_1[i]);
		free(size_mybound_2[i]);
		free(size_mybound_3[i]);
	}
	//free(center);
	free(prof);
	free(size_mybound_0);
	free(size_mybound_1);
	free(size_mybound_2);
	free(size_mybound_3);

	for (i=0;i<num_dm;i++)
		free(mybound[i]);
	free(mybound);
	
}



/******************************************************** Main Function **************************************************/
//for normalized data, there are five variables:
//cluster_ID
//population_center
//grid_clusterID
//grid_ID
//grid_Center

//the same five variables exist for the original data
//however, the three IDs (cluster_ID, grid_ID, grid_clusterID) don't change in either normalized or original data
//also, data + cluster_ID -> population_center
//data + grid_ID -> grid_Center

/* what is the final output */
//the final things we want are grid_Center in the original data and grid_clusterID //or population_center in the original data
//Sparse grids will not be considered when computing the two centroids (centroids of grids and centroids of clusters)

/*  what information should select_bin output */
//the size of all IDs are unknown to function main because we only consider the events in dense grids, and also the number of dense grids
//is unknown, therefore I must use a prescreening to output 
//how many bins I should use
//the number of dense grids
//total number of events in the dense grids

/* basic procedure of main function */ 
//1. read raw file and normalize the raw file
//2. select_bin
//3. allocate memory for eventID_To_denseventID, grid_clusterID, grid_ID, cluster_ID. 
//4. call select_dense and merge_grids with grid_clusterID, grid_ID, cluster_ID.
//5. release normalized data; allocate memory for grid_Center and population_center
//6. output grid_Center and population_center using ID2Center together with grid_clusterID //from IDs to centers

int main (int argc, char **argv)
{
  //inputs
  FILE *f_src; //source file pointer
 
  FILE *f_out; //coordinates
  FILE *f_cid; //population-ID of events
  FILE *f_ctr; //centroids of populations
  FILE *f_results; //coordinates file event and population column
  FILE *f_mfi; //added April 16, 2009 for mean fluorescence intensity
  FILE *f_parameters; //number of bins and density calculated by
                      //the algorithm. Used to update the database
  FILE *f_properties; //Properties file used by Image generation software
  FILE *f_spec; //user spec file for neg/pos separating coordinate on each dimenson

  FILE *f_profile_neg_pos; //Using the negative and positive to describe the cluster on each dimension
 
  FILE *f_majority_before; //proportion of a cluster inside the separating line on each dimension before dividing and merging
  //FILE *f_majority_after; //proportion of a cluster inside the separating line on each dimension in the final output

  FILE *f_range_before; //value range outside the separating line on each dimension before dividing and merging
  //FILE *f_range_after; //value range outside the separating line on each dimension in the final output

  
  char para_name_string[LINE_LEN];

  int time_ID=-1;
  int num_bin=0; //the bins I will use on each dimension
	
  int file_Len=0; //number of events		
  int num_dm=0;
  int num_clust=0;
  int num_dense_events=0;
  int num_dense_grids=0;
  int num_nonempty=0;
  int num_population=0;
  //int temp=0;

  //below are read from configuration file
  int i=0;
  int j=0;
  int max_num_pop=0;

  int den_t_event=0;
  int added_num_population=0;
  int added_num_pop=0;
  int ttt=0;

  int num_rows=0;
  int Constrained=0; //Constrained == 0 means it is not a constrained mode (i.e., original FLOCK without requring user input on negative positve separating line); otherwise it is

  int *grid_clusterID; //(dense)gridID to clusterID
  int *grid_ID; //(dense)eventID to gridID
  int *cluster_ID; //(dense)eventID to clusterID
  int *eventID_To_denseventID; //eventID_To_denseventID[i]=k means event i is in a dense grid and its ID within dense events is k
  int *all_grid_ID; //gridID for all events
  int *densegridID_To_gridevent;
  int *sorted_seq;
  int *dense_grid_reverse;
  int *population_ID; //denseeventID to populationID
  int *cluster_populationID; //clusterID to populationID
  int *grid_populationID; //gridID to populationID
  int *all_population_ID; //populationID of event

  int *constrained_d;
  int *constrained_line;

  int *size_c; //number of events in each cluster;
  int **size_f; //number of events outside each cluster on each dimension
  int **profile_neg_pos; //negative and positive profile based on user specification
  int **added_profile_neg_pos; //profile after extension

  //int *filtered_out;
  //int *filtered_p;

  int **position;
  int **dense_grid_seq;

  double constrain_proportion_t=0.4; //default value
  double constrain_range_t=0.25; //default value
  //double constrain_dist_t=0.0;

  double *interval;

  double *largest_value;  //for each dimension, the largest value
  double *smallest_value; //for each dimension, the smallest value

  double **percent_f; //percent_f records each population each dimension outside percentage
	
  double **population_center; //population centroids in the raw/original data
  double **cluster_center; //cluster centroids in the raw/original data

  double **input_data;
  double **normalized_data;

  double **range_f_max; //upper bound of the range on each dimension of each population
  double **range_f_min; //lower bound of the range on each dimension of each population

  double **range_f; //range of the outside region on each dimension of each population
  double **added_pop_center;
		
  int min = MAX_VALUE;
  int max = -MAX_VALUE;

  printf("Starting time:\t\t\t\t");
  fflush(stdout);
  system("/bin/date");
  /////////////////////////////////////////////////////////////

  if ((argc!=2) && (argc!=4) && (argc!=5) && (argc!=8))
  {
      //modified on Jul 23, 2010
	  //printf("usage:\n");
	  //printf("basic mode: flock data_file\n");
	  //printf("advanced mode: flock data_file num_bin density_index\n");       
      //exit(0);

	  fprintf(stderr,"Incorrect number of input parameters!\n"); //modified on July 23, 2010
	  fprintf(stderr,"basic mode: FLOCK_constrained data_file\n"); //modified on July 23, 2010
	  fprintf(stderr,"advanced mode 1: FLOCK_constrained data_file num_bin density_index\n"); //modified on July 23, 2010
	  fprintf(stderr,"advanced mode 2: FLOCK_constrained data_file num_bin density_index max_num_pop\n"); //added on Dec 16, 2010 for GenePattern
	  fprintf(stderr,"advanced mode 3: FLOCK_constrained data_file num_bins density_threshold max_num_clusters profile_spec_file pop_proportion_threshold pop_range_distance_threshold\n");
	  fprintf(stderr,"when using constraints, the number of constraints should be the same as number of dimensions, specified in the same order in the profile_spec_file\n");
	  abort();
  }	

  f_src=fopen(argv[1],"r");

  if (argc==2)
  {
	 max_num_pop=DEFAULT_NUM_POP; //default value, maximum 30 clusters
  }
  
  if (argc==4)
  {
	 num_bin=atoi(argv[2]);
	 printf("num_bin is %d\n",num_bin);

	 den_t_event=atoi(argv[3]);
	 printf("density_index is %d\n",den_t_event);

	 max_num_pop=DEFAULT_NUM_POP;

	 if (((num_bin<MIN_GRID) && (num_bin!=0)) || (num_bin>MAX_GRID))
	 { 
		fprintf(stderr,"Incorrect input range of number of bins, which should be larger than 5 and smaller than 31\n");
		abort();
	 }

	 if (((den_t_event<DEFAULT_DENSITY_T) && (den_t_event!=0)) || (den_t_event>MAX_D_T))
	 {
		fprintf(stderr,"Incorrect input range of density threshold, which should be larger than 2 and smaller than 51\n");
		abort();
	 }
  }

  if (argc==5)
  {
	  num_bin=atoi(argv[2]);
	  printf("num_bin is %d\n",num_bin);

	  den_t_event=atoi(argv[3]);
	  printf("density_index is %d\n",den_t_event);

	  max_num_pop=atoi(argv[4]);
	  printf("max number of clusters is %d\n",max_num_pop);

	  if (((num_bin<MIN_GRID) && (num_bin!=0)) || (num_bin>MAX_GRID))
	  { 
		fprintf(stderr,"Incorrect input range of number of bins, which should be larger than 5 and smaller than 31\n");
		abort();
	  }

	  if (((den_t_event<DEFAULT_DENSITY_T) && (den_t_event!=0)) || (den_t_event>MAX_D_T))
	  {
		fprintf(stderr,"Incorrect input range of density threshold, which should be larger than 2 and smaller than 51\n");
		abort();
	  }

	  if ((max_num_pop<LOWER_NUM_CLUSTER) || (max_num_pop>UPPER_NUM_CLUSTER))
	  {
		fprintf(stderr,"Incorrect input range of maximum number of populations, which should be larger than 1 and smaller than 1000\n");
		abort();
	  }
  }

  if (argc==8)
  {
	  Constrained=1;

	  num_bin=atoi(argv[2]);
	  printf("num_bin is %d\n",num_bin);

	  den_t_event=atoi(argv[3]);
	  printf("density_index is %d\n",den_t_event);

	  max_num_pop=atoi(argv[4]);
	  printf("max number of clusters is %d\n",max_num_pop);

	  f_spec=fopen(argv[5],"r");

	  constrain_proportion_t=atof(argv[6]);
	  constrain_range_t=atof(argv[7]);
	  //constrain_dist_t=atof(argv[8]);

	  get_spec_num(f_spec,&num_rows);
	  printf("Number of Dimensions to be profiled is %d\n",num_rows);

	  rewind(f_spec);

	  constrained_d=(int *)malloc(sizeof(int)*num_rows);
	  memset(constrained_d,0,sizeof(int)*num_rows);

	  constrained_line=(int *)malloc(sizeof(int)*num_rows);
	  memset(constrained_line,0,sizeof(int)*num_rows);

	  get_spec_info(f_spec,constrained_d,constrained_line,num_rows);

	  printf("The last row of the spec file is %d\t%d\n",constrained_d[num_rows-1],constrained_line[num_rows-1]);

	  fclose(f_spec);

	  if (((num_bin<MIN_GRID) && (num_bin!=0)) || (num_bin>MAX_GRID))
	  { 
		fprintf(stderr,"Incorrect input range of number of bins, which should be larger than 5 and smaller than 31\n");
		abort();
	  }

	  if (((den_t_event<DEFAULT_DENSITY_T) && (den_t_event!=0)) || (den_t_event>MAX_D_T))
	  {
		fprintf(stderr,"Incorrect input range of density threshold, which should be larger than 2 and smaller than 51\n");
		abort();
	  }

	  if ((max_num_pop<LOWER_NUM_CLUSTER) || (max_num_pop>UPPER_NUM_CLUSTER))
	  {
		fprintf(stderr,"Incorrect input range of maximum number of populations, which should be larger than 1 and smaller than 1000\n");
		abort();
	  } 
  }
    

  getfileinfo(f_src, &file_Len, &num_dm, para_name_string, &time_ID); //get the filelength, number of dimensions, and num/name of parameters

  /************************************************* Read the data *****************************************************/

  //if ((num_rows!=num_dm) && (Constrained==1))
  //{
	//  printf("every dimension needs to be specified in the constrain file, if not available, specify -1\n");
	//  abort();
  //}
	
  rewind(f_src); //reset the file pointer	

  input_data = (double **)malloc(sizeof(double*)*file_Len);
  memset(input_data,0,sizeof(double*)*file_Len);
  for (i=0;i<file_Len;i++)
  {
     input_data[i]=(double *)malloc(sizeof(double)*num_dm);
     memset(input_data[i],0,sizeof(double)*num_dm);
  }
	
  readsource(f_src, file_Len, num_dm, input_data, time_ID); //read the data;
	
  fclose(f_src);

  normalized_data=(double **)malloc(sizeof(double*)*file_Len);
  memset(normalized_data,0,sizeof(double*)*file_Len);
  for (i=0;i<file_Len;i++)
    {
      normalized_data[i]=(double *)malloc(sizeof(double)*num_dm);
      memset(normalized_data[i],0,sizeof(double)*num_dm);
    }
	
  largest_value=(double*)malloc(sizeof(double)*num_dm);
  memset(largest_value,0,sizeof(double)*num_dm);

  smallest_value=(double*)malloc(sizeof(double)*num_dm);
  memset(smallest_value,0,sizeof(double)*num_dm); 

  tran(input_data, file_Len, num_dm, NORM_METHOD, normalized_data,largest_value,smallest_value);
	

  position=(int **)malloc(sizeof(int*)*file_Len);
  memset(position,0,sizeof(int*)*file_Len);
  for (i=0;i<file_Len;i++)
    {
      position[i]=(int*)malloc(sizeof(int)*num_dm);
      memset(position[i],0,sizeof(int)*num_dm);
    }

  all_grid_ID=(int *)malloc(sizeof(int)*file_Len);
  memset(all_grid_ID,0,sizeof(int)*file_Len);

  sorted_seq=(int*)malloc(sizeof(int)*file_Len);
  memset(sorted_seq,0,sizeof(int)*file_Len);
	
  interval=(double*)malloc(sizeof(double)*num_dm);
  memset(interval,0,sizeof(double)*num_dm);

  /************************************************* select_bin *************************************************/
	
  if (num_bin>=1)  //num_bin has been selected by user
  	select_bin(normalized_data, file_Len, num_dm, MIN_GRID, MAX_GRID, position, sorted_seq, all_grid_ID, &num_nonempty, interval,num_bin);
  else  //num_bin has not been selected by user
  {
	num_bin=select_bin(normalized_data, file_Len, num_dm, MIN_GRID, MAX_GRID, position, sorted_seq, all_grid_ID, &num_nonempty, interval,num_bin);
	printf("selected bin number is %d\n",num_bin);    
  }
  printf("number of non-empty grids is %d\n",num_nonempty);
		
 

  /* Although we return sorted_seq from select_bin(), we don't use it for anything, except possibly diagnostics */
  free(sorted_seq);
	
	
  dense_grid_reverse=(int*)malloc(sizeof(int)*num_nonempty);
  memset(dense_grid_reverse,0,sizeof(int)*num_nonempty);	

  /************************************************* select_dense **********************************************/

  if (den_t_event>=1) //den_t_event must be larger or equal to 2 if the user wants to set it
	select_dense(file_Len, all_grid_ID, num_nonempty, &num_dense_grids, &num_dense_events, dense_grid_reverse, den_t_event, max_num_pop);
  else
  {
	den_t_event=select_dense(file_Len, all_grid_ID, num_nonempty, &num_dense_grids, &num_dense_events, dense_grid_reverse, den_t_event, max_num_pop);
	printf("automated selected density threshold is %d\n",den_t_event);
  }		

  printf("Number of dense grids is %d\n",num_dense_grids);

  densegridID_To_gridevent = (int *)malloc(sizeof(int)*num_dense_grids);
  memset(densegridID_To_gridevent,0,sizeof(int)*num_dense_grids);

  for (i=0;i<num_dense_grids;i++)
    densegridID_To_gridevent[i]=-1; //initialize all densegridID_To_gridevent values to -1
	

  eventID_To_denseventID=(int *)malloc(sizeof(int)*file_Len);
  memset(eventID_To_denseventID,0,sizeof(int)*file_Len);     //eventID_To_denseventID[i]=k means event i is in a dense grid and its ID within dense events is k


  grid_To_event(file_Len, dense_grid_reverse, all_grid_ID, eventID_To_denseventID, densegridID_To_gridevent);

	
  dense_grid_seq=(int **)malloc(sizeof(int*)*num_dense_grids);
  memset(dense_grid_seq,0,sizeof(int*)*num_dense_grids);
  for (i=0;i<num_dense_grids;i++)
    {
      dense_grid_seq[i]=(int *)malloc(sizeof(int)*num_dm);
      memset(dense_grid_seq[i],0,sizeof(int)*num_dm);
    }


  /* Look up the binned data values for each dense grid */
  generate_grid_seq(num_dm, num_dense_grids, densegridID_To_gridevent, position, dense_grid_seq); 	
	
	
  /************************************************* allocate memory *********************************************/
	
  grid_clusterID=(int *)malloc(sizeof(int)*num_dense_grids);
  memset(grid_clusterID,0,sizeof(int)*num_dense_grids);

  grid_ID=(int *)malloc(sizeof(int)*num_dense_events);
  memset(grid_ID,0,sizeof(int)*num_dense_events);

  cluster_ID=(int *)malloc(sizeof(int)*num_dense_events);
  memset(cluster_ID,0,sizeof(int)*num_dense_events);


  /*********************************************** merge_grids ***********************************************/
  //int merge_grids(int file_Len, int num_dm, int num_bin, int **position, int num_dense_grids, int *dense_grid_rank, int *dense_grid_reverse,
  //			 int **dense_grid_seq, int *eventID_To_denseventID, int *densegridID_To_gridevent, int *all_grid_ID,
  //			 int *cluster_ID, int *grid_ID, int *grid_clusterID)
	
  num_clust = merge_grids(normalized_data, interval, file_Len, num_dm, num_bin, position, num_dense_grids, dense_grid_reverse, dense_grid_seq, eventID_To_denseventID, densegridID_To_gridevent, all_grid_ID, cluster_ID, grid_ID, grid_clusterID);
	
  printf("computed number of groups is %d\n",num_clust);	

	
  /************************************** release unnecessary memory and allocate memory and compute centers **********************************/
	
	
  for (i=0;i<file_Len;i++)
    free(position[i]);
  free(position);

  for (i=0;i<num_dense_grids;i++)
    free(dense_grid_seq[i]);
  free(dense_grid_seq);

  free(dense_grid_reverse);
	
  free(densegridID_To_gridevent);
  free(all_grid_ID);
	
  // cluster_center ////////////////////////////////////////////////////////////////////////////////////////////////////////
	
  cluster_center=(double **)malloc(sizeof(double*)*num_clust);
  memset(cluster_center,0,sizeof(double*)*num_clust);
  for (i=0;i<num_clust;i++)
  {
     cluster_center[i]=(double*)malloc(sizeof(double)*num_dm);
     memset(cluster_center[i],0,sizeof(double)*num_dm);
  }
	
  ID2Center(normalized_data,file_Len,num_dm,eventID_To_denseventID,num_clust,cluster_ID,cluster_center); //produce the centers with normalized data
	
  //printf("pass the first ID2center\n");  //commented on July 23, 2010

  /*** population_ID and grid_populationID **/
	
  cluster_populationID=(int*)malloc(sizeof(int)*num_clust);
  memset(cluster_populationID,0,sizeof(int)*num_clust);

  grid_populationID=(int*)malloc(sizeof(int)*num_dense_grids);
  memset(grid_populationID,0,sizeof(int)*num_dense_grids);

  population_ID=(int*)malloc(sizeof(int)*num_dense_events);
  memset(population_ID,0,sizeof(int)*num_dense_events);

  num_population = merge_clusters(num_clust, num_dm, interval, cluster_center, cluster_populationID,max_num_pop);

  

  for (i=0;i<num_clust;i++)
    free(cluster_center[i]);
  free(cluster_center);

  free(interval);
	
  for (i=0;i<num_dense_grids;i++)
    {
      grid_populationID[i]=cluster_populationID[grid_clusterID[i]];
    }

  for (i=0;i<num_dense_events;i++)
    {
      population_ID[i]=cluster_populationID[cluster_ID[i]];
    }

  printf("computed number of populations is %d\n",num_population);

  

  // population_center /////////////////////////////////////////////////////////////////////////////////////////////////////

 
  population_center=(double **)malloc(sizeof(double*)*num_population);
  memset(population_center,0,sizeof(double*)*num_population);
  for (i=0;i<num_population;i++)
    {
      population_center[i]=(double*)malloc(sizeof(double)*num_dm);
      memset(population_center[i],0,sizeof(double)*num_dm);
    }

  
	
  ID2Center(normalized_data,file_Len,num_dm,eventID_To_denseventID,num_population,population_ID,population_center); //produce population centers with normalized data
	

  // show ////////////////////////////////////////////////////////////////////////////////
  all_population_ID=(int*)malloc(sizeof(int)*file_Len);
  memset(all_population_ID,0,sizeof(int)*file_Len);

  kmeans(normalized_data, num_population, KMEANS_TERM, file_Len, num_dm, all_population_ID, population_center);

  /*
    The constrained algorithm is:
	Step 1: Run FLOCK to identify all clusters in full dimensional space
	Step 2: For each FLOCK population 
	For each dimension
	Calcluate its proportion outside each dimension separatining_coordinate
	Calculate the value range of the outside part on each dimension
	Check whether the proportion outside the line is smaller than the proportion_threshold AND the range of the dimension is smaller than the range_distance_threshold
	if yes, check the next dimension
	otherwise, generate a new centroid (the centroid of the outside events) on this dimension; check the next dimension
	Step 3: run K-means again; Output the final result
	
	A new profile file with the negative/positive for each population on each user specified dimension will be generated.
	proportion and value range out of the user boundary before and after generation of new centroids will be output
 */

  if (Constrained==1) //if it is the constrained mode, the population ID needs to be changed based on the above algorithm
  {
	  
	  percent_f=(double **)malloc(sizeof(double*)*num_population);
	  memset(percent_f,0,sizeof(double*)*num_population);
	  for (i=0;i<num_population;i++)
	  {
		  percent_f[i]=(double*)malloc(sizeof(double)*num_rows);
		  memset(percent_f[i],0,sizeof(double)*num_rows);
	  }

	 

	  size_c=(int *)malloc(sizeof(int)*num_population);
	  memset(size_c,0,sizeof(int)*num_population);

	  size_f=(int **)malloc(sizeof(int*)*num_population);
	  memset(size_f,0,sizeof(int*)*num_population);
	  for (i=0;i<num_population;i++)
	  {
		  size_f[i]=(int *)malloc(sizeof(int)*num_rows);
		  memset(size_f[i],0,sizeof(int)*num_rows);
	  }

	  range_f_max=(double **)malloc(sizeof(double*)*num_population);
	  memset(range_f_max,0,sizeof(double*)*num_population);
	  for (i=0;i<num_population;i++)
	  {
		  range_f_max[i]=(double*)malloc(sizeof(double)*num_rows);
		  memset(range_f_max[i],0,sizeof(double)*num_rows);
	  }

	
	  range_f_min=(double **)malloc(sizeof(double*)*num_population);
	  memset(range_f_min,0,sizeof(double*)*num_population);
	  for (i=0;i<num_population;i++)
	  {
		  range_f_min[i]=(double*)malloc(sizeof(double)*num_rows);
		  memset(range_f_min[i],0,sizeof(double)*num_rows);
	  }

	  range_f=(double **)malloc(sizeof(double*)*num_population);
	  memset(range_f,0,sizeof(double*)*num_population);
	  for (i=0;i<num_population;i++)
	  {
		  range_f[i]=(double*)malloc(sizeof(double)*num_rows);
		  memset(range_f[i],0,sizeof(double)*num_rows);
	  }

	

	  profile_neg_pos=(int **)malloc(sizeof(int*)*num_population); //+1 means positive; -1 means negative
	  memset(profile_neg_pos,0,sizeof(int*)*num_population);
	  for (i=0;i<num_population;i++)
	  {
		  profile_neg_pos[i]=(int *)malloc(sizeof(int)*num_rows);
		  memset(profile_neg_pos[i],0,sizeof(int)*num_rows);
	  }

	  for (i=0;i<num_population;i++)
	  {
		  size_c[i]=0;
		  
		  for (j=0;j<num_rows;j++)
		  {
				size_f[i][j]=0;
				percent_f[i][j]=0.0;
				range_f_max[i][j]=-MAX_VALUE;
				range_f_min[i][j]=MAX_VALUE;
		  }		  
	  }

	 

	  for (i=0;i<file_Len;i++)
	  {
		  size_c[all_population_ID[i]]++;
	  
		  for (j=0;j<num_rows;j++)
		  {
			if (input_data[i][constrained_d[j]-1]>range_f_max[all_population_ID[i]][j])
				range_f_max[all_population_ID[i]][j]=input_data[i][constrained_d[j]-1];

			if (input_data[i][constrained_d[j]-1]<range_f_min[all_population_ID[i]][j])
				range_f_min[all_population_ID[i]][j]=input_data[i][constrained_d[j]-1];
		  }
	  }

	
	  
	  for (i=0;i<file_Len;i++)
		  for (j=0;j<num_rows;j++)
			  if (input_data[i][constrained_d[j]-1]>(double)constrained_line[j]) //on dimension j the event i is positive
				  size_f[all_population_ID[i]][j]++;
		
	 
		  
	  for (i=0;i<num_population;i++)
		for (j=0;j<num_rows;j++)
		{
			if (size_c[i]!=0)
				percent_f[i][j]=(double)size_f[i][j]/(double)size_c[i];
			else
				percent_f[i][j]=1;

			if (percent_f[i][j]>0.5)
			{
				percent_f[i][j]=1-percent_f[i][j];
				profile_neg_pos[i][j]=1; //positive on dimension constrained_d[j]-1
				
				//range_f[i][j]=((double)constrained_line[j]-range_f_min[i][j])/(largest_value[constrained_d[j]-1]-smallest_value[constrained_d[j]-1]);  //the smaller (even negative), the better
				range_f[i][j]=((double)constrained_line[j]-range_f_min[i][j])/(double)UPPER_VALUE;
			}
			else //negative on dimension constrained_d[j]-1
			{
				profile_neg_pos[i][j]=-1;

				//range_f[i][j]=(range_f_max[i][j]-(double)constrained_line[j])/(largest_value[constrained_d[j]-1]-smallest_value[constrained_d[j]-1]); //the smaller (even negative), the better
				range_f[i][j]=(range_f_max[i][j]-(double)constrained_line[j])/(double)UPPER_VALUE;
			}
		}


	 //now you have range_f and percent_f to check
	 added_num_pop=0;

	  for (i=0;i<num_population;i++)
	  {
		  for (j=0;j<num_rows;j++) //you only need to check the dimensions that have constraints
		  {
			  if ((range_f[i][j]<=constrain_range_t) && (percent_f[i][j]<=constrain_proportion_t)) //this dimension has no problem
				  continue;
			  else //create a new centroid on this dimension
			  	  added_num_pop++;			  
		  }
	  }

	  printf("there will be %d new populations added\n",added_num_pop);

	  added_num_population=num_population+added_num_pop;

	  added_pop_center=(double **)malloc(sizeof(double*)*added_num_population);
	  memset(added_pop_center,0,sizeof(double*)*added_num_population);
	  for (i=0;i<added_num_population;i++)
	  {
		  added_pop_center[i]=(double*)malloc(sizeof(double)*num_dm);
		  memset(added_pop_center[i],0,sizeof(double)*num_dm);
	  }

	

	  for (i=0;i<num_population;i++)
		  for (j=0;j<num_dm;j++)
			  added_pop_center[i][j]=population_center[i][j];

	

	  added_profile_neg_pos=(int **)malloc(sizeof(int*)*added_num_population); //+1 means positive; -1 means negative
	  memset(added_profile_neg_pos,0,sizeof(int*)*added_num_population);
	  for (i=0;i<added_num_population;i++)
	  {
		  added_profile_neg_pos[i]=(int *)malloc(sizeof(int)*num_rows);
		  memset(added_profile_neg_pos[i],0,sizeof(int)*num_rows);
	  }

	

	  for (i=0;i<num_population;i++)  //copy the profiles of the current cell populations
		  for (j=0;j<num_rows;j++)
			  added_profile_neg_pos[i][j]=profile_neg_pos[i][j];

	  added_num_pop=0;

	 

	  for (i=0;i<num_population;i++)
	  {
		  for (j=0;j<num_rows;j++) //note we require num_rows == num_dm, i.e., the user needs to specify constraints for all dimensions in the order from 0 to num_dm-1; constrained_d[j] == j+1 is always true
		  {
			  if ((range_f[i][j]<constrain_range_t) && (percent_f[i][j]<constrain_proportion_t)) //there is no need to separate cluster i on dimension j
				  continue;
			  else //create a new centroid on this dimension
			  {
				  if (range_f[i][j]<=0) //when range_f[i][j]<=0, the whole population is on one side of the boundary, the percent_f[i][j] should always be zero;
				  {
					  printf("Something is wrong. The range outside the specified boundary should be larger than zero\n");
					  exit(0);
				  }

				  for (ttt=0;ttt<num_dm;ttt++)
				  	  added_pop_center[num_population+added_num_pop][ttt]=added_pop_center[i][ttt];

				  for (ttt=0;ttt<num_rows;ttt++)
					  added_profile_neg_pos[num_population+added_num_pop][ttt]=added_profile_neg_pos[i][ttt];

				  if (profile_neg_pos[i][j]==1) //a positive population, then the added one should be a negative one
				  {
						//added_pop_center[num_population+added_num_pop][constrained_d[j]-1]=((range_f_min[i][j]-smallest_value[constrained_d[j]-1])/(largest_value[constrained_d[j]-1]-smallest_value[constrained_d[j]-1]))+(range_f[i][j]/2.0);
					  added_pop_center[num_population+added_num_pop][constrained_d[j]-1]=(range_f_min[i][j]/(double)UPPER_VALUE)+(range_f[i][j]/2.0);
						added_profile_neg_pos[num_population+added_num_pop][j]=-1;
				  }
				  else //a negative population, then the added one should be a positive one
				  {
						//added_pop_center[num_population+added_num_pop][constrained_d[j]-1]=((range_f_max[i][j]-smallest_value[constrained_d[j]-1])/(largest_value[constrained_d[j]-1]-smallest_value[constrained_d[j]-1]))-(range_f[i][j]/2.0);
					  added_pop_center[num_population+added_num_pop][constrained_d[j]-1]=(range_f_max[i][j]/(double)UPPER_VALUE)-(range_f[i][j]/2.0);
						added_profile_neg_pos[num_population+added_num_pop][j]=1;
				  }

				  added_num_pop++;
			  }
		  }
	  }

	

	  f_majority_before=fopen("proportion_out_before.txt","w");
	  f_range_before=fopen("range_out_before.txt","w");
	  f_profile_neg_pos=fopen("Profile_NegPos.txt","w");

	  fprintf(f_majority_before,"FLOCK_Population_ID_ProportionOut\t%s\n",para_name_string);
	  fprintf(f_range_before,"FLOCK_Population_ID_RangeOut\t%s\n",para_name_string);
	  fprintf(f_profile_neg_pos,"FLOCK_Population_ID_Profile\t%s\n",para_name_string);
	  //fprintf(f_profile_neg_pos,"%s\n",para_name_string);

	  //for (i=0;i<num_rows;i++)
	  //{
		//  fprintf(f_majority_before,"Proportion_outof_Dimension%d\t",constrained_d[i]);
		//  fprintf(f_range_before,"Range_outof_Dimension%d\t",constrained_d[i]);
		  
		  
		  //if (i==num_rows-1)
		  //	 fprintf(f_profile_neg_pos,"Profile_User_Dimension%d\n",constrained_d[i]);
		  //else
		  //	 fprintf(f_profile_neg_pos,"Profile_User_Dimension%d\t",constrained_d[i]);
	 // }
	
	  fprintf(f_majority_before,"Proportion_in_WholeFile\n");
	  fprintf(f_range_before,"Proportion_in_WholeFile\n");
	  //fprintf(f_profile_neg_pos,"Proportion_in_WholeFile\n");
	

	  for (i=0;i<num_population;i++)
	  {
	  	  fprintf(f_majority_before,"%d\t",i+1);
		  fprintf(f_range_before,"%d\t",i+1);
		 
		  for (j=0;j<num_rows;j++)
		  {
			  fprintf(f_majority_before,"%.4f\t",percent_f[i][j]);
			  fprintf(f_range_before,"%.4f\t",range_f[i][j]);
		  }
			  
		  fprintf(f_majority_before,"%.4f\n",(double)size_c[i]/(double)file_Len);
		  fprintf(f_range_before,"%.4f\n",(double)size_c[i]/(double)file_Len);
	  }



	  for (i=0;i<added_num_population;i++)
	  {
	  	  fprintf(f_profile_neg_pos,"%d\t",i+1);

		  for (j=0;j<num_rows;j++) //this actually requires num_rows==num_dm, i.e., the user needs to specify the same number of constraints in the same order as in the file header
		  {
			if (j==num_rows-1)
			{
				if (added_profile_neg_pos[i][j]==1)
					fprintf(f_profile_neg_pos,"+\n");
				else
					fprintf(f_profile_neg_pos,"-\n");
			}
			else
			{
				if (added_profile_neg_pos[i][j]==1)
					fprintf(f_profile_neg_pos,"+\t");
				else
					fprintf(f_profile_neg_pos,"-\t");
			}
		  }
			  
		  //fprintf(f_profile_neg_pos,"%.4f\n",(double)size_c[i]/(double)file_Len);
		  
	  }

	
	  
	  fclose(f_majority_before);
	  fclose(f_range_before);
	  fclose(f_profile_neg_pos);

	  free(size_c);
	  	
	
	  for (i=0;i<num_population;i++)
	  {
			free(population_center[i]);
		    free(size_f[i]);
			free(percent_f[i]);
			free(profile_neg_pos[i]);
			free(range_f_max[i]);
			free(range_f_min[i]);
			free(range_f[i]);
	  }
	  free(size_f);
	  free(percent_f);
	  free(profile_neg_pos);
	  free(range_f_max);
	  free(range_f_min);
	  free(range_f);
	  free(population_center);

	 

	  num_population=added_num_population;

	  population_center=(double **)malloc(sizeof(double*)*num_population);
	  memset(population_center,0,sizeof(double*)*num_population);
      for (i=0;i<num_population;i++)
      {
        population_center[i]=(double*)malloc(sizeof(double)*num_dm);
        memset(population_center[i],0,sizeof(double)*num_dm);
      }

	
	  for (i=0;i<num_population;i++)
		 for (j=0;j<num_dm;j++)
			 population_center[i][j]=added_pop_center[i][j];

	  for (i=0;i<num_population;i++)
	  {
		  free(added_pop_center[i]);
		  free(added_profile_neg_pos[i]);
	  }
	  free(added_pop_center);
	  free(added_profile_neg_pos);

	 printf("The final number of populations is %d\n",num_population);

	  kmeans(normalized_data, num_population, KMEANS_TERM, file_Len, num_dm, all_population_ID, population_center);

	  //num_population = merge_clusters(num_clust, num_dm, interval, cluster_center, cluster_populationID,max_num_pop);

	 /* cent_dist=(double **)malloc(sizeof(double*)*added_num_population); //+1 means positive; -1 means negative
	  memset(cent_dist,0,sizeof(double*)*added_num_population);
	  for (i=0;i<added_num_population;i++)
	  {
		  cent_dist[i]=(double *)malloc(sizeof(double)*added_num_population);
		  memset(cent_dist[i],0,sizeof(double)*added_num_population);
	  }

	  for (i=0;i<(added_num_population-1);i++)
		  for (j=i+1;j<added_num_population;j++)
		  {
			  cent_dist[i][j]=0.0;
			  for (ttt=0;ttt<num_dm;ttt++)
			  {
				  ddd=added_pop_center[i][ttt]-added_pop_center[j][ttt];

				  if (ddd<0)
					  ddd=-ddd;

				  if (ddd>cent_dist[i][j]) //the biggest dimension difference is used as the distance between two clusters
					  cent_dist[i][j]=ddd;
			  }
		  }


	  for (i=0;i<(added_num_population-1);i++)
		  for (j=i+1;j<added_num_population;j++)
		  {
			  to_be_merge[j]=0;
			  if (cent_dist[i][j]<constrain_dist_t) //the clusters i and j need to be merged
				to_be_merge[j]=1;
		  }
		  */
	  	  
  }
  

  show(input_data, all_population_ID, file_Len, num_population, num_dm, para_name_string, largest_value, smallest_value, constrained_d, constrained_line, num_rows);

  ID2Center_all(input_data,file_Len,num_dm,num_population,all_population_ID,population_center);
  

  f_cid=fopen("population_id.txt","w");
  f_ctr=fopen("population_center.txt","w");
  f_out=fopen("coordinates.txt","w");
  f_results=fopen("flock_results.txt","w");

/*
  f_parameters=fopen("parameters.txt","w");
  fprintf(f_parameters,"Number_of_Bins\t%d\n",num_bin);
  fprintf(f_parameters,"Density\t%f\n",aver_index);
  fclose(f_parameters);
*/

  for (i=0;i<file_Len;i++)
	fprintf(f_cid,"%d\n",all_population_ID[i]+1); //all_population_ID[i] changed to all_population_ID[i]+1 to start from 1 instead of 0: April 16, 2009

  /*
   * New to check for min/max to add to parameters.txt
   *
  */
  
  fprintf(f_out,"%s\n",para_name_string);
  fprintf(f_results,"%s\tEvent\tPopulation\n",para_name_string);
  for (i=0;i<file_Len;i++)
  {
	for (j=0;j<num_dm;j++)
	{
		if (input_data[i][j] < min) {
			min = (int)input_data[i][j];
		}
		if (input_data[i][j] > max) {
			max = (int)input_data[i][j];
		}
		if (j==num_dm-1)
		{
			fprintf(f_out,"%d\n",(int)input_data[i][j]);
			fprintf(f_results,"%d\t",(int)input_data[i][j]);
		}
		else
		{
			fprintf(f_out,"%d\t",(int)input_data[i][j]);
			fprintf(f_results,"%d\t",(int)input_data[i][j]);
		}
	}
	fprintf(f_results,"%d\t",i + 1);
	fprintf(f_results,"%d\n",all_population_ID[i]+1); //all_population_ID[i] changed to all_population_ID[i]+1 to start from 1 instead of 0: April 16, 2009
  }

  f_parameters=fopen("parameters.txt","w");
  fprintf(f_parameters,"Number_of_Bins\t%d\n",num_bin);
  fprintf(f_parameters,"Density\t%d\n",den_t_event);
  fprintf(f_parameters,"Min\t%d\n",min);
  fprintf(f_parameters,"Max\t%d\n",max);
  fclose(f_parameters);

  f_properties=fopen("fcs.properties","w");
  fprintf(f_properties,"Bins=%d\n",num_bin);
  fprintf(f_properties,"Density=%d\n",den_t_event);
  fprintf(f_properties,"Min=%d\n",min);
  fprintf(f_properties,"Max=%d\n",max);
  fprintf(f_properties,"Populations=%d\n",num_population);
  fprintf(f_properties,"Events=%d\n",file_Len);
  fprintf(f_properties,"Markers=%d\n",num_dm);
  fclose(f_properties);


  for (i=0;i<num_population;i++) {
	/* Add if we want to include population id in the output
	*/
	fprintf(f_ctr,"%d\t",i+1);  //i changed to i+1 to start from 1 instead of 0: April 16, 2009

	for (j=0;j<num_dm;j++) {
		if (j==num_dm-1)
			fprintf(f_ctr,"%.0f\n",population_center[i][j]);
		else
			fprintf(f_ctr,"%.0f\t",population_center[i][j]);
	}
  }

  	//added April 16, 2009
	f_mfi=fopen("MFI.txt","w");

	for (i=0;i<num_population;i++)
	{
		fprintf(f_mfi,"%d\t",i+1);

		for (j=0;j<num_dm;j++)
		{
			if (j==num_dm-1)
				fprintf(f_mfi,"%.0f\n",population_center[i][j]);
			else
				fprintf(f_mfi,"%.0f\t",population_center[i][j]);
		}
	}
	fclose(f_mfi);

	//ended April 16, 2009
			
  fclose(f_cid);
  fclose(f_ctr);
  fclose(f_out);
  fclose(f_results);


  for (i=0;i<num_population;i++)
  {
	free(population_center[i]);
  }
  free(population_center);
 
  free(largest_value);
  free(smallest_value);

  for (i=0;i<file_Len;i++)
    free(normalized_data[i]);
  free(normalized_data);	
	
  free(grid_populationID);

  free(cluster_populationID);
  free(grid_clusterID);
  free(cluster_ID);

  for (i=0;i<file_Len;i++)
    free(input_data[i]);
  free(input_data);

  free(grid_ID);
  free(population_ID);
  free(all_population_ID);
  free(eventID_To_denseventID);

  if (Constrained==1)
  {
	  free(constrained_d);
	  free(constrained_line);
  }
		
  ///////////////////////////////////////////////////////////
  printf("Ending time:\t\t\t\t");
  fflush(stdout);
  system("/bin/date");

  return 0;

}
