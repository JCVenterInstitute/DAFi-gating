/************************************************************************************************************************

	User-directed unsupervised identification of cell populations - for Galaxy platform

	Authors: Yu "Max" Qian, Ph.D., mqian@jcvi.org or qianyu_cs@yahoo.com, Ivan Chang, Ph.D., ichang@jcvi.org, and Bob Sinkovits, Ph.D., sinkovit@sdsc.edu

	Copyright: Author and J. Craig Venter Institute

	Usage: dafi_gating data_file gating_spec_file filter_spec_file first_pass_max_num_clusters second_pass_max_num_clusters num_threads seed

	data_file is a tab delimited file with compensated and transformed values.

	gating_spec_file: a 11-column tab delimited file, events inside will be kept
	1       1       4       20      70      5       55      0       0       0       0
	2       6       9       30      90      0       110     1       1       0       0
	3       10      19      100     150     80      140     2       0       1       1
	4       19      17      76      140     106     200     3       1       0       0
	5       19      17      76      140     55      105     3       1       0       0
	6       8       7       81      140     50      120     3       1       0       0
	7       8       7       20      80      25      90      3       1       0       0

	filter_spec file: a 11-column tab delimited file with the same format as gating_spec_file, but the events inside the gate will be removed instead of being kept:
	1       1       4       0       85      100     200     0       0       1	0

    first_pass_max_num_clusters is the maximum number of clusters for the first-pass k-means clustering

    second_pass_max_num_clusters is the maximum number of clusters for the additional pass k-means clustering (also re-normalized)

    num_threads is the number of threads available to partition the events for parallelization

    seed is the initialization number for random number generator for the cluster seeding
***********************************************************************************************************************/
#include "DAFi-gating_omp.h"

void get_spec_num(FILE *f_spec, int *num_rows) {
    int number_rows = 0;
    char line[LINE_LEN];

    line[0] = '\0';

    while (fgets(line, LINE_LEN, f_spec) != NULL) {
        number_rows++;
        line[0] = '\0';
    }

    printf("Number of filtering steps is %d\n", number_rows);
    *num_rows = number_rows;
}

//filtered_ID: 1 to num_rows
//filtered_d_x: 1 means the first dimension (number of dimensions starts from 1, i.e., FSC-A)
//filtered_d_x_low: a value between 0 to 200 (corresponding to 0 to 4095)
//filtered_parent: 1 means the parent of this population is the one at Row #1 (number of rows starts from 1)
//filtered_type: 0 means clustering (centers of clusters need to be inside the box); 1 means bisecting (all events in the box will be kept)
//filtered_output: 0 means not outputting the results (flock_results.txt and the filtered data (_filt.txt)); 1 means outputting both results for this population
//filtered_2nd_pass: 0 means using results from 1st pass clustering; 1 means to undergo second pass clustering - added 2/5/17 by Ivan
void get_spec_info(FILE *f_spec, int *filtered_ID, int *filtered_d_x, int *filtered_d_y, int *filtered_x_low,
                   int *filtered_x_upper, int *filtered_y_low, int *filtered_y_upper, int *filtered_parent,
                   int *filtered_type, int *filtered_output, int *filtered_2nd_pass, int num_rows) {
    int temp = 0;
    char line[LINE_LEN];

    line[0] = '\0';

    while ((fgets(line, LINE_LEN, f_spec) != NULL) && (temp < num_rows)) {
        sscanf(line, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", &filtered_ID[temp], &filtered_d_x[temp],
               &filtered_d_y[temp], &filtered_x_low[temp], &filtered_x_upper[temp], &filtered_y_low[temp],
               &filtered_y_upper[temp], &filtered_parent[temp], &filtered_type[temp], &filtered_output[temp], &filtered_2nd_pass[temp]);
        line[0] = '\0';
        temp++;
    }

    if (temp == num_rows)
        printf("The last row of the spec file is %d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", filtered_ID[temp - 1],
               filtered_d_x[temp - 1], filtered_d_y[temp - 1], filtered_x_low[temp - 1],
               filtered_x_upper[temp - 1], filtered_y_low[temp - 1], filtered_y_upper[temp - 1],
               filtered_parent[temp - 1], filtered_type[temp - 1], filtered_output[temp - 1], filtered_2nd_pass[temp - 1]);

}


/************* Read basic info of the source file ****************************/
void getfileinfo(FILE *f_src, int *file_Len, int *num_dm, char *name_string, int *time_ID) {
    char src[LINE_LEN];
    char current_name[64];
    char prv;

    int num_rows = 0;
    int num_columns = 0;
    int ch = '\n';
    int prev = '\n';
    int time_pos = 0;
    int i = 0;
    int j = 0;

    src[0] = '\0';
    fgets(src, LINE_LEN, f_src);

    if ((src[0] == 'F') && (src[1] == 'C') && (src[2] == 'S')) {
        fprintf(stderr, "the correct input format is a tab-delimited txt file, instead of FCS file.\n");
        abort();
    }

    name_string[0] = '\0';
    current_name[0] = '\0';
    prv = '\n';

    // skip space and tab characters
    while ((src[i] == ' ') || (src[i] == '\t'))
        i++;

    // repeat until the end of line is reached
    while ((src[i] != '\0') && (src[i] != '\n') && (src[i] != '\r')) {
        current_name[j] = src[i]; //each time only one character is copied

        if ((src[i] == '\t') && (prv != '\t')) //at the end of a complete word
        {
            current_name[j] = '\0';

            if (0 != strcmp(current_name, "Time")) {
                num_columns++; //num_columns does not inlcude the column of Time
                time_pos++;
                strcat(name_string, current_name);
                strcat(name_string, "\t");
            } else {
                *time_ID = time_pos;
            	}


            current_name[0] = '\0';
            j = 0;
        }

        if ((src[i] == '\t') && (prv == '\t')) //a duplicate tab or space, current_name needs to be reset
        {
            current_name[0] = '\0';
            j = 0;
        }

        if (src[i] != '\t') //a meaningful character, j can move forward
            j++;

        prv = src[i];
        i++; //in any case, i needs to move forward
    }

    if (prv != '\t') //the last word hasn't been retrieved because it doesn't end with a tab delimiter
    {
        current_name[j] = '\0';

        if (0 != strcmp(current_name, "Time")) {
            num_columns++;
            strcat(name_string, current_name);
            time_pos++;
        } else {
            *time_ID = time_pos;
        }


    }
    if (DEBUG == 1) {
        printf("time_ID is %d\n", *time_ID);
        printf("name_string is %s\n", name_string);
    }

    //start computing # of rows

    while ((ch = fgetc(f_src)) != EOF) {
        if (ch == '\n') {
            ++num_rows;
        }
        prev = ch;
    }
    if (prev != '\n')
        ++num_rows;

    //added on July 23, 2010
    if (num_rows < 50) {
        fprintf(stderr,
                "Number of events in the input file is too few and should not be processed!\n"); //modified on July 23, 2010
        abort();
    }

    *file_Len = num_rows;
    *num_dm = num_columns;

    printf("original file size is %d; number of dimensions is %d\n", *file_Len, *num_dm);
}


/************************************* Read the source file into uncomp_data **************************************/
void readsource(FILE *f_src, int file_Len, int num_dm, double **uncomp_data, int time_ID) {
    int time_pass = 0; //to mark whether the time_ID has been passed
    int index = 0;

    int i = 0;
    int j = 0;
    int t = 0;

    char src[LINE_LEN];
    char xc[LINE_LEN / 10];

    src[0] = '\0';
    fgets(src, LINE_LEN, f_src); //skip the first line about parameter names

    while (!feof(f_src) && (index < file_Len)) //index = 0, 1, ..., file_Len-1
    {
        src[0] = '\0';
        fgets(src, LINE_LEN, f_src);
        i = 0;
        time_pass = 0;

        if (time_ID == -1) {
            for (t = 0; t < num_dm; t++) //there is no time_ID
            {
                xc[0] = '\0';
                j = 0;
                while ((src[i] != '\0') && (src[i] != '\n') && (src[i] != ' ') && (src[i] != '\t')) {
                    xc[j] = src[i];
                    i++;
                    j++;
                }

                xc[j] = '\0';
                i++;

                uncomp_data[index][t] = atof(xc);
            }
        } else {
            for (t = 0; t <= num_dm; t++) //the time column needs to be skipped, so there are num_dm+1 columns
            {
                xc[0] = '\0';
                j = 0;
                while ((src[i] != '\0') && (src[i] != '\n') && (src[i] != ' ') && (src[i] != '\t')) {
                    xc[j] = src[i];
                    i++;
                    j++;
                }

                xc[j] = '\0';
                i++;

                if (t == time_ID) {
                    time_pass = 1;
                    continue;
                }

                if (time_pass)
                    uncomp_data[index][t - 1] = atof(xc);
                else
                    uncomp_data[index][t] = atof(xc);
            }
        }
        index++;
        //fprintf(fout_ID,"%s",src);
    } //end of while

    if (DEBUG == 1) {
        printf("the last line of the source data is:\n");
        for (j = 0; j < num_dm; j++)
            printf("%f ", uncomp_data[index - 1][j]);
        printf("\n");
    }
}


/**************************************** Normalization ******************************************/
void tran(double **orig_data, int file_Len, int num_dm, int norm_used, double **matrix_to_cluster) {
    int i = 0;
    int j = 0;

    double biggest = 0;
    double smallest = MAX_VALUE;

    double *aver; //average of each column
    double *std; //standard deviation of each column

    aver = (double *) malloc(sizeof(double) * file_Len);
    memset(aver, 0, sizeof(double) * file_Len);

    std = (double *) malloc(sizeof(double) * file_Len);
    memset(std, 0, sizeof(double) * file_Len);

    if (norm_used == 2) //z-score normalization
    {
        for (j = 0; j < num_dm; j++) {
            aver[j] = 0;
            for (i = 0; i < file_Len; i++)
                aver[j] = aver[j] + orig_data[i][j];
            aver[j] = aver[j] / (double) file_Len;

            std[j] = 0;
            for (i = 0; i < file_Len; i++)
                std[j] = std[j] + (orig_data[i][j] - aver[j]) * (orig_data[i][j] - aver[j]);
            std[j] = sqrt(std[j] / (double) file_Len);

            for (i = 0; i < file_Len; i++)
                matrix_to_cluster[i][j] = (orig_data[i][j] - aver[j]) / std[j];  //z-score normalization
        }
    }

    if (norm_used == 1) //0-1 min-max normalization
    {
        for (j = 0; j < num_dm; j++) {
            biggest = 0;
            smallest = MAX_VALUE;
            for (i = 0; i < file_Len; i++) {
                if (orig_data[i][j] > biggest)
                    biggest = orig_data[i][j];
                if (orig_data[i][j] < smallest)
                    smallest = orig_data[i][j];
            }

            for (i = 0; i < file_Len; i++) {
                if (biggest == smallest)
                    matrix_to_cluster[i][j] = 0.5;
                else
                    matrix_to_cluster[i][j] = (orig_data[i][j] - smallest) / (biggest - smallest);
            }
        }
    }

    if (norm_used == 3) //0-1 min-max normalization based on VALUE_RANGE (by default it is 4095)
    {
        for (j = 0; j < num_dm; j++) {
            for (i = 0; i < file_Len; i++)
                matrix_to_cluster[i][j] = (orig_data[i][j]) / (double) VALUE_RANGE;
        }
    }

    if (norm_used == 0) //no normalization
    {
        for (i = 0; i < file_Len; i++)
            for (j = 0; j < num_dm; j++)
                matrix_to_cluster[i][j] = orig_data[i][j];
    }

    free(aver);
    free(std);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Compute Population Center with all events
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
ID2Center_all(double **data_in, int file_Len, int num_dm, int num_clust, int *cluster_ID, double **population_center) {
    int i = 0;
    int j = 0;
    int ID = 0;
    int *size_c;


    size_c = (int *) malloc(sizeof(int) * num_clust);
    memset(size_c, 0, sizeof(int) * num_clust);

    for (i = 0; i < num_clust; i++)
        for (j = 0; j < num_dm; j++)
            population_center[i][j] = 0;

    for (i = 0; i < file_Len; i++) {
        ID = cluster_ID[i];

        if (ID == -1) {
            //commented on July 23, 2010
            //printf("ID==-1! in ID2Center_all\n");
            //exit(0);
            fprintf(stderr,
                    "Incorrect file format or input parameters (resulting in incorrect population IDs)\n"); //modified on July 23, 2010
            abort();
        }

        for (j = 0; j < num_dm; j++)
            population_center[ID][j] = population_center[ID][j] + data_in[i][j];

        size_c[ID]++;
    }


    for (i = 0; i < num_clust; i++) {
        for (j = 0; j < num_dm; j++)
            if (size_c[i] != 0)
                population_center[i][j] = (population_center[i][j] / (double) (size_c[i]));
            else
                population_center[i][j] = 0;
    }

    free(size_c);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Merge neighboring grids to clusters
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double
recursive_optimize(double **Matrix, int k, double terms, int file_Len, int num_dm, int *shortest_id, double **center,
                   int random_init, int num_threads, int rand_seed) {

  int i, j, t;
  int random, random1, random2;
  int times;
  int dist_used=0;     //0 is Euclidean and 1 is Pearson
  int real_Len;
  int skipped=0;

  int *num;  //num[i]=t means the ith cluster has t points

  double vvv=1.0; // the biggest variation;
  double distance=0.0;
  double xv=0.0;
  double variation=0.0;
  double sum_var;
  double mean_dx, mean_dy, dx, dy, sd_x, sd_y;
  double diff;
  double distortion=0;
  double shortest_distance;

  //double distancev[k];               //RSS
  double center_transpose[num_dm*k]; //RSS
  int counter;                       //RSS

  double *temp_center;

  double **sum;

  temp_center = (double *)malloc(sizeof(double)*num_dm);
  memset(temp_center, 0, sizeof(double)*num_dm);

  //RSS Added call to srand so that we can be assured of
  //RSS repeatable random numbers
    if (random_init != 0) {
        srand(rand_seed);
        for (i=0; i<k; i++) {
            for (j=0; j<num_dm; j++) {
                random = rand() % (file_Len - 0);
                center[i][j] = Matrix[random][j];
            }
        }
    }


  num = (int *)malloc(sizeof(int)*k);
  memset(num, 0, sizeof(int)*k);

  sum = (double **)malloc(sizeof(double*)*k);
  memset(sum, 0, sizeof(double*)*k);
  for (i=0; i<k; i++) {
    sum[i] = (double *)malloc(sizeof(double)*num_dm);
    memset(sum[i], 0, sizeof(double)*num_dm);
  }

  times    = 0;
  real_Len = 0;

  while ( (vvv>terms && terms<1) || (times<terms && terms>=1) ) {
    for (i=0; i<k; i++) {
      num[i] = 0;
      for (j=0; j<num_dm; j++)
	sum[i][j] = 0.0;
    }

    //RSS Transpose of centers into 1D array (used in optimized code)
    counter = 0;
    for (t=0; t<num_dm; t++) {
      for (j=0; j<k; j++) {
	center_transpose[counter] = center[j][t];
	counter++;
      }
    }

    #pragma omp parallel for reduction(+:real_Len)
    for (i=0; i<file_Len; i++) { // for each data point i, compute the distance between Matrix[i] and center[j]

        int skipped           = 0;
        int j, t;
        double shortest_distance = MAX_VALUE;
        double distancev[k];
        double diff;

        //RSS Euclidean distance optimized
        if (dist_used == 0) {
        	for (j=0; j<k; j++) {
        	  distancev[j] = 0.0;
        	}

        	int counter = 0;
        	for (t=0; t<num_dm; t++) {
                for (j=0; j<k; j++) {
                    diff = center_transpose[counter] - Matrix[i][t];
                    diff = diff * diff;
                    distancev[j] = distancev[j] + diff;
                    counter++;
                }
        	}

        	for (j=0; j<k; j++) {
        	  if (distancev[j] < shortest_distance && skipped == 0) {
        	    shortest_distance=distancev[j];
        	    shortest_id[i]=j;
        	  }
        	}
        }
          // End Euclidean distance

        // Pearson distance
        // WARNING - have not addressed public/private variables for Pearson distance
        if (dist_used != 0) {
        	for (j=0; j<k; j++) { // iterate over centers
        	  mean_dx = 0.0;
        	  mean_dy = 0.0;
        	  sum_var = 0.0;
        	  dx      = 0.0;
        	  dy      = 0.0;
        	  sd_x    = 0.0;
        	  sd_y    = 0.0;
        	  for (t=0; t<num_dm; t++) {
        	    mean_dx += center[j][t];
        	    mean_dy += Matrix[i][t];
        	  }
        	  mean_dx = mean_dx/(double)num_dm;
        	  mean_dy = mean_dy/(double)num_dm;

        	  for (t=0; t<num_dm; t++) {
        	    dx       = center[j][t] - mean_dx;
        	    dy       = Matrix[i][t] - mean_dy;
        	    sum_var += dx*dy;
        	    sd_x    += dx*dx;
        	    sd_y    += dy*dy;
        	  }
        	  if (sqrt(sd_x*sd_y) == 0)
        	    distance = 1.0;
        	  else
        	    distance = 1.0 - (sum_var/(sqrt(sd_x*sd_y))); // distance ranges from 0 to 2;
        	  //printf("distance=%f\n",distance);
        	}	//pearson correlation ends

        	if (distance < shortest_distance && skipped == 0) {
        	  shortest_distance = distance;
        	  shortest_id[i]    = j;
        	}
        }
          // End Pearson distance

          real_Len++;
          //num[shortest_id[i]]=num[shortest_id[i]]+1;
          //for (t=0;t<num_dm;t++) {
          //sum[shortest_id[i]][t]=sum[shortest_id[i]][t]+Matrix[i][t];
          //}
    } //end for i

    // Moved outside of parallel loop to avoid crtical regions
    for (i=0; i<file_Len; i++) {
      num[shortest_id[i]]=num[shortest_id[i]]+1;
      for (t=0;t<num_dm;t++) {
	    sum[shortest_id[i]][t]=sum[shortest_id[i]][t]+Matrix[i][t];
      }
    }

    /* recompute the centers */
    //compute_mean(group);
    vvv = 0.0;
    for (j=0; j<k; j++) {
      memcpy(temp_center, center[j], sizeof(double)*num_dm);
      variation = 0.0;
      if (num[j] != 0) {
	    for (t=0; t<num_dm; t++) {
	        center[j][t] = sum[j][t]/(double)num[j];
	        xv           = (temp_center[t]-center[j][t]);
	        variation    = variation + xv*xv;
	    }
      }

      if (variation > vvv)
	    vvv=variation; //vvv is the biggest variation among the k clusters;
    }

    //compute_variation;
    times++;
    //printf("Times recursive optimize ran: %d\n", times);
    //printf("Biggest variation: %d\n", vvv);
  } //end for while


  free(num);
   for (i=0; i<k; i++)
    free(sum[i]);
  free(sum);
  free(temp_center);

  return distortion;
}



/*************************** Show - produce profile.txt and percentage.txt *****************************/
void show(double **Matrix, int *cluster_id, int file_Len, int k, int num_dm, char *name_string, int pop) {
    int situ1 = 0;
    int situ2 = 0;

    int i = 0;
    int id = 0;
    int j = 0;
    int t = 0;

    int *size_c;


    int **size_mybound_1;
    int **size_mybound_2;
    int **size_mybound_3;
    int **size_mybound_0;

    double interval = 0.0;

    double *big;
    double *small;


    double **center;
    double **mybound;

    int **prof; //prof[i][j]=1 means population i is + at parameter j

    FILE *fpcnt_id; //proportion id
    //FILE *fcent_id; //center_id, i.e., centers of clusters within the original data
    FILE *fprof_id; //profile_id

    big = (double *) malloc(sizeof(double) * num_dm);
    memset(big, 0, sizeof(double) * num_dm);

    small = (double *) malloc(sizeof(double) * num_dm);
    memset(small, 0, sizeof(double) * num_dm);

    for (i = 0; i < num_dm; i++) {
        big[i] = 0.0;
        small[i] = (double) MAX_VALUE;
    }


    size_c = (int *) malloc(sizeof(int) * k);
    memset(size_c, 0, sizeof(int) * k);

    center = (double **) malloc(sizeof(double *) * k);
    memset(center, 0, sizeof(double *) * k);
    for (i = 0; i < k; i++) {
        center[i] = (double *) malloc(sizeof(double) * num_dm);
        memset(center[i], 0, sizeof(double) * num_dm);
    }

    mybound = (double **) malloc(sizeof(double *) * num_dm);
    memset(mybound, 0, sizeof(double *) * num_dm);
    for (i = 0; i < num_dm; i++) //there are 3 mybounds for 4 categories
    {
        mybound[i] = (double *) malloc(sizeof(double) * 3);
        memset(mybound[i], 0, sizeof(double) * 3);
    }

    prof = (int **) malloc(sizeof(int *) * k);
    memset(prof, 0, sizeof(int *) * k);
    for (i = 0; i < k; i++) {
        prof[i] = (int *) malloc(sizeof(int) * num_dm);
        memset(prof[i], 0, sizeof(int) * num_dm);
    }


    for (i = 0; i < file_Len; i++) {
        id = cluster_id[i];
        for (j = 0; j < num_dm; j++) {
            center[id][j] = center[id][j] + Matrix[i][j];
            if (big[j] < Matrix[i][j])
                big[j] = Matrix[i][j];
            if (small[j] > Matrix[i][j])
                small[j] = Matrix[i][j];
        }

        size_c[id]++;
    }

    for (i = 0; i < k; i++)
        for (j = 0; j < num_dm; j++) {
            if (size_c[i] != 0)
                center[i][j] = (center[i][j] / (double) (size_c[i]));
            else
                center[i][j] = 0;
        }

    for (j = 0; j < num_dm; j++) {
        interval = ((big[j] - small[j]) / 4.0);
        //printf("interval[%d] is %f\n",j,interval);
        for (i = 0; i < 3; i++)
            mybound[j][i] = small[j] + ((double) (i + 1) * interval);
    }


    size_mybound_0 = (int **) malloc(sizeof(int *) * k);
    memset(size_mybound_0, 0, sizeof(int *) * k);

    for (i = 0; i < k; i++) {
        size_mybound_0[i] = (int *) malloc(sizeof(int) * num_dm);
        memset(size_mybound_0[i], 0, sizeof(int) * num_dm);
    }

    size_mybound_1 = (int **) malloc(sizeof(int *) * k);
    memset(size_mybound_1, 0, sizeof(int *) * k);

    for (i = 0; i < k; i++) {
        size_mybound_1[i] = (int *) malloc(sizeof(int) * num_dm);
        memset(size_mybound_1[i], 0, sizeof(int) * num_dm);
    }

    size_mybound_2 = (int **) malloc(sizeof(int *) * k);
    memset(size_mybound_2, 0, sizeof(int *) * k);

    for (i = 0; i < k; i++) {
        size_mybound_2[i] = (int *) malloc(sizeof(int) * num_dm);
        memset(size_mybound_2[i], 0, sizeof(int) * num_dm);
    }

    size_mybound_3 = (int **) malloc(sizeof(int *) * k);
    memset(size_mybound_3, 0, sizeof(int *) * k);

    for (i = 0; i < k; i++) {
        size_mybound_3[i] = (int *) malloc(sizeof(int) * num_dm);
        memset(size_mybound_3[i], 0, sizeof(int) * num_dm);
    }

    for (i = 0; i < file_Len; i++)
        for (j = 0; j < num_dm; j++) {
            if (Matrix[i][j] < mybound[j][0])// && ((Matrix[i][j]-small[j])>0)) //the smallest values excluded
                size_mybound_0[cluster_id[i]][j]++;
            else {
                if (Matrix[i][j] < mybound[j][1])
                    size_mybound_1[cluster_id[i]][j]++;
                else {
                    if (Matrix[i][j] < mybound[j][2])
                        size_mybound_2[cluster_id[i]][j]++;
                    else
                        //if (Matrix[i][j]!=big[j]) //the biggest values excluded
                        size_mybound_3[cluster_id[i]][j]++;
                }

            }
        }

    char prof_id_name[LINE_LEN];

    snprintf(prof_id_name, sizeof(char) * LINE_LEN, "./pop%i/profile.txt", pop);

    fprof_id = fopen(prof_id_name, "w");
    fprintf(fprof_id, "Population_ID\t");
    fprintf(fprof_id, "%s\n", name_string);

    for (i = 0; i < k; i++) {
        fprintf(fprof_id, "%d\t", i + 1); //i changed to i+1 to start from 1 instead of 0: April 16, 2009
        for (j = 0; j < num_dm; j++) {

            if (size_mybound_0[i][j] > size_mybound_1[i][j])
                situ1 = 0;
            else
                situ1 = 1;
            if (size_mybound_2[i][j] > size_mybound_3[i][j])
                situ2 = 2;
            else
                situ2 = 3;

            if ((situ1 == 0) && (situ2 == 2)) {
                if (size_mybound_0[i][j] > size_mybound_2[i][j])
                    prof[i][j] = 0;
                else
                    prof[i][j] = 2;
            }
            if ((situ1 == 0) && (situ2 == 3)) {
                if (size_mybound_0[i][j] > size_mybound_3[i][j])
                    prof[i][j] = 0;
                else
                    prof[i][j] = 3;
            }
            if ((situ1 == 1) && (situ2 == 2)) {
                if (size_mybound_1[i][j] > size_mybound_2[i][j])
                    prof[i][j] = 1;
                else
                    prof[i][j] = 2;
            }
            if ((situ1 == 1) && (situ2 == 3)) {
                if (size_mybound_1[i][j] > size_mybound_3[i][j])
                    prof[i][j] = 1;
                else
                    prof[i][j] = 3;
            }

            //begin to output profile
            if (j == num_dm - 1) {
                if (prof[i][j] == 0)
                    fprintf(fprof_id, "1\n");
                if (prof[i][j] == 1)
                    fprintf(fprof_id, "2\n");
                if (prof[i][j] == 2)
                    fprintf(fprof_id, "3\n");
                if (prof[i][j] == 3)
                    fprintf(fprof_id, "4\n");
            } else {
                if (prof[i][j] == 0)
                    fprintf(fprof_id, "1\t");
                if (prof[i][j] == 1)
                    fprintf(fprof_id, "2\t");
                if (prof[i][j] == 2)
                    fprintf(fprof_id, "3\t");
                if (prof[i][j] == 3)
                    fprintf(fprof_id, "4\t");
            }
        }
    }
    fclose(fprof_id);

    ///////////////////////////////////////////////////////////
    char pcnt_id_name[LINE_LEN];

    snprintf(pcnt_id_name, sizeof(char) * LINE_LEN, "./pop%i/percentage.txt", pop);

    fpcnt_id = fopen(pcnt_id_name, "w");
    fprintf(fpcnt_id, "Population_ID\tPercentage\n");

    for (t = 0; t < k; t++) {
        fprintf(fpcnt_id, "%d\t%.4f\n", t + 1, (double) size_c[t] * 100.0 /
                                               (double) file_Len);    //t changed to t+1 to start from 1 instead of 0: April 16, 2009
    }
    fclose(fpcnt_id);

    free(big);
    free(small);
    free(size_c);

    for (i = 0; i < k; i++) {
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

    for (i = 0; i < num_dm; i++)
        free(mybound[i]);
    free(mybound);

}

/**************************************** Calculate number of filtered events  ******************************************/
int numOfEventsInGate(int currentGate, int len, int *event_gate, int **event_parent_gates, int *filtered_parent) {
    int j, t;
    int filtered_out;
    int sizeOfGate=0;

    for (j = 0; j < len; j++) {
            filtered_out = 0;
            event_gate[j] = 1; //1 means to remove

            t = currentGate;   //note that filtered_parent[0]=0;
            while (t >= 0) //the intersection of multiple gates needs to be calculated; iterates through the parent gates
            {
                if (event_parent_gates[j][t] > 0) {
                    filtered_out = 1;  //this event has been removed in some parent gate
                    break;
                }
                t = filtered_parent[t] - 1; //please note that the filtered_parent[x] returns 0 to num_rows-1, however, it needs to be changed to -1 to num_rows-2 to be used
            }

            if (filtered_out == 0) // neither removed by parent gates nor by the current gate
            {
                event_gate[j] = 0;    //0 means to keep
                sizeOfGate++;    //increment number of events kept in current gate
            }
    }
    return sizeOfGate;
}

/**************************** Assign sub populations from parent population ********************************************/
void assignPopulation(int currentPop, int len, int *event_gate, double **input_data, double ***pop_data, int **pop_ID_map){

    int i = 0, j;
        for (j = 0; j < len; j++) {
            if (event_gate[j] == 0)
            {
                pop_data[currentPop][i] = input_data[j];
                pop_ID_map[currentPop][i] = j;
                i++;
            }
        }
}


/******************************************************** Main Function **************************************************/

int main(int argc, char **argv) {
    //inputs
    FILE *f_src; //source file pointer

    FILE *f_out; //coordinates
    FILE *f_cid; //population-ID of events
    FILE *f_ctr; //centroids of populations
    FILE *f_cent; //centroids of k-means clusters
    FILE *f_results; //coordinates file event and population column
    FILE *f_mfi; //added April 16, 2009 for mean fluorescence intensity
    FILE *f_parameters; //number of bins and density calculated by
    //the algorithm. Used to update the database
    FILE *f_properties; //Properties file used by Image generation software
    FILE *f_spec; //user spec file for data filtering
    FILE *f_spec_2;
    FILE *f_filtered; //filtered result
    FILE *f_final_filtered; //user-selected filtered result
    FILE *f_majority; //proportion of a cluster inside the hyperegion in filtering mode

    FILE *f_filtered_percentage;
    FILE *f_filtered_events;
    FILE *f_filtered_MFI;

    char para_name_string[LINE_LEN];
    char file_name[LINE_LEN];
    char f_name[LINE_LEN];
    char f_selected_name[LINE_LEN];
    char f_selected_file_name[LINE_LEN];


    int file_Len = 0; //number of events from whole file
    int num_dm = 0;
    int time_ID = -1;
    int num_population = 0;
    int *num_sub_pop;
    //int temp=0;

    //below are read from configuration file
    int i = 0;
    int j = 0;
    int t = 0;
    int p = 0;
    int x = 0;
    int start_num_pop = 0; //Ivan: start_num_pop to denote starting number of population for adaptive clustering
    int max_num_pop = 0;
    int NUM_THREADS;
    int SEED;

    int num_rows = 0;
    int num_rows_2 = 0;

    int filtered_out = 0;
    int filtered_output_finished = 0; //indicating whether this is the first user-specified output population

    int *all_population_ID; //populationID of event (event -> cluster mapping)
    int **sub_population_ID; //Ivan: populationID of event (event -> cluster mapping) for next clustering run

    int **pop_ID_map; //Ivan: mapping of sub-predefined cell population ID to the whole file

    int *filtered_ID;
    int *filtered_d_x;
    int *filtered_d_y;
    int *filtered_x_low;
    int *filtered_x_upper;
    int *filtered_y_low;
    int *filtered_y_upper;
    int *filtered_parent;
    int *filtered_type;
    int *filtered_output;
    int *filtered_2nd_pass;

    int *filtered_ID_2;
    int *filtered_d_x_2;
    int *filtered_d_y_2;
    int *filtered_x_low_2;
    int *filtered_x_upper_2;
    int *filtered_y_low_2;
    int *filtered_y_upper_2;
    int *filtered_parent_2;
    int *filtered_type_2;
    int *filtered_output_2;
    int *filtered_2nd_pass_2;

    int *size_c; //number of events in each cluster;

    int *filtered_p;
    int *tmp_filtered_p;

    int *size_filtering;
    int *final_gate_ID;
    int **all_gate_ID;
    int *tmp_gate_ID;

    double *filtered_percentage;

    double **population_center; //population centroids in the raw/original data
    double ***sub_population_center; //Ivan: population centroids in the sub population data clusters

    double **input_data;
    double **normalized_data;

    double ***pop_data; //Ivan: events kept for each sub-population
    double ***norm_pop_data;

    double **gate_center;


    int min = 999999;
    int max = 0;

    int pppp = 0;

    int optimized = 0;


    printf("Starting time:\t\t\t\t");
    fflush(stdout);
    system("/bin/date");
    /////////////////////////////////////////////////////////////

    if ((argc != 4) && (argc != 5) && (argc != 6) && (argc != 7) && (argc != 8)) {
        fprintf(stderr, "Incorrect number of input parameters!\n"); //modified on July 23, 2010
        fprintf(stderr,
                "DAFi data_file filter_spec_file1 filter_spec_file2\n"); //filter_spec_file 1 is those that need to be kept; filter_spec_file2 is those that need to be removed
        fprintf(stderr, "DAFi data_file filter_spec_file1 filter_spec_file2 num_initial_clusters num_2ndpass_clusters num_threads rand_seed\n");
        abort();
    }

    file_name[0] = '\0';
    strcpy(file_name, argv[1]);
    f_src = fopen(argv[1], "r");

    f_name[0] = '\0';
    f_selected_name[0] = '\0';
    while (file_name[i] != '.') {
        f_name[i] = file_name[i];
        f_selected_name[i] = file_name[i];
        i++;
    }
    f_name[i] = '\0';
    f_selected_name[i] = '\0';

    printf("file name is %s\n", f_name);
    i = 0;

    if (argc == 4) {
        start_num_pop = DEFAULT_NUM_POP; //default value = starting 100 clusters
        max_num_pop = DEFAULT_NUM_POP;
        NUM_THREADS = 2;
        SEED = 2;
    }

    if (argc == 5) {
        start_num_pop = atoi(argv[4]);
        max_num_pop = DEFAULT_NUM_POP;
        NUM_THREADS = 2;
        SEED = 2;
        printf("starting number of clusters is %d\n", start_num_pop);
    }

    if (argc == 6) {
        start_num_pop = atoi(argv[4]);
        max_num_pop = atoi(argv[5]);
        NUM_THREADS = 2;
        SEED = 2;
        printf("starting number of clusters is %d\n", start_num_pop);
        printf("number of thread used is %d\n", NUM_THREADS);
    }

    if (argc == 7) {
        start_num_pop = atoi(argv[4]);
        max_num_pop = atoi(argv[5]);
        NUM_THREADS = atoi(argv[6]);
        SEED = 2;
        printf("starting number of clusters is %d\n", start_num_pop);
        printf("maximum number of clusters is %d\n", max_num_pop);
        printf("number of thread used is %d\n", NUM_THREADS);
    }

    if (argc == 8) {
        start_num_pop = atoi(argv[4]);
        max_num_pop = atoi(argv[5]);
        NUM_THREADS = atoi(argv[6]);
        SEED = atoi(argv[7]);
        printf("starting number of clusters is %d\n", start_num_pop);
        printf("maximum number of clusters is %d\n", max_num_pop);
        printf("number of thread used is %d\n", NUM_THREADS);
    }

    f_spec = fopen(argv[2], "r");
    f_spec_2 = fopen(argv[3], "r");


    get_spec_num(f_spec, &num_rows);
    get_spec_num(f_spec_2, &num_rows_2);

    rewind(f_spec);
    rewind(f_spec_2);

    filtered_ID = (int *) malloc(sizeof(int) * num_rows);
    memset(filtered_ID, 0, sizeof(int) * num_rows);

    filtered_d_x = (int *) malloc(sizeof(int) * num_rows);
    memset(filtered_d_x, 0, sizeof(int) * num_rows);

    filtered_d_y = (int *) malloc(sizeof(int) * num_rows);
    memset(filtered_d_y, 0, sizeof(int) * num_rows);

    filtered_x_low = (int *) malloc(sizeof(int) * num_rows);
    memset(filtered_x_low, 0, sizeof(int) * num_rows);

    filtered_x_upper = (int *) malloc(sizeof(int) * num_rows);
    memset(filtered_x_upper, 0, sizeof(int) * num_rows);

    filtered_y_low = (int *) malloc(sizeof(int) * num_rows);
    memset(filtered_y_low, 0, sizeof(int) * num_rows);

    filtered_y_upper = (int *) malloc(sizeof(int) * num_rows);
    memset(filtered_y_upper, 0, sizeof(int) * num_rows);

    filtered_parent = (int *) malloc(sizeof(int) * num_rows);
    memset(filtered_parent, 0, sizeof(int) * num_rows);

    filtered_type = (int *) malloc(sizeof(int) * num_rows);
    memset(filtered_type, 0, sizeof(int) * num_rows);

    filtered_output = (int *) malloc(sizeof(int) * num_rows);
    memset(filtered_output, 0, sizeof(int) * num_rows);

    filtered_2nd_pass = (int *) malloc(sizeof(int) * num_rows);
    memset(filtered_2nd_pass, 0, sizeof(int) * num_rows);

    /////////////////

    filtered_ID_2 = (int *) malloc(sizeof(int) * num_rows_2);
    memset(filtered_ID_2, 0, sizeof(int) * num_rows_2);

    filtered_d_x_2 = (int *) malloc(sizeof(int) * num_rows_2);
    memset(filtered_d_x_2, 0, sizeof(int) * num_rows_2);

    filtered_d_y_2 = (int *) malloc(sizeof(int) * num_rows_2);
    memset(filtered_d_y_2, 0, sizeof(int) * num_rows_2);

    filtered_x_low_2 = (int *) malloc(sizeof(int) * num_rows_2);
    memset(filtered_x_low, 0, sizeof(int) * num_rows_2);

    filtered_x_upper_2 = (int *) malloc(sizeof(int) * num_rows_2);
    memset(filtered_x_upper_2, 0, sizeof(int) * num_rows_2);

    filtered_y_low_2 = (int *) malloc(sizeof(int) * num_rows_2);
    memset(filtered_y_low_2, 0, sizeof(int) * num_rows_2);

    filtered_y_upper_2 = (int *) malloc(sizeof(int) * num_rows_2);
    memset(filtered_y_upper_2, 0, sizeof(int) * num_rows_2);

    filtered_parent_2 = (int *) malloc(sizeof(int) * num_rows_2);
    memset(filtered_parent_2, 0, sizeof(int) * num_rows_2);

    filtered_type_2 = (int *) malloc(sizeof(int) * num_rows_2);
    memset(filtered_type_2, 0, sizeof(int) * num_rows_2);

    filtered_output_2 = (int *) malloc(sizeof(int) * num_rows_2);
    memset(filtered_output_2, 0, sizeof(int) * num_rows_2);

    filtered_2nd_pass_2 = (int *) malloc(sizeof(int) * num_rows);
    memset(filtered_2nd_pass_2, 0, sizeof(int) * num_rows);

    get_spec_info(f_spec, filtered_ID, filtered_d_x, filtered_d_y, filtered_x_low, filtered_x_upper, filtered_y_low,
                  filtered_y_upper, filtered_parent, filtered_type, filtered_output, filtered_2nd_pass, num_rows);
    get_spec_info(f_spec_2, filtered_ID_2, filtered_d_x_2, filtered_d_y_2, filtered_x_low_2, filtered_x_upper_2,
                  filtered_y_low_2, filtered_y_upper_2, filtered_parent_2, filtered_type_2, filtered_output_2, filtered_2nd_pass_2,
                  num_rows_2);

    fclose(f_spec);
    fclose(f_spec_2);

    getfileinfo(f_src, &file_Len, &num_dm, para_name_string, &time_ID); //get the filelength, number of dimensions, and num/name of parameters

    i=0;
    while (para_name_string[i]!='\0')
	i++;

    if ((para_name_string[i-1]=='\t') || (para_name_string[i-1]=='\r') || (para_name_string[i-1]=='\n'))
	    para_name_string[i-1]='\0';

    i=0;

    /************************************************* Read the data *****************************************************/

    rewind(f_src); //reset the file pointer

    input_data = (double **) malloc(sizeof(double *) * file_Len);
    memset(input_data, 0, sizeof(double *) * file_Len);
    for (i = 0; i < file_Len; i++) {
        input_data[i] = (double *) malloc(sizeof(double) * num_dm);
        memset(input_data[i], 0, sizeof(double) * num_dm);
    }

    readsource(f_src, file_Len, num_dm, input_data, time_ID); //read the data;

    fclose(f_src);

    normalized_data = (double **) malloc(sizeof(double *) * file_Len);
    memset(normalized_data, 0, sizeof(double *) * file_Len);
    for (i = 0; i < file_Len; i++) {
        normalized_data[i] = (double *) malloc(sizeof(double) * num_dm);
        memset(normalized_data[i], 0, sizeof(double) * num_dm);
    }

    /**** Normalization of the raw input data ****/
    printf("Normalization of the raw input data...\n");
    tran(input_data, file_Len, num_dm, NORM_METHOD, normalized_data);

    /**** Memory allocation for starting clustering routine ****/
    //printf("Memory allocation for starting clustering routine...\n");

    num_population = start_num_pop; //Ivan: starting number of clusters

    population_center = (double **) malloc(sizeof(double *) * num_population); //centroids of the clusters
    memset(population_center, 0, sizeof(double *) * num_population);

    for (i = 0; i < num_population; i++) {
        population_center[i] = (double *) malloc(sizeof(double) * num_dm);
        memset(population_center[i], 0, sizeof(double) * num_dm);
    }

    all_population_ID = (int *) malloc(sizeof(int) * file_Len);
    memset(all_population_ID, 0, sizeof(int) * file_Len);


    /**** Memory allocation for initial filtering ****/
    //printf("Memory allocation for initial filtering ...\n");

    printf("num_population = %d\n", num_population);
    filtered_p = (int *) malloc(sizeof(int) * num_population); //for each cluster, whether it is in or out in the data kept
    memset(filtered_p, 0, sizeof(int) * num_population);

    printf("file_Len = %d\n", file_Len);
    printf("num_rows = %d\n", num_rows);
    all_gate_ID = (int **) malloc(sizeof(int *) * file_Len); //for all predefined cell populations, each constrained by only the 2D boundary, but their intersection will be saved into final_gate_ID
    memset(all_gate_ID, 0, sizeof(int *) * file_Len);
    for (i = 0; i < file_Len; i++) {
        all_gate_ID[i] = (int *) malloc(sizeof(int) * num_rows);
        memset(all_gate_ID[i], 0, sizeof(int) * num_rows);
    }
    //printf("Memory allocation for initial filtering/2...\n");


    size_filtering = (int *) malloc(sizeof(int) * num_rows); //the final size (number of events) of each predefined cell population
    memset(size_filtering, 0, sizeof(int) * num_rows);

    filtered_percentage = (double *) malloc(sizeof(double) * num_rows); //based on the parent population
    memset(filtered_percentage, 0, sizeof(double) * num_rows);

    final_gate_ID = (int *) malloc(sizeof(int) * file_Len); //the final membership ID of each predefined cell population
    memset(final_gate_ID, 0, sizeof(int) * file_Len);

    //printf("Memory allocation for initial filtering/3s ...\n");

    gate_center = (double **) malloc(sizeof(double *) * 2); //for centroids of the two clusters (the kept and the discarded), which change for every predefined gate/cell population
    memset(gate_center, 0, sizeof(double *) * 2);
    for (i = 0; i < 2; i++) {
        gate_center[i] = (double *) malloc(sizeof(double) * num_dm);
        memset(gate_center[i], 0, sizeof(double) * num_dm);
    }

    /**** Memory allocation for adaptive clustering routine ****/
    //printf("Memory allocation for adaptive clustering routine...\n");

    pop_data = (double ***) malloc(sizeof(double **) * num_rows);
    memset(pop_data, 0, sizeof(double **) * num_rows);

    norm_pop_data = (double ***) malloc(sizeof(double **) * num_rows);
    memset(norm_pop_data, 0, sizeof(double **) * num_rows);

    num_sub_pop = (int *) malloc(sizeof(int) * num_rows); //Ivan: final number of clusters for each predefined cell population
    memset(num_sub_pop, 0, sizeof(int) * num_rows);

    pop_ID_map = (int **) malloc(sizeof(int *) * num_rows);
    memset(pop_ID_map, 0, sizeof(int *) * num_rows);

    sub_population_center = (double ***) malloc(sizeof(double **) * num_rows);
    memset(sub_population_center, 0, sizeof(double **) * num_rows);

    sub_population_ID = (int **) malloc(sizeof(int *) * num_rows);
    memset(sub_population_ID, 0, sizeof(int *) * num_rows);

    //You need to define size_filtering[num_rows] for number of events in each predefined population
    //You need to define filtered_percentage[num_rows] for percentage of each predefined population i, based on size_filtering[filtered_parent[i]]

    //Algorithm:
    //For each predefined population
    //check the type of extraction (clusters or bisecting)
    //extract it based on the coordinates of the user-specified boundaries and the type of extraction
    //calculate and save its number of events
    //based on its parent, calculate its percentage
    //save its percentage
    //based on its output type, calculate flock_results.txt for visualization
    //complete the above procedure for each predefined population
    //output the number of events and percentage for all predefined populations


    f_filtered_percentage = fopen("DUnSup_pop_percentage.txt", "w");
    f_majority = fopen("ClusterFilter_proportion.txt", "w");
    f_filtered_events = fopen("DUnSup_pop_events.txt", "w");
    f_filtered_MFI = fopen("DUnSup_pop_MFI.txt", "w");

    fprintf(f_filtered_percentage, "Predefined_Population_ID\tPercentageBasedonParent\n");
    fprintf(f_filtered_events, "Predefined_Population_ID\tNumofEvents\n");

    printf("Starting adaptive clustering for each predefined population...\n");

    int parentPop;
    int parentSize;
    int lastClustering = -1;


    size_c = (int *) malloc(sizeof(int) * num_population); //size of each cluster
    memset(size_c, 0, sizeof(int) * num_population);
    for (i = 0; i < num_population; i++)
	    size_c[i] = 0;

    filtered_output_finished = 0;
    //for each predefined population i, check if further clustering is necessary (adaptive clustering), generate population_ID, based on population_ID and filtered_parent, calculate number of events and percentage and write into file
    for (i = 0; i < num_rows; i++)
    {

        printf("Predefined cell population #: %d\n", i);
        size_filtering[i] = 0;
        filtered_percentage[i] = 0.0;
        parentPop = filtered_parent[i]-1;
        parentSize = size_filtering[parentPop];

        if (filtered_type[i] == 0) //clustering
        {
            /*** First-run k-means clustering for the whole file ****/
            if (lastClustering==-1) {
              printf("Starting initial clustering for whole cell population...\n");
                recursive_optimize(normalized_data, num_population, TERMS, file_Len, num_dm, all_population_ID, population_center,
                       1, NUM_THREADS, SEED);

                char cent_name[LINE_LEN];

                snprintf(cent_name, sizeof(char) * LINE_LEN, "centroids%i.txt", i);

                f_cent = fopen(cent_name, "w");

                fprintf(f_cent, "%s\t", para_name_string);
                fprintf(f_cent, "Event\tPopulation\n");

                for (j = 0; j < num_population; j++) {
                    for (p = 0; p < num_dm; p++) {
                        fprintf(f_cent, "%d\t", (int) (population_center[j][p]*4096));
                    }
                    fprintf(f_cent, "%d\t", j + 1);
                    fprintf(f_cent, "%d\n", i + 1);
                }
                fclose(f_cent);

                for (j = 0; j < file_Len; j++) size_c[all_population_ID[j]]++;
                //The new centroids can be created based on the boundary of the polygon. Still a threshold on deciding which clusters to partition is needed

                lastClustering = parentPop+1;

            }

            for (j = 0; j < num_population; j++) //you need to reset all clusters
                filtered_p[j] = 0; //0 = default to keep cluster

            //Filtering based on defined hyper region
            for (j = 0; j < num_population; j++) {
                if ((population_center[j][filtered_d_x[i] - 1] < ((double) filtered_x_low[i] / (double) DAG_BIN)) ||
                    (population_center[j][filtered_d_x[i] - 1] > ((double) filtered_x_upper[i] / (double) DAG_BIN)) ||
                    (population_center[j][filtered_d_y[i] - 1] < ((double) filtered_y_low[i] / (double) DAG_BIN)) ||
                    (population_center[j][filtered_d_y[i] - 1] > ((double) filtered_y_upper[i] / (double) DAG_BIN)))
                    filtered_p[j] = 1; //1 = discard (outside of the defined hyper region)

                if (i < num_rows_2) //there exists the reverse filtering configuration for this population
                {
                    if ((population_center[j][filtered_d_x_2[i] - 1] >=
                         ((double) filtered_x_low_2[i] / (double) DAG_BIN)) &&
                        (population_center[j][filtered_d_x_2[i] - 1] <=
                         ((double) filtered_x_upper_2[i] / (double) DAG_BIN)) &&
                        (population_center[j][filtered_d_y_2[i] - 1] >=
                         ((double) filtered_y_low_2[i] / (double) DAG_BIN)) &&
                        (population_center[j][filtered_d_y_2[i] - 1] <=
                         ((double) filtered_y_upper_2[i] / (double) DAG_BIN)))
                        filtered_p[j] = 1; //1 = discard (inside the defined hyper region)
                }
            }

            for (j = 0; j < file_Len; j++) { //iterate through all events
                pppp = all_population_ID[j];

                if (filtered_p[pppp] > 0) //this cluster population needs to be discarded
                    all_gate_ID[j][i] = 1; //1 = discard current event j based on gate i
                else //this population needs to be kept
                    all_gate_ID[j][i] = 0; //0 = keep
            }

            if (filtered_2nd_pass[parentPop] != 0 && parentSize > 100)
            {

                if (filtered_2nd_pass[parentPop] == 2)
                    num_sub_pop[i] = max_num_pop;
                else if (filtered_2nd_pass[parentPop] == 1)
                    num_sub_pop[i] = start_num_pop;
                else
                    num_sub_pop[i] = max_num_pop;

                //recluster now if not done so earlier
                if(((parentPop+1)>lastClustering)||sub_population_center[parentPop]==0){
                    printf("Current number of k %d\n", num_sub_pop[i]);

                    sub_population_ID[parentPop] = (int *) malloc(sizeof(int) * parentSize);
                    memset(sub_population_ID[parentPop], 0, sizeof(int) * parentSize);

                    //printf("Memory allocation for temp population center...\n");
                    sub_population_center[parentPop] = (double **) malloc(sizeof(double *) * num_sub_pop[i]); //Ivan: based on the size of events in current predetermined cell population, create a data matrix for the cell population
                    memset(sub_population_center[parentPop], 0, sizeof(double *) * num_sub_pop[i]);
                    for(j = 0; j < num_sub_pop[i]; j++) {
                        sub_population_center[parentPop][j] = (double *) malloc(sizeof(double) * num_dm);
                        memset(sub_population_center[parentPop][j], 0, sizeof(double) * num_dm);
                    }

                    printf("Current cell population: %d\tParent cell population: %d\n", i+1, parentPop+1);
                    printf("Parent population size: %d\n", parentSize);

                    printf("Starting recursive_optimize...\n");
                    recursive_optimize(norm_pop_data[parentPop], num_sub_pop[i], TERMS, parentSize, num_dm, sub_population_ID[parentPop], sub_population_center[parentPop],
                          1, NUM_THREADS, SEED);
                    lastClustering = parentPop+1;

                    char cent_name[LINE_LEN];

                    snprintf(cent_name, sizeof(char) * LINE_LEN, "centroids%i.txt", parentPop+1);

                    f_cent = fopen(cent_name, "w");

                    fprintf(f_cent, "%s\t", para_name_string);

                    fprintf(f_cent, "Event\tPopulation\n");

                    for (j = 0; j < num_sub_pop[i]; j++) {
                        for (p = 0; p < num_dm; p++) {
                            fprintf(f_cent, "%d\t", (int) (sub_population_center[parentPop][j][p]*4096));
                        }
                        fprintf(f_cent, "%d\t", j + 1);
                        fprintf(f_cent, "%d\n", parentPop+1);
                    }

                    fclose(f_cent);
                }
                tmp_filtered_p = (int *) malloc(sizeof(int) * num_sub_pop[i]); //for each cluster, whether it is in or out in the data kept
                memset(tmp_filtered_p, 0, sizeof(int) * num_sub_pop[i]);


                printf("Starting filtering...\n");
                for (j = 0; j < num_sub_pop[i]; j++) {
                    if ((sub_population_center[parentPop][j][filtered_d_x[i] - 1] < ((double) filtered_x_low[i] / (double) DAG_BIN)) ||
                        (sub_population_center[parentPop][j][filtered_d_x[i] - 1] > ((double) filtered_x_upper[i] / (double) DAG_BIN)) ||
                        (sub_population_center[parentPop][j][filtered_d_y[i] - 1] < ((double) filtered_y_low[i] / (double) DAG_BIN)) ||
                        (sub_population_center[parentPop][j][filtered_d_y[i] - 1] > ((double) filtered_y_upper[i] / (double) DAG_BIN)))
                            tmp_filtered_p[j] = 1; //1 = discard (outside of the defined hyper region)

                    if (i < num_rows_2) //there exists the reverse filtering configuration for this population
                    {
                      if ((sub_population_center[parentPop][j][filtered_d_x_2[i] - 1] >=
                            ((double) filtered_x_low_2[i] / (double) DAG_BIN)) &&
                            (sub_population_center[parentPop][j][filtered_d_x_2[i] - 1] <=
                            ((double) filtered_x_upper_2[i] / (double) DAG_BIN)) &&
                            (sub_population_center[parentPop][j][filtered_d_y_2[i] - 1] >=
                            ((double) filtered_y_low_2[i] / (double) DAG_BIN)) &&
                            (sub_population_center[parentPop][j][filtered_d_y_2[i] - 1] <=
                            ((double) filtered_y_upper_2[i] / (double) DAG_BIN)))
                                tmp_filtered_p[j] = 1; //1 = discard (inside the defined hyper region)
                    }
                }
                for (j = 0; j < parentSize; j++) { //iterate through all events
                    pppp = sub_population_ID[parentPop][j];

                    if (tmp_filtered_p[pppp] > 0) //this cluster population needs to be discarded
                        all_gate_ID[pop_ID_map[parentPop][j]][i] = 1; //1 = discard current event j based on gate i
                    else //this population needs to be kept
                        all_gate_ID[pop_ID_map[parentPop][j]][i] = 0; //0 = keep
                }
                free(tmp_filtered_p);
            }

            //f_majority is only used in clustering mode, to show which clusters are in or out
            fprintf(f_majority, "------------------------------------------------------------\n");
            fprintf(f_majority, "For the predefined gate number %d:\n", i + 1);
            fprintf(f_majority, "FLOCK_Population_ID\tProportion_in_Sample\tFiltered=1orNot=0\n");
            for (j = 0; j < num_population; j++)
                fprintf(f_majority, "%d\t%.4f\t%d\n", j + 1, (double) size_c[j] / (double) file_Len, filtered_p[j]);
        }

        if (filtered_type[i] == 2) //ratio-based; FSC-A vs FSC-H or SSC-A vs SSC-H for singlets
        {
            for (j = 0; j < file_Len; j++) {
                if (((normalized_data[j][filtered_d_y[i] - 1] / normalized_data[j][filtered_d_x[i] - 1]) <
                     ((double) filtered_y_low[i] / (double) filtered_x_low[i])) ||
                    ((normalized_data[j][filtered_d_y[i] - 1] / normalized_data[j][filtered_d_x[i] - 1]) >
                     ((double) filtered_y_upper[i] / (double) filtered_x_upper[i])))
                    all_gate_ID[j][i] = 1;
                else {
                    if (i < num_rows_2) {
                        if ((normalized_data[j][filtered_d_x_2[i] - 1] >=
                             ((double) filtered_x_low_2[i] / (double) DAG_BIN)) &&
                            (normalized_data[j][filtered_d_x_2[i] - 1] <=
                             ((double) filtered_x_upper_2[i] / (double) DAG_BIN)) &&
                            (normalized_data[j][filtered_d_y_2[i] - 1] >=
                             ((double) filtered_y_low_2[i] / (double) DAG_BIN)) &&
                            (normalized_data[j][filtered_d_y_2[i] - 1] <=
                             ((double) filtered_y_upper_2[i] / (double) DAG_BIN)))
                            all_gate_ID[j][i] = 1;
                        else
                            all_gate_ID[j][i] = 0;
                    } else
                        all_gate_ID[j][i] = 0;
                }
            }
        }

        if (filtered_type[i] == 1) //bisecting
        {
            int tempcount = 0;
            for (j = 0; j < file_Len; j++) {
                if ((normalized_data[j][filtered_d_x[i] - 1] < ((double) filtered_x_low[i] / (double) DAG_BIN)) ||
                    (normalized_data[j][filtered_d_x[i] - 1] > ((double) filtered_x_upper[i] / (double) DAG_BIN)) ||
                    (normalized_data[j][filtered_d_y[i] - 1] < ((double) filtered_y_low[i] / (double) DAG_BIN)) ||
                    (normalized_data[j][filtered_d_y[i] - 1] > ((double) filtered_y_upper[i] / (double) DAG_BIN)))
                    all_gate_ID[j][i] = 1;
                else {
                    if (i < num_rows_2) {
                        if ((normalized_data[j][filtered_d_x_2[i] - 1] >=
                             ((double) filtered_x_low_2[i] / (double) DAG_BIN)) &&
                            (normalized_data[j][filtered_d_x_2[i] - 1] <=
                             ((double) filtered_x_upper_2[i] / (double) DAG_BIN)) &&
                            (normalized_data[j][filtered_d_y_2[i] - 1] >=
                             ((double) filtered_y_low_2[i] / (double) DAG_BIN)) &&
                            (normalized_data[j][filtered_d_y_2[i] - 1] <=
                             ((double) filtered_y_upper_2[i] / (double) DAG_BIN)))
                            all_gate_ID[j][i] = 1;
                        else {
                            all_gate_ID[j][i] = 0;
                            tempcount++;
                        }
                    } else {
                        all_gate_ID[j][i] = 0;
                        tempcount++;
                    }
                }
            }
        }

        size_filtering[i] = numOfEventsInGate(i, file_Len, final_gate_ID, all_gate_ID, filtered_parent); //Ivan: added function to calculate size of current predeterined cell population and assign gating of events
        printf("Final size of pop# %d is %d\n", i+1, size_filtering[i]);

        pop_data[i] = (double **) malloc(sizeof(double *) * size_filtering[i]); //Ivan: based on the size of events in current predetermined cell population, create a data matrix for the cell population
        memset(pop_data[i], 0, sizeof(double *) * size_filtering[i]);
        for(j = 0; j < size_filtering[i]; j++) {
            pop_data[i][j] = (double *) malloc(sizeof(double) * num_dm);
            memset(pop_data[i][j], 0, sizeof(double) * num_dm);
        }
        norm_pop_data[i] = (double **) malloc(sizeof(double *) * size_filtering[i]); //Ivan: based on the size of events in current predetermined cell population, create a data matrix for the cell population
        memset(norm_pop_data[i], 0, sizeof(double *) * size_filtering[i]);
        for(j = 0; j < size_filtering[i]; j++) {
            norm_pop_data[i][j] = (double *) malloc(sizeof(double) * num_dm);
            memset(norm_pop_data[i][j], 0, sizeof(double) * num_dm);
        }

        pop_ID_map[i] = (int *) malloc(sizeof(int) * size_filtering[i]);
        memset(pop_ID_map[i], 0, sizeof(int) * size_filtering[i]);

        //Assign population members' data, as well as a index mapping of the sub population to the original population
        //if (NORM_METHOD == 3 ||NORM_METHOD == 2) {
            assignPopulation(i, file_Len, final_gate_ID, normalized_data, pop_data, pop_ID_map);
        //}else {
            assignPopulation(i, file_Len, final_gate_ID, input_data, pop_data, pop_ID_map);
            tran(pop_data[i], size_filtering[i], num_dm, NORM_METHOD, norm_pop_data[i]);
        //}

        if (filtered_parent[i] == 0) //parent is the whole file
            filtered_percentage[i] = (double) size_filtering[i] * 100.0 / (double) file_Len;
        else {
            if ((filtered_parent[i] > 0) && (size_filtering[filtered_parent[i] - 1] !=
                                             0))  //note that the filtered_parent[x]-1 is the real row ID of the predefined population
                filtered_percentage[i] =
                        (double) size_filtering[i] * 100.0 / (double) size_filtering[filtered_parent[i] - 1];
        }

        fprintf(f_filtered_percentage, "%d\t%.4f\n", i + 1, filtered_percentage[i]);
        fprintf(f_filtered_events, "%d\t%d\n", i + 1, size_filtering[i]);
        ///////////////////////////////////////////////////////////////////////////

        /////////////////////MFI////////////////////////////////////////////////////
        //Compute Population Center of filtered and unfiltered events
        ID2Center_all(input_data, file_Len, num_dm, 2, final_gate_ID, gate_center);

        fprintf(f_filtered_MFI, "%d\t", i + 1);

        for (j = 0; j < num_dm; j++) {
            if (j == num_dm - 1)
                fprintf(f_filtered_MFI, "%.0f\n",
                        gate_center[0][j]); //only the population that is kept will need to have MFI
            else
                fprintf(f_filtered_MFI, "%.0f\t", gate_center[0][j]);
        }

        ///////////////////////////////////////////////////////////////////////////

        //Note that we need MFI from all predefined gates/populations; but this is different from the needs of the visualization purpose
        //that needs to have MFI from both the predefined gate and the rest of the cells. The second part will be output only when the filtered_output[i]==1.

        if (filtered_output[i] == 1)
        {

            ////////////////Initialize directory for each gating population//////////////
            int pop = i+1;
            struct stat st = {0};
            char popdir[32];

            snprintf(popdir, sizeof(char) * 32, "./pop%i", pop);

            if (stat(popdir, &st) == -1) {
                mkdir(popdir, 0700);
            }
            /////////////////////////////////////////////////////////////////////////////

            //////////////////Output filtered events for the selected population/////////
            snprintf(f_name, sizeof(char) * LINE_LEN, "./pop%i/_filtered.txt", pop);
            f_filtered = fopen(f_name, "w");
            fprintf(f_filtered, "%s\n", para_name_string);

            for (j = 0; j < file_Len; j++)
                if (final_gate_ID[j] == 0) {
                    for (p = 0; p < num_dm; p++)
                        if (p == num_dm - 1)
                            fprintf(f_filtered, "%d\n", (int) input_data[j][p]);
                        else
                            fprintf(f_filtered, "%d\t", (int) input_data[j][p]);
                }

            fclose(f_filtered);

            if (filtered_output_finished == 0)
            {
                f_selected_file_name[0]='\0';
                strcpy(f_selected_file_name,"./");
                strcat(f_selected_name, "_selected_filtered.txt");
                strcat(f_selected_file_name, f_selected_name);
                printf("user-selected filtered file name is %s\n",f_selected_file_name);
                f_final_filtered = fopen(f_selected_file_name,"w");
                fprintf(f_final_filtered, "%s\n", para_name_string);

            	for (j = 0; j < file_Len; j++)
                	if (final_gate_ID[j] == 0) {
                    	    for (p = 0; p < num_dm; p++)
                       	 	if (p == num_dm - 1)
                            	   fprintf(f_final_filtered, "%d\n", (int) input_data[j][p]);
                        	else
                            	   fprintf(f_final_filtered, "%d\t", (int) input_data[j][p]);
                    }

            fclose(f_final_filtered);
            filtered_output_finished = 1;
            }


            /////////////////////////////////////////////////////////////////////////////

            //////Output profile and filtered percentage for the selected population/////
            show(input_data, final_gate_ID, file_Len, 2, num_dm, para_name_string, pop);
            /////////////////////////////////////////////////////////////////////////////


            char cid_name[LINE_LEN];
            //char out_name[LINE_LEN];
            char results_name[LINE_LEN];

            snprintf(cid_name, sizeof(char) * LINE_LEN, "./pop%i/population_id.txt", pop);
            //snprintf(out_name, sizeof(char) * LINE_LEN, "coordinates.txt", pop);
            snprintf(results_name, sizeof(char) * LINE_LEN, "./pop%i/flock_results.txt", pop);

            f_cid = fopen(cid_name, "w");
            f_out = fopen("coordinates.txt", "w");
            f_results = fopen(results_name, "w");

            for (j = 0; j < file_Len; j++)
                fprintf(f_cid, "%d\n", final_gate_ID[j] + 1); //start from 1 instead of 0

            fprintf(f_out, "%s\n", para_name_string);
            fprintf(f_results, "%s\tEvent\tPopulation\n", para_name_string);
            for (j = 0; j < file_Len; j++) {
                for (p = 0; p < num_dm; p++) {
                    if (input_data[j][p] < min) {
                        min = (int) input_data[j][p];
                    }
                    if (input_data[j][p] > max) {
                        max = (int) input_data[j][p];
                    }
                    if (p == num_dm - 1) {
                        fprintf(f_out, "%d\n", (int) input_data[j][p]);
                        fprintf(f_results, "%d\t", (int) input_data[j][p]);
                    } else {
                        fprintf(f_out, "%d\t", (int) input_data[j][p]);
                        fprintf(f_results, "%d\t", (int) input_data[j][p]);
                    }
                }

                fprintf(f_results, "%d\t", j + 1);
                fprintf(f_results, "%d\n", final_gate_ID[j] + 1); //all_population_ID[i] changed to all_population_ID[i]+1 to start from 1 instead of 0: April 16, 2009
            }
            char parameters_name[LINE_LEN];
            char properties_name[LINE_LEN];

            snprintf(parameters_name, sizeof(char) * LINE_LEN, "./pop%i/parameters.txt", pop);
            snprintf(properties_name, sizeof(char) * LINE_LEN, "./pop%i/fcs.properties", pop);


            f_parameters = fopen(parameters_name, "w");
            fprintf(f_parameters, "Number_of_Bins\t0\n");
            fprintf(f_parameters, "Density\t0\n");
            fprintf(f_parameters, "Min\t%d\n", min);
            fprintf(f_parameters, "Max\t%d\n", max);
            fclose(f_parameters);

            f_properties = fopen(properties_name, "w");
            fprintf(f_properties, "Bins=0\n");
            fprintf(f_properties, "Density=0\n");
            fprintf(f_properties, "Min=%d\n", min);
            fprintf(f_properties, "Max=%d\n", max);
            fprintf(f_properties, "Populations=2\n");
            fprintf(f_properties, "Events=%d\n", file_Len);
            fprintf(f_properties, "Markers=%d\n", num_dm);
            fclose(f_properties);

            char mfi_name[LINE_LEN];
            char ctr_name[LINE_LEN];

            snprintf(mfi_name, sizeof(char) * LINE_LEN, "./pop%i/MFI.txt", pop);
            snprintf(ctr_name, sizeof(char) * LINE_LEN, "./pop%i/population_center.txt", pop);


            f_mfi = fopen(mfi_name, "w");
            f_ctr = fopen(ctr_name, "w");

            for (j = 0; j < 2; j++) {
                fprintf(f_mfi, "%d\t", j + 1);
                fprintf(f_ctr, "%d\t", j + 1);

                for (p = 0; p < num_dm; p++) {
                    if (p == num_dm - 1) {
                        fprintf(f_mfi, "%.0f\n", gate_center[j][p]);
                        fprintf(f_ctr, "%.0f\n", gate_center[j][p]);
                    } else {
                        fprintf(f_mfi, "%.0f\t", gate_center[j][p]);
                        fprintf(f_ctr, "%.0f\t", gate_center[j][p]);
                    }
                }
            }
            fclose(f_mfi);
            fclose(f_cid);
            fclose(f_ctr);
            fclose(f_out);
            fclose(f_results);

        } //end for if filtered_output[i]==1

    } //end for predefined population i


    f_results = fopen("DAFi_results_all.txt", "w");

    fprintf(f_results, "%s\t", para_name_string);
    for (i = 0; i < num_rows; i++) {
        fprintf(f_results, "pop%d\t", i + 1);
    }
    fprintf(f_results, "Event\tPopulation\n");

    for (j = 0; j < file_Len; j++) {
        for (p = 0; p < num_dm; p++) {
            if (input_data[j][p] < min) {
                min = (int) input_data[j][p];
            }
            if (input_data[j][p] > max) {
                max = (int) input_data[j][p];
            }
            if (p == num_dm - 1) {
                fprintf(f_results, "%d\t", (int) input_data[j][p]);
            } else {
                fprintf(f_results, "%d\t", (int) input_data[j][p]);
            }
        }

        i = 0;
        int highestPop = 0;
        for (i = 0; i < num_rows; i++) {
            t = i;   //note that filtered_parent[0]=0;
            while (t >= 0) //the intersection of multiple gates needs to be calculated; iterates through the parent gates
            {
                if (all_gate_ID[j][t] > 0) {
                    all_gate_ID[j][i] = 1;  //this event has been removed in some parent gate
                    break;
                }
                t = filtered_parent[t] - 1; //please note that the filtered_parent[x] returns 0 to num_rows-1, however, it needs to be changed to -1 to num_rows-2 to be used
            }

            fprintf(f_results, "%d\t", all_gate_ID[j][i]);
            if (all_gate_ID[j][i] == 0) {
                highestPop = i + 1;
            }
        }

        fprintf(f_results, "%d\t", j + 1);
        fprintf(f_results, "%d\n", highestPop); //Ivan 2-5-2017

    }
    fclose(f_results);

    fclose(f_filtered_percentage);
    fclose(f_filtered_events);
    fclose(f_filtered_MFI);
    fclose(f_majority);


    free(size_c);
    free(filtered_p);


    free(filtered_ID);
    free(filtered_d_x);
    free(filtered_d_y);
    free(filtered_x_low);
    free(filtered_x_upper);
    free(filtered_y_low);
    free(filtered_y_upper);
    free(filtered_parent);
    free(filtered_type);
    free(filtered_output);
    free(filtered_2nd_pass);

    free(filtered_ID_2);
    free(filtered_d_x_2);
    free(filtered_d_y_2);
    free(filtered_x_low_2);
    free(filtered_x_upper_2);
    free(filtered_y_low_2);
    free(filtered_y_upper_2);
    free(filtered_parent_2);
    free(filtered_type_2);
    free(filtered_output_2);
    free(filtered_2nd_pass_2);

    free(final_gate_ID);
    free(filtered_percentage);

    for (i = 0; i < file_Len; i++)
        free(all_gate_ID[i]);
    free(all_gate_ID);

    for (i = 0; i < 2; i++)
        free(gate_center[i]);
    free(gate_center);

    for (i = 0; i < num_population; i++)
        free(population_center[i]);
    free(population_center);

    for (i = 0; i < file_Len; i++)
        free(normalized_data[i]);
    free(normalized_data);

    for (i = 0; i < file_Len; i++)
        free(input_data[i]);
    free(input_data);

    free(all_population_ID);

    for (i = 0; i < num_rows; i++){
        free(sub_population_ID[i]);
        if (sub_population_center[i] != 0){
            for (j = 0; j < num_sub_pop[i]; j++) {
                free(sub_population_center[i][j]);
            }
        }
        for (j = 0; j < size_filtering[i]; j++) {
            //free(pop_data[i][j]);
            free(norm_pop_data[i][j]);
        }
        free(pop_data[i]);
        free(norm_pop_data[i]);
        free(sub_population_center[i]);
        free(pop_ID_map[i]);
    }
    free(sub_population_ID);
    free(sub_population_center);
    free(num_sub_pop);
    free(pop_data);
    free(norm_pop_data);
    free(pop_ID_map);
    free(size_filtering);


    ///////////////////////////////////////////////////////////
    printf("Ending time:\t\t\t\t");
    fflush(stdout);
    system("/bin/date");

    return 0;
}
