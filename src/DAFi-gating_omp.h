#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h> //added by Ivan
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <assert.h>

#define DEBUG 0
#define LINE_LEN 2048
#define FILE_NAME_LEN 128

#define MAX_VALUE 1000000000

#define NORM_METHOD 3 //2 if z-score; 0 if no normalization; 1 if min-max based on column range; 3 if min-max based on value_range (i.e., 4095)

#define DEFAULT_NUM_POP 200

#define DAG_BIN 200
#define VALUE_RANGE 4095
#define TERMS 100
//#define PROP_FILTER_T 0.3
