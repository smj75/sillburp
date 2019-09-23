
/* 
	Definition and parameter file for sillburp and lipburp.
	Statements common to both programs go in here.
 */

#ifndef _BURP_H
#define _BURP_H

#include <cfloat>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iomanip>
using namespace std;
#include <memory>
#include <string>
#include <vector>

#define CNULL ((char *)NULL)
#define DNULL ((double *)NULL)
#define FNULL ((float *)NULL)
#define INULL ((int *)NULL)
#define VNULL ((void *)NULL)

#define NR_END 1
#define FREE_ARG char*

#define GRAVITY (9.81)	/* in m/s2 */
#define PI (3.141592654)
#define PI_SQ (PI * PI)
#define TWOPI (PI * 2.0)
#define SQRT_TWO (sqrt(2.0))
#define SQRT_PI (sqrt(PI))
#define SQRT_TWOPI (sqrt(TWOPI))
#define SQRT_PIBYTWO (sqrt(PI/2))
#define RADIAN (PI / 180.0)
#define GMT_FLAG ">"
#define S_PER_YR (60 * 60 * 24 * 365.25)
#define S_PER_DAY (60 * 60 * 24)

/* For reading in from files */

#define MAX_FIELDS 		20
#define MAX_RECORDS 	1048576
#define GMT_IS_NAN		0	/* Returned by GMT_scanf routines when read fails */
#define GMT_IS_FLOAT	1	/* Generic (double) data type, no special format */
#define IO_MISMATCH		2
#define IO_EOF				4
#define is_dnan(x) ((x) != (x))


	static int io_n_bad_records;
	static int io_n_clean_rec;
	static int io_rec_no;
	static int io_status;
	static double *line_data; // POINTER TO STOREAGE FOR READING NUMBERS ON 1 LINE 
	static char io_current_record[BUFSIZ];

/* Symbols to indicate comment lines in input files */

const char comment_symbols[32] = {'C','c','#','%','/','>'};
const int n_comment_symbols = strlen(comment_symbols);
const bool skip_comment_symbols = true;

/* Subroutine and function definitions */

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

double **dmatrix(
	long nrl, 
	long nrh, 
	long ncl, 
	long nch);

double ***d3tensor(
	long nrl, 
	long nrh, 
	long ncl, 
	long nch, 
	long ndl, 
	long ndh);

double *dvector(
	long nl,
	long nh );

void free_dmatrix(
	double **m,
	long nrl,
	long nrh,
	long ncl,
	long nch	);

void free_dvector(
	double *v,
	long nl,
	long nh);

void free_ivector(
	int *v,
	long nl,
	long nh );

void get_user_parameters ( 
	int argc, 
	char **argv );

int	GMT_scanf_float (
	char *s, 
	double *val);

int GMT_strtok (
	const char *string, 
	const char *sep, 
	int *pos, 
	char *token);

int input_ascii_line (
	FILE *fp, 
	int *n, 
	double **ptr,
	int *n_rows	);

int *ivector(
	long nl,
	long nh );

int read_data_file (
	FILE *fp, 
	int *n_cols,
	int *n_rows,
	double **data, 
	double *time_end, 
	double *time_step,
	int n_header_rec	);

#endif /* _BURP_H */

