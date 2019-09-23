/* 
	Library of subroutines and functions common to sillburp and lipburp
 */

#include <burp.h>

/* Common variable for reading in files */


double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;
  m=(double **) malloc((unsigned int)((nrow+NR_END)*sizeof(double*)));
  if (!m) {
		cerr << "allocation failure 1 in dmatrix()" << "\n";
		exit (0);
	}
  m += NR_END;
  m -= nrl;
  m[nrl]=(double *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(double)));
  if (!m[nrl]) {
		cerr << "allocation failure 2 in dmatrix()" << "\n";
		exit (0);
	}
  m[nrl] += NR_END;
	m[nrl] -= ncl;
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  return m;
}


double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* Allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	double ***t;
	/* allocate pointers to pointers to rows */
	t=(double ***) malloc((size_t)((nrow+NR_END)*sizeof(double**)));
	if (!t) {
		cerr << " Allocation failure 1 in d3tensor()\n";
		exit(0);
	}
	t += NR_END;
	t -= nrl;
	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(double **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double*)));
	if (!t[nrl])  {
		cerr << " Allocation failure 2 in d3tensor()\n";
		exit(0);
	}
	t[nrl] += NR_END;
	t[nrl] -= ncl;
	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(double *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(double)));
	if (!t[nrl][ncl])  {
		cerr << " Allocation failure 3 in d3tensor()\n";
		exit(0);
	}
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;
	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}
	/* return pointer to array of pointers to rows */
	return t;
}


double *dvector(long nl, long nh)
{
	double *v;
	v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) cerr << "allocation failure in dvector()\n";
	return v-nl+NR_END;
}


void free_dvector(
	double *v,
	long nl,
	long nh )
{
        free((FREE_ARG) (v+nl-NR_END));
}



int read_data_file (FILE *fp, int *n_cols, int *n_rows, double **data, double *time_end, double *time_step, int n_header_rec)
{	

/* DATA INPUT ROUTINE BASED ON GMT ROUTINE
 * load_data_and_check_extreme FROM filter1d */
	
	int i, n_fields, n_expected_fields, n_extra_cols;
	double last_time, new_time;
	double *in;
	int t_col = 0; // Column for variable that must increase monotonically, in this case time

/* Temporary storeage for input line */

	line_data = dvector(0,MAX_FIELDS);

/* ASSUME FIRST LINE CONTAINS COPY OF COMMAND LINE THAT GENERATED THE FILE */

		char line[BUFSIZ];
		for ( int i=1; i<=n_header_rec; i++) fgets(line, BUFSIZ, fp);

/* Read in line and divide into fields */
  
	(*n_rows) = 1;
	n_expected_fields = (*n_cols);
	while ((n_fields = input_ascii_line (fp, &n_expected_fields, &in, n_rows)) >= 0 && !(io_status & IO_EOF)) { 

//    if ( (*n_rows) == 1 && prog_name==in[0] ) {
//			cerr << " File produced using command line" << io_current_record << "\n";
//			continue;
//		}
	
/* USE FIRST LINE OF DATA TO FIND OUT HOW MANY COLUMNS TO BE READ IN.
 * MAKE ARRAY OF APPROPRIATE SIZE */
	
		if ( (*n_rows) == 1 ) {
			n_expected_fields = n_fields;
			n_extra_cols = n_fields - (*n_cols);
			if ( n_extra_cols==0 ) {
				cerr << " First data record has " << n_fields << " fields, as expected\n";
			} else {
				cerr << " First data record has " << n_fields << " fields but " << (*n_cols) << " fields expected\n";
				cerr << " The " << n_extra_cols << " extra field(s) will be written out unchanged\n";
			}

/* SUBSEQUENT LINES MUST HAVE SAME NUMBER OF COLUMNS AS FIRST LINE */

		} else {
			if ( io_status & IO_MISMATCH ) {
					fprintf (stderr, " Mismatch between actual (%d) and expected (%d) fields near line %d\n", n_fields,  n_expected_fields, *n_rows);
					return (-1);
			}
		}
		for (i = 0; i < n_fields; i++) data[i][*n_rows] = in[i];

/* CHECK THAT time VALUE IS READABLE */

		if (!is_dnan(data[t_col][*n_rows]) ) {
			new_time = data[t_col][*n_rows];
			(*n_rows)++;
		}
	
/* MORE MEMORY REQUIRED FOR DATA ARRAY */
				
		if ( (*n_rows) == MAX_RECORDS ) {
			cerr << " Too many data records (lines): " << MAX_RECORDS << "\n";
			cerr << " Increase MAX_RECORDS in burp.h and recompile\n";
			return (-1);
		}

	} /* LOOP BACK TO READ IN NEXT LINE */

/* END-OF-FILE MARKER FOUND */

	fclose (fp);
	if ( *n_rows > 0 ) (*n_rows)--;
	cerr << " Read in " << *n_rows << " data records\n";

/* Time parameters */

	(*time_end) = data[0][*n_rows];
	(*time_step) = data[0][2] - data[0][1];

	return (0);
}



int input_ascii_line (FILE *fp, int *n, double **ptr, int *n_rows)
{

/* SIMPLIFIED FROM GMT ROUTINE GMT_ascii_input.
 * SKIP BLANK LINES AND C-SHELL COMMENT LINES WHICH START WITH #.
 * FIELDS MAY BE SEPARATED BY SPACES, TABS OR COMMAS.
 * RETURNS NUMBER OF ITEMS READ, OR -1 FOR EOF.
 * IF *n IS PASSED AS BUFSIZ IT IS RESET TO THE ACTUAL NUMBER OF FIELDS.
 */
/*
 * int input_ascii_line       MODIFIED FROM gmt_io.c
 * int GMT_strtok             COPIED DIRECTLY FROM gmt_support.c
 * int GMT_scanf_float        COPIED DIRECTLY FROM gmt_support.c
 */

	bool bad_record;
	bool comment_line = false;
	bool done = false;
	bool new_line;
	char line[BUFSIZ];
	char *p;
	char token[BUFSIZ];
	int col_no;
	int i;
	int len;
	int pos;
	int n_convert;
	double val;


	while ( !done ) {	/* BECOMES TRUE WHEN RECORD READ SUCCESSFULLY */
		
/* READ UNTIL WE GET A NON-BLANK, NON-COMMENT RECORD, OR REACH EOF */

		do {
			p = fgets(line, BUFSIZ, fp);
			new_line = (line[0] == '\n');
			if ( skip_comment_symbols ) {
				comment_line = false;
				for ( i=0; i<n_comment_symbols; i++ ) if (!comment_line) comment_line = (line[0]==comment_symbols[i]);
			}
		} while ( p && (new_line||comment_line) ) ;
//		while ((p = fgets (line, BUFSIZ, fp)) && (line[0] == '\n' || line[0] == '#' )) io_rec_no++;

// IF END OF FILE REACHED

		if ( !p ) {
			io_status = IO_EOF;
			if ( io_n_bad_records ) {	/* REPORT SUMMARY AND RESET */
				fprintf (stderr, " This file had %d records with invalid x and/or y values\n", io_n_bad_records);
				io_n_bad_records = io_n_clean_rec = 0;
			}
			return (-1);
		}

/* CHECK FOR DOS FORMAT AND CHOP OFF  
 * TRAILING WHITE SPACE AND COMMAS  */

		len = strlen (line);
		if (len >= (BUFSIZ-1)) {
			fprintf (stderr, " This file appears to be in DOS format - reformat with dos2unix\n");
			exit (EXIT_FAILURE);
		}
		for ( i=len-1; i>=0 && strchr(" \t,\r\n",(int)line[i]); i--);
		line[++i] = '\n';	
		line[++i] = '\0';	
		
/* Now have clean C string with \n\0 at end */

		bad_record = false;
		strcpy (io_current_record, line);
		line[i-1] = '\0';		/* CHOP OFF NEW LINE */
		col_no = pos = 0;
		while (!bad_record && col_no < *n && (GMT_strtok (line, " \t,", &pos, token))) {	/* EACH FIELD */
			if ((n_convert = GMT_scanf_float (token, &val)) == GMT_IS_NAN) {	/* GOT NAN OR FAILED TO DECODE */
				bad_record = true;
			} else {	/* SUCCESSFUL DECODE */
				line_data[col_no] = val;
			}
			col_no++;		/* NEXT FIELD */
			if (col_no >= MAX_FIELDS) {
				exit(0);
			}
		}

		if ( bad_record ) {
			io_n_bad_records++;
			if ( io_n_bad_records == 1 ) {	/* REPORT 1ST OCCURENCE */
				fprintf (stderr, " Encountered first invalid record near/at line %d\n", (*n_rows)+1);
			}
		} else done = true;

	}

	*ptr = line_data;
	io_status = (col_no == *n || *n == BUFSIZ) ? 0 : IO_MISMATCH;
	if (*n == BUFSIZ) *n = col_no;

	return (col_no);
}

int GMT_strtok (const char *string, const char *sep, int *pos, char *token)
{
	/* Reentrant replacement for strtok that uses no static variables.
	 * Breaks string into tokens separated by one of more separator
	 * characters (in sep).  Set *pos to 0 before first call.  Unlike
	 * strtok, always pass the original string as first argument.
	 * Returns 1 if it finds a token and 0 if no more tokens left.
	 * pos is updated and token is returned.  char *token must point
	 * to memory of length >= strlen (string).
	 */

	int i, j, string_len;
	
	string_len = strlen (string);
	
	/* Wind up *pos to first non-separating character: */
	while (string[*pos] && strchr (sep, (int)string[*pos])) (*pos)++;

	token[0] = '\0';	/* Initialize token to NULL in case we are at end */
	
	if (*pos >= string_len || string_len == 0) return 0;	/* Got NULL string or no more string left to search */

	/* Search for next non-separating character */
	for (i = *pos; string[i] && !strchr (sep, (int)string[i]); i++);
	
	/* Copy token */
	j = i - *pos;
	strncpy (token, &string[*pos], j);
	token[j] = 0;
	
	/* Wind up *pos to next non-separating character */
	while (string[i] && strchr (sep, (int)string[i])) i++;
	*pos = i;
	
	return 1;
}

int	GMT_scanf_float (char *s, double *val)
{
	/* Try to decode a value from s and store
	in val.  s should not have any special format
	(neither geographical, with suffixes or
	separating colons, nor calendar nor clock).
	However, D and d are permitted to map to e
	if this would result in a success.  This
	allows Fortran Double Precision to be readable.

	On success, return GMT_IS_FLOAT and store val.
	On failure, return GMT_IS_NAN and do not touch val.
	*/

	char	scopy[BUFSIZ], *p;
	double	x;
	int	j,k;

	x = strtod (s, &p);
	if (p[0] == 0) {
		/* Success (non-Fortran).  */
		*val = x;
		return (GMT_IS_FLOAT);
	}
	if (p[0] != 'D' && p[0] != 'd') return (GMT_IS_NAN);
	k = strlen(p);
	if (k == 1) {
		/* A string ending in e would be invalid  */
		return (GMT_IS_NAN);
	}
	/* Make a copy of s in scopy, mapping the d or D to an e:  */
	j = strlen(s);
	if (j > BUFSIZ) return (GMT_IS_NAN);
	j -= k;
	strncpy (scopy, s, (size_t)j );
	scopy[j] = 'e';
	strcpy (&scopy[j+1], &p[1]);
	x = strtod(scopy, &p);
	if (p[0] != 0) return (GMT_IS_NAN);
	*val = x;
	return (GMT_IS_FLOAT);
}

