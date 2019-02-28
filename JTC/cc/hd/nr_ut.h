/* header for routines "numerical recipes" to do with the alloc-

   cation and deallocation of memory to matrices - and some others.

   A typical use of these routines is :

	float **a;
	a=matrix(1,13,1,9);
		...
	a[3][6] = ...
	   ... + a[2][9]/3.0 ;
	someroutine(a, ... ) ;
	...
	free_matrix(a,1,13,1,9) ;

   This copy for Unix porting. 16 Feb. '98 on.   */

char nrerror(const char *error_text);


float  *vector(int nl,int nh);
int    *ivector(int nl,int nh);
double *dvector(int nl,int nh);



float  **matrix(int nrl,int nrh,int ncl,int nch);
double **dmatrix(int nrl,int nrh,int ncl,int nch);
int    **imatrix(int nrl,int nrh,int ncl,int nch);


float **submatrix(int a,int oldrl,int oldrh,int oldcl,int oldch,
		  int newrl,int newcl);



void free_vector(float *v,int nl,int nh);
void free_ivector(int *v,int nl,int nh);
void free_dvector(int *v,int nl,int nh);
void free_matrix(float **m,int nrl,int nrh,int ncl,int nch);
void free_dmatrix(double * * m, int nrl,int nrh,int ncl,int nch);
void free_imatrix(int    * * m, int nrl,int nrh,int ncl,int nch);

void free_submatrix(float * * b,int nrl,int nrh,int ncl,int nch);


/* I couldn't be bothered to type in the
   "float **convert_matrix(a,nrl,nrh,ncl,nch)"   */
