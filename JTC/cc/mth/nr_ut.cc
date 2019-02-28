/*           (Developped) copy of nr_util0.cpp for unix porting.         */
/*
   A typical use of these routines is :
	float **a;
	a=matrix(1,13,1,9);
		...
	a[3][6] = ...
	   ... + a[2][9]/3.0 ;
	someroutine(a, ... ) ;
	...
	free_matrix(a,1,13,1,9) ;
   note that the whole thing just asks for conversion to C++ classes ! */

#include <stdio.h>
#include <stdlib.h>     // in the Sparcs includes malloc, do not need
                        // ... #include <alloc.h>
#include "nr_ut.h"

char nrerror(const char error_text[])

// Numerical Recipes standard error handler .
{
 fprintf(stderr,"Numerical Recipes run-time error :\n");
 fprintf(stderr,"%s\n", error_text );
 fprintf(stderr,"...now exiting to system.\n") ;
 exit(-1) ;
 return(-1);  // line never reached
}

/* -------------------------------------------------------------------- */
float *vector(int nl,int nh)
/* allocates a float vector with range [nl,nh] */
{
 float *v ;

 v=(float *) malloc((unsigned) (nh-nl+1)*sizeof(float));
 if (!v) nrerror("allocation failure in vector()");
 return v-nl ;
}

/* -------------------------------------------------------------------- */
int *ivector(int nl,int nh)
/* allocates an int vector with range [nl,nh] */
{
 int *v ;

 v=(int *) malloc((unsigned) (nh-nl+1)*sizeof(int)) ;
 if (!v) nrerror("allocation failure in ivector()");
 return v-nl ;
}

/* ------------------------------------------------------------------------ */
double *dvector(int nl,int nh)
/* allocates a double vector with range [nl,nh] */
{
 double *v ;

 v=(double *) malloc((unsigned) (nh-nl+1)*sizeof(double)) ;
 if (!v) nrerror("allocation failure in dvector()");
 return v-nl ;
}

/* ------------------------------------------------------------------------ */
float **matrix(int nrl,int nrh,int ncl,int nch)
/* Allocates a float matrix with range [nrl...nrh][ncl...nch] */
{
 int i;
 float **m ;

 /* Allocate pointers to rows */
 m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*)) ;
 if (!m) nrerror("allocation failure 1 in matrix()") ;
 m -= nrl ;

 /*Allocate rows and set pointers to them */
 for (i=nrl; i<=nrh; i++) {
     m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float)) ;
     if (!m[i]) nrerror("allocation failure 2 in matrix()") ;
     m[i] -= ncl ;
 }
 /* return pointer to array of pointers to rows .*/
 return m ;
}

/* ------------------------------------------------------------------------ */
double **dmatrix(int nrl,int nrh,int ncl,int nch)
/* Allocates a double matrix with range [nrl...nrh][ncl...nch] */
{
 int i;
 double **m ;

 /* Allocate pointers to rows */
 m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*)) ;
 if (!m) nrerror("allocation failure 1 in dmatrix()") ;
 m -= nrl ;

 /*Allocate rows and set pointers to them */
 for (i=nrl; i<=nrh; i++) {
     m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double)) ;
     if (!m[i]) nrerror("allocation failure 2 in dmatrix()") ;
     m[i] -= ncl ;
 }
 /* return pointer to array of pointers to rows .*/
 return m ;
}

/* ----------------------------------------------------------------------- */
int **imatrix(int nrl,int nrh,int ncl,int nch)
/* Allocates an int matrix with range [nrl...nrh][ncl...nch] */
{
 int i, **m ;

 /* Allocate pointers to rows */
 m=(int **) malloc((unsigned) (nrh-nrl+1)*sizeof(int*)) ;
 if (!m) nrerror("allocation failure 1 in imatrix()") ;
 m -= nrl ;

 /*Allocate rows and set pointers to them */
 for (i=nrl; i<=nrh; i++) {
     m[i]=(int *) malloc((unsigned) (nch-ncl+1)*sizeof(int)) ;
     if (!m[i]) nrerror("allocation failure 2 in dmatrix()") ;
     m[i] -= ncl ;
 }
 /* return pointer to array of pointers to rows .*/
 return m ;
}

/* ------------------------------------------------------------------------ */
float **submatrix(float **a,int oldrl,int oldrh,int oldcl,int oldch,
			    int newrl,int newcl)
/* Returns a submatrix with range
      [newrl .. newrl+(oldrh-oldrl)][newcl..newcl_(oldch-oldcl)]
   pointing to the existing matrix range a[oldrl..oldrh][oldcl..oldch] . */
{
 int i,j ;
 float **m ;

 /* Allocate pointers to rows . */
 m=(float **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(float*)) ;
 if (!m) nrerror("allocation failure in submatrix()");
 m -= newrl ;

 /* Set pointers to rows . */
 for(i=oldrl,j=newrl ; i<=oldrh; i++,j++) m[j]=a[i]+oldcl-newcl ;

 /* Return pointer to array of pointers to rows. */
 return m ;
}

/* ----------------------------------------------------------------------- */
void free_vector(float *v,int nl,int nh)
/*  Frees a float vector allocated by vector() */
{
 free((char*) (v+nl)) ;
}
  
/* ----------------------------------------------------------------------- */
void free_ivector(int *v,int nl,int nh)
/*  Frees an int vector allocated by ivector() */
{
 free((char*) (v+nl)) ;
}
  
/* ----------------------------------------------------------------------- */
void free_dvector(double *v,int nl,int nh)
/*  Frees a double vector allocated by dvector() */
{
 free((char*) (v+nl)) ;
}
  
/* ----------------------------------------------------------------------- */
void free_matrix(float **m,int nrl,int nrh,int ncl,int nch)
/* Frees a matrix allocated with matrix() */
{
 int i;

 for (i=nrh; i>=nrl; i--) free((char*) (m[i]+ncl)) ;
 free((char*) (m+nrl)) ;
}

/* ----------------------------------------------------------------------- */
void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
/* Frees a double matrix allocated with dmatrix() */
{
 int i;

 for (i=nrh; i>=nrl; i--) free((char*) (m[i]+ncl)) ;
 free((char*) (m+nrl)) ;
}


/* ----------------------------------------------------------------------- */
void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch)
/* Frees an int matrix allocated with imatrix() */
{
 int i;

 for (i=nrh; i>=nrl; i--) free((char*) (m[i]+ncl)) ;
 free((char*) (m+nrl)) ;
}

/* ----------------------------------------------------------------------- */
void free_submatrix(float **b,int nrl,int nrh,int ncl,int nch)
/* Frees submatrix allocated by submatrix() */
{
 free((char*) (b+nrl)) ;
}

/* ------------------------------------------------------------------------ */
/* I couldn't be bothered to type in the
   "float **convert_matrix(a,nrl,nrh,ncl,nch)"   */










































