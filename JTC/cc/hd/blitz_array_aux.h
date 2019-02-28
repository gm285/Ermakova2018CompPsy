/* ******  Auxiliary functions to do things with blitz++ arrays  *********

     1.x - general purpose C functions
     2.x - Interpolation classes.
*/ 

#ifndef BLITZ_ARRY_AUX_H
#define BLITZ_ARRY_AUX_H 

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>    // input & output streams
#include <iomanip>    // i-o  parametrized manipulators
#include <sstream>    // in-memory string manipulations for i/o, new standard
#include <string>
#include <cstdlib>
#include <time.h>
using namespace std; 

#include <blitz/array.h>
#include <gsl/gsl_statistics_double.h>

BZ_USING_NAMESPACE(blitz)
using namespace blitz::tensor;

/* 1.x general purpose C functions ------------------------------------------ */

// Sorting & related
// output 2-D arrays to disk in csv format, e.g. for gnuplot or spreadsheet. mode=1
// opens w. ios::out, otherwise opens w. ios::app. Useful in loops where first time
// round may want to open w. ios::out, then append.
void write2DArray_csv(const char *fname, const Array<double, 2>& arr, int mode=1, char sepC = '\t') ;
void write2DArray_csv(const char *fname, const Array<float, 2>& arr,  int mode=1, char sepC = '\t') ;
void write2DArray_csv(ofstream &output, const Array<float, 2>& arr, char sepC = '\t') ;
void write2DArray_csv(const char *fname, const Array<int, 2>& arr,    int mode=1, char sepC = '\t') ;
// very similar, for Python format (first add some preliminaries,
// and at end close it properly) :
void write2DArPrel_py(const char *fname, const char *arname, int mode=1);
void write2DArEnd_py(const char *fname, const char *dtype);
// write just the array, without its name / prelim or end bits:
void write2DArray_py(const char *fname, const Array<int, 2>& arr,    int mode=1, char sepC = ',') ;
// Together with array name and dtype='i' :
void write2DArray_py(const char *fname, const char *arname, const Array<int, 2>& arr, int mode=1, char sepC = ',') ;


// Take a 2D double array arr(ind1st ... ind1fin, ind2st ... ind2fin) and values for
// the underlying coordinates val1st, val1fin, val2st, val2fin, and write a csv
// array suitable for plotting with gnuplot splot as a 3-D surface plot. It is assumed
// that the underlying coordinate values are equally spaced.
void write_splot_csv(  const char *fname, const Array<double, 2>& arr
		     , double val1st, double val1fin, double val2st, double val2fin);
void write_splot_csv(  const char *fname, const Array<float, 2>& arr
		     , double val1st, double val1fin, double val2st, double val2fin);

// 1b simple stats etc 
double mean1D(const Array<double, 1>& arr);
double mean1D(const Array<float, 1>& arr);

double  var1D(const Array<float, 1>& arr, double *mean = NULL);

// Simple Simpson for integrn. REM arr needs even no of intervals,
// i.e. odd no. of interval endpoints.
double Simpson1DInt(const Array<double, 1>& arr, double start, double end);
double Simpson1DInt(const Array<float, 1>& arr, double start, double end);

// Treat the array as containing a pdf, and calculate its basic stats:
double pdfSimpsMean1D(const Array<double, 1>& arr, double start, double end);
double pdfSimpsMean1D(const Array<float, 1>& arr, double start, double end);
// The following returns the mean if the last entry isn't NULL
double pdfSimpsVar1D(const Array<double, 1>& arr, double start, double end
                     , double *mean = NULL);
double pdfSimpsVar1D(const Array<float, 1>& arr, double start, double end
                     , double *mean = NULL);

// ----------------------------------------------------------------------------
// 2.x  - Interpolation -------------------------------------------------------

// 2.1 
class Interpolation1D
{
 protected:
  double avoid_div_by_0_a ;  // a tiny additive constant ? 8.e-307
  double avoid_div_by_0_m ;  // a multiplic. constant just > 1
  // double NaNdouble ;          Carries 'Not a Number' value for error work.
  int n; // determines order of the approximation function. Default n=4 ...
         // e.g. order of approximating polynomial would be n-1.
  Array<double, 1> xa;
  Array<double, 1> ya;  // used in PolynInt etc, as is n above.
  Array<double, 1> *X, *Y; // pointers to ext. arrays w. data -
                           // 1-D for here.

  // if input array is equally spaced, much faster to find where we are 
  // using the following :
  double xStart, xStep, xEnd;
  int NDat;  // Number of input data (Y) points.
  int xType; // regular -> 0; irregular (as per xa) -> 1

 public:
  Interpolation1D();

  int orderPar(){ return n ; };
  void setOrderPar(int newOrdPar);
  // sometimes we may need to set ya values directly. Note no checking ! :
  void set_ya_el(double new_ya, int ind){ ya(ind)=new_ya ; };
  void setRegInp_x(double xInStart, double xInStep) {
    xType=0; xStart=xInStart;  xStep=xInStep; } ;
  void setInp_y(Array<double, 1> *yIn ){ Y=yIn;  NDat=Y->extent(firstDim); };
  void setRegInp(double xInStart, double xInStep,  Array<double, 1> *yIn);
  void setGenInp(Array<double, 1> *xIn,  Array<double, 1> *yIn);
  // The following prepares xa and ya for use with PolynInt... etc.
  int prepXbox(double x);   // also return first index of Y-window
  void prepXYbox(double x); // as per prepXbox, but also fill in ya.
  double xa_val(int ind){ return xa(ind); };
  double X_val(int ind){ return xStart+(ind-1)*xStep; } ;

  //    *********** For Accuracy (cv. p. 88 NR in C) ******************

  // much like NR polint, except main output y is returned directly 
  // and the error only returned if needed:
  double PolynInt0(double x, double *err=NULL);
  // similarly, like NR ratint (Rational function Bulirsch-Stoer Interpolation) :
  double RationInt0(double x, double *err=NULL);
  // and one that tries Burlish-Stoer and falls back to Polynomial
  // if the Rational one really fails (esp. has pole) :
  double RationInt_w_Poly_rescue0(double x, double *err=NULL);

  // Full thing, using full input arrays -  both for x and y,  
  // or using equally spaced abscissa points, acc. to curr. val. of xType :
  double PolynInt(double x, double *err=NULL);
  double RationInt(double x, double *err=NULL);

};
// ----------------------------------------------------------------------------
class Interpolation2D : public Interpolation1D
{
 protected:
  Array<double, 1> Y1,Y2;   // Not very clever, hopefully will do.

  Interpolation1D intX1;    // For the two indep. variables
  Interpolation1D intX2;

  Array<double, 2> y2a; // 2-D cell with dep. var.
  Array<double, 2> *Y2DP;  // 2-D ext. data w. dep. var.
  
  // auxiliaries for regular arrays:
  int NDatTot;   // NDat1*NDat2

 public: 
  Interpolation2D();
  double extY(int r, int c){ return (*Y2DP)(r,c); };
  void setOrderPars(int newOrd1, int newOrd2);
  void setRegInp2D_x(  double xIn1Start, double xIn1Step
		      ,double xIn2Start, double xIn2Step );
  void setInp2D_y(Array<double, 2> *yIn2D) ;
  void setRegInp2D( double xIn1Start, double xIn1Step
                   ,double xIn2Start, double xIn2Step, Array<double, 2> *yIn2D);
  void prepX1X2Ybox(double x1, double x2);
  double X1_val(int ind){ return intX1.X_val(ind); };
  double X2_val(int ind){ return intX2.X_val(ind); };

   
   //    *********** For Accuracy (cv. p. 105 NR in C) ******************

  // Given  intX1.xa[1..m] and intX1.xa[1..n] as independent variables, and a 
  // matrix of function values y2a[1..m][1..n], tabulated at the grid
  // points defined by  intX1.xa and  intX2.xa; and given values x1 and x2 of the 
  // independent variables, find  interpolated function value , and an
  // accuracy indication err (based only on the interpolation in the X1 direction,
  // however !).
  //  1a. Based on rational function interpolation :
  double RationInt(double x1, double x2, double *err=NULL);
  //  1b.    same, but using externally given ext. provided 2-D ext. data 
  //         w. dep. var. directly - MUST BE CONSISTENT W. Y1, Y2 ETC !!
  double RationInt(double x1, double x2, const Array<double, 2>& y2dp, double *err=NULL);
  //  2   Based on polynomial interpol:
  //  2b.    same, but using externally given ext. provided 2-D ext. data 
  //         w. dep. var. directly - MUST BE CONSISTENT W. Y1, Y2 ETC !!
  double PolynInt(double x1, double x2, const Array<double, 2>& y2dp, double *err=NULL);
  //  3   Based on rational fn, but with 'rescue' if the rational interpol.
  //      fails, esp. if it hits pole
  //  3a. - here the 'rescue' is based on Polynomial fn. interpoln:
  //      Uses externally given ext. 2-D ext. CONSISTENT W. Y1, Y2 ETC as above !!
  double RationInt_w_Poly_rescue(  double x1, double x2
                                 , const Array<double, 2>& y2dp, double *err=NULL);



};
// ----------------------------------------------------------------------------
// dont't forget the }; at the end of classes !

#endif    // BLITZ_ARRY_AUX_H 

//                     End of file

