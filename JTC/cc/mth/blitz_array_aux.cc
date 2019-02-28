// methods for blitz_array_aux.h
 
#include "blitz_array_aux.h"
#include "maths.h"
#include "utility.h"

// 1.x general purpose c functions
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
void write2DArray_csv(const char *fname, const Array<double, 2> &arr, int mode,
                      char sepC) {
  ofstream outp;
  if (mode == 1)
    outp.open(fname, ios::out);
  else
    outp.open(fname, ios::app);

  for (int i1 = arr.lbound(firstDim); i1 <= arr.ubound(firstDim); i1++) {
    for (int i2 = arr.lbound(secondDim); i2 <= arr.ubound(secondDim); i2++) {
      outp << arr(i1, i2) << sepC;
    }
    outp << endl;
  }

  outp.close();
}

// Sorry - I've just copy-pasted the following two:
// ----------------------------------------------------------------------------
void  write2DArray_csv(const char *fname, const Array<float, 2>& arr, int mode, char sepC)
{
  // Error code goes here

  ofstream outp ;  
  if (mode==1)  outp.open(fname, ios::out);
  else outp.open(fname, ios::app);

  write2DArray_csv(outp, arr, sepC);
  outp.close();
}

void write2DArray_csv(ofstream &outp, const Array<float, 2> &arr, char sepC) {
  for (int i1=arr.lbound(firstDim); i1<=arr.ubound(firstDim); i1++) {
     for  (int i2=arr.lbound(secondDim); i2<=arr.ubound(secondDim); i2++) 
       outp<<arr(i1,i2)<<sepC  ;
     outp << endl ;
  }
}

// ----------------------------------------------------------------------------
void  write2DArPrel_py(const char *fname, const char *arname, int mode)
{
  ofstream outp ;  
  if (mode==1)  outp.open(fname, ios::out);
  else outp.open(fname, ios::app);

  outp<<"import numpy as np\nfrom numpy import *\n\n";
  outp << arname <<" = array(";
}
// ----------------------------------------------------------------------------
void write2DArray_py(  const char *fname, const char *arname, const Array<int, 2>& arr
                     , int mode, char sepC)
{
  write2DArPrel_py(fname, arname,mode);
  write2DArray_py(fname,arr,2); // open again in append mode !
  write2DArEnd_py(fname,"i");
}
// ----------------------------------------------------------------------------
void  write2DArEnd_py(const char *fname, const char *dtype)
{
  ofstream outp ;  
  // only append - we dont' put end-bits at the start of a file !
  outp.open(fname, ios::app);

  outp <<", dtype='"<<dtype<<"' )\n";
}
// ----------------------------------------------------------------------------
void  write2DArray_py(const char *fname, const Array<int, 2>& arr, int mode, char sepC)
{
  // Error code goes here

  ofstream outp ;  
  if (mode==1)  outp.open(fname, ios::out);
  else outp.open(fname, ios::app);
  int i1, i2;

  outp<<'['; // whole array opening
  for ( i1=arr.lbound(firstDim); i1<=arr.ubound(firstDim); i1++) {
    outp<<'['; // line opening
    for  ( i2=arr.lbound(secondDim); i2<arr.ubound(secondDim); i2++) 
      outp<<arr(i1,i2)<<sepC  ;  // for all except last element in row !
    outp<<arr(i1,arr.ubound(secondDim))<<']'; // for last element in row only.
    if (i1<arr.ubound(firstDim)) // if not last row ...
      outp<<",\n"; // ... string to prepare for next row.
    else // for last row
      outp<<']'; // whole array closing.
  }

  outp.close();
}
// ----------------------------------------------------------------------------
void  write2DArray_csv(const char *fname, const Array<int, 2>& arr, int mode, char sepC)
{
  // Error code goes here

  ofstream outp ;  
  if (mode==1)  outp.open(fname, ios::out);
  else outp.open(fname, ios::app);

  int i1, i2;
  for ( i1=arr.lbound(firstDim); i1<=arr.ubound(firstDim); i1++) {
     for  ( i2=arr.lbound(secondDim); i2<=arr.ubound(secondDim); i2++) 
       outp<<arr(i1,i2)<<sepC  ;
     outp << endl ;
  }

  outp.close();
}// ----------------------------------------------------------------------------
void write_splot_csv(  const char *fname, const Array<float, 2>& arr
		     , double val1st, double val1fin, double val2st, double val2fin)
{
  int i1, i2, n1, n2;    
  n1=arr.extent(firstDim);   n2=arr.extent(secondDim);
  if ((n1<2) || (n2<2) || (val1st==val1fin) || (val2st==val2fin)) {
    cerr<<"\n arr.extent(firstDim)="<<n1<<" arr.extent(secondDim)="<<n2<<'\n'; 
    err_msg(" \nwrite_splot_csv: invalid extents or ranges");
  }
  char   sepC='\t'; 
  double var1, var2 ;                var1=val1st;       var2=val2st; 
  double step1, step2;       
  step1= (val1fin-val1st)/(n1-1);    step2= (val2fin-val2st)/(n2-1);  

  ofstream outp ;  outp.open(fname, ios::out);

  //  for ( i1=arr.lbound(firstDim); i1<=arr.ubound(firstDim); i1++) {
  //   for  ( i2=arr.lbound(secondDim); i2<=arr.ubound(secondDim); i2++) 
  for ( i1=1; i1<arr.lbound(firstDim)+n1; i1++) {
    for  ( i2=1; i2<arr.lbound(secondDim)+n2; i2++) 
       outp<<val1st+(i1-1)*step1<<sepC<<val2st+(i2-1)*step2<<sepC<<arr(i1,i2)<<endl  ;
     outp << endl ;
  }

  outp.close();
}
// ----------------------------------------------------------------------------
void write_splot_csv(  const char *fname, const Array<double, 2>& arr
		     , double val1st, double val1fin, double val2st, double val2fin)
{
  int i1, i2, n1, n2;    
  n1=arr.extent(firstDim);   n2=arr.extent(secondDim);
  if ((n1<2) || (n2<2) || (val1st==val1fin) || (val2st==val2fin)){
    cerr<<"\n arr.extent(firstDim)="<<n1<<" arr.extent(secondDim)="<<n2<<'\n'; 
    err_msg(" \nwrite_splot_csv: invalid extents or ranges");
  }
  char   sepC='\t'; 
  double var1, var2 ;                var1=val1st;       var2=val2st; 
  double step1, step2;       
  step1= (val1fin-val1st)/(n1-1);    step2= (val2fin-val2st)/(n2-1);  

  ofstream outp ;  outp.open(fname, ios::out);

  //  for ( i1=arr.lbound(firstDim); i1<=arr.ubound(firstDim); i1++) {
  //   for  ( i2=arr.lbound(secondDim); i2<=arr.ubound(secondDim); i2++) 
  for ( i1=1; i1<arr.lbound(firstDim)+n1; i1++) {
    for  ( i2=1; i2<arr.lbound(secondDim)+n2; i2++) 
       outp<<val1st+(i1-1)*step1<<sepC<<val2st+(i2-1)*step2<<sepC<<arr(i1,i2)<<endl  ;
     outp << endl ;
  }

  outp.close();
}
// ----------------------------------------------------------------------------
double  mean1D(const Array<double, 1>& arr)
{
  double sum = 0;
  int i1;
  for ( i1=arr.lbound(firstDim); i1<=arr.ubound(firstDim); i1++) 
    sum +=arr(i1);

  return sum/arr.extent(firstDim);
}
// ----------------------------------------------------------------------------
double  mean1D(const Array<float, 1>& arr)
{
  double sum = 0;
  int i1;
  for ( i1=arr.lbound(firstDim); i1<=arr.ubound(firstDim); i1++) 
    sum +=arr(i1);

  return sum/arr.extent(firstDim);
}
// ----------------------------------------------------------------------------
double  var1D(const Array<float, 1>& arr, double *mean)
  // DO NOT use 'more efficient' sum(sq_x) - sq_sum(x) 
  //  - error of difference of two large quantities is large !
{
  if ( arr.extent(firstDim) <3 )
    err_msg("please don't try to calculate varinaces with < 3 data pts !");
  double sumSqDif = 0;
  double avg; 
  if (mean != NULL)  // :-) don't need to recalculate ...
    avg = *mean;
  else
    avg = mean1D(arr) ;  // :-( need to recalc ...

  int i1;
  for ( i1=arr.lbound(firstDim); i1<=arr.ubound(firstDim); i1++) 
    sumSqDif += (avg-arr(i1))*(avg-arr(i1)) ;

  return sumSqDif/(arr.extent(firstDim)-1) ;
}
// ----------------------------------------------------------------------------
double Simpson1DInt(const Array<double, 1>& arr, double start, double end)
{
  div_t temp = div( arr.extent(firstDim), 2) ;
  if ( temp.rem != 1)
    err_msg("Simpson's rule only w. even no. of intervals (odd No. of endpoints)");
  
  double h, Sum0, Sum1, Sum2;
  h = (end-start)/(arr.extent(firstDim)-1) ;
  Sum0 = arr(arr.lbound(firstDim))+arr(arr.ubound(firstDim)) ;
  int i;
  Sum1 = Sum2 = 0.0;
  for ( i=arr.lbound(firstDim)+1; i<arr.ubound(firstDim); i +=2)
    Sum1 += arr(i);
  for ( i=arr.lbound(firstDim)+2; i<arr.ubound(firstDim); i +=2)
    Sum2 += arr(i);

  // cout <<endl<< Sum0<<'\t'<<Sum1<<'\t'<<Sum2<< endl;
  // cout << (h/3)*(Sum0 + 4*Sum1 +2*Sum2);
  return (h/3)*(Sum0 + 4*Sum1 +2*Sum2);
}
// ----------------------------------------------------------------------------
double Simpson1DInt(const Array<float, 1>& arr, double start, double end)
{
  div_t temp = div( arr.extent(firstDim), 2) ;
  if ( temp.rem != 1)
    err_msg("Simpson's rule only w. even no. of intervals (odd No. of endpoints)");
  
  double h, Sum0, Sum1, Sum2;
  h = (end-start)/(arr.extent(firstDim)-1) ;
  Sum0 = arr(arr.lbound(firstDim))+arr(arr.ubound(firstDim)) ;
  int i;
  Sum1 = Sum2 = 0.0;
  for ( i=arr.lbound(firstDim)+1; i<arr.ubound(firstDim); i +=2)
    Sum1 += arr(i);
  for ( i=arr.lbound(firstDim)+2; i<arr.ubound(firstDim); i +=2)
    Sum2 += arr(i);

  // cout <<endl<< Sum0<<'\t'<<Sum1<<'\t'<<Sum2<< endl;
  // cout << (h/3)*(Sum0 + 4*Sum1 +2*Sum2);
  return (h/3)*(Sum0 + 4*Sum1 +2*Sum2);
}
// ----------------------------------------------------------------------------
double pdfSimpsMean1D(const Array<double, 1>& arr, double start, double end)
{
  div_t temp = div( arr.extent(firstDim), 2) ;
  if ( temp.rem != 1)
    err_msg("Simpson's rule only w. even no. of intervals (odd No. of endpoints)");
  
  double x, h, Sum0, Sum1, Sum2;
  h = (end-start)/(arr.extent(firstDim)-1) ;
  Sum0 = start*arr(arr.lbound(firstDim)) + end*arr(arr.ubound(firstDim)) ;
  int i;
  Sum1 = Sum2 = 0.0;
  for ( i=arr.lbound(firstDim)+1; i<arr.ubound(firstDim); i +=2) {
    x = start+ (i-arr.lbound(firstDim))*h ;
    Sum1 += x*arr(i);
  }
  for ( i=arr.lbound(firstDim)+2; i<arr.ubound(firstDim); i +=2) {
    x = start+ (i-arr.lbound(firstDim))*h ;
    Sum2 += x*arr(i);
  }

  // cout <<endl<< Sum0<<'\t'<<Sum1<<'\t'<<Sum2<< endl;
  // cout << (h/3)*(Sum0 + 4*Sum1 +2*Sum2);
  return (h/3)*(Sum0 + 4*Sum1 +2*Sum2);
}
// ----------------------------------------------------------------------------
double pdfSimpsVar1D(const Array<double, 1>& arr, double start, double end, double *mean)
  // DO NOT use 'more efficient' sum(sq_x) - sq_sum(x) 
  //  - error of difference of two large quantities is large !
{
  double avg;
  if (mean != NULL) 
    avg = *mean;
  else {
    avg = pdfSimpsMean1D(arr, start, end);
    *mean = avg;
  }

  double x, h, Sum0, Sum1, Sum2;
  h = (end-start)/(arr.extent(firstDim)-1) ;
  Sum0 = (start-avg)*(start-avg)*arr(arr.lbound(firstDim)) + (end-avg)*(end-avg)*arr(arr.ubound(firstDim)) ;
  int i;
  Sum1 = Sum2 = 0.0;
  for ( i=arr.lbound(firstDim)+1; i<arr.ubound(firstDim); i +=2) {
    x = start+ (i-arr.lbound(firstDim))*h ;
    Sum1 += (x-avg)*(x-avg)*arr(i);
  }
  for ( i=arr.lbound(firstDim)+2; i<arr.ubound(firstDim); i +=2) {
    x = start+ (i-arr.lbound(firstDim))*h ;
    Sum2 += (x-avg)*(x-avg)*arr(i);
  }

  // cout <<endl<< Sum0<<'\t'<<Sum1<<'\t'<<Sum2<< endl;
  return (h/3)*(Sum0 + 4*Sum1 +2*Sum2) ;
}
// ----------------------------------------------------------------------------
double pdfSimpsMean1D(const Array<float, 1>& arr, double start, double end)
{
  div_t temp = div( arr.extent(firstDim), 2) ;
  if ( temp.rem != 1)
    err_msg("Simpson's rule only w. even no. of intervals (odd No. of endpoints)");
  
  double x, h, Sum0, Sum1, Sum2;
  h = (end-start)/(arr.extent(firstDim)-1) ;
  Sum0 = start*arr(arr.lbound(firstDim)) + end*arr(arr.ubound(firstDim)) ;
  int i;
  Sum1 = Sum2 = 0.0;
  for ( i=arr.lbound(firstDim)+1; i<arr.ubound(firstDim); i +=2) {
    x = start+ (i-arr.lbound(firstDim))*h ;
    Sum1 += x*arr(i);
  }
  for ( i=arr.lbound(firstDim)+2; i<arr.ubound(firstDim); i +=2) {
    x = start+ (i-arr.lbound(firstDim))*h ;
    Sum2 += x*arr(i);
  }

  // cout <<endl<< Sum0<<'\t'<<Sum1<<'\t'<<Sum2<< endl;
  // cout << (h/3)*(Sum0 + 4*Sum1 +2*Sum2);
  return (h/3)*(Sum0 + 4*Sum1 +2*Sum2);
}
// ----------------------------------------------------------------------------
double pdfSimpsVar1D(const Array<float, 1>& arr, double start, double end, double *mean)
  // DO NOT use 'more efficient' sum(sq_x) - sq_sum(x) 
  //  - error of difference of two large quantities is large !
{
  double avg;
  if (mean != NULL) 
    avg = *mean;
  else {
    avg = pdfSimpsMean1D(arr, start, end);
    *mean = avg;
  }

  double x, h, Sum0, Sum1, Sum2;
  h = (end-start)/(arr.extent(firstDim)-1) ;
  Sum0 = (start-avg)*(start-avg)*arr(arr.lbound(firstDim)) + (end-avg)*(end-avg)*arr(arr.ubound(firstDim)) ;
  int i;
  Sum1 = Sum2 = 0.0;
  for ( i=arr.lbound(firstDim)+1; i<arr.ubound(firstDim); i +=2) {
    x = start+ (i-arr.lbound(firstDim))*h ;
    Sum1 += (x-avg)*(x-avg)*arr(i);
  }
  for ( i=arr.lbound(firstDim)+2; i<arr.ubound(firstDim); i +=2) {
    x = start+ (i-arr.lbound(firstDim))*h ;
    Sum2 += (x-avg)*(x-avg)*arr(i);
  }

  // cout <<endl<< Sum0<<'\t'<<Sum1<<'\t'<<Sum2<< endl;
  return (h/3)*(Sum0 + 4*Sum1 +2*Sum2) ;
}
// ----------------------------------------------------------------------------
// 2.x  - Interpolation -------------------------------------------------------

// 2.1 
Interpolation1D::Interpolation1D()
{
  // Defaults:
  avoid_div_by_0_a = 1.0e-15;        // additively different from 0
  avoid_div_by_0_m = 1.0 + 1.0e-10 ; // multiplies x (!=0) to something else.
  xType = 0;  // ie default is regularly spaced indep. var. values.
  n = 4;  // n-1 is order of the approximating polynomial in PolynInt etc.
  xa.resize(Range(1,n)); ya.resize(Range(1,n)); 
  xa=0;
  ya=0;                      
  // unsigned long long raw = 0x7ff0000000000000 ;
  // NaNdouble = *( double* )&raw ;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
void Interpolation1D::setOrderPar(int newOrdPar)
{
 n=newOrdPar; 
 xa.resize(Range(1,n)); 
 ya.resize(Range(1,n));
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
void Interpolation1D::setRegInp(    double xInStart, double xInStep
                                ,  Array<double, 1> *yIn)
{
 xStart=xInStart; 
 xStep=xInStep; 
 Y=yIn; 
 xType=0;  
 NDat=Y->extent(firstDim); 
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
int  Interpolation1D::prepXbox(double x)
{
  double aux;
  int i0Y = Y->lbound(firstDim)-1;  // 1st index of ordinate minus one
  int nw=0;      // in the end, xa(1) will be = xStart+nw*xStep
  int n2 = n/2 ; // integer division rounds twd. 0 i.e. down as n and 2 are > 0
  if (xType == 0) {
    aux = (x-xStart)/xStep;
    if (aux > 0) {                    // i.e. xStep moves from xStart towards x
      nw = (int) (floor(aux) + 0.1) ; // no. of whole intervals to x
      if (nw > NDat-n )  // i.e. x is within or past the last n-tuplet of data pts.
	nw = NDat - n; 
      else if  (nw >= n2)  // the usual case
	nw -= n2-1 ; // so we'll center on x, approx.
      else // if nw < n/2
	nw = 0;
    }
    else if (aux < 0) { // i.e. xStep moves from xStart away from x
      nw = 0;
    }
    for (int i = 1; i<=n ; i++) {
      xa(i) = xStart + (nw+i-1)*xStep ;     // this is the key result !
      // ya(i) = (*Y)( i0Y + i + nw);
    }
  }
  else if (xType == 1)
    err_msg("Sadly Interpolation1D::prepXYbox not ready for irreg. spaced x inputs ...");
  else 
    err_msg("in Interpolation1D::prepXYbox, funny xType encountered.");

  // debug:   cout<<"\n xStart: "<<xStart<<"  xStep: "<<xStep<<endl; 
  // cout<<" ya() dim and elements:"<<ya<<"\n xa() dim and elements:"<<xa<<endl;

  return i0Y+1+nw ;  // 
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
void  Interpolation1D::prepXYbox(double x)
{
  int nyoffset = prepXbox(x)-1 ;
  for (int i = 1; i<=n ; i++) 
      ya(i) = (*Y)( nyoffset + i);
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
double  Interpolation1D::PolynInt0(double x, double *err )
{
  //  Various obsolete declarations in the original:
  // void polint(xa,ya,n,x,y,dy)
  // float xa[],ya[],x,*y,*dy;  int n;

  // now n, xa and ya are held by the class and the latter two are blitz arrays.

  int i,m,ns=1;
  double den,dif,dift,ho,hp,w,   y, dy ; // working y and dy moved inside fn.
  double *c,*d; 
  dy = 6.66E+66;         // this is nonsense which should make no difference,
                         // just to keep compiler from moaning.  
  dif=fabs(x-xa(1));
  c= new double[n+1];    d= new double[n+1];
  // Find the index of the closest table entry: 
  for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa(i))) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya(i); // ... and initialise the tableau of c's and d's 
    d[i]=ya(i);
  }
  y=ya(ns--);   // this is the initial approximation to y 
  for (m=1;m<n;m++) {   // for each column of the tableau, loop over the 
    for (i=1;i<=n-m;i++) {       // current c's and d's and update them.
      ho=xa(i)-x;
      hp=xa(i+m)-x;
      w=c[i+1]-d[i];
      // Error if two input xa's are identical to within roundoff: 
      if ( (den=ho-hp) == 0.0) 
	err_msg("in Interpolation1D::PolynInt, xa(i)-x == xa(i+m)-x");
      den=w/den;
      d[i]=hp*den;   // c's and d's updated 
      c[i]=ho*den;
    }
    y += (dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
    /* After each colum i n the tableau is completed, we decide which correction, 
       c or d, we want to add to our accumulating value of y, i.e. which path to take 
       through the tableau - forking up or down. We do this in such a way as to take 
       the most 'straight line' route through the tableau to its apex, updating ns 
       accordingly to keep track of where we are. This route keeps the partial 
       approximations centered (in so far as possible) on the target x. The last dy 
       added is thus the error indication. */
  }
  delete c;  delete d; 
  if (err != NULL) *err = dy; 
  return y;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
double  Interpolation1D::RationInt0(double x, double *err)
{
  // obsolete bits from original: 
  // void ratint(xa,ya,n,x,y,dy)
  //float xa[],ya[],x,*y,*dy; int n;

  // now n, xa and ya are held by the class.
  
  int m,i,ns=1;
  double  w,t,hh,h,dd,  y, dy ; //  working y and dy moved inside fn.
  double *c,*d;                 // ,*vector(); void nrerror(),free_vector();
  
  c= new double[n+1];    d= new double[n+1];
  dy = 6.66E+66;         // this is nonsense which should make no difference,
                         // just to keep compiler from moaning.
  hh=fabs(x-xa(1));
  for (i=1;i<=n;i++) {
    h=fabs(x-xa(i));
    if (h == 0.0) {
      y=ya(i);
      dy=0.0;
      delete d; delete c; 
      if (err != NULL) *err = dy;    return y;
    } 
    else if (h < hh) {
      ns=i;
      hh=h;
    }
    c[i]=ya(i);
    d[i]=ya(i) + avoid_div_by_0_a ;
    // (d[i]==0.0) ?  d[i]=avoid_div_by_0_a : d[i] *= avoid_div_by_0_m ;
     // perturbation to prevent (relatively uncommon)  zero-over-zero  condition.
  }                               
  y=ya(ns--);
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      w=c[i+1]-d[i];
      h=xa(i+m)-x;
      t=(xa(i)-x)*d[i]/h;
      dd=t-c[i+1];
      if (dd == 0.0) {
	cerr <<"in Interpolation1D::RationInt0, interpolating fn. has pole at x=";
	cerr <<x<<";\n xa : "<<xa<<"\n ya : "<<ya<<"\n ... returning NaN: \n\n";
	if (err != NULL) *err = NAN ;     // REM maths.h has :
	return NAN ;                      // #define  NAN   __builtin_nan("")
      } 
      dd=w/dd;
      d[i]=c[i+1]*dd;
      c[i]=t*dd;
    }
    y += (dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  delete d; delete c; 
  if (err != NULL) *err = dy;   return y;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
double  Interpolation1D::RationInt(double x, double *err)
{
  prepXYbox(x);  // this prepares the n-point little array xa
  return RationInt0(x,err);
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
double  Interpolation1D:: RationInt_w_Poly_rescue0(double x, double *err)
{
  double aux;
  aux = RationInt0( x, err) ;
  // debugs: cerr<<"\n Rational interp. output: "; cerr<<'\n'<<aux<<'\t'; 
  if ( isnan( aux) ) {
    cerr<<"\n You may have to worry a bit: RationInt returned NaN at x="<<x;
    cerr<<"\n   RationInt_w_Poly_rescue0 proceeded to use polynomial extrapolation.\n";
    aux = PolynInt0( x, err) ;
    // debug:  cerr<<aux;
    } // end if aux is not a number, ie RationInt0 failed / hit pole
  return aux;   // err should be already taken care of by Ration/PolynInt0
}
// ----------------------------------------------------------------------------
Interpolation2D::Interpolation2D() :
  intX1(), intX2()
{
  y2a.resize(Range(1,intX1.orderPar()),Range(1,intX2.orderPar())) ;
  y2a = 0;
 
  // the following isn't very clever, but will do:
  Y1.resize(Range(1,3));  Y2.resize(Range(1,3));  Y1=0; Y2=0;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
void Interpolation2D::setRegInp2D_x(  double xIn1Start, double xIn1Step
                                    , double xIn2Start, double xIn2Step )
{
  intX1.setRegInp_x(xIn1Start,xIn1Step);
  intX2.setRegInp_x(xIn2Start,xIn2Step);
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
void Interpolation2D::setInp2D_y( Array<double, 2> *yIn2D)
{
  Y1.resize(Range(1, yIn2D->extent(firstDim)));
  Y1 = (*yIn2D)(Range::all(),1); // second index at 1 - hopefully correct !
  Y2.resize(Range(1, yIn2D->extent(secondDim)));
  Y2 = (*yIn2D)(1,Range::all()); // first index at 1 - hopefully correct !

  Y2DP = yIn2D;

  intX1.setInp_y( &Y1 );   // a bit of a waste of space, but we need some
  intX2.setInp_y( &Y2 );   // of this stuff for the class members.
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
void Interpolation2D::setRegInp2D( double xIn1Start, double xIn1Step
                   ,double xIn2Start, double xIn2Step, Array<double, 2> *yIn2D)
{
  Y1.resize(Range(1, yIn2D->extent(firstDim)));
  Y1 = (*yIn2D)(Range::all(),1); // second index at 1 - hopefully correct !
  Y2.resize(Range(1, yIn2D->extent(secondDim)));
  Y2 = (*yIn2D)(1,Range::all()); // first index at 1 - hopefully correct !

  intX1.setRegInp(xIn1Start,xIn1Step,&Y1);
  intX2.setRegInp(xIn2Start,xIn2Step,&Y2);

  Y2DP = yIn2D;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
void Interpolation2D::setOrderPars(int newOrd1, int newOrd2)
{
 intX1.setOrderPar(newOrd1); 
 intX2.setOrderPar(newOrd2); 
 y2a.resize(Range(1,newOrd1),Range(1,newOrd2));
 y2a = 0;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
void Interpolation2D::prepX1X2Ybox(double x1, double x2)
{
  int n1=intX1.orderPar(); int n2=intX2.orderPar();
  int ny1offset = intX1.prepXbox(x1)-1;
  int ny2offset = intX2.prepXbox(x2)-1;
  double aux;

  // cout<<(*Y);

  for (int i1=1; i1<=n1 ; i1++) 
    for (int i2=1; i2<=n2 ; i2++) {
      aux = (*Y2DP)(ny1offset+i1, ny2offset+i2);
      y2a(i1,i2) = aux ;
    }
 
  cout<<"\nx1 to est.: "<<x1<<"\nx2 to est.: "<<x2;
  cout<<"\nFirst y1 index: "<< ny1offset+1<<"\tfirst y2 index: "<<ny2offset+1<<'\n';
  cout<<"\n y2a : "<<y2a;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
double Interpolation2D::RationInt(double x1, double x2, double *err)
{ 
  int k,j;
  double err1, aux;
  int n1   = intX1.orderPar();        int n2 = intX2.orderPar();
  // starting offsets: 
  int yso1= intX1.prepXbox(x1)-1 ;    int yso2 = intX2.prepXbox(x2)-1;

  // Original, obsolete: void polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
  // float x1a[],x2a[],**ya,x1,x2,*y,*dy;

  for (j=1;j<=n1;j++) {

    // first, fill the ya vector along the X2 direction with the
    // right row segment from the external dep. variable matrix:
    for (k=1; k<=n2; k++) 
      // aux =  (*Y2DP)(yso1+j,yso2+k);
      intX2.set_ya_el((*Y2DP)(yso1+j,yso2+k), k) ;     // aux, k) ;
    
    // intX2.set_ya_el( (*Y2DP)(yso1+j,yso2+k), k) ;

    // Now fill the appropriate element of the ya vector along
    // the X1 direction with the interpolation based on the X2 dir.:
    aux = intX2.RationInt0( x2 ); // no track of error here.
    intX1.set_ya_el(aux, j);


  }
  aux = intX1.RationInt0( x1, &err1 );

  if (err !=NULL) *err = err1;
  return aux;

}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// REST OF CLASS MUST BE SET UP CONSISTENT W. y2dp - NOT CHECKED HERE FOR SPEED !
double Interpolation2D::RationInt(  double x1, double x2
                                  , const Array<double, 2>& y2dp, double *err)
{ 
  int k,j;
  double err1, aux;
  int n1   = intX1.orderPar();        int n2 = intX2.orderPar();
  // starting offsets: 
  int yso1= intX1.prepXbox(x1)-1 ;    int yso2 = intX2.prepXbox(x2)-1;

  for (j=1;j<=n1;j++) {

    // first, fill the ya vector along the X2 direction with the
    // right row segment from the external dep. variable matrix:
    for (k=1; k<=n2; k++) 
      intX2.set_ya_el( y2dp(yso1+j,yso2+k), k) ;    
    

    // Now fill the appropriate element of the ya vector along
    // the X1 direction with the interpolation based on the X2 dir.:
    // aux = intX2.RationInt0( x2 ); // no track of error here.
    intX1.set_ya_el( intX2.RationInt0( x2 ), j);

  }
  aux = intX1.RationInt0( x1, &err1 );

  if (err !=NULL) *err = err1;
  return aux;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// REST OF CLASS MUST BE SET UP CONSISTENT W. y2dp - NOT CHECKED HERE FOR SPEED !
double Interpolation2D::PolynInt(  double x1, double x2
                                  , const Array<double, 2>& y2dp, double *err)
{ 
  int k,j;
  double err1, aux;
  int n1   = intX1.orderPar();        int n2 = intX2.orderPar();
  // starting offsets: 
  int yso1= intX1.prepXbox(x1)-1 ;    int yso2 = intX2.prepXbox(x2)-1;

  for (j=1;j<=n1;j++) {

    // first, fill the ya vector along the X2 direction with the
    // right row segment from the external dep. variable matrix:
    for (k=1; k<=n2; k++) 
      intX2.set_ya_el( y2dp(yso1+j,yso2+k), k) ;    
    

    // Now fill the appropriate element of the ya vector along
    // the X1 direction with the interpolation based on the X2 dir.:
    // aux = intX2.RationInt0( x2 ); // no track of error here.
    intX1.set_ya_el( intX2.PolynInt0( x2 ), j);

  }
  aux = intX1.PolynInt0( x1, &err1 );

  if (err !=NULL) *err = err1;
  return aux;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// REST OF CLASS MUST BE SET UP CONSISTENT W. y2dp - NOT CHECKED HERE FOR SPEED !
double Interpolation2D::RationInt_w_Poly_rescue (  double x1, double x2
                                           , const Array<double, 2>& y2dp, double *err)
{ 
  int k,j;                        double err1, aux;
  int n1   = intX1.orderPar();    int n2 = intX2.orderPar();
  // starting offsets: 
  int yso1= intX1.prepXbox(x1)-1 ;    int yso2 = intX2.prepXbox(x2)-1;

  for (j=1;j<=n1;j++) {
    // first, fill the ya vector along the X2 direction with the
    // right row segment from the external dep. variable matrix:
    for (k=1; k<=n2; k++) 
      intX2.set_ya_el( y2dp(yso1+j,yso2+k), k) ;    
    
    // Now fill the appropriate element of the ya vector along
    // the X1 direction with the interpolation based on the X2 dir.:
    intX1.set_ya_el( intX2.RationInt_w_Poly_rescue0( x2 ), j);
  }
  aux = intX1.RationInt_w_Poly_rescue0( x1, &err1 );

  if (err !=NULL) *err = err1;
  return aux;
}
// ----------------------------------------------------------------------------
//                     End of file
