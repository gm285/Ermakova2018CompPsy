/*          ******  Maths functions, Unix version  *********
 
     1.x - Functions - e.g. y = mx + c
     2.x - Piecewise linear fns.
     3.x - Misc. c fns : rounding, max & min of arrays, trapezium integrn.,
           rms-diff of 2 fns, some random-based e.g. bin_dec, Kroneker delta ...
     4.x - 'Iter 'ation class etc ...
     5.x - Ramps
     6.x - Sigmoids (but WJF sigmoids are in WJF1.H) ...
     7.x - Fitting functions to data
     8.x - Some standard dynamical systems .1 Silnikov, 
           .2 harmonic oscillators HO (incl. SHO->DHO->FDHO)
     9.x - Physics derived functions such as Gibbs
    10.x - Receiver Operating Characteristic (ROC) ...
*/ 

#ifndef MATHS_H
#define MATHS_H 

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>   // for debugging etc.
#include <cstring>    // these last two are needed from GCC 4.3
#include <memory>

#include <blitz/array.h>
BZ_USING_NAMESPACE(blitz)

#ifndef PI
#define PI 3.14159265458
#endif
#ifndef SMALLFLTSIGN
#define SMALLFLTSIGN 0.00001
#endif

/* following was const enum Boolean { False=0, True=1 } ;  Linux C++ has
   its own, 'bool' type which is 'true' or 'false', Swan book p. 103 */
typedef enum { False, True } Boolean ;   

typedef int      *int_p ;      typedef int      *  *intArr_2p ;
typedef char     *char_p ;     typedef char     *  *charArr_2p ;
typedef float    *float_p ;    typedef float    *  *floatArr_2p ;
typedef double   *double_p ;   typedef double   *  *doubleArr_2p ;

// 3.x general purpose C functions ------------------------------------------

// a. General error (derived from nrerror)
int err_msg( const char *error_text) ;

/* 3b.floatArrayMAX & 3b floatArrayMIN find the maximum or minimum element
      of an array of AT LEAST TWO reals . They assume the data is found
      between variable[firstIndex] & variable[ArPosn] (as ArPosn initially
      given) . Similarly for ...AVG, for the mean etc. 
*/
float floatArrayMAX(float variable[], int ArPosn, int firstIndex=0) ;
float floatArrayMIN(float variable[], int ArPosn, int firstIndex=0) ;
float floatArrayAVG(float variable[], int ArPosn, int firstIndex=0) ;
float floatArraySD( float variable[], int ArPosn, int firstIndex=0) ;
// Next a lttle unusual: *Dev != NULL, rms dev. from value pointed to by Dev. 
float floatArrayDev(float variable[], int ArPosn, int firstIndex=0
                    , float *Dev = NULL  ) ;

// 3c simple lin. algebra - det. of 1-3x1-3 matrix etc :                      
float Det3( float **a ) ;     // a from [1][1] DANGEROUS - NO CHEKING AT ALL
float Det3(Array<float, 2> A);  // much safer !

// 3d 
int roundfl( float x ) ;          // Rounds 'float' to NEAREST int.          
long int rounddb( double x ) ;    // Rounds 'double' to NEAREST long int.
float  FAbs( float x) ;           // my version of fabs
double DAbs( double x) ;          //  ... similar.

// 3e lin. interpolate:  
// i. one point : "return *y+(xNew-*x)*(*(y+1)-*y)/(*(x+1)-*x) "
float LinInter(const float *y, const float *x, const float xNew);

// ii. Get interpolate vect. var. array into regular interval array.
//     put interpolated values in array **yInt, starting from [0][0].
//     yDat[0][indStart ... indEnd]  points to the parametre / indep. var., 
//     yDat[1 to depDatDim] to the dep. vars / ordinates.
//     yDat[0] must be in ascending order. Nx is the no. of
//     points needed in the output. Function searches yDat[0] for xStart &
//     xEnd containing intervals.
int  ArrLinInter(float **yInt,
                 float **yDat, 
                 const float xStart, const float xEnd, const int   Nx,  
                 const int   indEnd, const int   indStart = 0,
		 const int   depDatDim  = 1 ) ;

/* 3e Integration based C functions */

// i.  Trapezoidal integration of tabulated fn.
double TrapInt(  const float *y, const float *x
	       , const float xStart, const float xEnd
	       , int   lastInd, int firstInd = 1 );

// ii. mean sqrt of  [U(T) - V(T)]^2 i.e. rms diff of 2 fns. 
// USE ONLY *U, *Tu etc going from 1 to Nu and are in ascending order
// of Tu and Tv values !
double rms_tser(const float *U,     const float *Tu, int Nu,
                const float *V,     const float *Tv, int Nv,
                const float tStart, const float tEnd         ) ;

/* 3f - Range adapting / checking */ 
/* i. to save typing in range checks. Is x outside the interval (a,b) ? */ 
inline char outOf(double x, double a, double b) {return ((x-a)*(x-b) > 0) ; };

/* Adapting angles to ranges */
/* ii. take an angle in rad, which can be +ve or -ve, reduce it from 0 to 2*PI
       and optionally return the no. of full rotns. in the angle : */ 
float phi0to2PI( const float phi, int *kappa = NULL ) ;
/* iii. similar to above (and uses it), but rets within -PI to PI : */ 
float phiMinusPItoPI( const float phi, int *kappa = NULL ) ;

/* 3g - Random-based functions ASSUME srand has been run, e.g. by ... */ 
/*      ...  srand( ( time( NULL ) % 604800) ) ;                      */
/* i. bin_ded returns 1 (accept) or 0 (reject) decision based on inputed  */
/*    probability P:                                                      */
char bin_dec(double P) ; 
/* ii. shorthand for producing random numbers in interval (0,1):          */
inline double rand01(){ return 1.0*rand()/RAND_MAX; };

/* 3h - Kroneker delta (for clarity - doesn't save much coding ! */ 
inline int  KronD(int a, int b){ return (a==b) ? 1:0 ; } ;

/* for storing a filename type string of the form filename045.csv . Not 
   actuall a math fn. but put here for convenience. Number index 0 to
   999, user has to make sure there's memory for returned string. Returns
   length actually used.  */
int fnstring(char *wholename, char *stem, int ind, char *ext);

// 4.x Data Structures for plotting etc.  -----------------------------------
class BasicData {
protected:
  int alpha, omega;    // first & last index
public:
  int Num(){return (omega-alpha+1) ;};
  int firstInd() {return alpha;};
  int lastInd() {return omega;};
};

// 4a. ----------------------------------------------------------------------
// 4.b Iteration classes ----------------------------------------------------- 
class Iter  // iterate over a bunch of parametres.                             
            // 1st param refers to  innermost in loop.                         
   /* to use with while loop :   Iter paramIter( .... );                       */
   /* while(paramIter.Do()) { .... paramIter.Step(); }                         */
{
  char *auxStr;       /* will store 6-dig accurate stuff etc */ 
  float **parPtr ;    // points to par. to change each time 
  float **iterNo ;    //    "    " inputted iter numbers to take
  float **endPtr;     
  float *start;      // within class stored starting points
  float *step ;      //   "     "    stored step sizes
  int   parNo ;
  int   doAgain ;
  double  epsilon ;  // this to compare float 'equality'
  unsigned long itersDone ;    // set to 0 by Set(..., augmented by Step()
  int *stepIndex ;   // monitors where we are, in terms of step count, for
                     // each parametre. Usu. set to 1 1 1 ... initially.
  unsigned long totN ;         // total N. of iterations (product of 
                                    // each loop) cf TotN()

 public:
   
  Iter( int   newParamNo ) ; // Only save space for some var.
  ~Iter();
  /* 'Set(...' is required before any use of the class. It sets the pointers,
     but also the start[i] internal vars acc. to **newparam vals, the epsi-
     lon, sets doAgain = 1 and calls TotN()  */ 
  void Set( float **newparam,  float **newEnd  , float **newIterNo );  
  
  int Do(){ return doAgain ;} ;
  int Step();  /* returns whether more iterns left to do or not */ 
  float Start(int i){ return start[i] ; };
  float End(int i){ return *endPtr[i] ; };
  int ParNo(){ return parNo ; } ;
  unsigned long ItersDone() { return itersDone ; } ;
  int AtStep(int parIndex) ;  // Where iteration is w.r.t. param parIndex
  unsigned long TotN();  /* find & ret. total N. of iter. & store in totN */
  float FractDone() { return 1.0*itersDone / totN ; } ;
  /* Step, then FractDone ; slower than separately, but in one package. Rets.
     whether to proceed, fraction done and also no. of iterns. done : */ 
  int StepFD(float *fractDone = NULL, unsigned long  *itDone = NULL ) ;
  /* Output  the No. of the itern. UNDER WAY, and the params. to the error 
     device => ParsToErr is used INSIDE LOOP, BEFORE stepping fn. Can use
     a header, outside iter. loop, such as 
     cerr<<"\nn/Tot.N\tPar1\tPar2\tPar3\tPar4\n" ;                          */ 
  void ParsToErr() ;
  /* Value of specific parametre : */
  float valPar( int i) { return *parPtr[i] ; } ;  
  /* a par 'v. near' specific value; Used in selecting iterns : */
  int nearParVal(int i, float val) 
    { return (DAbs(val - *parPtr[i]) < epsilon) ? 1 : 0 ; } ;
} ;

// 1.1 Functions e.g. y = mx + c ----------------------------------------------
class Lin_fn_float 	      // simple linear fn. for 'float's
{
    float Slope ;
    float Constant ;

  public:
    void SetSlope(float NewSlope){Slope = NewSlope;} ;
    void SetConstant(float NewConst){Constant = NewConst;}  ;
    float GetSlope(){return Slope ;} ;
    float GetConstant(){ return Constant ;} ;
    float Val(float x){	return Slope*x+Constant; } ;
};

// 2.x Piecewise linear fns. --------------------------------------------------
class LinMapBasics
// code common to all piecewise linear functions, which are defined below
//   each one separately for speed 
{
 protected:
    Lin_fn_float *Piece ;       // Array of linear pieces
    unsigned int PieceNumber ; 	// Number of pieces that make up the function
    float *Lim ;                // Limit Points defining piecewise ranges

 public:
    // Constructors differ --> leave for derived classes
    // Destructor common :
    ~LinMapBasics() { delete Piece; delete Lim; };

    // Seting various params. Note piece no. start from 0 not 1 ! :
     void SetLimit(unsigned int No, float NewVal){
      (No <= PieceNumber) ?  Lim[No] = NewVal:err_msg("Lim No too big" );};
     void SetPiecePars(int No, float M, float C){
      Piece[No].SetSlope(M); Piece[No].SetConstant(C);};

     void SetPcByEndpts(unsigned int PcNo, float y1, float y2);
	// calculates and sets slope & const. by SEGMENT ENDPTS.
	// It doesn't alter the limits!

     void Set_pts(float **data);        //sets all PARAMETRES: slopes etc.
 
     int PieceNo(){ return PieceNumber ;} ;
     float MapContent(int Row, int Col) ;
	//Row 0->limits,1->slopes,2->constants ; Col's-># of pts.
}; // end  of Linear basics class declaration

// 2.1 ---------------------------------------------------------------------
class OneLinPcMap : public LinMapBasics
{
 public:
    // Constructor
    OneLinPcMap() { PieceNumber=1;
      Piece = new Lin_fn_float; Lim = new float[2] ;};
      // note labels of pieces start from 0 not 1
    // Destructor as in base class
     float Val1Pc(float X){ return  (X < Lim[0] || X > Lim[1]) ?
     err_msg("out of range in Val1Pc(X)" ) : Piece[0].Val(X) ;   } ;
};

// 2.2 ---------------------------------------------------------------------
class TwoLinPcMap : public LinMapBasics        //explanations in OneLinPcMap
{
 public:
    TwoLinPcMap() { PieceNumber=2;
      Piece = new Lin_fn_float[2]; Lim = new float[3] ;};
    float Val2Pc(float X){ return  (X < Lim[0] || X > Lim[2]) ?
     err_msg("out of range in Val2Pc(X)" ) :
     ((X<Lim[1]) ?  Piece[0].Val(X) : Piece[1].Val(X) ) ;} ;
};
// 2.3 ---------------------------------------------------------------------
class ThreeLinPcMap : public LinMapBasics      //explanations in OneLinPcMap
{
 public:
    ThreeLinPcMap() { PieceNumber=3;
      Piece = new Lin_fn_float[3]; Lim = new float[4] ;};
    float Val3Pc(float X){ return  (X < Lim[0] || X > Lim[3]) ?
     err_msg("out of range in Val3Pc(X)" ) : ((X<Lim[1]) ?  Piece[0].Val(X) :
      (X<Lim[2]) ? Piece[1].Val(X):(Piece[2].Val(X)) ) ;} ;
};
// 2.4 ----------------------------------------------------------------------
class FourLinPcMap : public LinMapBasics        //explanations in OneLinPcMap
{
 public:
    FourLinPcMap() { PieceNumber=4;
      Piece = new Lin_fn_float[4]; Lim = new float[5] ;};
    float Val4Pc(float X){ return  (X < Lim[0] || X > Lim[4]) ?
     err_msg("out of range in Val4Pc(X)" ) : ( (X<Lim[2]) ?
      ( (X<Lim[1]) ? Piece[0].Val(X):Piece[1].Val(X) ) :
      ( (X<Lim[3]) ? Piece[2].Val(X):Piece[3].Val(X) ) ) ;} ;
};
// 2.5 ----------------------------------------------------------------------
class FiveLinPcMap : public LinMapBasics	//explanations in OneLinPcMap
{
 public:
    FiveLinPcMap() { PieceNumber=5;
      Piece = new Lin_fn_float[5]; Lim = new float[6] ;};
    float Val5Pc(float X){ return  (X < Lim[0] || X > Lim[5]) ?
     err_msg("out of range in Val5Pc(X)" ) : ((X<Lim[2]) ? ((X<Lim[1]) ?
     Piece[0].Val(X):Piece[1].Val(X)):((X<Lim[3]) ? Piece[2].Val(X):
     ((X<Lim[4]) ? Piece[3].Val(X):Piece[4].Val(X)))) ;} ;
};

// 6.x C++ fns.: Starting with 6.1  (non-WJF) Sigmoids -----------------------
class plainSig {
 public:
  float Val0( float x )  { return 1.0/(1+exp(-x)); };

  // Heavyside step function as limiting case of sigmoid :
  int Heavyside(float x) { return ((x>0) ? 1 : 0); };
};
//  6.2   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
class LinSig : public OneLinPcMap {
 protected:
  float lVal, rVal;  // values for x-> -inf and +inf . Def. 0 and 1.

 public:
  LinSig();
  float LinV(float x) { if (x >= Lim[1]) return Val1Pc(Lim[1]) ;
	  if (x> Lim[0]) return Val1Pc(x) ;
          else         return Val1Pc(Lim[0]) ;              } ;
};
//  6.3   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
class Sig : public plainSig {             // logistic sigmoid with parametres
 protected:
  float midPt;                            // point of inflection of sigma
  float m, c, a, b ;                      // auxiliary
 public :
  // infl is abcissa of point of inflection, riseT is 1/max.slope.
  Sig( float infl, float riseT, float baseY, float ceilingY )  {
    m= 4/riseT ; c=-4*(midPt=infl)/riseT ;  a=-baseY+ceilingY ; b= baseY ;
   } ;

  // fully parametrised logistic sigmoid
  float Val( float X ) { return a*Val0(m*X+c)+b ; };

  // parametrized Heavyside step function
  float Step(float x) { return ((x>midPt) ? a+b : b);};
};
// 5.1 Ramps  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class Ramp {   // smoothly ramping a value up to a plateau and back again
 protected:
   float t1, t6, Dt, ylo, yhi; //= tstart,tend,ramp-time,resting val,hi val.
   float t2 ,t3, t4, t5, a   ; //  auxiliary times and parametre
   float W,  ymid, yAmp      ; //  aux. params. for Sinusoidal ramp.
		  
/*           <-Dt->                        <-Dt->  = tramp
  
     yhi	  .------------------------.
  		 /t3                      t4\
     ymid	/                            \
  	       /t2                          t5\
     ylo  ___./                                \.________
  	    t1=tstart                            t6=tend  
*/       

 public:
   // Ramp();
   void SetParam(float tstart, float tend,    float tramp,
   			       float yhigh=1, float ylow=0);
   // Val. if Ramp= product of 2 logistic sigmoids, centred at t2 and t5 :
   float lgVal( float t) {
    return ylo + (yhi-ylo)/((1+exp(4*(t2-t)/Dt))*(1+exp(4*(t-t5)/Dt))) ;   
   } ;
   //Each Ramp end made out of 2 parabolic pieces:
   float ParabVal(float t) {
    if (t<t4) {
       if (t>t3) {
          return yhi;
         }
       else if (t<t1) {
          return ylo;
         }
       else if (t<t2) {
          return ylo + a*(t-t1)*(t-t1) ;
         }
       else return yhi - a*(t3-t)*(t3-t) ;
      }
    else {
       if (t>t6) {
          return ylo;
         }
       else if (t<t5) {
          return yhi - a*(t-t4)*(t-t4) ;
         }
       else return ylo + a*(t6-t)*(t6-t) ;
      }
   };
   //Each ramp end made out of 1 sinusoidal piece:
   float SinVal(float t) {
    if (t<t4) {
       if (t>t3) {
          return yhi;
         }
       else if (t<t1) {
          return ylo;
         }
       else return ymid + yAmp*sin(W*(t-t2));
      }
    else {
       if (t>t6) {
          return ylo;
         }
       else return ymid + yAmp*sin(W*(t5-t)) ;
      }

   };
   // Various params :
   float midX1()  { return t2 ; } ;
   float midX2()  { return t5 ; } ;
   float riseDX() { return Dt ; } ;
};
/* 7.1 Fitting a parabola to 3 points -------------------------------------- */
class Quad3Fit {
  Array<float, 2> A, B ;
  float D        ;
  float a, b, c  ;

 public:
  Quad3Fit()  ;
  ~Quad3Fit() ;

  float Set(  float x1,  float y1     // float to return Det B, used
	    , float x2,  float y2     // as denom. in soln. of eqns.
	    , float x3,  float y3 ) ; // For user error testing.

  float Coef2()   ;      /* a ; rets error on D == 0.0       */
  float Coef1()   ;      /* b ; rets error on D == 0.0       */
  float Coef0()   ;      /* c ; rets error on D == 0.0       */
  float OrdStat() ;      /* ret ordinate of stationary point */
  float AbcStat() ;      /* "   abcissa  "     "       "     */

};

// 8.x Some dynamical systems  ----------------------------------------------
// 8.1 Silnikov w. deflt chaotic params , optional            ---------------
//     offset & time scaling and variable space.              ---------------
// See test_Siln_w_Y_YDotSer.cc for usage. Advice is to always keep variables
// in order, and if need to store x an xdot side by side to arrange the ap-
// priate storage pointers, esp. w.  Y_YDotSeries::SelectForStorage 

class Silnikov {

 protected:
  double xOffset, yOffset, zOffset ;
  double a, b, c;
  double timeScale ;

 public:
  float *s ;      // init. cond & current state can be stored here by ODEs
                  // it's public: USE WITH CARE !
  Silnikov()  {                                              // 1. dflt, for arrays etc.
    s = new float[5] ;     a=0.65; b=0.2217; c=0.55; 
    timeScale=1.0;         xOffset=yOffset=zOffset=0.0; } ;
  Silnikov( const Silnikov & rval) ;                         // 2. copy ctor

  ~Silnikov() { if(s) delete s ;  } ;

  // Setting (or re-setting to deflt) the basic parametres :
  void SetBasics( double newA=0.65, double newB=0.2217, double newC=0.55)
    { a = newA; b= newB; c=newC; } ;

  // Setting auxilliary params. Deflt. vals are such that they APPROXIMATELY
  // zero the corresponding means :
  void SetAux(  double newTScale,     double newXOff=0.67
	      , double newYOff = 0.0, double newZOff=0.0  ) ;

  // Setting the values of variables further down ...

  // fast fns. to ret. the basic rates, without auxilliary params. See
  // below - for rate of xdot here is simply z !
  double XDot(double x, double y, double z) { return y; } ;
  double YDot(double x, double y, double z) { return z; } ;
  double ZDot(double x, double y, double z) { return (x*(a-b*x-a*x*x)-y-c*z); };

  // rates With Auxilliary params. NB if we need an x, xdot pair of time series,
  // and x(t) is given by ODE solving the 1st var, then the t series of its ti-
  // me der. will be v(t)= d/dt( x(t) )= XDot..(), and the RATE of this 
  // (that can be directly used by ODE routines) will be 
  // dv/dt  = d/dt(timeScale*(y+yOffset)) = timeScale*YDot...(...) =
  //        = timeScale^2 *(z+zOffset)
  double XDotWAux(double x, double y, double z) { 
    return timeScale*(y+yOffset); } ;
  double YDotWAux(double x, double y, double z) { 
    return timeScale*(z+zOffset); } ;
  double ZDotWAux(double x, double y, double z) {
    return timeScale*((x+xOffset)*(a -b*(x+xOffset) 
		       - a*(x+xOffset)*(x+xOffset))      
	               - (y+yOffset) - c*(z+zOffset)) ; } ;
  
  // Some applications need Y - YDot pairs so can use this instead of y
  // to solve for t series of derivative of first var. Or could simply 
  // just calculate it from solved t series for y algebraically. 
  double XD_DotWAux(double x, double y, double z) {  // t-der of above
    return timeScale*timeScale*(z+zOffset); } ;

  // Setting the values of the variables. The fourth one is set automatically,
  //   as it' s just the deriv. of the first.
  void SetVar(double x= 0.764095, double y=0.173443, double z=0.118232 ) 
   {     // dflt. val. nr. dflt. attractor
    s[1]=x; s[2]=y; s[3]=z;
    s[4]=y;                     } ;  /* see above - simple xDot equals y.     */
  void SetVarWAux(double x= 0.19613, double y= 0.287184, double z= -0.185299) 
    {     /* dflt values are nr. approx-zeroed-mean-x attr.                   */
    s[1]=x; s[2]=y; s[3]=z;                             
    s[4]=XDotWAux(x, y, z);     } ;  
  void RandomSVar() ;               /* ought to be improved ...               */

  int EqN() { return 3; } ;   /* how trivial can OOP get ?? Yet there's space */
  int VarN(){ return 4; } ;   /* for a redundant 4th var in the array s[].    */
};

/* 8.2 Harmonic oscillator class hierarchy   -------------------------------- */
/*                                           -------------------------------- */

class HO { /* ~ inert*xddot+resis*xdot+resto*x = 0 (or = force*cos(W*t + phi) )
               qdot = p, qddot = pdot = -2*b*p - W0SQ*q + F*cos(W*t + phi)
               Note p is the 'velocity', not the momentum  ! This is i. !     */
 protected:
  float W0SQ, b, F;                           /* as for calculation           */
  float inert, resis, resto, force, W, phi ;  /* as set up                    */

 public:
  float *x ;  /* x[1] is q, x[2] is p = qDot. dflt init cond q=0, p=1.        */

  HO() ;                               /* 1. dflt, for arrays etc.            */
  HO( const HO &rval );                /* 2. copy ctor                        */
  ~HO() { if (x) delete x; } ;

  void SetDHOPar( float new_resto, float new_resis=1, float new_inert=1)  ;
  void SetForcPar(float newforce, float newphi=PI/2, float newW=1 ) ;
  void SetVar(float q, float p) { *(x+1)=q; *(x+2)=p; } ;
  
  /* Get Parametres :                                                         */
  void GetDHOPar(float *restor,      float *resist=NULL, float *inerti =NULL, 
                 float *forcin=NULL, float *Wforci=NULL, float *phiforc=NULL ) ;

  /* Non - forced damped HO rates of change :                                 */
  float QDotNF(float q, float p) { return p; } ;    /* as by construction     */
  float PDotNF(float q, float p) { return (-2*b*p -W0SQ*q) ;} ; 

  // Forced damped HO rates of change :
  float QDot(float q, float p, float t) { return p; } ;   // as by construction
  float PDot(float q, float p, float t) { 
    return (-2*b*p -W0SQ*q + F*cos(W*t + phi)) ;} ; 


  /* Resonant fcy Wres, decay time Tdec1 (overdamped - also Tdec2)
      Returns underdamped->-1, overdamped or crit. damped -> 1               */
  int OutPar( float *Wres=NULL, float *Tdec1=NULL, float *Tdec2=NULL ) ;
  /* forced steady state amplitude and phase (formulae in maths.cc) :        */
  void ForcedAmpPh( float *amp, float *phase) ;

  int EqN() { return 2; } ;   /* how trivial can OOP get ??                  */
};

/* 9.1 Gibbs distribution   ------------------------------------------------- */
/*                                           -------------------------------- */
class Gibbs  { /*      */
 protected:
    double T;
    double Len0();  // Sum_j[ exp(E_j/T) ]  NOTE NO K_Boltzman !
    Array<double, 1> E ;    // ('handle' to) ext. vector w. 'energies'
 
 public:
    void SetEPtr(Array<double, 1> energy_vals){ E=energy_vals; };
    void SetT(double temp){ T=temp; };
    // this is the Py0 = exp(E_state/T) / Sum_j[ exp(E_j/T) ] NOTE NO K_Boltzman !
    double Py0( int state) { return exp((E(state))/T) / Len0(); };
};

/*  10.x - Receiver Operating Characteristic (ROC) ... ---------------------- */
/*                                           -------------------------------- */
class ROC  { /*      */
 protected:
  float m ;        // general purpose parameter
 public:
  float Sn, Sp;
  float FNC, FPC, TNC, TPC ; // various types of cost.

  void set_m(double new_m){ m = new_m ; } ; 
  void set_costs(double fnc, double fpc, double tpc, double tnc) {
    FNC=fnc; FPC=fpc; TPC=tpc; TNC=tnc; };
  float mVal(){return m;} ;
  float currSn(){return Sn ;};  // Current sensitivity, as stored.
  float currSp(){return Sp ;};  // Current specificity, as stored.
};

/* ROC made of one hyperbolic piece ----------------------------------------  */
class HypROC : public ROC {
 protected:
  double a,b,c ; // basic params
 public:
  HypROC();      // default: m=4.0
  void SetPar(double New_M); // set a,b,c,  also checks New_M > 1
  // the following don't change the Sn and Sp stored in the base class:
  double SpFnSn(double SnVal) { return c+1/(a*SnVal+b) ; } ;
  double SlopeFnSn(double SnVal) {  
    double aux= a*SnVal+b ;      return -a/(aux*aux) ; } ;
  double CostMinSn(double P) {
    return (1/a)*(-b + sqrt( ((a*P)/(1-P))*(FNC/FPC) ) ) ; } ;
  // These also change the values stored in the base class (nb floats !) :
  float SetSpFnSn(float SnVal) { return (Sp=(c+1/(a*(Sn=SnVal)+b))) ; } ;
  float SetCostMinSn(float P) {
    Sn=(1/a)*(-b+sqrt(((a*P)/(1-P))*(FNC/FPC))); Sp=c+1/(a*Sn+b); return Sn; };

};

// ----------------------------------------------------------------------------


#endif    // MATHS_H

//                     End of file

