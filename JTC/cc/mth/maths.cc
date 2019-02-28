// MATHS.CC : Maths routines - Code for functions in MATHS.H

#include "maths.h"     // which contains the various constants etc.

// A simple error messaging method: --------------------------------------
int err_msg( const char *error_text)
{
 fprintf(stderr,"Run-time error :\n");
 fprintf(stderr,"%s\n", error_text );
 fprintf(stderr," :-( sadly we have to exit. \n") ;
 exit(-1) ;
 return -1; // never reached - for consistency.
}
// ----------------------------------------------------------------------------
// 4.2 Iteration classes
Iter::Iter( int   newParamNo )
{
 parNo = newParamNo;  
 start  = new float[parNo+1] ;   // to begin from start[1] 
 step   = new float[parNo+1] ;  
 auxStr = new char [12] ;        // will store 6-dig accurate stuff etc   
 stepIndex = new int[parNo+1] ;  

 doAgain = 0 ;             // ie. requires Set function before use of class.
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Iter::~Iter()
{
  delete start; delete step ; delete stepIndex ;

} 
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int Iter::Step()
{
  if (totN == 0)  return (doAgain = 0) ;

  int i;
  for(i=1; i<=parNo; i++) {
   if( *iterNo[i] > 1) {
     if( DAbs(*endPtr[i] - *parPtr[i]) > epsilon ) {  // i.e if looping at this 
       *parPtr[i] += step[i] ;                        // level incomplete
       itersDone += 1;        // augment overall stepping monitor
       stepIndex[i] += 1;     // augment local stepping monitor
       return (doAgain=i) ;   // rough index of where we are
     }
     // else done at this level: reset this and allow i-> i+1 :
     *parPtr[i] = start[i] ;    stepIndex[i] = 1 ;
   } // if *stepNo[i] == 1 just go to next level.
  }

  // if all levels done, set doAgain to False :
  itersDone = totN ;  // hopefully covers 0, 1, >1
  return (doAgain = 0) ;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int Iter::StepFD( float *fractDone, unsigned long  *itDone ) 
{
  int aux = Step() ;
  if (fractDone) {
    *fractDone = FractDone() ;
    if (itDone)
      *itDone = itersDone ;
  }
  return aux ;
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void Iter::Set(  float **newparam, float **newEnd, float **newIterNo)
{
 parPtr = newparam;       endPtr = newEnd  ;           iterNo= newIterNo ;

 int i; 
 double aux ;
 for(i=1; i<=parNo; i++) {    // store starting points and step sizes
   start[i] = *parPtr[i] ;
   if ( (*iterNo[i]=FAbs(roundfl(*iterNo[i]))) > 1 )
    // steps are one less in number than iterations :  
    step[i] = (*endPtr[i] - start[i])/(*iterNo[i]-1) ;
   else
    step[i] = 0 ;             // Nonsense value for flagging etc.
 }
 
 epsilon = 666 ;  // will retain this if no stepping to do at all, 
                  // i.e. if 0 or 1 iterations to do overall.
 for(i=1; (*iterNo[i] <= 1 && i<parNo); i++) ;
 if (*iterNo[i] > 1) epsilon = DAbs(step[i]/4.0) ;
 while( i<parNo ) {   // loop to further adjust epsilon if nec.
   i++ ;
   if( *iterNo[i] > 1 && (aux=DAbs(step[i]/4.0)) < epsilon )
     epsilon = aux ;
 }

 itersDone = 0 ;
 TotN() ;       // find & ret. total No. of iterations & store in totN
 // No steps to do if the tot. N. of steps to do is 0, i.e. here if any
 // of the *stepNo[i] is 0. 
 if (totN > 0) { // the usual situation, where some iteration will be done !
   for (i=0; i<=parNo; i++)  stepIndex[i] = 1 ; // even sets the redundant i=0 
                                                // member to 1 
   doAgain = 1;
 } 
 else { // i.e. if totN signifies we don't iterate
   for (i=0; i<=parNo; i++)  stepIndex[i] = 0 ;   // cf. above
   doAgain =   0 ;         
 }

}
// ---------------------------------------------------------------------------
unsigned long int Iter::TotN()  // find, set & return total No. of iteration 
{  
  totN = 1;
  for (int i=1; i<=parNo ; i++ ) 
    totN *= rounddb( *iterNo[i] )  ; 
  return totN ;
}
// ---------------------------------------------------------------------------
void Iter::ParsToErr()
{
  int i, j ; 
  float x;
  static char totLen = 9 ; /* each no. will occupy this many spaces, 
			      including the gap. NB auxStr is a char [12] */
  // cerr<<"\n" ;
  for (i=1; i<=parNo; i++) {
    x = *parPtr[i] ;
    /* | | <= 0.0001 is auto-cov to 1e-04 ; | | >= 10^6 to 1e+06
       and e.g. "-3.56e-06' ''\0'" needs  10 chars, writes 9 spaces. */
    if (DAbs(x) <= 0.0001  ||  DAbs(x) >= 1000000 ) 
      gcvt( *parPtr[i], 3, auxStr ) ;
    else                                  // i.e. if has no e*XX  
      gcvt( *parPtr[i], 6, auxStr ) ;
    for (j = strlen( auxStr ); j<= totLen; j++ ) 
      auxStr[j] = ' ' ;
    auxStr[totLen+1] = '\0' ;
    cerr<< auxStr ;
  }
  cerr<<itersDone+1<<'/' <<totN<<'\t' ;
}
// ---------------------------------------------------------------------------
int Iter::AtStep( int parIndex )  // Where iteration is w.r.t. param parIndex
{  
  if ((parIndex < 1) || (parIndex > parNo)) 
      err_msg("Iter::AtStep : parIndex out of range in fn. AtStep " ) ;
  return  stepIndex[ parIndex ] ;
}

// ---------------------------------------------------------------------------
// 1. Code for piecewise linear maps.

void LinMapBasics::Set_pts(float **data)
  // note that the matrix used must start from 0 not 1 though 
  // derived with the routine in nr_ut.cpp
{
 int i;
 for(i=PieceNumber-1; i>=0; i--)
   {
    SetLimit(i,data[0][i]);
    SetPiecePars(i,data[1][i],data[2][i]);
   }
 SetLimit(PieceNumber,data[0][PieceNumber]) ; // last lim[] is a loner !
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void LinMapBasics::SetPcByEndpts(unsigned int PcNo, float y1, float y2)
{
 float m ;
 if (PcNo>PieceNumber) err_msg("LinMapBasics::SetPcByEndpts :\nPiece number out of range in SetPcByEndpts");

 Piece[PcNo].SetSlope( m=(y2-y1)/(Lim[PcNo+1]-Lim[PcNo]) );

 Piece[PcNo].SetConstant( y1-m*Lim[PcNo] );
}
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
float LinMapBasics::MapContent(int Row, int Col)
{
  float content=-1; // just to keep compiler from spewing up warnings.
  if ( (Col*(Col-PieceNumber) > 0) || ( (Row-2)*Row > 0) )
    err_msg("LinMapBasics::MapContent :\nRow or Col. wrong in MapContent(Row,Col)");
  else if (Row == 0)  content = Lim[Col];
  else if (Row == 2)  content = Piece[Col].GetConstant() ;
  else if (Row == 1)  content = Piece[Col].GetSlope() ;
  else
    err_msg("LinMapBasics::MapContent :\nMistake in looking up lin. map parametre");

  return content;
}
//--------------------------------------------------------------------------
// 6.2 code for linear sigmoid
LinSig::LinSig() : OneLinPcMap() {    
    SetLimit(0, 0.0); SetLimit(1, 1.0);
    SetPiecePars(0,1.0,0.0); 
    lVal = Piece[0].Val(Lim[0]); rVal =  Piece[0].Val(Lim[1]);           
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

//--------------------------------------------------------------------------
// 2. code for ramping function
void Ramp::SetParam(float tstart, float tend,    float tramp,
				  float yhigh,   float ylow    )
{
 if (tend-tstart<2*tramp) err_msg("Ramp::SetParam :\nParabRamp constr.: ramp too short");

 t1=tstart; t6=tend;  Dt=tramp;

 yhi=yhigh; 		ylo=ylow;
 ymid=(ylo+yhi)/2;      yAmp = (yhi-ylo)/2;

 t2=t1+Dt/2; t3=t1+Dt; t4=t6-Dt; t5=t6-Dt/2;

 a = 2*(yhi-ylo)/(Dt*Dt) ;

 W = PI/Dt ;
}

//--------------------------------------------------------------------------
// 7.1 code for Quadratic fit to 3 points function
Quad3Fit::Quad3Fit()
{
  A.resize(Range(1,3),Range(1,3)); B.resize(Range(1,3),Range(1,3)); 
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Quad3Fit::~Quad3Fit()
{
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
float Quad3Fit::Set(  float x1, float y1, float x2,  float y2
		    , float x3, float y3 )
{
 int i,j ;

 A(1,1)=x1*x1;		A(1,2)=x1; 	     //Basic
 A(2,1)=x2*x2;		A(2,2)=x2;	     //coef.
 A(3,1)=x3*x3;		A(3,2)=x3; 
 A(1,3)=A(2,3)=A(3,3)= 1.0 ;
 for(i=1; i<=3; i++)
   for(j=1; j<=3; j++) B(i,j)= A(i,j) ;  //B is copy.
 if( (D = Det3(A)) == 0.0 )
   return D ;                                //D as error flag.

 A(1,1)=y1;	  A(2,1)=y2;	 A(3,1)=y3;
 a = Det3(A)/D ;
 for(i=1; i<=3; i++)
    for(j=1; j<=3; j++)
       A(i,j)= B(i,j) ;                   // restore A.
 A(1,2)=y1; 	  A(2,2)=y2; 	 A(3,2)=y3;
 b = Det3(A)/D ;
 for(i=1; i<=3; i++)
    for(j=1; j<=3; j++)
       A(i,j)= B(i,j) ;                   // restore A.
 A(1,3)=y1; 	 A(2,3)=y2; 	 A(3,3)=y3;
 c = Det3(A)/D ;

 return D;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
float Quad3Fit::Coef2()
{
 if (D != 0.0) return a;
 else  err_msg("Quad3Fit::Coef2 :\n Det=zero in Coef2() " ) ; 
 return -1    ;  			// never reached.
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
float Quad3Fit::Coef1()
{
 if (D != 0.0) return b;
 else err_msg("Quad3Fit::Coef1 :\n Det=zero in Coef1() ") ;
 return -1    ;  			// never reached.
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
float Quad3Fit::Coef0()
{
 if (D != 0.0) return c;
 else err_msg("Quad3Fit::Coef0 :\n Det=zero in Coef0() " ) ;
 return -1    ;  			// never reached.
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
float Quad3Fit::OrdStat()
{
 if (D != 0.0) return (c - b*b/(4*a)) ;
 else err_msg( "Quad3Fit::OrdStat :\n Det=zero in OrdStat() ") ;
 return -1    ;  			// never reached.
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
float Quad3Fit::AbcStat()
{
 if (D != 0.0) return ( -.5 * b/a ) ;
 else err_msg("Quad3Fit::AbcStat :\n Det=zero in AbcStat() ") ;
 return -1    ;  			// never reached.
}
//--------------------------------------------------------------------------
// 8.1 code for Silnikov systems :
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// 2. copy ctor
Silnikov::Silnikov( const Silnikov & rval)
{
  s = new float[4] ;     
  a=rval.a;                       b =rval.b;                 c =rval.c; 
  xOffset  =rval.xOffset;   yOffset =rval.yOffset;    zOffset  =rval.zOffset; 
  timeScale=rval.timeScale;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void Silnikov::SetAux(  double newTScale, double newXOff
                      , double newYOff, double newZOff)
{
  timeScale = newTScale; 
  xOffset = newXOff;      yOffset = newYOff;   zOffset = newZOff ;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void Silnikov::RandomSVar()
{
  s[1]=SMALLFLTSIGN + /* noiseInDisp */ .25*rand()/RAND_MAX ;
  for (int i=2; i<4; i++)
    s[i] = 0.1+SMALLFLTSIGN + /* noiseInDisp */ .25*rand()/RAND_MAX ;
  s[4] = timeScale*(s[2]+yOffset);  // as per def. of XDotWAux
}

//--------------------------------------------------------------------------
// 8.2 code for linear oscillatory systems and derivatives :
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
HO::HO()
{ 
 x= new float[3];  
 *(x+2)= 1. ;           *(x+1)= 0. ;       //'kicked' init cond
 
 inert = 1. ;  resis = 1. ;   resto = 1. ; 
 force = 0. ;  W     = 2. ;   phi   = 0. ;                 
 
 W0SQ = resto/inert;   F= force/inert;  b= .5*resis/inert ;

}  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// 2. copy ctor
HO::HO(const HO &rval )
{ 
 x= new float[3];  
 *(x+2)= rval.x[2] ;        *(x+1)= rval.x[1] ;       // state ;
 
 inert = rval.inert ;  resis = rval.resis ;   resto = rval.resto ; 
 force = rval.force ;     W  = rval.W ;         phi = rval.phi ;       
 W0SQ  = rval.W0SQ;       F  = rval.F ;           b = rval.b ;

}  

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void HO::SetDHOPar( float new_resto, float new_resis, float new_inert )
{ 
  if ((inert = new_inert) < 0)
    err_msg("HO::SetDHOPar :\n inert<0 in SetDHOPar" ) ; 
  if (( resis = new_resis) < 0)
    err_msg("HO::SetDHOPar :\n resis<0 in SetDHOPar" ) ;  
 if ((  resto = new_resto) < 0)
   err_msg("HO::SetDHOPar :\n resto<0 in SetDHOPar" ) ;  
  W0SQ = resto/inert;     F = force/inert;      b = .5*resis/inert ;
}  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void HO::SetForcPar( float newforce, float newphi, float newW )
{ 
 force = newforce ;  phi = newphi ;   W = newW ; 
  F= force/inert;
}  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void HO::GetDHOPar(float *restor, float *resist, float *inerti, 
                   float *forcin, float *Wforci, float *phiforc )
{ 
 if (restor)  *restor  = resto ;
 if (resist)  *resist  = resis ;
 if (inerti)  *inerti  = inert ;
 if (forcin)  *forcin  = force ;
 if (Wforci)  *Wforci  = W     ;
 if (phiforc) *phiforc = phi   ;
}  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int HO::OutPar( float *Wres, float *Tdec1, float *Tdec2 )
//    Returns underdamped->-1, overdamped or crit. damped-> 1
{
  double discSQ = b*b - W0SQ ;
  if (discSQ < 0) {                  // ie oscill. soln. exists
    if (Tdec1) *Tdec1 = 1/b ;
    if (Tdec2) *Tdec2 = 0.0 ;
    if (Wres ) *Wres  = sqrt( -1*discSQ ) ;
    return -1;
  }
  else { 
    if (Tdec1) *Tdec1 = 1.0/(b + sqrt( discSQ )) ;
    if (Tdec2) *Tdec2 = 1.0/(b - sqrt( discSQ )) ;
    if (Wres ) *Wres  = 0.0 ;
    return 1;
  }
  
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void HO::ForcedAmpPh( float *amp, float *phase)
{ 
  if (!amp || !phase)
   err_msg("HO::ForcedAmpPh :\n ForcedAmpPh asked to write to NULL" ) ; 
  
  float delWSQ = W0SQ-W*W ;
  *amp   =  F / sqrt( delWSQ*delWSQ + 4*b*b*W*W) ;
  *phase =  phi - atan2( delWSQ, -2*b*W ) ;
}  

//--------------------------------------------------------------------------

//                   ***** GEN. PURPOSE C FUNCTIONS *****                  

// ---------------------------------------------------------------------------
// Rounds to NEAREST integer

int roundfl( float x )
 {
  if( x <= -32768.5 || x > 32767.5 ) err_msg("int roundfl(x): x out of range.");
  const float midpoint = 0.5 ;

  if ((x>0) && ( x-((int) x) > midpoint))  return ((int) x)+1 ;
  else if ((x<0) && ( ((int) x)-x >= midpoint))  return ((int) x)-1 ;
  else return (int) x ;
 }
// ---------------------------------------------------------------------------
// Rounds to NEAREST LONG integer

long int rounddb( double x )
 {
  if( x <= -2147483648.5 || x > 2147483647.5 ) 
    err_msg("long int rounddb(x): x out of range.");
  const float midpoint = 0.5 ;

  if ((x>0) && ( x-((int) x) > midpoint))  return ((long int) x)+1 ;
  else if ((x<0) && ( ((int) x)-x >= midpoint))  return ((long int) x)-1 ;
  else return (long int) x ;
 }
// ---------------------------------------------------------------------------
float FAbs( float x)
{
  return ( x>0 ? x : -x ) ;
}
// ---------------------------------------------------------------------------
double DAbs( double x)
{
  return ( x>0 ? x : -x ) ;
}
// ---------------------------------------------------------------------------
// floatArrayMAX & floatArrayMIN find the maximum or minimum
// element of an array of AT LEAST TWO reals . They assume the 
// data is found between variable[firstIndex] & variable[ArPosn]
// (as ArPosn initially given) . 

float floatArrayMAX(float variable[], int ArPosn, int firstIndex)
 {
  float largestSoFar ;

  for  ( largestSoFar = variable[ArPosn] ; ArPosn >= firstIndex ; --ArPosn)
   if (variable[ArPosn] > largestSoFar ) largestSoFar = variable[ArPosn] ;

  return largestSoFar  ;
 }
// --------------------------------------------------------------------------
float floatArrayMIN(float variable[], int ArPosn, int firstIndex)
 {
  float leastSoFar ;

  for  ( leastSoFar = variable[ArPosn] ; ArPosn >= firstIndex ; --ArPosn)
   if (variable[ArPosn] < leastSoFar ) leastSoFar = variable[ArPosn] ;

  return leastSoFar  ;
 }
// --------------------------------------------------------------------------
float floatArrayAVG(float variable[], int ArPosn, int firstIndex)
 {
  int i ;
  float Sum = 0.0 ;
  int N = ArPosn - firstIndex + 1 ;
  for (i=firstIndex; i<= ArPosn ; i++ ) Sum += variable[i] ;

  return Sum/N ;
 }
// --------------------------------------------------------------------------
float floatArraySD(float variable[], int ArPosn, int firstIndex)
 {
  int i ;
  float SumSQ = 0.0 ;
  float mean  = floatArrayAVG( variable, ArPosn, firstIndex ) ;
  int N_1     = ArPosn - firstIndex ;
  for (i=firstIndex; i<= ArPosn ; i++ ) SumSQ += pow( (mean - variable[i]), 2) ;

  return sqrt( SumSQ/N_1 ) ;
 }
// --------------------------------------------------------------------------
float floatArrayDev(float variable[], int ArPosn, int firstIndex
                    , float *Dev )
 {
  int i ;   float SumSQ = 0.0 ;   int N_1 = ArPosn - firstIndex ;
  float ref = (Dev != NULL) ? (*Dev) : (floatArrayAVG( variable, ArPosn, firstIndex )) ;
  for (i=firstIndex; i<= ArPosn ; i++ ) SumSQ += pow( (ref - variable[i]), 2) ;
  return sqrt( SumSQ/N_1 ) ;
 }
// -------------------------------------------------------------------------
float LinInter(const float *y, const float *x, const float xNew)
 {
  return *y+(xNew-*x)*(*(y+1)-*y)/(*(x+1)-*x) ;
 }
// -------------------------------------------------------------------------
int  ArrLinInter(float **yInt,
                 float **yDat, 
                 const float xStart, const float xEnd, const int   Nx,  
                 const int   indEnd, const int   indStart, 
		 const int   depDatDim )
 {
   // yInt[0][0] is at the start of the output array. yDat[0] points to the 
   // parametre / indep. var., yDat[1 to depDatDim] to the dep. vars / ordinates. 
   // yDat[0] must be in ascending order.
   
   if  (yDat[0][indStart]>xStart || yDat[0][indEnd]<xEnd || xStart>xEnd ) 
     err_msg("ArrLinInter input boundaries" );

   int N = 0;           // tracks output array posn.
   int Nd = indStart;   // tracks data (input) array posn.
   float X=xStart;    float dX = (xEnd - xStart)/(Nx-1) ;

   while (N<Nx) { // rem Nx is the inputed no. of points needed
     while ( *(*yDat+Nd) <X) Nd++ ;  // jump over first ordinate point to use
     Nd -= 1;                           // ... and return to just below it !
     for (int i=0; i<depDatDim; i++)
       yInt[i][N] = LinInter( *(yDat+i+1)+Nd, *yDat+Nd, X) ;
     X += dX ;  
     N += 1 ;
   }
   return 0;
 }
// -------------------------------------------------------------------------
double TrapInt(  const float *y, const float *x
               , const float xStart, const float xEnd
	       , int   lastInd, int firstInd )
 {
  double Sum, x1, x2, y1, y2 ; 
  x1 = xStart; 
  int i=firstInd;

  if( (x[firstInd]-xStart)*(x[lastInd]-xStart) > 0  ||
      (x[firstInd]-xEnd  )*(x[lastInd]-xEnd  ) > 0      )
    err_msg("TrapInt: x boundaries");
  if( ( x[lastInd] - x[firstInd] ) <= 0  ||
      ( xEnd - xStart <=0        )           ) //Can't cope with backwards...
    err_msg("TrapInt: ind. Var. order") ;      //... indep. variable.

  while( (x[i]-xStart)*(x[i+1]-xStart) > 0 ) i++ ;
  y1 = LinInter( y+i, x+i, x1 );

  Sum = (x2=x[i+1])*y1-x1*(y1+(y2=y[i+1])) ;
  i++ ;
  x1 =  x2 ;
  y1 =  y2 ;

  while( x[i+1] < xEnd) {
     Sum += (x2=x[i+1])*y1-x1*(y2=y[i+1]) ;
     x1 =  x2 ;
     y1 =  y2 ;
     i++ ;
  }
  y2   = LinInter( y+i, x+i, (x2=xEnd) ) ;
  Sum += x2*y1+x2*y2-x1*y2 ;

  return Sum/2 ;
 }  

// ------------------------------------------------------------
// Simple determinant of 3x3 matrix using C array (dangerous!)
float Det3( float **a )
{
 return   a[1][1] * ( a[2][2]*a[3][3] -  a[2][3]*a[3][2])
	+ a[1][2] * ( a[2][3]*a[3][1] -  a[2][1]*a[3][3])
	+ a[1][3] * ( a[2][1]*a[3][2] -  a[2][2]*a[3][1])  ;
}
// ----------------------------------------------------------------------------
// Simple determinant of 3x3 matrix using blitz++ arrays (better):
float Det3(Array<float, 2> A)
{
  if ((A.extent(firstDim) !=3) || (A.extent(secondDim) !=3))
    err_msg("Det3(Array<float, 2> A) :\n A has to be 3x3 ...");
  int ro, co; // row offset and column offset.
  ro = A.lbound(firstDim)-1;  co=A.lbound(secondDim)-1;
  
  return A(ro+1,co+1) * ( A(ro+2,co+2)*A(ro+3,co+3) - A(ro+2,co+3)*A(ro+3,co+2))
       + A(ro+1,co+2) * ( A(ro+2,co+3)*A(ro+3,co+1) - A(ro+2,co+1)*A(ro+3,co+3))
       + A(ro+1,co+3) * ( A(ro+2,co+1)*A(ro+3,co+2) - A(ro+2,co+2)*A(ro+3,co+1));
}
// ------------------------------------------------------------
// 3e ii. mean sqrt of [U(T) - V(T)]^2 over tStart to tEnd
double rms_tser(const float *U,     const float *Tu, int Nu,
                const float *V,     const float *Tv, int Nv,
                const float tStart, const float tEnd        )
{
 if (   (( tStart < Tu[1]) || ( tStart < Tv[1] )) 
     || (( tEnd   > Tu[Nu])|| ( tEnd   > Tv[Nv])) )
   err_msg( "rms_tser: problem with boundaries" ) ;
  
 double   SqDiffInt;            SqDiffInt = 0.0 ;
 double   xlo, xhi ;            xlo = tStart ;
 int      Nlo, Nhi, Mlo, Mhi ;  Nlo = Nhi = Mlo = Mhi = 1 ;

 while ( Tu[Nhi] < tStart ) Nhi++ ;  Nlo = Nhi-1 ;
 while ( Tv[Mhi] < tStart ) Mhi++ ;  Mlo = Mhi-1 ;

 xhi = (Tu[Nhi]  >  Tv[Mhi]) ?  Tv[Mhi] : Tu[Nhi] ;

 while (xhi < tEnd) {
   SqDiffInt += (xhi-xlo) * pow(
                    LinInter(U+Nlo, Tu+Nlo, (xhi+xlo)/2.0)
		   -LinInter(V+Mlo, Tv+Mlo, (xhi+xlo)/2.0), 2) ;

   xlo = xhi ; // first boundary of next interval considered
   
   (xhi < Tu[Nhi]) ? Mlo=Mhi++ : Nlo=Nhi++ ;
   // equivalent to: if xhi equals Tv[Mhi], advance the Ms etc.

   xhi = ( (Tv[Mhi] < Tu[Nhi]) ? Tv[Mhi] : Tu[Nhi] ) ;
   // set xhi to the lowest of the Tu or Tv current upper vals.

 } // at exit of this loop, (xlo,xhi) must cont. tEnd
 
 xhi = tEnd ;
 SqDiffInt += (xhi-xlo) * pow(
                    LinInter(U+Nlo, Tu+Nlo, (xhi+xlo)/2.0)
		   -LinInter(V+Mlo, Tv+Mlo, (xhi+xlo)/2.0), 2) ;
 
 return sqrt( SqDiffInt/(tEnd-tStart) ) ;
}     

// ------------------------------------------------------------
// 3f ii. take an angle in rad, which can be +ve or -ve, reduce it from 0 to 
//        2*PI and optionally return the no. of full rotns. in the angle :
float phi0to2PI( const float phi, int *kappa )
{
  int k=0 ;
  float theta = phi ;
  double  pi2 = 2*PI  ;
  if (theta >= pi2) {
    while (theta >= pi2) {
      theta -= pi2;
      k += 1 ;
    }
  }
  if (theta < 0.0)
    k = 1 ;  // so that first, duff reduction below gets it O.K. 
    while (theta < 0.0 ) {
    theta += pi2 ;
    k -= 1 ;
  }
  if (kappa) *kappa = k ;
  return theta ;
}
// ------------------------------------------------------------
/* iii. similar to above, but ret. within -PI to PI : */ 
float phiMinusPItoPI( const float phi, int *kappa /* = NULL */ ) {
  int k=0 ;
  float theta = phi0to2PI( phi, &k) ;
  if (theta > PI) theta -= 2*PI ;
  if (kappa) *kappa = k;
  return theta ;
}

// ------------------------------------------------------------

/* 3g - Random-based functions ASSUME srand has been run, e.g. by ... */ 
/*      ...  srand( ( time( NULL ) % 604800) ) ;                      */

/* i. bin_ded returns 1 (accept) or 0 (reject) decision based on inputed  */
/*    probability P:                                                      */
char bin_dec(double P) {
    if (P*(P-1) >  0)
	err_msg("bin_dec(double P) :  needs P from 0 to 1 !");
    else
	return (1.0*rand()/RAND_MAX <= P) ? 1:0 ;
    return -1 ;   // never reached
}
// ------------------------------------------------------------

/* 3i - filename string utility */
int fnstring(char *wholename, char *stem, int ind, char *ext)
{
  if ((ind<0) || (ind>999)) err_msg("fnstring: ind <0 or >999 is no good");
  if (ind<10)
    return sprintf(wholename, "%s00%d.%s",stem,ind,ext );
  else if (ind<100)
    return sprintf(wholename, "%s0%d.%s",stem,ind,ext );
  else 
    return sprintf(wholename, "%s%d.%s",stem,ind,ext );
}

/* 9.1 Gibbs distribution   ------------------------------------------------- */
/*                                           -------------------------------- */
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
double Gibbs::Len0()
{
    double Sum_j=0.0; 
    for(int i=E.lbound(firstDim); i<=E.ubound(firstDim); i++) 
      Sum_j += exp( (E(i))/T ) ;

   return Sum_j ;
}
//  10.x - Receiver Operating Characteristic (ROC) ... ----------------------
//                                           --------------------------------

// ROC made of one hyperbolic piece ----------------------------------------  
HypROC::HypROC(): ROC()  {
  SetPar(4.0) ;     // a good-ish test as default.
  set_costs(10.0,2.0,0.0,0.0); // TPC zero isn't very realistic ... 
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
void HypROC::SetPar(double New_M)
{
  if (New_M <= 1.0)
    err_msg("ROC::SetPar :\n HypROC par. m must be > 1.0" ) ;
  m = New_M ;
  b = 1-m*m ;
  a = b*b/(m*m) ;
  c = -m*m/b ; 
  
}

//  ------------------- End of file ! -------------------------




