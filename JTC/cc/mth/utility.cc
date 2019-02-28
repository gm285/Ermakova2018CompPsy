//   ====================< UTILITY.CPP >========================
//   * File containing the basic utilities                     *
//   * Description: Chapter 8                                  *
//   * Scientific C++ Building Numerical Libraries             *
//   *                the Object-Oriented Way                  *
//   * G. Buzzi-Ferraris, Addison-Wesley (1993)                *
//   ===========================================================

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>

#include "utility.h"

//   ===========================================================
//   ======================   GLOBAL   =========================
//   ===========================================================
FILE *bzzFileOut = stdout;
FILE *bzzFileIn = stdin;
char bzzYesNo = 'Y';
const float MACH_EPS = MachEps();     

const float SMALLFLTSIG = 0.00001 ; // Near 0 w.r.t sig. places; for comparisons

const float BIG_FLOAT =                        3.4e+38;
const float TINY_FLOAT =                       3.4e-38;
const double BIG_DOUBLE =                      1.7e+308;
const double TINY_DOUBLE =                     8.e-307;
const long double BIG_LONG_DOUBLE =            3.4e+4932;
const long double TINY_LONG_DOUBLE =           3.4e-4932;
const char BIG_CHAR =                          127;
const char TINY_CHAR =                         -128;
const unsigned char BIG_UNSIGNED_CHAR =        255;
const unsigned char TINY_UNSIGNED_CHAR =       0;
const int BIG_INT =                            32767;
const int TINY_INT =                           -32768;
const unsigned int BIG_UNSIGNED_INT =          65535;
const unsigned int TINY_UNSIGNED_INT =         0;
const long int BIG_LONG =                      2147483647;
const long int TINY_LONG =                     -2147483648;
const unsigned long int BIG_UNSIGNED_LONG  =   4294967295;
const unsigned long int TINY_UNSIGNED_LONG =   0;

//   ===========================================================
//   ======================   ERRORS   =========================
//   ===========================================================


/*
//========== Error Type ==============
const char *const ERR_RANGE=
      "\nIndex out of range\n";
const char *const ERR_SPACE=
      "\nNo room in heap\n";
const char *const ERR_OPEN_FILE=
      "\nCan't open FILE\n";
const char *const ERR_CLOSE_FILE=
      "\nCan't close FILE\n";
const char *const ERR_READING_FILE=
      "\nCan't read FILE\n";
const char *const ERR_WRITING_FILE=
      "\nCan't write FILE\n";
const char *const ERR_CHECK_DIMENSION=
      "\nDimension check failure\n";
const char *const ERR_OUT_OF_SCOPE=
      "\nOut of Scope\n";
const char *const ERR_FACTORIZED=
      "\nFactorized Matrix\n";
const char *const ERR_IMPLEMENTATION=
      "\nFunction not implemented\n";

//============ Function Type =============
const char *const ERR_FUNCTION="\nFunction: ";
const char *const ERR_CONSTRUCTOR="\nConstructor: ";
const char *const ERR_OPERATOR="\nOperator: ";

*/

//========== Error Type ==============
const char *const ERR_RANGE=
      "\nRange!\n";
const char *const ERR_SPACE=
      "\nHeap!\n";
const char *const ERR_OPEN_FILE=
      "\nOpen FILE\n";
const char *const ERR_CLOSE_FILE=
      "\nClose FILE\n";
const char *const ERR_READING_FILE=
      "\nRead FILE\n";
const char *const ERR_WRITING_FILE=
      "\nWrite FILE\n";
const char *const ERR_CHECK_DIMENSION=
      "\nDim.chk.failure\n";
const char *const ERR_OUT_OF_SCOPE=
      "\nScope!\n";
const char *const ERR_FACTORIZED=
      "\nFact.Mtx\n";
const char *const ERR_IMPLEMENTATION=
      "\nFn implementn\n";

//============ Function Type =============
const char *const ERR_FUNCTION="\nFn: ";
const char *const ERR_CONSTRUCTOR="\nConst:";
const char *const ERR_OPERATOR="\nOpr:";

//   ===========================================================
//   ====================   FUNCTIONS   ========================
//   ===========================================================

//   **********************< MachEps >**************************
//   * Purpose: Calculating the precision of the machine       *
//   * Description: It initialises the const MACH_EPS which    *
//   *              can be used  anywhere in the program       *
//   * Example: if( x < MACH_EPS)                              *
//   ***********************************************************
float MachEps(void)
   {
   float macheps = 1.,eps = 2.;
   while(eps != 1.)
      {
      macheps /= 2.;
      eps = 1. + macheps; // in simple rounded form
      }
   return macheps*2.;
   }

//   **********************< Message >**************************
//   * Purpose: Sending messages to the FILE bzzFileOut        *
//   * Description: For handling in the same manner as printf. *
//   *              It does not interrupt the program          *
//   *              To disable: bzzYesNo = 'N'                 *
//   *              To enable: bzzYesNO = 'Y' (default)        *
//   * Example: Message("The value of i = %d",i);              *
//   ***********************************************************
void Message(const char *myFormat,...)
{
if(bzzYesNo != 'Y')return;
va_list argPointer;
va_start(argPointer,myFormat);
vfprintf(bzzFileOut,myFormat,argPointer);
va_end(argPointer);
}

//   ***********************< Error >***************************
//   * Purpose: Sending error messages to the file stderr      *
//   * Description: As with Message, but it interrupts         *
//   *              the program                                *
//   * Example: if(x < 1) Error("x = %f",x);                   *
//   ***********************************************************
void Error(const char *myFormat,...)
{
va_list argPointer;
va_start(argPointer,myFormat);
vfprintf(stderr,myFormat,argPointer);
va_end(argPointer);
exit(1);
}

//   **********************< Control >**************************
//   * Purpose: For controlling overflow                       *
//   * Description: It is used both for float and for double   *
//   * Example: float x = Control(y);                          *
//   ***********************************************************
float Control(double value)
   {
   if(value > BIG_FLOAT)return BIG_FLOAT;
   else if(value < -BIG_FLOAT)return -BIG_FLOAT;
   return value;
   }

double Control(long double value)
   {
   if(value > BIG_DOUBLE)return BIG_DOUBLE;
   else if(value < -BIG_DOUBLE)return -BIG_DOUBLE;
   return value;
   }

//   ***********************< Swap >****************************
//   * Purpose: Having a unique function for Swap              *
//   * Description: Overload for float, int, double,           *
//   *              *float, Vector, Matrix,  etc.              *
//   * Example: Swap(&x,&y);                                   *
//   ***********************************************************
// Swap float
void Swap(float *x,float *y)
    {float temp = *x; *x = *y; *y = temp;}

// Swap double
void Swap(double *x,double *y)
    {double temp = *x; *x = *y; *y = temp;}

// Swap int
void Swap(int *x,int *y)
    {int temp = *x; *x = *y; *y = temp;}

// Swap float pointers
// Swap(&p1,&p2);
void Swap(float **x,float **y)
   {float *temp = *x; *x = *y; *y = temp;}

// Swap double pointers
// Swap(&p1,&p2);
void Swap(double **x,double **y)
   {double *temp = *x; *x = *y; *y = temp;}

// Swap int pointers
// Swap(&p1,&p2);
void Swap(int **x,int **y)
   {int *temp = *x; *x = *y; *y = temp;}

//   ************************< Max >****************************
//   * Purpose: Max of an array                                *
//   * Description: It returns the Max of an array.            *
//   *              If the index im is put in the argument     *
//   *              it also returns the position.              *
//   * Example: y = Max(10,x,&im);                             *
//   ***********************************************************
float Max(int n,float *x,int *im)
   {
   if(n < 0) return x[0];
   float temp = x[0];
   if(im != 0) *im = 0;
   for(int i = 1;i < n;i++)
     if(temp < x[i])
       { temp = x[i]; if(im != 0) *im = i; }
   return temp;
   }

//   **********************< MaxAbs >***************************
//   * Purpose: Max Abs of an array                            *
//   * Example: y = MaxAbs(10,x,&im);                          *
//   ***********************************************************
float MaxAbs(int n,float *x,int *im)
   {
   float temp = Abs(x[0]);
   if(n < 0) return temp;
   if(im != 0) *im = 0;
   for(int i = 1;i < n;i++)
     if(temp < Abs(x[i]))
       { temp = Abs(x[i]); if(im != 0) *im = i; }
   return temp;
   }

//   ************************< Min >****************************
//   * Purpose: Min of an array                                *
//   * Example: y = Min(10,x,&im);                             *
//   ***********************************************************
float Min(int n,float *x,int *im)
   {
   if(n < 0) return x[0];
   float temp = x[0];
   if(im != 0) *im = 0;
   for(int i = 1;i < n;i++)
     if(temp > x[i])
       { temp = x[i]; if(im != 0) *im = i; }
   return temp;
   }

//   **********************< MinAbs >***************************
//   * Purpose: Min abs of an array                            *
//   * Example: y = MinAbs(10,x,&im);                           *
//   ***********************************************************
float MinAbs(int n,float *x,int *im)
   {
   float temp = Abs(x[0]);
   if(n < 0) return temp;
   if(im != 0) *im = 0;
   for(int i = 1;i < n;i++)
     if(temp > Abs(x[i]))
       { temp = Abs(x[i]); if(im != 0) *im = i; }
   return temp;
   }

//   ************************< Sum >****************************
//   * Purpose: The sum of two vectors                         *
//   * Description: The sum of two arrays with control         *
//   * Example: Sum(n,x,y,z); z = x + y;                       *
//   ***********************************************************
void Sum(int n,float *lval,float *rval,float *result)
   {
   double sum;
   for(int i=0;i < n;i++)
       {
       sum = (*lval++) + (*rval++);
       *result++ = Control(sum);
       }
   }

//   **************************< Sum >**************************
//   * Purpose: The sum of two vectors                         *
//   * Description: The sum substitutes the first vector       *
//   * Example: Sum(n,x,y) x = x + y                           *
//   ***********************************************************
void Sum(int n,float *lvalAndResult,float *rval)
   {
   double sum;
   for(int i=0;i < n;i++)
       {
       sum = (*lvalAndResult) + (*rval++);
       *lvalAndResult++ = Control(sum);
       }
   }

//   **************************< Sum >**************************
//   * Purpose: The sum of two equal vectors                   *
//   * Description: The sum substitutes the vector             *
//   * Example: Sum(n,x) x = x + x                             *
//   ***********************************************************
void Sum(int n,float *lvalRvalAndResult)
   {
   double sum;
   for(int i=0;i < n;i++)
       {
       sum = (*lvalRvalAndResult) + (*lvalRvalAndResult);
       *lvalRvalAndResult++ = Control(sum);
       }
   }

//   *********************< Difference >************************
//   * Purpose: The difference between two vectors             *
//   * Description: The difference between two arrays          *
//   *              with control                               *
//   * Example: Difference(n,x,y,z); z = x - y;                *
//   ***********************************************************
void Difference(int n,float *lval,float *rval,float *result)
   {
   double diff;
   for(int i=0;i < n;i++)
       {
       diff = (*lval++) - (*rval++);
       *result++ = Control(diff);
       }
   }

//   *********************< Difference >************************
//   * Purpose: The difference between two vectors             *
//   * Description: The result substitutes                     *
//   *              the first vector                           *
//   * Example: Difference(n,x,y) x = x - y                    *
//   ***********************************************************
void Difference(int n,float *lvalAndResult,float *rval)
   {
   double diff;
   for(int i=0;i < n;i++)
       {
       diff = (*lvalAndResult) - (*rval++);
       *lvalAndResult++ = Control(diff);
       }
   }

//   ******************< Difference >***************************
//   * Purpose: The difference between two vectors             *
//   * Description: The result substitutes                     *
//   *              the second vector                          *
//   * Example: Difference(n,x,y,1) y = x - y                  *
//   ***********************************************************
void Difference(int n,float *lval,float *rvalAndResult,int)
   {
   double diff;
   for(int i=0;i < n;i++)
       {
       diff = (*lval++) - (*rvalAndResult);
       *rvalAndResult++ = Control(diff);
       }
   }

//   ************************< Dot >****************************
//   * Purpose: The scalar (dot) product of two vectors        *
//   * Description: The sum of the products of the coefficients*
//   * Example: float dot = Dot(n,x,y);                        *
//   ***********************************************************
float Dot(int n,float *lval,float *rval)
   {
   double result = 0.;
   for(int i=0;i < n;i++)
       result += (*lval++) * (*rval++);
   return Control(result);
   }

//   **********************< Product >**************************
//   *Purpose: The product of a float with a vector            *
//   *Description: It multiplies n terms of the vector         *
//   *             v by c, resulting in r                      *
//   *Example: Product(n,c,v,r)                                *
//   ***********************************************************
void Product(int n,float lval,float *rval,float *result)
   {
   for(int i=0;i < n;i++)
       *result++ = Control(lval * (*rval++));
   }

//   **********************< Product>***************************
//   *Purpose: The product of a float with a vector            *
//   *Description: It multiplies n terms of the vector         *
//   *             v by c, substituting in v                   *
//   *Example: Product(n,c,v)                                  *
//   ***********************************************************
void Product(int n,float lval,float *rvalAndResult)
   {
   for(int i=0;i < n;i++)
       *rvalAndResult++ = 
         Control(lval * (*rvalAndResult));
   }

//   *********************< Division >**************************
//   *Purpose: Dividing a vector by a float                    *
//   *Description: It divides n terms of the vector v by c,    *
//   *               resulting in r                            *
//   *Example: Division(n,v,c,r)                               *
//   ***********************************************************
void Division(int n,float *lval,float rval,float *result)
   {
   if(rval == 0.)rval = TINY_FLOAT;
   for(int i=0;i < n;i++)
       *result++ = Control((*lval++)/rval);
   }

//   **********************< Division >*************************
//   *Purpose: Dividing a vector by a float                    *
//   *Description: It divides n terms of the vector v by c,    *
//   *             substituting in v                           *
//   *Example: Division(n,c,v)                                 *
//   ***********************************************************
void Division(int n,float *lval,float rval)
   {
   if(rval == 0.)rval = TINY_FLOAT;
   for(int i=0;i < n;i++)
       *lval++ = Control((*lval)/rval);
   }

//   ********************< SqrtSumSqr >*************************
//   * Purpose: The Euclidean Norm calculation for a float     *
//   * Description: Chapter 8                                  *
//   * Example: float x[100]; x[0]=.....                       *
//   *          float norm = SqrtSumSqr(n,x);                  *
//   ***********************************************************
float SqrtSumSqr(int n,float *x)
   {
   if(n <= 0)
      Error("%s%sSqrtSumSqr",ERR_RANGE,ERR_FUNCTION);
   double norm = 0.,aux;
   for(int i = 0;i < n;i++)
      {
      aux = x[i];
      norm += aux*aux;
      }
      norm=sqrt(norm);
      if(norm > BIG_FLOAT)norm = BIG_FLOAT;
   return norm;
   }

//   ********************< SqrtSumSqr >*************************
//   * Purpose: The Euclidean Norm for a double                *
//   * Description: Chapter 8                                  *
//   * Example: double x[100]; x[0]=.....                      *
//   *          norm = SqrtSumSqr(n,x);                        *
//   ***********************************************************
double SqrtSumSqr(int n,double *x)
   {
   if(n <= 0)
      Error("%s%sSqrtSumSqr",ERR_RANGE,ERR_FUNCTION);
   double aux, xmax = 0.,xmin = BIG_DOUBLE;
   for(int j = 0;j < n;j++)
      {
      aux = Abs(x[j]);
      if(xmax < aux)xmax = aux;
      if(xmin > aux)xmin = aux;
      }
   if(xmax == 0.)return xmax;
   if (xmin == 0.)xmin = TINY_DOUBLE;
   long double longaux = 
      (long double)xmax/(long double)xmin;
   aux = sqrt(BIG_DOUBLE/((double)n));
   // to avoid the problems of
   if(xmax < aux &&                            // overflow
      xmax > TINY_DOUBLE/MACH_EPS &&       // small numbers
      longaux < 1./MACH_EPS)                     // sort
      {
      double norm = 0.;  // without problems: double
      for(int i = 0;i < n;i++)
        {
        aux = x[i];
        norm += aux*aux;
        }
      return sqrt(norm);
      }
   else  // if there are problems it works in long double
      { 
      long double norm = 0.;
      for(int i = 0;i < n;i++)
        {
        longaux = x[i];
        norm += longaux*longaux;
        }
      if(norm < BIG_DOUBLE && norm > TINY_DOUBLE)
         return sqrt(norm);
      longaux = (long double)xmax*(long double)n;
      norm /= longaux; // avoids overflow
      norm /= longaux;
      norm = longaux*sqrt(norm); // renormalises
           // avoids overflow
      if(norm > BIG_DOUBLE) norm = BIG_DOUBLE; 
      return norm;
      }
   }

//   ***********************< Sort >****************************
//   * Purpose: To sort an array of floats                     *
//   * Description: Chapter 5                                  *
//   * Example:                                                *
//   *        float x[5]={3.,2.,5.,1.,4.};                     *
//   *        Sort(5,x);                                       *
//   ***********************************************************
void Sort(int n,float *x)
   {
   if(n <= 1)return;
   int node,i,j,k,ik,jk;
   for(node = 1;node < n;node++)
      {
      i = node;
      j = ((i+1)/2)-1;
      while(i != 0 && x[j] <= x[i])
        {
        Swap(x+j,x+i);
        i=j;
        j=((i+1)/2)-1;
        }
      }
   for(i = n-1;i >= 1;i--)
      {
      Swap(x+i,x);
      k = i-1;
      ik = 0;
      jk = 1;
      if(k >= 2 && x[2] > x[1])jk=2;
      while(jk <= k && x[jk] > x[ik])
        {
        Swap(x+jk,x+ik);
        ik = jk;
        jk = (2*(ik+1))-1;
        if(jk+1 <= k)
        if(x[jk+1] > x[jk])jk++;
        }
      }
   }

/* 'TempFile'below commented out as has various problems (see comments 
   in utility.h ...

//   ===========================================================
//   =====================   TempFile   ========================
//   *          Class for creating temporary files             *
//   ***********************************************************

int TempFile::countTempFile = 0;

//   ********************< NewFileName >************************
//   * Purpose: Initialisation when the default                *
//   *          constructor has been used.                     *
//   * Description: Chapter 10                                 *
//   * Example: TempFile ff[10];                               *
//   *               ff[0].NewFileName("D:");                  *
//   ***********************************************************
void TempFile::NewFileName(char *directoryTemp)
   {
   nameTempFile = new char [strlen(directoryTemp)+10];
   strcpy(nameTempFile,directoryTemp);
   countTempFile++;
   strcat(nameTempFile,"T");
   char buff[5];
   itoa(countTempFile,buff,10);  // NON-PORTABLE : NOT ANSI !!!
   strcat(nameTempFile,buff);
   strcat(nameTempFile,".TMP");
   }

*/ 

// -------------------- end of utility.cc -------------------

