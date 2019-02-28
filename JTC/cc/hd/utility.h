/*   =================< Adapted from UTILITY.HPP >==============
    Ch. 8 of Scientific C++. 
    Has v. useful ut.s for temp. files which don't work yet bec.
    of an itoa at l. 569 of utility.cc plus the dos path stuff
    over there ... - COMMENTED OUT !
    ===========================================================
*/

#ifndef UTILITY_H
#define UTILITY_H

// for inline functions :
#include <math.h> 

//   ===========================================================
//   ======================   GLOBAL   =========================
//   ===========================================================
extern FILE *bzzFileOut;
extern FILE *bzzFileIn;

extern char bzzYesNo;
extern const float MACH_EPS;

extern const float SMALLFLTSIG;  // Near 0 w.r.t sig. places; for comparisons

extern const float BIG_FLOAT;
extern const float TINY_FLOAT;
extern const double BIG_DOUBLE;
extern const double TINY_DOUBLE;
extern const long double BIG_LONG_DOUBLE;
extern const long double TINY_LONG_DOUBLE;
extern const char BIG_CHAR;
extern const char TINY_CHAR;
extern const unsigned char BIG_UNSIGNED_CHAR;
extern const unsigned char TINY_UNSIGNED_CHAR;
extern const int BIG_INT;
extern const int TINY_INT;
extern const unsigned int BIG_UNSIGNED_INT;
extern const unsigned int TINY_UNSIGNED_INT;
extern const long int BIG_LONG;
extern const long int TINY_LONG;
extern const unsigned long int BIG_UNSIGNED_LONG;
extern const unsigned long int TINY_UNSIGNED_LONG;

//   ===========================================================
//   ======================   ERRORS   =========================
//   ===========================================================

//========== Error Type ==============
extern const char *const ERR_RANGE;
extern const char *const ERR_SPACE;
extern const char *const ERR_OPEN_FILE;
extern const char *const ERR_CLOSE_FILE;
extern const char *const ERR_READING_FILE;
extern const char *const ERR_WRITING_FILE;
extern const char *const ERR_CHECK_DIMENSION;
extern const char *const ERR_OUT_OF_SCOPE;
extern const char *const ERR_FACTORIZED;
extern const char *const ERR_IMPLEMENTATION;

//============ Function Type =============
extern const char *const ERR_FUNCTION;
extern const char *const ERR_CONSTRUCTOR;
extern const char *const ERR_OPERATOR;

//   ===========================================================
//   =================  Inline  Functions  =====================
//   ===========================================================

//   ************************< Abs >****************************
//   * Purpose: To have a unique function for absolute value   *
//   * Description: Overload for float, int, double            *
//   * Example: float x = Abs(y);                              *
//   ***********************************************************
inline float Abs(float a)
   {return fabs(a);}
inline double Abs(double a)
   {return fabs(a);}
inline int Abs(int i)
   {return abs(i);}

//   *************< Max, Min, MaxAbs, MinAbs >******************
//   * Purpose: Maximum and minimum for a pair of numbers      *
//   * Description: Functions overload float, int, double      *
//   * Example: x=Max(a,b);i=Min(l,k); x=MaxAbs(a,b);          *
//   ***********************************************************
inline float Max(float a,float b)
   {return (a > b ? a : b);}
inline float Min(float a,float b)
   {return (a < b ? a : b);}
inline double Max(double a,double b)
   {return (a > b ? a : b);}
inline double Min(double a,double b)
   {return (a < b ? a : b);}
inline int Max(int l,int k)
   {return (l > k ? l : k);}
inline int Min(int l,int k)
   {return (l < k ? l : k);}

inline float MaxAbs(float a,float b)
   {return (fabs(a) > fabs(b) ? fabs(a) : fabs(b));}
inline float MinAbs(float a,float b)
   {return (fabs(a) < fabs(b) ? fabs(a) : fabs(b));}
inline double MaxAbs(double a,double b)
   {return (fabs(a) > fabs(b) ? fabs(a) : fabs(b));}
inline double MinAbs(double a,double b)
   {return (fabs(a) < fabs(b) ? fabs(a) : fabs(b));}
inline int MaxAbs(int l,int k)
   {return (abs(l) > abs(k) ? abs(l) : abs(k));}
inline int MinAbs(int l,int k)
   {return (abs(l) < abs(k) ? abs(l) : abs(k));}

//   ===========================================================
//   ====================   PROTOTYPES   =======================
//   ===========================================================
float MachEps(void);

void Message(const char *myFormat,...);
void Error(const char *myFormat,...);

float Control(double value);
double Control(long double value);

void Swap(float *x,float *y);
void Swap(double *x,double *y);
void Swap(int *x,int *y);
void Swap(float **x,float **y);
void Swap(double **x,double **y);
void Swap(int **x,int **y);

float Max(int n,float *x,int *imax = 0);
float MaxAbs(int n,float *x,int *imax = 0);
float Min(int n,float *x,int *imin = 0);
float MinAbs(int n,float *x,int *imin = 0);

void Sum(int n,float *lval,float *rval,float *result);
void Sum(int n,float *lvalAndResult,float *rval);
void Sum(int n,float *lvalRvalAndResult);
void Difference(int n,float *lval,float *rval,float *result);
void Difference(int n,float *lvalAndResult,float *rval);
void Difference(int n,float *lval,float *rvalAndResult,int);
float Dot(int n,float *lval,float *rval);
void Product(int n,float lval,float *rval,float *result);
void Product(int n,float lval,float *rvalAndResult);
void Division(int n,float *lval,float rval,float *result);
void Division(int n,float *lvalAndResult,float rval);

float SqrtSumSqr(int n,float *x);
double SqrtSumSqr(int n,double *x);

void Sort(int n,float *x);

#endif // UTILITY_H

/* 

//   *********************< TempFile >**************************
//   * Class for creating temporary files                      *
//   ***********************************************************

#ifndef TEMPFILE_H
#define TEMPFILE_H


//   ===========================================================
//   ==================   class TempFile   =====================
//   ===========================================================

class TempFile
   {
private:
   static int countTempFile;
   char *nameTempFile;
public:
   // default constructor TempFile ff[10];
   TempFile(void){nameTempFile = 0;}

   // constructor TempFile f1("D:"),f2("C:\\TEMP\\");
   TempFile(char *directoryTemp)
      {NewFileName(directoryTemp);}

   // destructor
   ~TempFile(void){delete nameTempFile;}

   // provides the unique name f1.FileName
   char *FileName(void){return nameTempFile;}

   // for using after default ff[0].NewFileName("D:");
   void NewFileName(char *directoryTemp);
   };

#endif // TEMPFILE_H

*/

// ---------- end of utility.h -----------------
