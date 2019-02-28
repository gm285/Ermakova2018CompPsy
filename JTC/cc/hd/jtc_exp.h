// File with experimental setups for jumping to conclusions / 
// beads-in-a-jar and  related systems. 


#ifndef JTC_EXP_H
#define JTC_EXP_H

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

#include "maths.h"
#include "utility.h"

#include "jtc.h"
#include "blitz_array_aux.h"

#include <blitz/array.h>
#include <gsl/gsl_statistics_double.h>

BZ_USING_NAMESPACE(blitz)

/* A.a. Summary data from Beads and related experiments :
     a.  
*/
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
/* B.b. from Bentall_et_al_06 ('Welcome study'-
                             - The Phenomenology and Cognitive Structure...) 
		 and also the Mackintosh - El-Deredy variant                  */

class RBWelcomeDat {
 protected:
    static const char *const ERROR;
    int sN;   // auxiliary for No. of sequences presented
    void CheckInd(int ind){ 
	if (ind<1 || ind>20) Error("%s%s \n",RBWelcomeDat::ERROR,"Index off range" ); };

 public:
  int SamplMaxN() const { return 20; };
  // allow free access - we trust our users ! This is the standard bead / word
  // sequences used in Young & Bentall 97 and various others, and its inverse.
  int **seq;  // seqN entries. Each one a float array for unknown reasons...
  // These are the number-of-yellow/good as fn. of no. sampled for Welcome:
  int **ng_tot;  // will point to the running totals below

    RBWelcomeDat();
    ~RBWelcomeDat();

  void setSequences(int seqN, const string *seqstr);

//  ---------------   Methods of general applicability in this class ---
    // General Lookups :
    int seqN() const { return sN; }; // Either for Wellcome or ME variants

 //  .------------------------------------------------------------------.
 //  |             Methods concerning Welcome I - RPB data              |
 //  |__________________________________________________________________|

  private:
    void fillSequence(int *target, int start, int length, const string& str) const;
};
/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
/* C.a. Classes with the 'Beads' machinery that can do more */
/*      First, can also do just iteration:                  */
class Beads_Exp_I : public Beads
// See attenmod.h -> class Exp_Dual_At_K_II : public DualAtK_II ...
{
 protected:
  static const char *const ERROR;
  // The following are used by RefreshIter: 
  float **iterParPtr, **iterIterNo, **iterEndPtr ;
  Iter iter ;
  // auxiliary stuff for file output :  
  static const int fN  = 21;  // auxiliary for No. of P components / output files.

  int seqNum, seqLen;  // local copies of number and length of  draw-sequences
                       //  that the experiment uses - def. 20 (-> P00 to P20 !), 3
  int  retNum; // [local copy of] the number of different patterns of returns that
               // the expt. uses - e.g. 3 for ME1 that has zero, low and high. 

 public:
  float CGB_e, CGB_N;        // endpoint and number of CGB 
  float kT_e,   kT_N;  
  float CS_e,   CS_N;
  float rCW, rCW_e, rCW_N;  // This is the ratio CGB/CBG, putatively similar
                            // to false-neg.cost/false-pos.cost 
  float q_e,     q_N;       // in case we need to iterate the P(g|G)
  float PoG_e, PoG_N;       // in case we need to iterate the PoG.

  float CDPar1, CDPar1_e, CDPar1_N;  // to handle & iterate contr. distr. param. 1
  float CDPar2, CDPar2_e, CDPar2_N;  // ... and 2.

  Array<double, 4> PVecs;   // 4-D array: param1 x param2 x sequence_Num x Probability_vector
                            //       e.g:  CS_N  x  CGB_N x  3 (for Wellcome) x fN
                            // If trials have different 
     // BayesianLik3Beads BayesianLik3Words ContrLik3Beads ContrLik3Words BestParBayesBeads1 2 3 BestBayesianParamWords1 2 3 ControlParBeads1 2 ContolParWords1 2
  //  OptL3B    OptPar1B          OptPar2B    OptPar3B  \
  //  OptL3W    OptPar1W          OptPar2W    OptPar3W  \
  //  ContrL3B    ContrPar1B                            \
  //  ContrL3W    ContrPar1W   
  // for file output, plan is to add the pt. numbers,
  // [the control distribution best-likelihoods and corresponding params]
  //  and the two derivative ratios,  Opt/Contr_L3B and Opt/Contr_L3W, 
  // as additional cols.

  Interpolation2D IP2;

  Beads_Exp_I();
  ~Beads_Exp_I();

  // paramIter access fns. Initial 'Set'ting in  ctor  
  void RefreshIter() ; // for after each menu-dep input               
  int DoIter() { return iter.Do() ; } ;
  int StepIter() { return iter.Step() ; } ;
  float StartIter(int j) { return iter.Start(j) ; } ;
  float EndIter(int j)   { return iter.End(j)   ;};
  // determine total no. of iterations to do; also assigns result to totRunN 
  unsigned long TotIterN() { return iter.TotN() ; } ;
  unsigned long ItersDone() { return iter.ItersDone() ; } ;
  int AtStepIter(int k){ return iter.AtStep(k); };
  float FractRunsDone(){ return iter.FractDone() ; } ;

  // resetting etc.
  // PreparePVecs1() etc. use existing params. Inline for clarity rather than speed. 
  void PreparePVecs1(){
    PVecs.resize(Range(1,roundfl(*iterIterNo[1])),Range(1,roundfl(*iterIterNo[2])),Range(1,seqNum),Range(0,seqLen));
    PVecs = 0; } ;

  // ** Methods for writing results (mostly 2D arrays) to disk **

  // Another version, similar to above, but designed so that first iter par 
  // is kT and second is CS :
  int WriteAllPd2d_2(  float *seqAr, int seqArEnd=20, int seqArStart=1
                     , char sepChar = '\t', int check=1 );
  // Method for writing 21-D Bayesian-der. Prob. vectors to 4-D array of the form
  // PVecs.resize(Range(1,roundfl(GB.CS_N)), Range(1,roundfl(GB.CGB_N)), Range(1,RBWelc.seqN()), Range(0,20))
  // validZero var says if decision at 0 draws is allowed. 
  // Version using RBWelcomeDat :
  int StoreAllPVecs(const RBWelcomeDat *RBWptr, int validZero=1 ) ; 

  // Control distributions & Py products :
  // 1. Poisson, except stage 20 which mops up the remainder of all the rest.
  //    If param is <0, it uses CDPar1 as the lambda (or lambda*t) of the Poisson 
  //    distr. 
  double CPDistr1(int d2d, double param = -666 ) ;
  // 2.   Sequential Probability Ratio Test using CDPar1
  //   a. Directly calculate the error-free log-likelihood-ratio for a particular 
  //      step of particular sequence. 
  double SPRT(int step, int seq, const RBWelcomeDat *RBWptr);

  int seqN() const { return seqNum; };
};

// Don't forget the semilcolon at the end of class declarations !
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

#endif    // JTC_EXP_H

// eof
