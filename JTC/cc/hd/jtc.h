/* jtc.h - Classes for models of jumping-to-conclusions, and possibly 
           getting stuck there ! Starting with beads-in-a-jar-guessing 
	   simulations. 

  The jtc and related (jtc_exp ...) methods make use of the GNU 
  Scientific Libray and the Array methods of blitz++. Also moves from
  using floats to doubles, mostly.

*/                   

#ifndef JTC_H
#define JTC_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>    //input & output streams
#include <iomanip>    // i-o  parametrized manipulators
#include <sstream>    // in-memory string manipulations for i/o, new standard
#include <string>
#include <cstdlib>
#include <time.h>
using namespace std; 

#include <blitz/array.h>
#include <gsl/gsl_statistics_double.h>

BZ_USING_NAMESPACE(blitz)

// A.a.  Simplest bead-jar-guessing RL model.
//
// -----------------------------------------------------------------------------
// A.a. Simplest beads-jar-guessing 'RL' model, with 20 draws.

class Beads
{
 protected:
  static const char *const ERROR;
  int Glab, Blab, Slab;    // useful labels, 1,-1,0
  int N ;   // No. no. of samples. MATRICES BELOW ARE N+1 x N+1  !

  inline void CheckS(int s); //make sure s is less than or equal to N
  int meN;  // for Makintosh - El-Deredy; is <= N.
 
 public:

  Beads(int maxSampl);
  ~Beads();            // trivial ...

  // Parameters that users may set:  
  float q;           // Proportion of g beads in jar G (the higher), P(g|G)
  float PoG;         // Initial (prior) prob. that cause is jar G
  float CBG,CGB,CS;  // Cost coefficients for calling B when G is true, calling
                     // G when B is true and for asking for a further Sample. 
                     // Have -ve values to denote cost/punishment.
                     // In Makintosh - El-Deredy sims CS is the deduction per step,
                     //    0 -5 -10, but CGB and CGB are the subject's 'own' or 
                     //    'internal' assumed returns.
  float R0;          // Initial (external) 'prize', for Makintosh - El-Deredy sims.
  float meCS;        // Makintosh - El-Deredy loss-per-extra-sample (eg 0,-5,-10)
  float kT;          // Used in some action-selection Softmaxes as noise param.

                     // In what follows, P is initialized to zero but the rest are
                     //    initially set to all elements = -666.0 (i.e. 'invalid')
  Array<double, 2> PG, V;       // Probabilities and Values for all states (ng,S-ng)
  Array<double, 2> QG, QB, QS;  // Action values for each action available
  Array<double, 2> AS;          // Advantage values for Sample action. Note Sampling is not
	                   // an option for last step - so given 'symbolic' value of -666.0	
  Array<double, 2> MS, MG, MB ;  // Policy (probability of selection) vals for Actions,
                                  // 'decide Sample', 'decide G',...

  // look-up methods:
  int MaxSamplN(){ return N; };
  int meMaxSamplN(){ return meN; };

  // Setting non-public variables:
  void Set_meN(int new_meN); // will create error if attempt to set to greater than N
  // in version ME2, set meCS and R0 acc. to CS, depending on whether
  // in zero (cond=1) , low or high-cost condition cond :
  void SetME2_R0_meCS_acc_to_CS(int cond) ;

  // post. probs. that  G (or B) are the  causes
  // given ng g's drawn out of S samples :
  float GProb(int ng, int S){
    return 1/(1+ (PoG/(1-PoG))*pow( q/(1-q) , S-2*ng ) ) ; }
  float BProb(int ng, int S){ return 1-GProb(ng,S); } ;  //as above for cause=B.
  // log-likelihood ratio of G over b, if ng out of S draws were g:
  // i.e. ln( P(d1,d2,...,dS | G) / P(d1,d2,...,dS | B)
  float LLR_GoB(int ng, int S){  return (S-2*ng)*log((1-q)/q) ; } ;
  // unit of accumulation of llr for 'consonant' decision, i.e. G if last 
  // draw is g, or B if last draw is b :
  double consDeltaLLR(){ return log(q/(1-q)) ; };

  //  For action-value calculation. First, the very last step, where no further
  // sampling can take place.
  float QlastB(int ng){ return CBG*GProb(ng,N) ; };
  float QlastG(int ng){ return CGB*BProb(ng,N) ; };
  // Makintosh - El-Deredy schema is in principle a bit more complicated,
  // but can boils down to the same as often R0+N*CS=0 (in first version
  // of simulations). 
  // The following code allows freedom from this constraint but uses CS 
  // instead of meCS (a quick-and-dirty initial shortcut ...) 
  float QlastBME0(int ng){ return CBG*GProb(ng,N)+BProb(ng,N)*(R0+N*CS) ; };
  float QlastGME0(int ng){ return CGB*BProb(ng,N)+GProb(ng,N)*(R0+N*CS) ; };
  // Less dirty version 1. Note meN-1 as first bead is given 'free'.
  float QlastBME1(int ng){ return CBG*GProb(ng,meN)+BProb(ng,N)*(R0+(meN-1)*meCS) ; };
  float QlastGME1(int ng){ return CGB*BProb(ng,meN)+GProb(ng,N)*(R0+(meN-1)*meCS) ; };
  // ME2, 'Clean' but like ME0. Again meN-1 as first bead is given 'free' and 
  // in addition CS is chosen from distro and scales and R0 and meCS, which are used
  // here, MUST BE DERIVED FROM CS AND BE READY FOR THIS !
  float QlastBME2(int ng){ return QlastBME1(ng); };
  float QlastGME2(int ng){ return QlastGME1(ng); };


  // Calculate softmax probability of all 3 actions, based on advantages based 
  // on VALID Qs, so as to avoid 'inf' and 'nan' exponential calculations.
  // Also return label (as above) for preferred action.
  int PQSoft(int ng, int s, float *pPG,  float *pPB, float *pPS ); 

  // Actually I don't think the following is that useful - not coded yet:
  void OptVasQmaxBGS(int ng, int S); // V(ng,S) for optimal policy is 
                                     // the max. Q from that state

  // Explicit calculation of Q (Action) values, given all that is needed. Version 1
  // has cost just proportional to CBG, CGB and CS :
  float QSVal(int ng, int S, float nextV_if_g, float nextV_if_b);
  float QBVal(int ng, int S){ return CBG*GProb(ng, S); } ;
  float QGVal(int ng, int S){ return CGB*BProb(ng, S); } ;

  // Makintosh - El-Deredy experiments:
  // Explicit calculation of Q values, given all that is needed: Version 0 'dirty':
  float QSValME0(int ng, int S, float nextV_if_g, float nextV_if_b);
  float QBValME0(int ng, int S){ return GProb(ng, S)*CBG+BProb(ng,S)*(R0+S*CS); } ;
  float QGValME0(int ng, int S){ return BProb(ng, S)*CGB+GProb(ng,S)*(R0+S*CS); } ;
  // Explicit calculation of Q values, given all that is needed: Version 1 'cleaner' :
  float QSValME1(int ng, int S, float nextV_if_g, float nextV_if_b);
  float QBValME1(int ng, int S){ return (S>0) ? CBG*GProb(ng,S)+(R0+(S-1)*meCS)*BProb(ng,S) : CBG*PoG    +R0*(1-PoG); } ;
  float QGValME1(int ng, int S){ return (S>0) ? CGB*BProb(ng,S)+(R0+(S-1)*meCS)*GProb(ng,S) : CGB*(1-PoG)+R0*PoG    ; } ;
  // Version 2: NB CS needs to be used as the population variable, chosen from 
  // distro and scales while CGB & CBG are fixed.
  // R0 and meCS, which are used here, MUST BE DERIVED FROM CS AND BE READY FOR THIS !
  float QSValME2(int ng, int S, float nextV_if_g, float nextV_if_b){ return QSValME1(ng, S, nextV_if_g, nextV_if_b );} ;
  float QBValME2(int ng, int S){ QBValME1(ng,S); } ;
  float QGValME2(int ng, int S){ QGValME1(ng,S); } ;


  // MOST IMPORTANT METHODS: FILL IN THE ARRAYS ACC. TO CURRENT VALUES OF
  // PARAMETRES. First, noise-free case. Note that deciding before any
  // draws is an option - no variable to disable this, as in SoftM_... 
  void CalcPG_V_QG_QB_QS_AS();
  // then, SOFTMAX-NOISY case. It disallows declaring before
  // any draws if zeroDrawChoice = 0
  void SoftM_PG_V_QG_QB_QS_AS_MS(int zeroDrawChoice = 0);
  //  SOFTMAX-NOISY case FOR Sequential Prob. Ratio Test. It disallows 
  // declaring before any draws if zeroDrawChoice = 0 and fills in  the
  // PG and MS matrices meaningfully. Use the QG and QG to store the 
  // log-likelihood-ratios.  The rest are set to invalid flag -666.
  void SPRT_Soft_PG_and_MS(  double phi_H, double L_over_H=-1.0
                           , int zeroDrawChoice = 0);
  // Makintosh - El-Deredy scheme: MDP approach, 
  // SOFTMAX-NOISY case. It disallows declaring before
  // any draws if zeroDrawChoice = 0.
  // 0. Quick and dirty:
  void SoftME0_M_PG_V_QG_QB_QS_AS_MS(int zeroDrawChoice = 0);
  // 1. Charlotte's first experiment, cleaner:
  void SoftME1_M_PG_V_QG_QB_QS_AS_MS(int zeroDrawChoice = 0);
  // 2. Charlotte's first experiment, using ME2 versions. NEED TO MAKE SURE
  //    THAT R0 and meCS ARE CORRECTLY DERIVED FROM CS BEFORE USING THIS:
  void SoftME2_M_PG_V_QG_QB_QS_AS_MS(int zeroDraws){ SoftME1_M_PG_V_QG_QB_QS_AS_MS(zeroDraws); };


  // THE FOLLOWING *ASSUME* THAT Q VALUES HAVE BEEN ALL FILLED IN VALIDLY !
  // First, return the label of the choice w. the greatest action val. at
  // draw S. Does various checks to ensure correct use. May be slower !
  int DeciQmax(int sampl, float *seqAr, int seqArEnd, int seqArStart=1);
  // Return the sample number where QS ceases to be the prefered option, i.e.
  // how many draws-to-declare. Tailored to classes as e.g. Garety91 (rl_exp.h)
  int DeciDrawQmax(float *seqAr, int seqArEnd, int seqArStart=1);
  //  Another variation that returns the prefered action, like DeciQmax,
  //  but here ng and s are given, not sequence identity and no_of_draws.
  //  Also fills in, deterministically, the MS, MG and MB element in question.
  //  Again, assumes that Q's have been filled in correctly.
  int DeterMQmax(int ng, int s);

  // THE FOLLOWING ONLY NEEDS MS TO BE FILLED IN CORRECTLY.
  // Return the probability that the decision to declare, rather than sample
  // again, is taken at draw dNum GIVEN THAT we start from draw 0 - i.e.
  // sample has been depleted ! Assumes MS contains up-to-date policy Py's.
  float CumulPDec(int dNum, float *seqAr, int seqArEnd, int seqArStart=1);

  // Write to a file one_seq.csv the draw num,
  // sequence, PG, QG, QB, QB, DeciQmax like a matrix, then the params and
  // comments in a bit of text.
  void WriteSeqDeciPGQVPar(float *seqAr, int seqArEnd, int seqArStart=1); 
};

// -----------------------------------------------------------------------------
#endif  // JTC_H
