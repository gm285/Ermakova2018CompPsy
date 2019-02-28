// methods for jtc_exp.h, experiments that relate to reinforcement-learning 
//  systems.   
/* Some default val. in C quotes */
                
#include "jtc_exp.h"
#include "nr_ut.h"
  
const char *const RBWelcomeDat::ERROR="  >>---> RBWelcomeDat error : " ;
const char *const Beads_Exp_I::ERROR="  >>---> Beads_Exp_I error : " ;

/* = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = */
/* B.b. from Bentall_et_al_06 ('Welcome study'-
                             - The Phenomenology and Cognitive Structure...) */

RBWelcomeDat::RBWelcomeDat()
{
}

void RBWelcomeDat::setSequences(int seqN, const string *seqstr) {
 cerr << "setSequences(" << seqN << ")\n";
 sN = seqN;
 seq = new int* [sN + 1];
 for (int i = 0; i < 10; i++) {
     seq[i+1] = ivector(1, 20);
     fillSequence(seq[i+1], 1, 20, seqstr[i]);
 }

 // Now set up number-of-yellow/good/red as fn. of no. sampled :
 ng_tot = new int* [sN + 1];
 // Indexed from 1, for historical reasons
 for (int i = 1; i <= sN; i++) {
     ng_tot[i] = ivector(0, 20);
     ng_tot[i][0] = 0;
     for (int j = 1; j <= 20; j++) {
         ng_tot[i][j] = ng_tot[i][j-1] + seq[i][j];
     }
 }
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
RBWelcomeDat::~RBWelcomeDat()
{
 
  for (int i = 1; i <= seqN(); i++) {
    free_ivector(seq[i], 1, 20);
    seq[i] = NULL;
    free_ivector(ng_tot[i], 0, 20);
    ng_tot[i] = NULL;
  }
  delete seq;
  delete ng_tot;
}

void RBWelcomeDat::fillSequence(
        int *target, int start, int length, const string& str) const
{
    cerr << "str: " << str << "\n";
    for (int i = 0; i < length; ++i) {
        if (!str[i]) {
            err_msg("Not enough characters in str");
            return;
        }
        if (str[i] == '1') {
            target[i + start] = 1;
        } else if (str[i] == '0') {
            target[i + start] = 0;
        } else {
            err_msg("Invalid character in str");
            return;
        }
    }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
Beads_Exp_I::Beads_Exp_I() :
  Beads(20), iter(6), IP2()  {   // NB iter.parNo assignement hidden here ! 

  // First, complicated stuff to do with all the iteration parametres :    
  // (file stuff -> look up or future !)                                            
  int iparN = iter.ParNo() ;
  iterParPtr = new float* [iparN+1];   iterIterNo = new float* [iparN+1];
  iterEndPtr = new float* [iparN+1];

  // Items Costs (cf. menu1 in os2t6m.cc) :
  //   kT :range     ... and No. of pts. for iteration :
  kT_e  = 10.01 ;            kT_N      = 40. ;  
  iterParPtr[1]= &kT;      iterEndPtr[1]= &kT_e;     iterIterNo[1]= &kT_N;

  //   CS :range     ... and No. of pts. for iteration :             
  CS_e  = -0.001 ;         CS_N      = 80. ;  
  iterParPtr[2]= &CS;      iterEndPtr[2]= &CS_e;     iterIterNo[2]= &CS_N;

  // in beads05rbw.cc only the two above are actually iterated !
  CGB_e = -5.001;          CGB_N = 1.0;
  iterParPtr[3]= &CGB;     iterEndPtr[3]= &CGB_e;    iterIterNo[3]= &CGB_N;

  // ( potentially !)  the rest:
  q_e  = .60 ;            q_N      = 1. ;  
  iterParPtr[4]= &q;      iterEndPtr[4]= &q_e;     iterIterNo[4]= &q_N;

  PoG_e  = 0.05;           PoG_N      = 1. ;  
  iterParPtr[5]= &PoG;     iterEndPtr[5]= &PoG_e;     iterIterNo[5]= &PoG_N;

  // to iterate contr. distr. param. 1; esp. SPRT Threshold_hi 
  CDPar1 = 0.0;  CDPar1_e = 30.0;  CDPar1_N=150. ; 
  iterParPtr[6]= &CDPar1;     iterEndPtr[6]= &CDPar1_e;    iterIterNo[6]= &CDPar1_N;

  // REM MINUS  rCW, -rCW, used as the ratio of lower to higher thresh. in SPRT,
  // and should have -rCW < 1.
  // rCW :               range:        ... and No. of pts. for iteration :             
  rCW = 1. ;             rCW_e=0.5 ;               rCW_N = 1. ;

  // ... and  contr. distr. param. 2:
  CDPar2 = -1.0; CDPar2_e = -.5;   CDPar2_N=1. ; 

  // Strings to help w. names of output files:
  // initialise default sequence var. & related storage array. Note 0-20 component !
  seqNum = 10;      seqLen=20;
  PVecs.resize(Range(1,roundfl(*iterIterNo[1])),Range(1,roundfl(*iterIterNo[2])),Range(1,seqNum),Range(0,seqLen));
  PVecs = 0;

  // Prepare the Itepolation object :
  IP2.setOrderPars(4,4);  // use 4x4 set of data for interpolating a point.

  // for Mackintosh - El-Deredy stuff
  retNum = 3;   // [local copy of] the number of different patterns of returns that
                // the expt. uses - e.g. 3 for ME1 that has zero, low and high. 

} // End of ctor for Beads_Exp_I
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
Beads_Exp_I::~Beads_Exp_I()
{
  delete iterParPtr;  delete iterIterNo; delete iterEndPtr;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
void Beads_Exp_I::RefreshIter()
{ 
  CBG = CGB/rCW ; // no need for CDPar2 = -rCW*CDPar1 or such like as yet.
  iter.Set(iterParPtr,iterEndPtr,iterIterNo); 

} // End of RefreshIter()
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
int Beads_Exp_I::StoreAllPVecs(const RBWelcomeDat *RBWptr, int validZero)
{ 
  // Rem:  PVecs is 4-D array: param1 x param2 x sequence_Num x Probability_vector
  //       e.g. for main costed bayesian:  kT_N  x  CS_N x  3 (for Wellcome) x fN     
  // if key extends don't match, abort. Should be OK if PreparePVecs1 was called.
  if( (RBWptr->seqN()) != PVecs.extent(thirdDim))
      Error("%s%s \n",Beads_Exp_I::ERROR,"RBWptr->seqN() != PVecs.extent(thirdDim) in StoreAllPVecs.");
  if ((RBWptr->SamplMaxN()) != (fN-1))
      Error("%s%s \n",Beads_Exp_I::ERROR,"max. draw number mismatch in StoreAllPVecs.");
  if(   (roundfl(*iterIterNo[1]) != PVecs.extent(firstDim)) 
     || (roundfl(*iterIterNo[2]) != PVecs.extent(secondDim))) 
    Error("%s%s \n",Beads_Exp_I::ERROR,"itern. number doesn't match PVecs extent in StoreAllPVecs."); 
    /* cout << "*iterIterNo[1]: "<< *iterIterNo[1] << "PVecs.extent(firstDim)" <<  PVecs.extent(firstDim) << endl ;
    cout << "*iterIterNo[6]: "<< *iterIterNo[6] << "PVecs.extent(secondDim)" <<  PVecs.extent(secondDim) << endl ; */

  int seq_n, s, initDraw, endDraw, ng;                      // auxiliaries
  double  prod;

  initDraw = PVecs.lbound(fourthDim);  endDraw=PVecs.ubound(fourthDim); 

  for (seq_n = 1; seq_n <= (RBWptr->seqN()); seq_n++) {
 
    // now start calculating. First P(d2d=0), by hand :
    if (validZero){  // If subject allowed to declare at 0 draws ...
      PVecs(AtStepIter(1),AtStepIter(2), seq_n,initDraw) = 1-MS(1,1) ;
      prod = MS(1,1);  // zero-draws 'sample again' Py
    }
    else {
      PVecs(AtStepIter(1),AtStepIter(2), seq_n,initDraw) = 0.0;
      prod = 1.0 ;  // zero-draws OBLIGATORY 'sample again' Py !
    }
    
    // now P(d2d=1) to  P(d2d=SamplMaxN()) :
    for (s=initDraw+1; s<=endDraw; s++) {
      ng = RBWptr->ng_tot[seq_n][s];
      // the following should work for last step too:
      PVecs(AtStepIter(1),AtStepIter(2), seq_n, s) = prod*(1-MS(ng+1,s-ng+1)); 
      prod *= MS(ng+1,s-ng+1);               // ready for next step
    }
    
  } // end loop over sequences,  for (seq_n = 1; seq_n <=...

  return 0; 
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
double Beads_Exp_I::CPDistr1(int d2d, double param)
{  
  double fact, PSum, lt;
  int i,j;
  if (param<0.0) lt=CDPar1 ;
  else lt = param;

  if (d2d<1) return 0.0;
  else if (d2d<(fN-1)){    //e.g. if d2d<20 for usual setup
    fact=1.0;
    // note that in both step that follow, we use d2d-1 - Poisson distro
    // is displaced by 1 as constrained to have P(d2d=0) = 0
    for (i=1; i<d2d; i++) fact*=i;
    return pow(lt,(d2d-1))*exp(-lt)/fact ;
  }
  else if (d2d==(fN-1)) {
    PSum = 0;
    for (j=1; j<(fN-1); j++) { 
      // first calculate the pmf value of the Poisson at each of first
      // 20 states (actually bother with 19 only, as stage 0 is constrained to have P=0)
      fact=1.0;
      // again next two lines are for f((j-1),lt)
      for (i=1; i<j; i++) fact*=i;
      PSum += pow(lt,(j-1))*exp(-lt)/fact ; 
    }
    return 1.0 - PSum ;
  }
  else {  // i.e. invalid argument has been given to fn. 
    cout<< "d2d: "<<d2d<< endl;
    Error("%s%s \n",Beads_Exp_I::ERROR,"invalid d2d argument in CPDistr1");
  }

}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
double Beads_Exp_I::SPRT(int step, int seq, const RBWelcomeDat *RBWptr)
{
  if (step<0 || step>20) err_msg("Beads_Exp_I::SPRT needs step 0 to 20");
  if (seq < 1 || seq > seqN()) err_msg("Beads_Exp_I::SPRT needs seq 1 to 3");
  return LLR_GoB(RBWptr->ng_tot[seq][step], step);
}
//=============================================================================
// end of File jtc_exp.cc
