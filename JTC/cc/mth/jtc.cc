// jtc.cc - methods for Classes for models of jumping-to-conclusions
//   and other paranoid phenomena.  


#include "jtc.h"
#include "maths.h"
#include "utility.h"

const char *const Beads::ERROR="  >>---> Beads error : " ;

//------------------------------------------------------------------------------
/// A.a.  Simplest beads-jar-guessing 'RL' model
Beads::Beads(int maxSampl) : PG(), V(), QG(), QB(), QS(), AS(), MS(), MG(), MB()
{
  int i,j;
  // Structural param. values :
  Glab=1; Blab=-1; Slab=0;
  meN = N = maxSampl ;

  // Default param. values :
  q   = 0.60 ;                  // As per Bentall et al Wellcome study
  PoG = 0.50;
  CBG = -100.0 ; 
  CGB = -100.0 ;
  CS  = -5.001 ;
  kT  = 0.01;

  R0  = 0.0 ;  // as per no-cost (or benefit) NCC or 'normal beads'

  // Now the main objects (matrices). Allocate space for the arrays according to N:

  PG.resize(Range(1,N+1), Range(1,N+1));  // These will all be set to zero ...
  V.resize(Range(1,N+1),  Range(1,N+1)); 
  QG.resize(Range(1,N+1), Range(1,N+1)); 
  QB.resize(Range(1,N+1), Range(1,N+1)); 
  QS.resize(Range(1,N+1), Range(1,N+1));

  PG = 0;  V = 0;   QG = 0;   QB = 0; 
  
  AS.resize(Range(1,N+1), Range(1,N+1));  // ... but these set to 'invalid' value:	  
  for (i=1; i<=N+1; i++)  for (j=1; j<=N+1; j++) { AS(i,j)=-666; QS(i,j)=0; } 

  MG.resize(Range(1,N+1), Range(1,N+1));   // These will all be set to zero ...
  MB.resize(Range(1,N+1), Range(1,N+1));
  MS.resize(Range(1,N+1), Range(1,N+1));

  MS = 0;   MG = 0;   MB = 0; 

}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// hmmm ..
Beads::~Beads()
{
}  
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
inline void Beads::CheckS( int s)
{ 
    if (s<0 || s>N)
	Error("%s%s \n", Beads::ERROR,"Sampling index out of range." ); 
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void Beads::Set_meN( int new_meN)
{ 
  if (new_meN<1 || new_meN>N)
	Error("%s%s \n", Beads::ERROR,"Attempt to set invalid meN" ); 
  meN = new_meN; 
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
void Beads:: SetME2_R0_meCS_acc_to_CS(int cond) 
{ 
  if (cond<1 || cond>3) 
    err_msg("SetME2_R0_meCS_acc_to_CS: cond must be bet. 1 and 3");    
  meCS=(cond-1)*CS;      
  R0= -meN*meCS;    
  //debug: cerr<<"\nCS: "<<CS<<" R0: "<<R0<<"  meCS:"<<meCS<<endl;  
} 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
int Beads::PQSoft(int ng, int s, float *pPG,  float *pPB, float *pPS )
{
 CheckS(s);
 int auxQmax, labQmax;
 double qg, qb, qs, qmax, exp_qs_term ;   exp_qs_term = 0.0;  // non-trivial !!
 double pg, pb, ps, denom ;               pg=pb=ps=denom=0.0;  // various auxil

 // qg and qb always exist, and matrices start from 1 not 0:
 qg = QG(ng+1,s-ng+1);                     qb = QB(ng+1,s-ng+1); 

 if (qg >= qb) {
   labQmax = Glab;  qmax = qg; 
 }
 else {
   labQmax = Blab;  qmax = qb; 
 }
 
 if (s<N) {                     // i.e. if not on the last step
   qs = QS(ng+1,s-ng+1);
   if (qs > qmax) {
     labQmax = Slab; qmax = qs;  exp_qs_term = 1.0;
   }
   else // non-trivial exponential qs term only if non-max qs:
     exp_qs_term = exp((qs-qmax)/kT)  ;  
 }

 // now calculate the dreaded denominator of the softmax !
 if (labQmax == Glab) {
   denom = 1 + exp((qb-qg)/kT) + exp_qs_term ; // calc in double !
   pg = 1/denom ;                    // conv. to float !
   pb = (exp((qb-qg)/kT))/denom ;
   ps = exp_qs_term/denom ;
 }
 else if  (labQmax == Blab) {
   denom = 1 + exp((qg-qb)/kT) + exp_qs_term ;
   pg = (exp((qg-qb)/kT))/denom ;                    // conv. to float !
   pb = 1/denom ;
   ps = exp_qs_term/denom ;
 }
 else { // i.e. if labQmax == Slab
   denom = 1 + exp((qg-qs)/kT) + exp((qb-qs)/kT) ;
   pg =  (exp((qg-qs)/kT))/denom ;                 // conv. to float !
   pb =  (exp((qb-qs)/kT))/denom ;
   ps = 1/denom ;
 }

 *pPG  = pg;     *pPB = pb;   *pPS = ps;
 return labQmax ;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
float Beads::QSVal(int ng, int S, float nextV_if_g, float nextV_if_b)
{
 CheckS(S);
 // In the following, forbid ng=N as then Sampling not an option.
 if (ng>=N) Error("%s%s \n", Beads::ERROR,"ng is >= N in QSVal." );

 return CS+ GProb(ng,S)*(q*nextV_if_g + (1-q)*nextV_if_b) + BProb(ng,S)*((1-q)*nextV_if_g + q*nextV_if_b) ;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
float Beads::QSValME0(int ng, int S, float nextV_if_g, float nextV_if_b)
{
 CheckS(S);
 // In the following, forbid ng=N as then Sampling not an option.
 if (ng>=N) Error("%s%s \n", Beads::ERROR,"ng is >= N in QSValME0." );

 // In the ME versions, CS is absorbed in the 'prizes' available in the 
 // subsequent states rather than added separately:
 return GProb(ng,S)*(q*nextV_if_g + (1-q)*nextV_if_b) + BProb(ng,S)*((1-q)*nextV_if_g + q*nextV_if_b) ;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
float Beads::QSValME1(int ng, int S, float nextV_if_g, float nextV_if_b)
{
 CheckS(S);
 // In the following, forbid ng=meN as then Sampling not an option.
 if (ng>=meN) Error("%s%s \n", Beads::ERROR,"ng is >=meN in QSValME1." );

 // In the ME versions, CS is absorbed in the 'prizes' available in the 
 // subsequent states rather than added separately:
 return GProb(ng,S)*(q*nextV_if_g + (1-q)*nextV_if_b) + BProb(ng,S)*((1-q)*nextV_if_g + q*nextV_if_b) ;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // MOST IMPORTANT METHOD 1: FILL IN THE ARRAYS ACC. TO CURRENT VALUES OF
  // PARAMETRES, in the noise-free case:
void Beads::CalcPG_V_QG_QB_QS_AS()
{
  // NOT USED
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // MOST IMPORTANT METHOD 2: FILL IN THE ARRAYS ACC. TO CURRENT VALUES OF
  // PARAMETRES, in the SOFTMAX-NOISY case:
void Beads::SoftM_PG_V_QG_QB_QS_AS_MS(int zeroDrawChoice)
{
  int i,s;
  float aux1, aux2, aux3;
  float pG, pB, pS ;

  // First the last step, where there is no 'sample' choice:
  s=MaxSamplN();
  for (i=s; i>=0; i--) {  //all possible numbers of g draws
    PG(i+1,s-i+1) =GProb(i,s); // Row i+1 has i x g. Note
    QG(i+1,s-i+1) =aux1 =QlastG(i); //  +1 as indices OF
    QB(i+1,s-i+1) =aux2 =QlastB(i); // MATRICES from 1 not 0.
    PQSoft(i,s, &pG, &pB, &pS); 
    MS(i+1,s-i+1) = pS;
    V(i+1,s-i+1)  = pG*aux1 + pB*aux2 ;   
  }
  for (s=MaxSamplN()-1; s>=0; s--) { // remaining sampling steps, 
    // backwards.
    for (i=s; i>=0; i--) {       //all possible numbers of g draws
      PG(i+1,s-i+1) =GProb(i,s);  // Row i+1 has i x  g. 
      QG(i+1,s-i+1) =aux1 =QGVal(i,s); //+1 as indices OF
      QB(i+1,s-i+1) =aux2 =QBVal(i,s); //MATRICES from 1
      // Rem only non-last steps have valid QS (and AS) : 
      QS(i+1,s-i+1) = aux3 = QSVal(i,s,V(i+2,s-i+1),V(i+1,s-i+2));
      // Now find the greatest of the action-values and set the state Value:
      PQSoft(i,s, &pG, &pB, &pS); 
      MS(i+1,s-i+1) = pS;
      V(i+1,s-i+1)  = pG*aux1 + pB*aux2 + pS*aux3 ; 
      
      // Rem only non-last steps have valid QS (and AS) : 
      AS(i+1,s-i+1) = QS(i+1,s-i+1) - V(i+1,s-i+1) ;
      
    } // end loop over possible g - b combinations
  }  // end loop over sampling steps

  // Now correct the zero-draw elements if appropriate:
  if ( zeroDrawChoice == 0) {
    MS(1,1) = 1; // force to always sample
    // QS(1,1) doesn't need recalculation !
    V(1,1) = QS(1,1); // as S selected w. Py=1
    AS(1,1) = 0;
    QG(1,1) = QB(1,1) =-6.6E+6;  // Nonsense 'no go' values
  }
/*
    // DEBUG LINES
    cout<<"\n "<<" s:"<<s<<" ng:"<<i; // <<" pG:"<<pG<<" pB:"<<pB;   
    cout<<" PG:"<<PG(i+1,s-i+1)<<" QG:"<<QG(i+1,s-i+1)<<" QB:"<<QB(i+1,s-i+1);
    cout<<" QS:"<<QS(i+1,s-i+1);
    cout<<"\n"; 
*/
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // MOST IMPORTANT METHOD 3: FILL IN THE ARRAYS ACC. to Wael El-Deredy &
  // Charlotte Makintosh scheme ('quick and dirty' initial version) and 
  // TO CURRENT VALUES OF PARAMETRES, in the SOFTMAX-NOISY case:
void Beads::SoftME0_M_PG_V_QG_QB_QS_AS_MS(int zeroDrawChoice)
{
  // NOT USED
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // MOST IMPORTANT METHOD 4: FILL IN THE ARRAYS ACC. to Wael El-Deredy &
  // Charlotte Makintosh scheme (first 'clean' version - 2nd overall) and 
  // TO CURRENT VALUES OF PARAMETRES, in the SOFTMAX-NOISY case:
void Beads::SoftME1_M_PG_V_QG_QB_QS_AS_MS(int zeroDrawChoice)
{
  // NOT USED
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //  SOFTMAX-NOISY case FOR Sequential Prob. Ratio Test. It disallows 
  // declaring before any draws if zeroDrawChoice = 0 and fills in  the
  // PG and MS matrices meaningfully. Use the QG and QG to store the 
  // log-likelihood-ratios.  The rest are set to invalid flag -666.
void Beads:: SPRT_Soft_PG_and_MS(double phi_H, double L_over_H, int zeroDrawChoice)
{
  // NOT USED
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// return the label of the choice w. the greatest action val. at draw S,
// this time given number of green balls and total (i.e. not the sequence),
// but also fill the (deterministic) choice of best actions MG, MB, MS.
// if act. vals. equal, return 0 (for sample) or 1 (for say G). No randomness!
// Does various checks to ensure correct use as above.
int Beads::DeterMQmax(int ng, int s)
{
  if ( ng<0 || ng>N)
    Error("%s%s \n", Beads::ERROR,"ng out of range in DeterMQmax" );
  if ( s > N )
    Error("%s%s \n", Beads::ERROR,"s out of range in DeterMQmax" );

  float Qs, Qb, Qg ;             // auxiliaries for comparisons  
  int k=ng+1;  int l=s-ng+1;     // matrix indices that turn up a lot

  Qb=QB(k,l); 
  Qg=QG(k,l);
  Qs=QS(k,l); // these also hold for ng=0

  if ((Qs>=Qb && Qs>=Qg) && (ng<N)) {
    MS(k,l)=1.0;   /* MG(k,l)=0.0;   MB(k,l)=0.0;  // deterministic choice */
    return Slab ;                               // i.e. label for 'sample'
  }
  else                   // now Qs is excluded !
      if(Qg>=Qb) {
	MS(k,l)=0.0;  
	/* MG(k,l)=1.0;  
	MB(k,l)=0.0;  // determ. choice of G  */
	return Glab;
      }
      else {
	MS(k,l)=0.0;   /* MG(k,l)=0.0;   MB(k,l)=1.0;  // determ. choice of B */
	return Blab;
      }
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// return the label of the choice w. the greatest action val. at draw S.
// if act. vals. equal, return 0 (for sample) or 1 (for say G). No randomness!
// Does various checks to ensure correct use. May be slower !
int Beads::DeciQmax(int sampl, float *seqAr, int seqArEnd, int seqArStart)
{
  if ( sampl<0 || sampl>(seqArEnd-seqArStart+1))
    Error("%s%s \n", Beads::ERROR,"sampl num out of range in DeciQmax" );
  if ( seqArEnd-seqArStart+1 > N )
    Error("%s%s \n", Beads::ERROR,"sequence passed to Beads::DeciQmax too long." );
  int   s;
  int   ng = 0;      // important, will be retained if sampl=0.
  float nG = 0.0;    // also for counting, but float on purpose !
  float Qs, Qb, Qg ; // auxiliaries for comparisons  

  if (sampl > 0) {
    for (s=seqArStart; s<=sampl+seqArStart-1; s++)
      nG += *(seqAr+s);
    ng = roundfl(nG);
  }

  Qb=QB(ng+1,sampl-ng+1); 
  Qg=QG(ng+1,sampl-ng+1);
  Qs=QS(ng+1,sampl-ng+1); // these also hold for sampl=0.

  if ((Qs>=Qb && Qs>=Qg) && (sampl<seqArEnd-seqArStart+1)) 
    return Slab ;        // i.e. label for 'sample'
  else                   // now Qs is excluded !
    if(Qg>=Qb)
      return Glab;
    else
      return Blab;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // return the sample number where QS ceases to be the prefered option, i.e.
  // how many draws-to-declare. Tailored to classes as e.g. Garety91 (rl_exp.h)
  // *ASSUMES* THAT Q VALUES HAVE BEEN ALL FILLED IN VALIDLY !
int Beads::DeciDrawQmax(float *seqAr, int seqArEnd, int seqArStart)
{
  if ( seqArEnd-seqArStart+1 > N )
    Error("%s%s \n", Beads::ERROR,"sequence passed to Beads::DeciDrawQmax too long." );

  int   s=0;        int   ng = 0;     int   nb = 0;
  float Qs, Qg, Qb;

  while (s<seqArEnd-seqArStart+1) {
    
    // first, keep track of how many G and how many B :
    if (s>0) {    
      ng += roundfl( *(seqAr+seqArStart+s-1) );
      nb = s-ng ;
    }
    else
      ng = nb = 0;

    // Now look up Q's
    Qs=QS(ng+1,nb+1); // these also hold for sampl=0.
    Qb=QB(ng+1,nb+1); 
    Qg=QG(ng+1,nb+1);

    s+=1 ;
    if (Qs<=Qb || Qs<=Qg)  // i.e. Qs not greatest -> report 'declare here'
      return s-1;          // Qs were recorded for previous val. of s !
 
  } // end while looking for states where further sampling is an option.

  // ... so now we must be where further sampling is not an option -> return
  // 'error' or 'invalid' value as flag:
  return -666 ;

}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Return the probability that the decision to declare, rather than sample
  // again, is taken at draw dNum GIVEN THAT we start from draw 0 - i.e.
  // sample has been depleted !
float Beads::CumulPDec(int dNum, float *seqAr, int seqArEnd, int seqArStart)
{
  if ( dNum <0 || dNum>(seqArEnd-seqArStart+1))
    Error("%s%s \n", Beads::ERROR,"sampl num out of range in CumulPDec" );
  if ( seqArEnd-seqArStart+1 > N )
    Error("%s%s \n", Beads::ERROR,"sequence passed to Beads::CumulPDec too long." );

  int   s;
  int   ng = 0;      // important, will be retained if sampl=0.
  float nG = 0.0;    // also for counting, but float on purpose !
  float aux = 0.0 ;
  double cumP = 1;  //  'cumulative' probability. Important !

  if (dNum > 0) {
    cumP *= MS(1,1);  // Py of going past draw 0 !
    for (s=1; s<=dNum-1; s++) {          // up to draw BEFORE LAST, i.e.
      nG += *(seqAr+seqArStart-1+s);     // ... possibly not at all !
      ng = roundfl(nG);    // the number of G balls up to this step.
      cumP *= MS(ng+1,s-ng+1);
    }
    ng += roundfl( *(seqAr+seqArStart-1+dNum) );  // grand total ng, if dNum>=1.
  }
 
  // last step multiplicative factor is  different, and also works for dN==ng==0,
  // so is outside loop above. Should also work for dN==N, where the MS term is 0.
  cumP *= 1-MS(ng+1,dNum-ng+1); 
  
  aux = cumP;  // just to force conversion to float.
  return aux ;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#define outfsTb13 outfs<<'\t';outfs.width(13);outfs<<

void Beads::WriteSeqDeciPGQVPar(float *seqAr, int seqArEnd, int seqArStart) 
{
  // NOT USED
}
#undef outfsTb13

//------------------------------------------------------------------------------

// eof rl.cc
