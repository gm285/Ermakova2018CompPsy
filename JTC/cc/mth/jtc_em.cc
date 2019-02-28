// methods for jtc_em.h, EM and other methods to deal with rbw data etc.
/* Some default val. in C quotes  */  

#include "jtc_em.h"

const char *const Beads_Exp_II::ERROR="  >>---> Beads_Exp_II error : " ;

// FIRST, NON-CPP STUFF (Plain C) to link w. gsl Monte - - - - - - - - - - - - - 
double gammaPDF(double q, double Gmean, double Gvar) 
{
  // negative Gmean -> use mirror-image distro :
  if (Gmean<0) { Gmean *= -1.0 ;  q *= -1.0; }
  double Th =  Gvar/Gmean ;
  double K  =  Gmean/Th ;
  return gsl_ran_gamma_pdf (q, K, Th);
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
MCIntegr::MCIntegr()
{
  // set up the gsl random number generation machinery:
  gsl_rng_env_setup ();
  T = gsl_rng_default;     r = gsl_rng_alloc (T);

  integrand = NULL ;
  parNum = 0;
  dim = calls = 0;    // unless set up separately, nothing works !

}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
void MCIntegr::SetUpMCI(  double (*new_integrand)(double *q, size_t dim, void *param)
                        , size_t new_dim
                        , size_t new_calls )
{
  dim = new_dim;
  calls = new_calls ;
  // parArr.resize(newParArr.shape());
  // parArr = newParArr ;

  integrand = new_integrand ;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
double  MCIntegr::vegasIntegrate(  double *xl, double *xu, int auxCalls )
// Uses an initial warm-up run of calls/20 function calls to prepare, or 'warm up', the grid.
// This is followed by a main run with five iterations of calls/5 function calls. 
// The chi-squared per degree of freedom for the five iterations are checked for consistency 
// with 1, and the run is repeated if the results have not converged.
{
  int new_calls; 
  int do_again;
  new_calls = ( (auxCalls==0) ? calls : auxCalls ) ;

  gsl_monte_function G = { integrand, dim , 0 };
  {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (dim );   
    gsl_monte_vegas_integrate (&G, xl, xu, dim, new_calls/20, r, s,
                               &cur_result, &cur_error);
    //cout << "\nVegas warm-up result:"<<cur_result<<"\t error:"<< cur_error; 
    //cout <<"\nconverging...\n"; 
    do
      {
        gsl_monte_vegas_integrate (&G, xl, xu, dim, new_calls/5, r, s,
                                   &cur_result, &cur_error);
        // printf ("result = % .6f sigma = % .6f "
        //   "chisq/dof = %.1f\n", cur_result, cur_error, (cur_chisq=s->chisq));
        do_again =  (fabs (s->chisq - 1.0) > 0.5) ? 1 : 0;
        if (do_again) {
          cerr <<"Vegas converging (res="<<cur_result<<", err="<<cur_error;
          cerr <<"  chisq/dof ="<<(cur_chisq=s->chisq)<<" ) \n";
        }
      }
    while  (do_again);   
    // cout << "\nVegas final result:"<<cur_result<<"\t error:"<< cur_error;    
    gsl_monte_vegas_free (s);
  }
  return cur_result;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
Beads_Exp_II::Beads_Exp_II() : mcInt() 
 // (int dimSize1, int dimSize2, int dimSize3=1, int wrkPDim=2)
{
  int i;
  ApproxInf = 6.0 ;  // ignore contribution of Gamma pdf over ApproxInf * SD.
  maxMicroParNum = 3 ;
  wrkDim = 2;
  min_curr_Par = 0.0001 ;             // Will it manage this ? ! was 0.001 OK.
  xl = new double[maxMicroParNum] ;
  xu = new double[maxMicroParNum] ;
  for (i=0; i<maxMicroParNum; i++)  xl[i]= xu[i]=0.0 ;

  // the following are defaults, not nonsense flags.
  curr_pt = 66;  curr_draw_type = 1.0; 
  // the following are meant to be 'not unreasonable' for pt. 81 of Wellcome I :
  currPoidPar_m = 3.15; currPoidPar_v = 0.25 ;
  // currCSm = -1.77  ;    currCSv = 0.16 ; 
  currCSm = -1.77  ;    currCSv = 0.16 ;  currTm  = 2.00   ;    currTv  = 0.01 ;
  currCGBm = -100.0;    currCGBv = 49.0 ; currCBGm = -100.0;    currCBGv = 49.0 ;

  AuxPV.resize(Range(1,seqNum),Range(0,seqLen));
  PVecSlice.resize(Range(1,roundfl(kT_N)), Range(1,roundfl(CS_N))) ;
} 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
Beads_Exp_II::~Beads_Exp_II()
{
 if (xl) delete xl; 
 if (xu) delete xu;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
void Beads_Exp_II::ResetPDumD_L_etc2(int new_subjN){ 
  PDumD_L.resize(Range(1,new_subjN)); 
  PDumD_L=0;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
void Beads_Exp_II::fillDumDatBays(RBWelcomeDat *RBWptr, int subjN, long int rndSeed)
{ 
  int i,j,k,l ;
  int datN = RBWptr->seqN();
  float *SumP; SumP = new float[21];
  float auxP, PPar ;
  // ofstream dumDat_out ;   dumDat_out.open(fname, ios::out);
  double Psum, Pcur;
  double ThCS, K_CS, ThT, K_T;  //,ThCGB, C_CGB  if nec...

  ThCS = currCSv/currCSm;     K_CS=currCSm/ThCS;
  ThT  = currTv/currTm;       K_T =currTm/ThT;    
  // ThCGB = currCGBv/currCGBm;     C_CGB=currCGBm/ThCGB;   // if nec... 

  DumDat.resize(Range(1,subjN), Range(1,datN+3));
  DumDat = 0;

  if( rndSeed == 0)  // if no specific seed given ...
    rndSeed = ( time( NULL ) % 604800) ;
  // various bizarre bits for GSL random no generators 
  const gsl_rng_type *T;    gsl_rng *r; T = gsl_rng_default; 
  r = gsl_rng_alloc (T);    gsl_rng_set (r,  rndSeed );

  for (k=1; k<=subjN; k++) {
    // for each subject, set the parametres of the base classes to the
    // desired values so that the base classes methods can be used in a mo:
    CS = gsl_ran_gamma(r, K_CS, ThCS) ; // each subject has its own parametres
    kT = gsl_ran_gamma(r, K_T,  ThT ) ; 
    // cout <<'\n'<< k <<'\t'<< CS ;
    DumDat(k,1) = k; DumDat(k,2) = CS;  // other params after generated data ...

    SoftM_PG_V_QG_QB_QS_AS_MS();     // update all the probability tables
    StoreAuxPV(RBWptr,0); 
    // write2DArray_csv("PVecs.csv",PVecs(1,1,Range::all(),Range::all())); // debug

    for (i=1; i<= datN; i++) {
      Psum =  Pcur = 0.0 ;
      for (j=1; j<=20; j++) {  // each trial AND each subject has own mdf

        Pcur =  AuxPV(i,j);     Psum += Pcur;      SumP[j] = Psum;
      }
      l = 1;
      auxP = gsl_ran_flat( r, 0, 1) ;
      while ( (auxP > SumP[l]) && (l<20) ) l+= 1 ;
      // cout<<'\t'<<l ;
      DumDat(k,2+i) = l ;  // The generated data !
    }    
    DumDat(k,datN+3) = kT ;  // rest of parametres, for interest etc ...

  } // end loop over subjects

  // cout <<'\n'<< DumDat <<'\n' ;
  delete SumP;
  gsl_rng_free (r); 
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
void Beads_Exp_II::SetApproxInf(int model_type)
{
  int i ; double sd;
  for (i=0; i<maxMicroParNum; i++) xl[i] = 0.0;
  if (model_type == 0) {
    xu[0] = currPoidPar_m + ApproxInf*sqrt( currPoidPar_v ) ;
    xu[1] = 1.0;   // so that this dim. can be used as dummy.
    xu[2] = 0.0; 
  }
  else if ((model_type == 1) || (model_type == 4))  { // basic Costed Bayesian model OR ...
                                                      // ... ME_expt1-model_ME2
    // boundaries for kT :
    if (currTm<0)
      Error("%s%s \n",Beads_Exp_II::ERROR,"Can't cope with negative noise param.");
    sd = sqrt(currTv) ;           xu[0] = currTm  + ApproxInf*sd ;
    xl[0] = ((currTm-ApproxInf*sd < min_curr_Par) ? min_curr_Par : currTm-ApproxInf*sd ) ;

    // boundaries for CS :
    if (currCSm >= 0) {
      sd = sqrt(currCSv) ;        xu[1] = currCSm + ApproxInf*sd ;
      xl[1] = ((currCSm-ApproxInf*sd < min_curr_Par ) ? min_curr_Par : currCSm-ApproxInf*sd ) ;
    }
    else {  // i.e. if currCSm < 0
      sd = sqrt(currCSv) ;        xl[1] = currCSm - ApproxInf*sd ;
      xu[1] = ((currCSm+ApproxInf*sd > -min_curr_Par ) ? -min_curr_Par : currCSm+ApproxInf*sd ) ;
    }

    // boundaries for CGB :
    if (currCGBm >= 0) {
      sd = sqrt(currCGBv) ;       xu[2] = currCGBm + ApproxInf*sd ;
      xl[2] = ((currCGBm-ApproxInf*sd < min_curr_Par) ? min_curr_Par : currCGBm-ApproxInf*sd ) ;
    }
    else {
      sd = sqrt(currCGBv) ;       xl[2] = currCGBm - ApproxInf*sd ;
      xu[2] = ((currCGBm+ApproxInf*sd > -min_curr_Par) ? -min_curr_Par : currCGBm+ApproxInf*sd ) ;
    }
  } // end if  (model_type == 1),  basic Costed Bayesian model, or ME2

  else if (model_type == 2) {  // i.e. SPRT model
    min_curr_Par = 0.005;  // Huge !
    // boundaries for kT, which is used as the SD/st. error of the increments :
    if (currTm<0) // *1.001<min_curr_Par) 
      Error("%s%s \n",Beads_Exp_II::ERROR,"Can't cope with such low/-ve currTm (increment SD).");
    sd = sqrt(currTv) ;
    xu[0] = currTm  + ApproxInf*sd ;
    xl[0] = ((currTm-ApproxInf*sd < min_curr_Par) ? min_curr_Par : currTm-ApproxInf*sd ) ;

    // boundaries for Threshold_hi, i.e. CDPar1 :
    // Hi Threshold in SPRT is derived from a positive log-lik-ratio, so currHiTh_m must be >=0,
    // so no need for 'if (currHiTh_m>= 0), but may need the following:
    if ( currHiTh_m < min_curr_Par)
      cerr << "\nWARNING -  currHiTh_m ("<< currHiTh_m<<")has fallen below min_curr_Par ("<<min_curr_Par<<")\n";
    sd = sqrt(currTh_v) ;
    xu[1] = currHiTh_m + ApproxInf*sd ;
    xl[1] = ((currHiTh_m-ApproxInf*sd < min_curr_Par) ? min_curr_Par : currHiTh_m-ApproxInf*sd ) ;
    // } if for different ranges of currHiTh_m not needed ...

  } // end if  (model_type == 2), i.e. SPRT model

  else if (model_type == 3) { //  == 3, i.e.  basic Mackintosh_El-Deredy model ME1
    // boundaries for kT :
    if (currTm<0)
      Error("%s%s \n",Beads_Exp_II::ERROR,"Can't cope with negative noise param.");
    sd = sqrt(currTv) ;
    xu[0] = currTm  + ApproxInf*sd ;
    xl[0] = ((currTm-ApproxInf*sd < min_curr_Par) ? min_curr_Par : currTm-ApproxInf*sd ) ;

    // boundaries for CGB (== CGB, set when needed in other methods) :
    if (currCGBm >= 0) {
      sd = sqrt(currCGBv) ;
      xu[1] = currCGBm + ApproxInf*sd ;
      xl[1] = ((currCGBm-ApproxInf*sd < min_curr_Par) ? min_curr_Par : currCGBm-ApproxInf*sd ) ;
    }
    else {
      sd = sqrt(currCGBv) ;
      xl[1] = currCGBm - ApproxInf*sd ;
      xu[1] = ((currCGBm+ApproxInf*sd > -min_curr_Par) ? -min_curr_Par : currCGBm+ApproxInf*sd ) ;
    }

    // debug: cout<<"\n sd = sqrt(currCGBv) = "<<sd<<"\t CGB from "<< xl[1]<<" to "<< xu[1]<<'\n';
    
    //  CS is not really to be integrated over in basic ME version, 
    //  so sd used above here is essentially set to zero. Expect currCSm also to be zero.
    xl[2] = currCSm  ;        xu[2] = currCSm  ;

  } // end if  (model_type == 3), i.e.  basic Mackintosh_El-Deredy model ME1


  else
    Error("%s%s \n",Beads_Exp_II::ERROR,"SetApproxInf not ready for inputted model type.");    

}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// Rem model_types: 0 for control, 1 for the 'real', e.g.  MDP, 2 for SPRT,
//                  3 for simple Mackintosh_Elderedy ME1 ...
void Beads_Exp_II::SetCurrAndBounds( const double *newTm, const double *newTv, const double *newV2m
                                   , const double *newV2v,const double *newCWm,const double *newCWv
                                   , const double *newPParm, const double *newPParv
                                   , int new_model_type )
{
  // All model types are allowed a noise term:
  if (newTm) currTm = *newTm;         if (newTv) currTv = *newTv; 
  // All are allowed 'internal cost of getting it wrong' terms too, tho
  // in basic Bayesian and SPRT this CW will be a fixed reference, whereas
  // in ME1 it will be central:
  if (newCWm) currCGBm = *newCWm;     if (newCWv) currCGBv = *newCWv; 
  // debug: cout<<"\n currCGBm = *newCWm = "<<currCGBm<<" currCGBv = *newCWv = "<<currCGBv<<'\n';

  // V2m can be used for threshold or for sampling costs in 
  // SPRT or basic Bayesian:
  if (new_model_type == 1) {
    if (newV2m) currCSm = *newV2m;     if (newV2v) currCSv = *newV2v; 
  }
  else if (new_model_type == 2) {
    if (newV2m) currHiTh_m = *newV2m;  if (newV2v) currTh_v = *newV2v; 
  }

  // Poissonoid 'control / testing' model:
  if (newPParm) currPoidPar_m = *newPParm;  if (newPParv) currPoidPar_v = *newPParv;    

  SetApproxInf(new_model_type);
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
int Beads_Exp_II::StoreAuxPV(const RBWelcomeDat *RBWptr, int validZero )
{ 
  // Rem:  AuxPV is 2-D array: sequence_Num x Probability_vector
  //                              e.g:  3 (for Wellcome) x fN
  // if key extends don't match, abort:
  if( (RBWptr->seqN()) != AuxPV.extent(firstDim))
    Error("%s seqN(%d) != extent(firstDim) (%d) in StoreAuxPV.",
          Beads_Exp_II::ERROR, RBWptr->seqN(), AuxPV.extent(firstDim));
  if ((RBWptr->SamplMaxN()) != (fN-1))
      Error("%s%s \n",Beads_Exp_II::ERROR,"max. draw number mismatch in StoreAuxPV.");

  int seq_n, s, initDraw, endDraw, ng;                      // auxiliaries
  double  prod;

  initDraw = AuxPV.lbound(secondDim);  endDraw=AuxPV.ubound(secondDim); 

  for (seq_n = 1; seq_n <= (RBWptr->seqN()); seq_n++) {
    // First P(d2d=0), by hand :
    if (validZero){  // If subject allowed to declare at 0 draws ...
      AuxPV(seq_n,initDraw) = 1-MS(1,1) ;
      prod = MS(1,1);  // zero-draws 'sample again' Py
    }
    else {
      AuxPV(seq_n,initDraw) = 0.0;
      prod = 1.0 ;  // zero-draws OBLIGATORY 'sample again' Py !
    } 
    // now P(d2d=1) to  P(d2d=SamplMaxN()) :
    for (s=initDraw+1; s<=endDraw; s++) {
      ng = lround(RBWptr->ng_tot[seq_n][s]) ;
      // the following should work for last step too:
      AuxPV(seq_n, s) = prod*(1-MS(ng+1,s-ng+1)); 
      prod *= MS(ng+1,s-ng+1);               // ready for next step
    }
  }
  return 0; 
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
double Beads_Exp_II::update_PDumD_L(double (*new_PDF)(double *q, size_t dim, void *param),
                                  int mcCalls, int model_type, int verbose)
{
  // int model_type = 1; initially was just for Bayesian model ...
  int i;
  double aux0;
  if (verbose) cout<<"\n Now working in Beads_Exp_II::update_PDumD_L :";
  sum_ln_PDumD_L = 0.0; 
  SetApproxInf(model_type); 
  for (i=1; i<=DumDat.rows(); i++) {
     curr_pt = i;
     mcInt.SetUpMCI(new_PDF, wrkDim) ; // this hopefully uses the above curr_pt !
     
     if ((verbose) && (i==1)) {
       if (model_type==1)
         cout <<"\n Parm. for update_PDumD_L: CS_mean: "<<currCSm<<" ; CS_var: "<< currCSv;
       else if (model_type==2)
         cout <<"\n Parm. for update_PDumD_L: Thresh_mean: "<<currHiTh_m<<" ; Thresh_var: "<< currTh_v ;
       else if (model_type==3)
         cout <<"\n Parm. for update_PDumD_L: CGB_mean: "<<currCGBm<<" ; CGB_var: "<< currCGBv ;
       // common to model types: error / uncertainty / temperature parameter :
       cout <<"; kT_mean : "<<currTm<<" ; kT_var: "<< currTv<<endl;
       // sumRes(nIter+2,1) = 0; sumRes(nIter+2,2)=GB.currCSm; sumRes(nIter+2,3)=GB.currCSv; 
       // sumRes(nIter+2,4) = GB.currTm;  sumRes(nIter+2,5) = GB.currTv; 
       cout <<"Number of pts (total) = "<< DumDat.rows()<<endl;
       cout<<"Pt.N :\tP[D ; L] :\tln(P[D ; L]) :"<<endl; // header for what follows later.
     } // end initial verbose display
     aux0 = mcInt.vegasIntegrate(xl, xu, mcCalls);
     PDumD_L(i) = aux0;
     sum_ln_PDumD_L += log(aux0);
     if (verbose) 
       cout<<i<<'\t'<<aux0<<'\t'<<log(aux0)<<endl; 
  } // end loop over i, i.e. rows of DumDat.
  // sumRes(nIter+2,6)=sum_ln_PDumD_L;
  return sum_ln_PDumD_L;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
double Beads_Exp_II::estCurrDumTm(double (*new_PDF)(double *q, size_t dim, void *param),
                                  int mcCalls, int model_type, int verbose)
{
  // int model_type = 1; was for Bayesian model, at least initially ...
  int i;  double aux0, tempTm ;   tempTm = 0.0; 
  if (verbose) cout<<"\n Running Beads_Exp_II::updateCurrDumTm (RELIES ON PdumD_L):\n";
  SetApproxInf(model_type); 
  for (i=1; i<=DumDat.rows(); i++) {
     curr_pt = i;
     mcInt.SetUpMCI(new_PDF, wrkDim) ; // this hopefully uses the above curr_pt !    
     aux0 = mcInt.vegasIntegrate(xl, xu, mcCalls);
     tempTm += aux0/PdumD_L(i);
     if (verbose)
       cout<<"Tm est ("<<i<<")="<<aux0/PdumD_L(i)<<"     ... done : "<<i<<"/"<<DumDat.rows()<<endl; // debug
   }
 
  return tempTm/DumDat.rows();
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
double Beads_Exp_II::estCurrDumCSm(double (*new_PDF)(double *q, size_t dim, void *param),
                                  int mcCalls, int verbose)
{
  int model_type = 1; // for Bayesian model, at least initially ...
  int i;  double aux0, tempCSm ;   tempCSm = 0.0; 
  if (verbose) cout<<"\n Running Beads_Exp_II::updateCurrDumCSm (RELIES ON PdumD_L):\n";
  SetApproxInf(model_type); 
  for (i=1; i<=DumDat.rows(); i++) {
     curr_pt = i;
     mcInt.SetUpMCI(new_PDF, wrkDim) ; // this hopefully uses the above curr_pt !    
     aux0 = mcInt.vegasIntegrate(xl, xu, mcCalls);
     tempCSm += aux0/PdumD_L(i);
     if (verbose)
       cout<<"CSm est ("<<i<<")="<<aux0/PdumD_L(i)<<"     ... done : "<<i<<"/"<<DumDat.rows()<<endl; // debug
   }
 
  return tempCSm/DumDat.rows();
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
double Beads_Exp_II::estCurrDum_QlnT(double (*new_PDF)(double *q, size_t dim, void *param),
                                       int mcCalls, int model_type, int verbose)
{
  // int model_type = 1; was for Bayesian model, at least initially ...
  int i;  double aux0, tempIntQlnT ;   tempIntQlnT = 0.0; 
  if (verbose) cout<<"\n Running Beads_Exp_II::estCurrDum_QlnT (RELIES ON PdumD_L and currXX):\n";
  SetApproxInf(model_type); 
  for (i=1; i<=DumDat.rows(); i++) {
     curr_pt = i;
     // the following hopefully uses the above curr_pt; new_PDF SHOULD make use of currXX  !
     mcInt.SetUpMCI(new_PDF, wrkDim) ;   
     aux0 = mcInt.vegasIntegrate(xl, xu, mcCalls);
     tempIntQlnT += aux0/PdumD_L(i);
     if (verbose)
       cout<<"QlnT est ("<<i<<")="<<aux0/PdumD_L(i)<<"     ... done : "<<i<<"/"<<DumDat.rows()<<endl; // debug
   }
 
  return tempIntQlnT/DumDat.rows();
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
double Beads_Exp_II::estCurrDum_QlnCS(double (*new_PDF)(double *q, size_t dim, void *param),
                                  int mcCalls, int verbose)
{
  int model_type = 1; // for Bayesian model, at least initially ...
  int i;  double aux0, tempIntQlnCS ;   tempIntQlnCS = 0.0; 
  if (verbose) cout<<"\n Running Beads_Exp_II::estCurrDum_QlnCS (RELIES ON PdumD_L AND, implicitly, currXX):\n";
  SetApproxInf(model_type); 
  for (i=1; i<=DumDat.rows(); i++) {
     curr_pt = i;
     // the following hopefully uses the above curr_pt; new_PDF SHOULD make use of currXX  !
     mcInt.SetUpMCI(new_PDF, wrkDim) ;   
     aux0 = mcInt.vegasIntegrate(xl, xu, mcCalls);
     tempIntQlnCS += aux0/PdumD_L(i);
     if (verbose)
       cout<<"QlnCS est ("<<i<<")="<<aux0/PdumD_L(i)<<"     ... done : "<<i<<"/"<<DumDat.rows()<<endl; // debug
   }
 
  return tempIntQlnCS/DumDat.rows();
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
double Beads_Exp_II::estCurrDum_QlnCGB(double (*new_PDF)(double *q, size_t dim, void *param),
                                  int mcCalls, int verbose)
{
  int model_type = 3; // for Mackintosh - El-Deredy Bayesian model, at least initially ...
  int i;  double aux0, tempIntQlnCGB ;   tempIntQlnCGB = 0.0; 
  if (verbose) cout<<"\n Running Beads_Exp_II::estCurrDum_QlnCGB (RELIES ON PdumD_L AND, implicitly, currXX):\n";
  SetApproxInf(model_type); 
  for (i=1; i<=DumDat.rows(); i++) {
     curr_pt = i;
     // the following hopefully uses the above curr_pt; new_PDF SHOULD make use of currXX  !
     mcInt.SetUpMCI(new_PDF, wrkDim) ;   
     aux0 = mcInt.vegasIntegrate(xl, xu, mcCalls);
     tempIntQlnCGB += aux0/PdumD_L(i);
     if (verbose)
       cout<<"QlnCGB est ("<<i<<")="<<aux0/PdumD_L(i)<<"     ... done : "<<i<<"/"<<DumDat.rows()<<endl; // debug
   }
 
  return tempIntQlnCGB/DumDat.rows();
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
double Beads_Exp_II::estCurrDum_QlnTh(double (*new_PDF)(double *q, size_t dim, void *param),
                                  int mcCalls, int verbose)
{
  int model_type = 2; // for SPRT model only.
  int i;  double aux0, tempIntQlnTh ;   tempIntQlnTh = 0.0; 
  if (verbose) cout<<"\n Running Beads_Exp_II::estCurrDum_QlnTh (RELIES ON PdumD_L AND, implicitly, currXX):\n";
  SetApproxInf(model_type); 
  for (i=1; i<=DumDat.rows(); i++) {
     curr_pt = i;
     // the following hopefully uses the above curr_pt; new_PDF SHOULD make use of currXX  !
     mcInt.SetUpMCI(new_PDF, wrkDim) ;   
     aux0 = mcInt.vegasIntegrate(xl, xu, mcCalls);
     tempIntQlnTh += aux0/PdumD_L(i);
     if (verbose)
       cout<<"QlnTh est ("<<i<<")="<<aux0/PdumD_L(i)<<"     ... done : "<<i<<"/"<<DumDat.rows()<<endl; // debug
   }
 
  return tempIntQlnTh/DumDat.rows();
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
double Beads_Exp_II::estCurrDumCSv(double (*new_PDF)(double *q, size_t dim, void *param),
                                  int mcCalls, int verbose)
{
  int model_type = 1; // for Bayesian model, at least initially ...
  int i;  double aux0, tempCSv ;   tempCSv = 0.0; 
  if (verbose) cout<<"\n Running Beads_Exp_II::updateCurrDumCSv (RELIES ON PdumD_L AND, implicitly, currCSm):\n";
  SetApproxInf(model_type); 
  for (i=1; i<=DumDat.rows(); i++) {
     curr_pt = i;
     // the following hopefully uses the above curr_pt; new_PDF SHOULD make use of currCSm  !
     mcInt.SetUpMCI(new_PDF, wrkDim) ;   
     aux0 = mcInt.vegasIntegrate(xl, xu, mcCalls);
     tempCSv += aux0/PdumD_L(i);
     if (verbose)
       cout<<"CSv est ("<<i<<")="<<aux0/PdumD_L(i)<<"     ... done : "<<i<<"/"<<DumDat.rows()<<endl; // debug
   }
 
  return tempCSv/DumDat.rows();
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
void Beads_Exp_II::prepIntBasics(int mod_type)
{
  // Rem:  PVecs is 4-D array: param1 x param2 x sequence_Num x Probability_vector
  //       e.g. for SPRT: kT_N x CDPar1_N (usu. in iterPar 6) x 3 x fN
  SetApproxInf(mod_type); 
  int N1, N2;               N1 = lround(kT_N);     // for most !
  if (mod_type == 1)        N2 = lround(CS_N);     // for basic costed Bayesian
  else if (mod_type == 2)   N2 = lround(CDPar1_N); // for SPRT
  else if (mod_type == 3)   N2 = lround(CGB_N);    // for ME1 Bayesian 
  else if (mod_type == 4)   N2 = lround(CS_N);     // for ME2 Bayesian MODEL (exp ME1)
  else err_msg("Beads_Exp_II::prepIntBasics : Not ready for mod_type value");
  IP2.setRegInp2D_x(xl[0],(xu[0]-xl[0])/(N1-1),xl[1],(xu[1]-xl[1])/(N2-1) ) ; 
  // The following contain no data, but are needed:
  PVecSlice.resize(Range(1,N1), Range(1,N2)) ;
  IP2.setInp2D_y( &PVecSlice ); 
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
void Beads_Exp_II::prepPVecSlice(int d2d, int seq_num, int mod_type)
{
  if (d2d*(d2d-20) > 0) err_msg("d2d off range in Beads_Exp_II::prepPVecSlice");
  if (mod_type != 3) // ME1 typ e models have special third dimension !
    { if ((seq_num-1)*(seq_num - seqNum) > 0) 
        err_msg("seq_num off range in Beads_Exp_II::prepPVecSlice");  }
  else
    { if ((seq_num-1)*(seq_num - seqNum*retNum) > 0) 
        err_msg("seq_num off range in Beads_Exp_II::prepPVecSlice (mod_type=3)");  }


  int N1, N2;               N1 = lround(kT_N);  
  if (mod_type == 1)        N2 = lround(CS_N);    // for basic costed Bayesian
  else if (mod_type == 2)   N2 = lround(CDPar1_N);   // for SPRT
  else if (mod_type == 3)   N2 = lround(CGB_N);      // for ME1 Bayesian
  else err_msg("Beads_Exp_II::prepPVecSlice : Not ready for mod_type value");

  IP2.setRegInp2D_x(xl[0],(xu[0]-xl[0])/(N1-1),xl[1],(xu[1]-xl[1])/(N2-1) ) ; 
  PVecSlice.resize(Range(1,N1), Range(1,N2)) ;
  PVecSlice=PVecs( Range::all(),Range::all(),seq_num, d2d ) ;
  IP2.setInp2D_y( &PVecSlice ); 
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
double Beads_Exp_II::xVarEst(  double IntQx, double IntQlnx, int NRiterNum
                             , double *frLastCorr, int verbose)
{
  double s = log(IntQx) - IntQlnx ;
  double b = s-3.0;
  double k0, kappa, dk;
  // Initial approximation:
  kappa = k0 = ( sqrt(b*b+24*s) - b)/(s*12) ;;

  // Newton-Raphson :
  for (int i=1; i<= NRiterNum; i++) {
    dk = (log(kappa)- gsl_sf_psi(kappa) - s)/((1/kappa) - gsl_sf_psi_1(kappa)) ;
    kappa -= dk ;
  }

  // Diagnostics :
  if (verbose) {
    cerr<<"\n xVarEst:\n IntQx="<<IntQx<<"\tIntQlnx="<<IntQlnx;
    cerr<<"    s="<<s    <<"\t      b="<<b<<"\t init. kappa"<<k0;
    cerr<<"\tfinal kappa="<<kappa<<" Last correction to kappa:"<<dk<<endl;
  }

  // Output:
  if (frLastCorr != NULL) *frLastCorr = dk/kappa ;  // last corrn. as fraction of Kappa
  return IntQx*IntQx/kappa ;  // return the variance, not kappa itself.
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
//=============================================================================
// end of File rl_em.cc
