/* File with methods to perform Expectation-minimisation and
   such like, inheriting the experiment stuff in jtc_exp.h

*/

#ifndef JTC_EM_H
#define JTC_EM_H

#include "jtc_exp.h"

// For multiple integration :
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
     
// for Gamma distro, incl. dummulative dist fns. :
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

// for special functions, esp digamma, trigamma:
#include <gsl/gsl_sf.h>


/* D.a. ---------------------------- first simple C declarations ----------- */

  // pdf of a Gamma distribution at q given its mean Gmean and its variance Gvar.
  // To rely on gsl necessities being already initialised in a ctor ...
  double gammaPDF(double q, double Gmean, double Gvar) ;

// Basic classes assembling bits for MC integration: --------------------------
class MCIntegr {
 protected:
    Array<double, 1> parArr; 
    int parNum;  // will be  from parArr size (i.e. rows() ...)
    // various bizarre bits that GSL random no generators etc. need:
    const gsl_rng_type *T;         gsl_rng *r;
    
    int dim, calls; 
    double cur_result, cur_error, cur_chisq;  // outputs of gsl, incl. chi square check
    double (*integrand)(double *q, size_t dim, void *param); // fn. itself HAS TO be 
                                                             // given externally.

 public:
    MCIntegr(); 
    // ~MCIntegr();
    
    void SetUpMCI(  double (*new_integrand)(double *q, size_t dim, void *param)
		   , size_t new_dim
		   , size_t new_calls = 100000 ) ;
    double vegasIntegrate(  double *xl  // calls the integrand with param = NULL !
			  , double *xu, int auxCalls = 10000 ) ;   
};

/* Classes : -------------------------------------------------------------   */

class Beads_Exp_II : public Beads_Exp_I
// See jtc_exp.h
{
 protected:
  static const char *const ERROR;

  int currModelType; //  model type, 0 for control, 1 for the 'real', basic 
                     // Bayesian, one, 2 for SPRT,3 for basic Mackintosh_El-Deredy ...
  int maxMicroParNum ; // THIS IS 3, BECAUSE Bayesian-based model can have a max. 
                       // of 3 micro-params, CS,CW and T, also covering poissonoid 
                       // control (1). SOME METHODS RETURN ERROR IF IT'S NOT 3 !
  int wrkDim;  // if 1, we are only working with 1 param - can go up to maxMicroParNum

  // for Monte-Carlo integration methods
  float ApproxInf; // can ignore contribution of Gamma pdf over this many times its SD.
  // for storage ... :
  double min_curr_Par;    // Don't start scannig from zero, quite !
  Array<double, 2> AuxPV; // AuxPV: sequence_Num x Probability_vector e.g:  3
                          //  (for Wellcome) x fN - usu. holds Py(d2d=drawNum)
  Array<double, 2> PVecSlice;  // Another auxiliary, where to copy a slice of PVec.

  Array<double, 1> PDumD_L;  // will have range 1 to subject_number - for DumDat work.
  double sum_ln_PDumD_L;     // (sum of ln of elements of above)

 public:
  // ctor initialises various bizarre bits that GSL random no generators etc. need, i.e.
  // const gsl_rng_type *T; gsl_rng *r; T = gsl_rng_default; r = gsl_rng_alloc (T);  :
  Beads_Exp_II();
  ~Beads_Exp_II();

  MCIntegr mcInt;      // Integration engine

  // copies of 'parameters' that will be needed for estimation of various probabili-
  // ty densities etc., esp. incl. recognition distro :
  //   Data point currently considered, by participant number and draw type. It 
  //   has slightly different meaning in different types of expt., Welcome vs. ME1
  //   etc :
  float curr_pt, curr_gp, curr_draw_type;   // integers, really, but float_ed to use w. menus.

  //   'macro' parameters currentl considered: for control poissonoid and SPRT
  //   and for Bayesian distr.:
  double currPoidPar_m, currPoidPar_v ;
  double currHiTh_m, currTh_v ;    // also to use currTm and currTv for 'noise'
  double currCSm,  currCSv,  currTm,   currTv ;
  double currCGBm, currCGBv, currCBGm, currCBGv ;
  //  Boundaries for integrations over hypercube :
  double *xl ;  // lower limit  NOTE THIS GOES FROM ZERO INDEX: ZERO *NOT* IGNORED
  double *xu ;  // upper limit  NOTE THIS GOES FROM ZERO INDEX: ZERO *NOT* IGNORED

  // Look-ups etc:
  int WDim(){ return wrkDim; };
  double PdumD_L(int i){ return PDumD_L(i); } ;
  double sum_ln_PdumD_L(){ return sum_ln_PDumD_L; };
  int currModType(){ return currModelType; };
  // Adjustments:
  void setCurrModType(double newModT){ currModelType=lround(newModT); };
  // The following don't just copy data from array MS etc. into AuxPV, but it also
  // converts from P[Sample if in State] to P[Sample AND State-is-reached].
  // First, for default, Welcome study:
  int StoreAuxPV(const RBWelcomeDat *RBWptr, int validZero=0 );

  double cProd3P_2(int d1, int d2, int d3, double L=-666.0) { 
    return  CPDistr1(d1,L)*CPDistr1(d2,L)*CPDistr1(d3,L); } ;

  // Various bits for dummy data (can be not so dummy at all, but this one, DumDat,
  // contains mainly decisions-to-draw data.
  Array<float,2> DumDat ;   // has ranges starting from 1.
 
  // the following produces output formatted as 
  // ptNum   CS    data1    data2 ...  data_datN   kT   , with datN usu. = 3
  // this does NOT allow for decision to be made at zero draws.
  void fillDumDatBays(RBWelcomeDat *RBWptr, int subjN=30, long int rndSeed = 0);  

  // Set ... :
  // set xl to zero and xu to respective mean+ApproxInf*SD acc. to the model type,
  // 0 for control, 1 for the 'real', basic Bayesian, one, 2 for SPRT,
  // 3 for basic Mackintosh_El-Deredy ...
  void SetApproxInf(int model_type = 0);
  // Set curr.. variables, then run SetApproxInf. model_types as above. In the case
  // of more complicated parameter structures, arrays can be passed with the pointers
  // e.g. *newCWx can be (CGBx, CBGx) :
  void SetCurrAndBounds( const double *newTm,  const double *newTv, const double *newV2m, const double *newV2v
                        ,const double *newCWm=NULL, const double *newCWv=NULL , const double *newPParm=NULL
                        ,const double *newPParv=NULL, int new_model_type = 1 );
  // reset no. of subjects & what depends on it:
  void ResetAndSubjN(int new_subjN, int datN=3){
   PDumD_L.resize(Range(1,new_subjN)); 
   PDumD_L=0;
   DumDat.resize(Range(1,new_subjN),Range(1,6));
   DumDat = 0; } ;
  void ResetPDumD_L_etc2(int new_subjN); // resize & zero PDumD_L 

  // laborious stuff ! 1st, update PD_L for data in DumDat, and return the
  // updated sum_ln_PDumD_L. Modeltype=1 is Bayesian, 2 is SPRT:
  double update_PDumD_L( double (*new_PDF)(double *q, size_t dim, void *param),
                        int mcCalls=10000, int model_type=1, int verbose= 1);
  // Estimate (return new mean but don't update) currTm and currCSm :
  // RELY ON UPDATED PdumD_L !
  double estCurrDumTm( double (*new_PDF)(double *q, size_t dim, void *param),
                        int mcCalls=10000, int model_type=1, int verbose = 1);
  double estCurrDumCSm( double (*new_PDF)(double *q, size_t dim, void *param),
                        int mcCalls=10000, int verbose = 1);

  // Estimate (return vals. but don't update) on way to currTv and currCSv :
  // RELY ON UPDATED PdumD_L AND ON UPDATED currTm AND currCSm !
  double estCurrDum_QlnT(double (*new_PDF)(double *q, size_t dim, void *param),
                          int mcCalls=10000, int model_type=1, int verbose = 1);
  double estCurrDum_QlnCS(double (*new_PDF)(double *q, size_t dim, void *param),
                          int mcCalls=10000, int verbose = 1);
  double estCurrDum_QlnTh(double (*new_PDF)(double *q, size_t dim, void *param),
                          int mcCalls=10000, int verbose = 1);
  double estCurrDum_QlnCGB(double (*new_PDF)(double *q, size_t dim, void *param),
                          int mcCalls=10000, int verbose = 1); // for ME1 +

  // THE FOLLOWING TWO ARE SUBOPTIMAL W.R.T. EM THEORY !
  double estCurrDumTv( double (*new_PDF)(double *q, size_t dim, void *param),
                        int mcCalls=10000, int verbose = 1);
  double estCurrDumCSv( double (*new_PDF)(double *q, size_t dim, void *param),
                        int mcCalls=10000, int verbose = 1);

  // Interpolation machinery :
  // These two prepare PVecSlice and call setRegInp2D_x and setInp2D_y :
  void prepIntBasics(int mod_type);
  void prepPVecSlice(int d2d, int seq_num, int mod_type=1);	
  double wrkP(int trial, int step_to_decl){ return AuxPV(trial, step_to_decl); };
  double wrkP2(int trial, int step_to_decl){  // Costed Bayesian version
    return IP2.PolynInt(kT,CS, PVecs(Range::all(),Range::all(),trial,step_to_decl)); };
  // M phase related methods:
  // 1. Use averaged integration results, obtained separately and including 
  //    new estimate of mean of gamma distr., to derive new estimate of variance
  //    of gamma distr. Param. NRiterNum is number of Newton-Raphson iterations 
  //    used in final phase of calculation; frLastCorr is a rough error estimate
  //    (see code) :
  double xVarEst(  double IntQx, double IntQlnx, int NRiterNum=10
                 ,  double *frLastCorr=NULL, int verbose = 1);
};
// Don't forget the semilcolon at the end of class declarations !
// = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 

#endif    // JTC_EM_H

// eof
