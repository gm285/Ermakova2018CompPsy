/*  Like 06, but version for Words or Beads.

    Using interpolation for faster MC integrations to implement
    better M phase formulae, as of Apr 08.
  
    Again two 'jars', B and G, which have a higher proportion of b and g
    type of beads. There are actions choose B (chB), choose G (chG) and sample (S).
    Action values - Q(X,Ng,N) = value of taking action x when there have been N
    draws so far, from which Ng are of type g.

Basic procedure: 
1. Estimate 'ideal observer' probabilities re. Jar-as-Cause prob.s
2. Take into account possible costs and derive Q values for chosing
   actions 
3. Fit with one set of results, and apply to another ...

Version 05rbw uses data from Wellcome 07 Richard Bentall et al study,
and concentrates on exploration of 
noise-param. vs. cost-of-slowness
slice of the parametre space.

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>   // input & output streams
#include <iomanip>   // i-o  parametrized manipulators
#include <sstream>     // in-memory string manipulations for i/o, new standard
#include <string.h>
#include <time.h>
#include <getopt.h> // For getopt_long, to read command line options
#include <fenv.h> 

#include "jtc.h"       // Beads etc. setup classes.
#include "jtc_exp.h"   // Experimental procedures/setups based on jtc.h
#include "jtc_em.h"    // EM and such like operating on stuff in jtc_exp.h;
                       // contains random distro, cummul distro & multiple integrn. stuff.
#include "maths.h"   

#include "blitz_array_aux.h"


// ******************** Global Variables & Objects *****************************

// model Objects & Variables ---------------------------------------------------
Beads_Exp_II GB  ;         // General Beads with iteration class
RBWelcomeDat RBWelc ;      // provides sequence & exp. data from Wellcome study 


// Various model Parameters  ----------------------------------------------------
int groupN =1;  // would be 3 for RB expts, 4 for Fear & Hayley expt. simulation
int deprLab=0;                                          // For RPB simulations.
int cntrLab=1;    int paranLab =2;                      // for all
int ocdLab=3;     int mixedLab =4;                      // for F&H etc.
float nIter = 10;                         // Default no. of EM iterations
float mcalls = 60000;                     // default calls to Vegas for integr
float expType = 1;    // 1 for beads, 2 for words
bool approximateFunctions = true;
bool outputStandardFiles = false;

// ************************ In-file Functions **********************************
// versions using interpolation:
static double baysPDat_if_macropPDF2( double *q, size_t dim, void *param); 

// to  be used to calculate means and variances of kT and CS:
static double Tm_baysPDat_if_macropPDF2( double *q, size_t dim, void *param);
static double QlnT_baysPDat_if_macropPDF2( double *q, size_t dim, void *param);
static double CSm_baysPDat_if_macropPDF2( double *q, size_t dim, void *param);
static double QlnCS_baysPDat_if_macropPDF2( double *q, size_t dim, void *param);
static Array<int, 1> getDraws(const Array<float, 2> &, int);
static double getWrkPProd(const Array<int, 1> &, Beads_Exp_II &);
static double getWrkP2Prod(const Array<int, 1> &, Beads_Exp_II &);
static int loadDataFromFile(const string& dataFile);


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~           MAIN         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char **argv)
{
  #ifndef __APPLE__
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  #endif
  float rCSm, rCSv, rTm, rTv;             // 'reference' storage values.
  // could have 'expType' for Beads, Words or combined, but not no ! 

  // Parametrest to stay const (not all needed ...):
  GB.CGB_N= GB.kT_N = GB.CS_N= GB.rCW_N= GB.q_N= GB.PoG_N= GB.CDPar1_N= 1.0;
  GB.CGB =GB.CBG = -100.0 ;
  int ptN =33;    // 69; // 36;

  string dataFile("");
  Array<double, 2> sumRes(Range(0,roundfl(nIter)),Range(1,8)); sumRes =0; 
  GB.ResetAndSubjN(ptN);

  // Initialize  :
  GB.currCSm = rCSm = -1.5 ;      GB.currCSv = rCSv = 17.0  ;  
  GB.currTm = rTm = 4.75 ;        GB.currTv = rTv = 11.0  ;
  GB.kT_N = 80;                    GB.CS_N  = 120. ;

 int flag;
 static struct option long_options[] = {
     // No flag parameter since there's a short option.
     {"iterations", required_argument, NULL, 'i'},
     {"data_file", required_argument, NULL, 'f'},
     {"approximate", no_argument, &flag, 'a'},
     {"noapproximate", no_argument, &flag, 'A'},
     {"ref_csm", required_argument, &flag, 'c'},
     {"ref_csv", required_argument, &flag, 'C'},
     {"ref_tm", required_argument, &flag, 't'},
     {"ref_tv", required_argument, &flag, 'T'},
     {"output_standard_files", no_argument, &flag, 's'},
     {NULL, 0, NULL, 0},
 };
 int option_index = 0;
 int c;
 while ((c = getopt_long(argc, argv, "i:f:", long_options, &option_index)) != -1) {
     switch (c) {
         case 0:
            // Long option
            switch (flag) {
                case 'a':
                    approximateFunctions = true;
                    break;
                case 'A':
                    approximateFunctions = false;
                    break;
                case 'c':
                    std::istringstream(optarg) >> rCSm;
                    break;
                case 'C':
                    std::istringstream(optarg) >> rCSv;
                    break;
                case 't':
                    std::istringstream(optarg) >> rTm;
                    break;
                case 'T':
                    std::istringstream(optarg) >> rTv;
                    break;
                case 's':
                    outputStandardFiles = true;
                    break;
                default:
                    cerr << "Unexpected long option: "
                            << long_options[option_index].name << '\n';
                    exit(1);
            }
            break;
        case 'i':
           std::istringstream(optarg) >> nIter;
           break;
        case 'f':
           dataFile.assign(optarg);
           break;
        default:
           cerr << "Invalid option: " << c << '\n';
           exit(1);
     }
 }
 cout << "iterations: " << nIter << '\n';
 cout << "dataFile: " << dataFile << '\n';
 cout << "approximate: " << approximateFunctions << '\n';
 cout << "ref_csm: " << rCSm << '\n';
 cout << "ref_csv: " << rCSv << '\n';
 cout << "ref_tm: " << rTm << '\n';
 cout << "ref_tv: " << rTv << '\n';

 // convert menu item re. calls of MC routine etc. to integers
 long n = 1000*(roundfl(mcalls/1000.0));

 // Find / load appropriate data set(s) :
 ptN = loadDataFromFile(dataFile);

 // set various bits influenced  by menu:				   
 GB.ResetPDumD_L_etc2(ptN);
 sumRes.resize(Range(0,roundfl(nIter)),Range(1,8));   sumRes =0;    GB.ResetPDumD_L_etc2(ptN);
 // set starting guess :
 GB.currCSm = rCSm ;       GB.currCSv = rCSv  ;   GB.currTm = rTm ;     GB.currTv = rTv  ; 
 GB.SetApproxInf(1);  
 
 if (outputStandardFiles) {
   write2DArray_csv("data_used.csv",GB.DumDat);
 }
 cout << GB.DumDat ;
 
 cout <<"\nPriors: CS_mean: "<<GB.currCSm<<" ; CS_var: "<< GB.currCSv; 
 cout <<"; kT_mean : "<<GB.currTm<<" ; kT_var: "<< GB.currTv<<endl<<endl;
 sumRes(0,1) = 0; sumRes(0,2) = GB.currCSm; sumRes(0,3) = GB.currCSv; 
 sumRes(0,4) = GB.currTm;  sumRes(0,5) = GB.currTv; 
 sumRes(0,7) = sqrt(GB.currCSv);  sumRes(0,8) = sqrt(GB.currTv); 
 
 // ******************** Main loop *****************
 // , setting up GB.curr_pt too, which will tell the in-file fn. where we are
 
 for (int j=1; j<=roundfl(nIter); j++) { 
   
   double tempCSm, tempCSv, tempTm, tempTv;  //temporary values.
   tempCSm= tempCSv= tempTm= tempTv= 0.0;
   cout<<endl<<"Cycle "<<j<<" : "<<"w. bounds for kT and CS: "<<GB.xl[0] ;
   cout<<"  "<<GB.xu[0]<<"  "<<GB.xl[1]<<"  "<<GB.xu[1]<<" ..."<<endl;
   sumRes(j,1)=j; 
   
   // Now that we have integration bounaries, prepare the table of probabilities:
   // hard hat area ...
   // Parametrest to iterate ( kT_N already set by menu ): 
   GB.kT=GB.xl[0];  GB.kT_e=GB.xu[0];  GB.CS=GB.xl[1] ;  GB.CS_e=GB.xu[1];   


   cout<<"\n First set up the pmf on the  grid ... ";
   GB.RefreshIter(); 
   GB.PreparePVecs1();       
   while ( GB.DoIter() ) {
     GB.SoftM_PG_V_QG_QB_QS_AS_MS(); 
     GB.StoreAllPVecs(&RBWelc,0); // THIS HAS TO BE CALLED BEFORE e.g. Write2Prod3P !
     GB.StepIter() ;     // take step in parametre space
   } // end loop while ( GB.DoIter() ) ...
   GB.prepIntBasics(1); 


   cout<<"... then find Sum ln PD_L : "<<GB.update_PDumD_L( baysPDat_if_macropPDF2,n)<<endl ;
   sumRes(j-1,6)=GB.sum_ln_PdumD_L(); 

   // Parameter estimation. Do not update the GB.currXX values, as are used by 
   // the integration functions (e.g. the gammas in the integrations)
   cout<<"... now the macroparameters : \n";
   cout<<"New Tm estimate : "<<(sumRes(j,4)=tempTm=GB.estCurrDumTm( Tm_baysPDat_if_macropPDF2,n,1))<<endl ;

   double IntQlnT = GB.estCurrDum_QlnT(QlnT_baysPDat_if_macropPDF2,n,1) ;
   tempTv  = GB.xVarEst( tempTm, IntQlnT, 10);
   cout<<"New Tv estimate : "<<(sumRes(j,5) =tempTv)<<endl ;
   sumRes(j,8) =sqrt(tempTv);

   cout<<"New CSm estimate : "<<(sumRes(j,2)=tempCSm=GB.estCurrDumCSm(CSm_baysPDat_if_macropPDF2,n))<<endl ;

   double IntQlnCS = GB.estCurrDum_QlnCS(QlnCS_baysPDat_if_macropPDF2,n,0) ;
   // Here we need to use fabs(tempCSm) explicitly, as derivation of xVarEst is based on 
   // gamma distribution with positive domain :
   tempCSv  = GB.xVarEst( fabs(tempCSm), IntQlnCS, 10);
   cout<<"New CSv estimate : "<<(sumRes(j,3) =tempCSv)<<endl ;
   sumRes(j,7) =sqrt(tempCSv);

   // cout<<" intr. new CSv est: "<<(sumRes(j,3)=tempCSv=GB.estCurrDumCSv(CSv_baysPDat_if_macropPDF2,n,0))<<endl ;

   GB.currCSm=tempCSm;  GB.currCSv=tempCSv ; GB.currTm =tempTm;   GB.currTv=tempTv  ;
   
   // Determine sum_ln_PdumD_L for final estimate of parameters:
   if (j==roundfl(nIter)) {
     cout<<"Final Sum ln PD_L : "<<GB.update_PDumD_L( baysPDat_if_macropPDF2,n)<<endl ;
     sumRes(j,6)=GB.sum_ln_PdumD_L();
   }

   if (outputStandardFiles) {
     write2DArray_csv("summ_res.csv",sumRes); // overwrite with new results.
   }
 
 } // end loop over j (nIter)

 cout <<"Final sumRes: "<< sumRes ;
}  // end of main()

//------------------------ In-file Functions -------------------------------------

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// 2a. fn to integrate w. MC for denom. of Bayesian model recog. distr.
//  using Interpolated results.
static double baysPDat_if_macropPDF2( double *q, size_t dim, void *param)
{
  // use of *param not needed as fn. has access to global objects 

  double  aux1, aux2, aux3; 
  Array<int, 1> draws = getDraws(GB.DumDat, roundfl(GB.curr_pt));

  GB.kT = q[0];  GB.CS = q[1]; // NEEDED by wrkP2 ...

  if (approximateFunctions) {
      aux1 =  getWrkP2Prod(draws, GB);
  } else {
      GB.SoftM_PG_V_QG_QB_QS_AS_MS();
      GB.StoreAuxPV(&RBWelc);
      aux1 = getWrkPProd(draws, GB);
  }
  aux2 =  gammaPDF(q[1], GB.currCSm, GB.currCSv) ;  // P[micropar_1|macropar_1]  
  aux3 =  gammaPDF(q[0], GB.currTm, GB.currTv) ; 

  return aux1*aux2*aux3;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// 3a. Similar, to integrate w. MC for mean of kT using interpolation
// NOTE USE OF GLOBAL OBJECT MEMBERS AS PARAMS, *param redundant in this approach !
// RELIES ON currTm etc. TO BE UP TO DATE !
static double Tm_baysPDat_if_macropPDF2( double *q, size_t dim, void *param)
{
  double  aux1, aux2;
  Array<int, 1> draws = getDraws(GB.DumDat, roundfl(GB.curr_pt));

  GB.CS = q[1];   GB.kT = q[0];  // NEEDED by wrkP2 ...

  if (approximateFunctions) {
      aux1 =  getWrkP2Prod(draws, GB);
  } else {
      GB.SoftM_PG_V_QG_QB_QS_AS_MS();
      GB.StoreAuxPV(&RBWelc);
      aux1 = getWrkPProd(draws, GB);
  }
  aux2 = gammaPDF(q[1],GB.currCSm, GB.currCSv)*gammaPDF(q[0],GB.currTm,GB.currTv);

  return q[0]*aux1*aux2;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// 4b. Similar, to integrate w. MC for mean of CS, w. interp.
static double CSm_baysPDat_if_macropPDF2( double *q, size_t dim, void *param)
{
  double  aux1, aux2;
  Array<int, 1> draws = getDraws(GB.DumDat, roundfl(GB.curr_pt));
  GB.CS = q[1];   GB.kT = q[0];  

  if (approximateFunctions) {
      aux1 =  getWrkP2Prod(draws, GB);
  } else {
      GB.SoftM_PG_V_QG_QB_QS_AS_MS();
      GB.StoreAuxPV(&RBWelc);
      aux1 = getWrkPProd(draws, GB);
  }
  aux2 =  gammaPDF(q[1],GB.currCSm, GB.currCSv)*gammaPDF(q[0],GB.currTm,GB.currTv) ; 
  return q[1]*aux1*aux2;
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// 5. Similar, for MC int. for var of kT; RELIES ON currTm etc. TO BE UP TO DATE !
static double QlnT_baysPDat_if_macropPDF2( double *q, size_t dim, void *param)
{
  Array<int, 1> draws = getDraws(GB.DumDat, roundfl(GB.curr_pt));
  GB.kT = q[0]; GB.CS=q[1]; 

  double aux1;
  if (approximateFunctions) {
      aux1 =  getWrkP2Prod(draws, GB);
  } else {
      GB.SoftM_PG_V_QG_QB_QS_AS_MS();
      GB.StoreAuxPV(&RBWelc);
      aux1 = getWrkPProd(draws, GB);
  }
  return log(q[0]) * aux1
        * gammaPDF(q[1],GB.currCSm, GB.currCSv)
        * gammaPDF(q[0],GB.currTm,GB.currTv);
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// 6. Similar, for MC int. for var of CS; RELIES ON currCSm etc. TO BE UP TO DATE !
static double QlnCS_baysPDat_if_macropPDF2( double *q, size_t dim, void *param)
{
  Array<int, 1> draws = getDraws(GB.DumDat, roundfl(GB.curr_pt));
  GB.kT = q[0]; GB.CS=q[1]; 

  // in the following, note that CS is usually -ve, as a cost, so fabs'd !
  double aux = log(fabs(q[1]));
  double aux1;
  if (approximateFunctions) {
      aux1 =  getWrkP2Prod(draws, GB);
  } else {
      GB.SoftM_PG_V_QG_QB_QS_AS_MS();
      GB.StoreAuxPV(&RBWelc);
      aux1 = getWrkPProd(draws, GB);
  }
  return aux * aux1
        * gammaPDF(q[1],GB.currCSm, GB.currCSv)
        * gammaPDF(q[0],GB.currTm,GB.currTv);
}

static Array<int, 1> getDraws(const Array<float, 2> &DumDat, int row) {
  // Casting float values into int requires creating the array explicitly
  Array<int, 1> result(GB.seqN());
  int firstIndex = 3;
  result = DumDat(row, Range(firstIndex, firstIndex + result.size() - 1));
  return result;
}

static double getWrkPProd(const Array<int, 1> &draws, Beads_Exp_II &GB) {
  double wrkPProd = 1.0;
  for (int i = 0; i < draws.size(); ++i) {
    wrkPProd *= GB.wrkP(i + 1, draws(i));
  }
  return wrkPProd;
}

static double getWrkP2Prod(const Array<int, 1> &draws, Beads_Exp_II &GB) {
  double wrkP2Prod = 1.0;
  for (int i = 0; i < draws.size(); ++i) {
    wrkP2Prod *= GB.wrkP2(i + 1, draws(i));
  }
  return wrkP2Prod;
}

static int loadDataFromFile(const string& dataFile) {
 ifstream inpStr(dataFile.c_str()); 
 if (inpStr.bad()) { cerr<<"Arr. read error... /n"; exit(1); }
 int sN;
 inpStr >> sN;
 string *seqs = new string[sN];
 for (int i = 0; i < sN;) {
     std::getline(inpStr, seqs[i]);
     if (seqs[i].length()) {
         cerr << "seq[i]: " << seqs[i] << "\n";
         i++;
     }
 }
 RBWelc.setSequences(sN, seqs);
 // TODO(nimrod): delete individual strings
 delete[] seqs;

 int rows = 0;
 inpStr >> rows;
 if (rows <= 0) {
     cerr << "Couldn't read number of rows from " << dataFile << "\n";
     exit(1);
 }
 GB.DumDat.resize(Range(1,rows),Range(1,6)); GB.DumDat=0;
 inpStr >> GB.DumDat;  // has to be 2D !
 inpStr.close();
 return GB.DumDat.rows();
}
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

// ------------------------- end of file -----------------------------------

