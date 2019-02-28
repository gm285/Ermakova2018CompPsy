/*
 * Generates fake data for the given parameters.
 */

#include <fenv.h>
#include <getopt.h>  // For getopt_long, to read command line options
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fstream>  // input & output streams
#include <iomanip>  // i-o  parametrized manipulators
#include <sstream>  // in-memory string manipulations for i/o, new standard

#include "jtc.h"      // Beads etc. setup classes.
#include "jtc_em.h"   // EM and such like operating on stuff in jtc_exp.h;
#include "jtc_exp.h"  // Experimental procedures/setups based on jtc.h
// contains random distro, cummul distro & multiple integrn. stuff.
#include "maths.h"

#include "blitz_array_aux.h"

static int loadDataFromFile(const string& sequencesFile, int* ptN);
static void writeSequences(ofstream&, const RBWelcomeDat&);
static void writeParticipantData(ofstream&, const Array<float, 2>&);
// ******************** Global Variables & Objects *****************************

// model Objects & Variables ---------------------------------------------------
Beads_Exp_II GB;      // General Beads with iteration class
RBWelcomeDat RBWelc;  // provides sequence & exp. data from Wellcome study

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~           MAIN         ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main(int argc, char** argv) {
#ifndef __APPLE__
  feenableexcept(FE_INVALID | FE_OVERFLOW);
#endif
  // Parametrest to stay const (not all needed ...):
  GB.CGB_N = GB.kT_N = GB.CS_N = GB.rCW_N = GB.q_N = GB.PoG_N = GB.CDPar1_N =
      1.0;
  GB.CGB = GB.CBG = -100.0;

  // Initialize  :
  GB.currCSm = -1.5;
  GB.currCSv = 17.0;
  GB.currTm = 4.75;
  GB.currTv = 11.0;
  GB.kT_N = 80;
  GB.CS_N = 120.;

  string sequencesFile("");
  string outputFile("");
  // Use constant random seed for reproducibility
  int seed = 1337;
  double rTm = 3.0;
  double rTv = 5.0;
  double rCSm = -2.0;
  double rCSv = 17;

  int flag;
  static struct option long_options[] = {
      // No flag parameter since there's a short option.
      {"csm", required_argument, &flag, 'c'},
      {"csv", required_argument, &flag, 'C'},
      {"tm", required_argument, &flag, 't'},
      {"tv", required_argument, &flag, 'T'},
      {"seed", required_argument, &flag, 's'},
      {NULL, 0, NULL, 0},
  };
  int option_index = 0;
  int c;
  while ((c = getopt_long(argc, argv, "f:o:", long_options, &option_index)) !=
         -1) {
    switch (c) {
      case 0:
        // Long option
        switch (flag) {
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
            std::istringstream(optarg) >> seed;
            break;
          default:
            cerr << "Unexpected long option: "
                 << long_options[option_index].name << '\n';
            exit(1);
        }
        break;
      case 'f':
        sequencesFile.assign(optarg);
        break;
      case 'o':
        outputFile.assign(optarg);
        break;
      default:
        cerr << "Invalid option: " << c << '\n';
        exit(1);
    }
  }
  cerr << "csm: " << rCSm << '\n';
  cerr << "csv: " << rCSv << '\n';
  cerr << "tm: " << rTm << '\n';
  cerr << "tv: " << rTv << '\n';
  cerr << "seed: " << seed << '\n';
  cerr << "sequencesFile: " << sequencesFile << '\n';
  cerr << "outputFile: " << outputFile << '\n';

  if (outputFile.empty()) {
    // TODO(nimrod): Just spit out to stdout in that case
    cerr << "No output file specified.\n";
    exit(1);
  }

  int ptN;
  // Load sequences from real data experiment.
  // We could potentially generate new ones, but they were constant in our
  // experiment.
  loadDataFromFile(sequencesFile, &ptN);

  Array<double, 2> sumRes(Range(0, roundfl(1)), Range(1, 8));
  sumRes = 0;
  GB.ResetPDumD_L_etc2(ptN);
  sumRes.resize(Range(0, roundfl(1)), Range(1, 8));
  sumRes = 0;
  GB.ResetPDumD_L_etc2(ptN);

  // Create pseudodata and make appropriate adjustments:

  double aux4, aux5, aux6, aux7;
  aux4 = rTm;
  aux5 = rTv;
  aux6 = rCSm;
  aux7 = rCSv;

  GB.SetCurrAndBounds(&aux4, &aux5, &aux6, &aux7);
  // NB fillDum... automatically resizes appropriately.
  // Simulated d2d data:
  GB.fillDumDatBays(&RBWelc, ptN, seed);
  // Simulated decisions taken:
  for (int j = 1; j <= ptN; j++) GB.DumDat(j, 1) = 1 + 0.001 * j;
  ofstream output;
  output.open(outputFile.c_str(), ios::out);
  writeSequences(output, RBWelc);
  writeParticipantData(output, GB.DumDat);
  // write2DArray_csv(output, GB.DumDat);
  output.close();

}  // end of main()
static int loadDataFromFile(const string& sequencesFile, int *ptN) {
  ifstream inpStr(sequencesFile.c_str());
  if (inpStr.bad()) {
    cerr << "Arr. read error... /n";
    exit(1);
  }
  int sN;
  inpStr >> sN;
  string* seqs = new string[sN];
  for (int i = 0; i < sN;) {
    std::getline(inpStr, seqs[i]);
    if (seqs[i].length()) {
      cerr << "seq[i]: " << seqs[i] << "\n";
      i++;
    }
  }
  RBWelc.setSequences(sN, seqs);
  inpStr >> *ptN;
  // TODO(nimrod): delete individual strings
  delete[] seqs;
  inpStr.close();
}

static void writeSequences(
    ofstream &output, const RBWelcomeDat &rbwelc) {
  int sN = rbwelc.seqN();
  output << sN << endl;

  // One-based. :(
  for (int i = 1; i < sN + 1; i++) {
    for (int j = 1; j < rbwelc.SamplMaxN() + 1; j++) {
      output << (int) rbwelc.seq[i][j];
    }
    output << endl;
  }
}
static void writeParticipantData(ofstream &output, const Array<float, 2> &data) {
  int participants = data.ubound(firstDim) - data.lbound(firstDim) + 1;
  // There are 3 extra items which aren't part of the draws to decision output
  int sequences = data.ubound(secondDim) - data.lbound(secondDim) - 3 + 1;
  // write number of participants
  output << participants << endl;
  // write dimensions of array
  output << "(1," << participants << ") x (1," << sequences + 2 << ")" << endl;
  // write array with entries
  // participant id, another id, dtd * 20
  output << "[ ";
  for (int row = data.lbound(firstDim); row <= data.ubound(firstDim); row++) {
    output << " " << row << " " << 1000 + row;
    for (int col = data.lbound(secondDim) + 2; col <= data.ubound(secondDim) - 1; col++) {
      output << " " << (int) (data(row, col));
    }
    output << endl;
  }

  output << "]" << endl;
}
