#ifndef _PHENOANALYZER_H_
#define _PHENOANALYZER_H_

#include "TF2.h"
#include "TH1F.h"
#include "TChain.h"

#include <iostream>
#include <utility>
#include <vector>

#include "TROOT.h"
#include "TSystem.h"

#include "TString.h"

#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TApplication.h"
#include "DelphesFunctions.h"
#include "TDirectory.h"
#include "TFile.h"

using namespace std;

class PhenoAnalysis {
public :
   PhenoAnalysis(TChain&, TFile*, TDirectory* dir[], int nDir);
   ~PhenoAnalysis();
   void crateHistoMasps (int);
   std::map<unsigned int, TH1*> _hmap_jet_pT;
   std::map<unsigned int, TH1*> _hmap_jet_eta;

private :

};

#endif
