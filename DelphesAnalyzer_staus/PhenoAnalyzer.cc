////////////////////////////////////////////////////////////////
//
//                                                            //
// Author: Andrés Flórez, Universidad de los Andes, Colombia  //
//                                                            //
////////////////////////////////////////////////////////////////

#include <iostream>
//Hola  Comment the following lines if you do not need ROOT or HepMC
#include "ROOTFunctions.h"
#include "PhenoAnalyzer.h"

int main(int argc, char *argv[]) {

  TApplication app("App",&argc, argv);
  TChain chain("Delphes");
  chain.Add("output_delphes.root");
  TFile * HistoOutputFile = new TFile("HistoOutputFile.root", "RECREATE");
  int nDir = 5;
  TDirectory *theDirectory[nDir];
  theDirectory[0] = HistoOutputFile->mkdir("After_cut_1");
  theDirectory[1] = HistoOutputFile->mkdir("After_cut_2");
  theDirectory[2] = HistoOutputFile->mkdir("After_cut_3");
  theDirectory[3] = HistoOutputFile->mkdir("After_cut_4");
  theDirectory[4] = HistoOutputFile->mkdir("After_cut_5");
  PhenoAnalysis BSM_analysis(chain, HistoOutputFile, theDirectory, nDir);

}

using namespace std;
PhenoAnalysis::PhenoAnalysis(TChain& chain, TFile* theFile, TDirectory *cdDir[], int nDir)
{
   crateHistoMasps(nDir);

   ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
   Long64_t numberOfEntries = treeReader->GetEntries();
   
   TClonesArray *branchJet = treeReader->UseBranch("Jet");
   TClonesArray *branchElectron = treeReader->UseBranch("Electron");
   TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");  
 
   MissingET *METpointer; 

   numberOfEntries = 10000;
   for(Int_t entry = 0; entry < numberOfEntries; ++entry)
     {
       int pass_cuts[nDir] = {0};
       TLorentzVector Jet_leading_vec(0., 0., 0., 0.);
       treeReader->ReadEntry(entry);
       METpointer = (MissingET*) branchMissingET->At(0);
       Double_t MET = METpointer->MET;
   
       if(branchJet->GetEntries() > 0)
	 {
           double jet_highest_pt = 0.;
           for (int j = 0; j < branchJet->GetEntriesFast(); j++)
             {
               Jet *jet = (Jet*) branchJet->At(j);
               if (jet->PT > jet_highest_pt){jet_highest_pt = jet->PT;}
             }
	   for (int j = 0; j < branchJet->GetEntriesFast(); j++)
	     {
	       Jet *jet = (Jet*) branchJet->At(j);
               if (jet->PT == jet_highest_pt){
                 // Jet energy is not correct! YOU NEED TO CALCULATE THE P OF THE JET 
                 // USING THE JET ETA and PT
                 cout <<"JET PT "<<jet->PT<<endl;
		 double p = sqrt(pow(exp(jet->Eta)*jet->PT/2, 2)-pow(jet->PT,2));
                 double jet_energy = sqrt(pow(p, 2) + pow(jet->Mass, 2));
                 Jet_leading_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
                 if(jet->PT > 30.0 ){pass_cuts[0] = 1;}
                 if((jet->PT > 30.0) && (abs(jet->Eta) < 5.0) ){pass_cuts[1] = 1;}
                 if(jet->PT > 100.0 ){pass_cuts[2] = 1;}
                 if((jet->PT > 100.0) && (abs(jet->Eta) < 5.0) ){pass_cuts[3] = 1;}
               }
	     }
            if((Jet_leading_vec.Pt() > 100.0) && (abs(Jet_leading_vec.Eta()) < 5.0)  && (MET > 30.)){pass_cuts[4] = 1;} 
	 }
       for (int i = 0; i < nDir; i++){
         if ( pass_cuts[i] == 1){
            _hmap_jet_pT[i]->Fill(Jet_leading_vec.Pt());
            _hmap_jet_eta[i]->Fill(Jet_leading_vec.Eta());
         }
       }  
       Electron *elec1, *elec2;
     }

     theFile->cd();
     for (int d = 0; d < nDir; d++)
       {
         cdDir[d]->cd();
         _hmap_jet_pT[d]->Write();
         _hmap_jet_eta[d]->Write();
       }
     theFile->Close();
   
}

PhenoAnalysis::~PhenoAnalysis()
{
  // do anything here that needs to be done at desctruction time
}


void PhenoAnalysis::crateHistoMasps (int directories)
{
   for (int i = 0; i < directories; i++)
     {
       _hmap_jet_pT[i]     = new TH1F("jet_pT", "j p_{T}", 300, 0., 300.);
       _hmap_jet_eta[i]    = new TH1F("jet_eta", "j #eta", 50, -2.5, 2.5);
     }
}
