////////////////////////////////////////////////////////////////
//                                                            //
// Author: Andrés Flórez, Universidad de los Andes, Colombia  //
//                                                            //
////////////////////////////////////////////////////////////////


#include <iostream>
#include "ROOTFunctions.h"
#include "PhenoAnalyzer.h"

int main(int argc, char *argv[]) {

  //TApplication app("App",&argc, argv);
  TChain chain("Delphes");
  chain.Add(argv[1]);
  TFile * HistoOutputFile = new TFile(argv[2], "RECREATE");
  int nDir = 9;
  TDirectory *theDirectory[nDir];
  theDirectory[0]  = HistoOutputFile->mkdir("After_tau_pt");
  theDirectory[1]  = HistoOutputFile->mkdir("After_b_jet_veto");
  theDirectory[2]  = HistoOutputFile->mkdir("After_MET");
  theDirectory[3]  = HistoOutputFile->mkdir("After_Ntau_1");
  theDirectory[4]  = HistoOutputFile->mkdir("After_Ntau_2");
  theDirectory[5]  = HistoOutputFile->mkdir("After_OnlyJets_pt");
  theDirectory[6]  = HistoOutputFile->mkdir("After_OnlyJets_eta");
  theDirectory[7] = HistoOutputFile->mkdir("After_Ntau_1_pass_VBF");
  theDirectory[8] = HistoOutputFile->mkdir("After_Ntau_2_pass_VBF");
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

   //numberOfEntries = 10000;
   for(Int_t entry = 0; entry < numberOfEntries; ++entry)
     {
       int pass_cuts[nDir] = {0};
       TLorentzVector Jet_leading_vec(0., 0., 0., 0.);
       TLorentzVector Tau1_vec (0., 0., 0., 0.);
       TLorentzVector Tau2_vec (0., 0., 0., 0.);
       vector<int> tau_index;
       bool is_b_jet = false;
       treeReader->ReadEntry(entry);
       METpointer = (MissingET*) branchMissingET->At(0);
       Double_t MET = METpointer->MET;
   
       if(branchJet->GetEntries() > 0)
	 {
           // For Jets
           double jet_highest_pt = 0.;
           int lead_ljet_index = 10000;

           if (branchJet->GetEntriesFast() < 3) continue;

           for (int j = 0; j < branchJet->GetEntriesFast(); j++)
             {
               Jet *jet = (Jet*) branchJet->At(j);
               if ((jet->PT > 20.0) && (jet->BTag == 1)){is_b_jet = true;}
               if (jet->TauTag == 1){
                 if ((abs(jet->Eta) < 2.4)){
                   tau_index.push_back(j);
                 }
               }
               if (jet->TauTag == 1) continue;
               if (jet->PT > jet_highest_pt){
                  jet_highest_pt = jet->PT; 
                  lead_ljet_index = j;
               }
             }
	   for (int j = 0; j < branchJet->GetEntriesFast(); j++)
	     {
	       Jet *jet = (Jet*) branchJet->At(j);
               // Jet energy is not correct! YOU NEED TO CALCULATE THE P OF THE JET 
               // USING THE JET ETA and PT
               double theta =TMath::ATan(TMath::Exp(-jet->Eta));  //Variable de angulo
               double cos_theta=TMath::Cos(2*theta);
               double Jet_p=jet->PT/cos_theta;
               double jet_energy = sqrt(pow(Jet_p, 2) + pow(jet->Mass, 2));
               if (lead_ljet_index == j){
                 Jet_leading_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
               }
	     }
           
            if (tau_index.size() > 0){
              for (int i = 0; i < tau_index.size(); i++)
                {
                  int tau_i = tau_index.at(i);
                  Jet *jet = (Jet*) branchJet->At(tau_i); 
                  double theta =TMath::ATan(TMath::Exp(-jet->Eta));  //Variable de angulo
                  double cos_theta=TMath::Cos(2*theta);
                  double Jet_p=jet->PT/cos_theta;
                  double jet_energy = sqrt(pow(Jet_p, 2) + pow(jet->Mass, 2));
                  if (tau_index.size() == 1){
                    Tau1_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
                  }
                  if (tau_index.size() == 2){
                    if(i == 0) Tau1_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
                    if(i == 1) Tau2_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
                  }
                }
            }
            // Apply central event selection criteria
            bool pass_lead_jet_cuts = false;

            if ((Tau1_vec.Pt() > 10.0) || (Tau2_vec.Pt() > 10.0)){
               pass_cuts[0] = 1; 
            }
            if ((pass_cuts[0] == 1) && (!is_b_jet)){pass_cuts[1] = 1;}
            if ((pass_cuts[1] == 1) && (MET > 10.)){pass_cuts[2] = 1;}
            if ((pass_cuts[2] == 1) && (tau_index.size() == 1)){pass_cuts[3] = 1;}
            if ((pass_cuts[2] == 1) && (tau_index.size() == 2)){pass_cuts[4] = 1;}

            // Apply VBF selections
            if (Jet_leading_vec.Pt() > 30.0 ){pass_cuts[5] = 1;} 
            if ((pass_cuts[5] == 1) && (abs(Jet_leading_vec.Eta()) < 5.0)) {pass_cuts[6] = 1;}
            if ((pass_cuts[3] == 1) && (pass_cuts[6] == 1)){pass_cuts[7] = 1;}
            if ((pass_cuts[4] == 1) && (pass_cuts[6] == 1)){pass_cuts[8] = 1;}

	 }
       for (int i = 0; i < nDir; i++){
         if ( pass_cuts[i] == 1){
if (i==7){cout<<"PASSING 7"<<endl;}
            _hmap_lead_jet_pT[i]->Fill(Jet_leading_vec.Pt());
            _hmap_lead_jet_eta[i]->Fill(Jet_leading_vec.Eta());
            _hmap_lead_jet_phi[i]->Fill(Jet_leading_vec.Phi());
            if(Tau1_vec.Pt() > 10.0){
              _hmap_tau1_pT[i]->Fill(Tau1_vec.Pt());
              _hmap_tau1_eta[i]->Fill(Tau1_vec.Eta());
              _hmap_tau1_phi[i]->Fill(Tau1_vec.Phi());
            }
            if(Tau2_vec.Pt() > 10.0){
              _hmap_tau2_pT[i]->Fill(Tau2_vec.Pt());
              _hmap_tau2_eta[i]->Fill(Tau2_vec.Eta());
              _hmap_tau2_phi[i]->Fill(Tau2_vec.Phi());
            }
         }
       }  
     }

     theFile->cd();
     for (int d = 0; d < nDir; d++)
       {
         cdDir[d]->cd();
         _hmap_lead_jet_pT[d]->Write();
         _hmap_lead_jet_eta[d]->Write();
         _hmap_lead_jet_phi[d]->Write();
         _hmap_tau1_pT[d]->Write();
         _hmap_tau1_eta[d]->Write();
         _hmap_tau1_phi[d]->Write();
         _hmap_tau2_pT[d]->Write();
         _hmap_tau2_eta[d]->Write();
         _hmap_tau2_phi[d]->Write();
       }
     theFile->Close();
   
}

PhenoAnalysis::~PhenoAnalysis()
{
  // do anything here that needs to be done at desctruction time
}

bool PhenoAnalysis::overlapingObjects(double eta1, double eta2, double phi1, double phi2, double dR)
{

  double dEta = eta1 - eta2;
  double dPhi = phi1 = phi2;
  double DR = sqrt(pow(dEta, 2.) + pow(dPhi, 2.));
  bool pass = false;
  if (DR > dR){pass = true;}

  return pass;
}

void PhenoAnalysis::crateHistoMasps (int directories)
{
   for (int i = 0; i < directories; i++)
     {
       _hmap_lead_jet_pT[i]     = new TH1F("jet_lead_pT",      "j p_{T}", 1000, 0., 1000.);
       _hmap_lead_jet_eta[i]    = new TH1F("jet_lead_eta",     "j #eta", 100, -5.0, 5.0);
       _hmap_lead_jet_phi[i]    = new TH1F("jet_lead_phi",     "j #phi", 72, -3.6, 3.6); 
       _hmap_tau1_pT[i]         = new TH1F("tau1_pT",          "p_{T}(#tau_{1})", 300, 0., 300.);
       _hmap_tau1_eta[i]        = new TH1F("tau1_eta",         "#eta(#tau_{1})", 50, -2.5, 2.5);
       _hmap_tau1_phi[i]        = new TH1F("tau1_phi",         "#phi(#tau_{1})", 72, -3.6, 3.6);
       _hmap_tau2_pT[i]         = new TH1F("tau2_pT",          "p_{T}(#tau_{2})", 300, 0., 300.);
       _hmap_tau2_eta[i]        = new TH1F("tau2_eta",         "#eta(#tau_{2})", 50, -2.5, 2.5);
       _hmap_tau2_phi[i]        = new TH1F("tau2_phi",         "#phi(#tau_{2})", 72, -3.6, 3.6);
     }
}
