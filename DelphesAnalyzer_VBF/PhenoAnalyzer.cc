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
  int nDir = 14;
  TDirectory *theDirectory[nDir];
  theDirectory[0]  = HistoOutputFile->mkdir("After_tau_pt");
  theDirectory[1]  = HistoOutputFile->mkdir("After_b_jet_veto");
  theDirectory[2]  = HistoOutputFile->mkdir("After_MET");
  theDirectory[3]  = HistoOutputFile->mkdir("After_Ntau_1");
  theDirectory[4]  = HistoOutputFile->mkdir("After_Ntau_2");
  theDirectory[5]  = HistoOutputFile->mkdir("After_Ntau_3");
  theDirectory[6]  = HistoOutputFile->mkdir("After_OnlyJets_pt");
  theDirectory[7]  = HistoOutputFile->mkdir("After_OnlyJets_eta");
  theDirectory[8]  = HistoOutputFile->mkdir("After_OnlyJets_backToback");
  theDirectory[9]  = HistoOutputFile->mkdir("After_OnlyJets_deltaEta");
  theDirectory[10] = HistoOutputFile->mkdir("After_OnlyJets_mjj");
  theDirectory[11] = HistoOutputFile->mkdir("After_Ntau_1_pass_VBF");
  theDirectory[12] = HistoOutputFile->mkdir("After_Ntau_2_pass_VBF");
  theDirectory[13] = HistoOutputFile->mkdir("After_Ntau_3_pass_VBF");
  PhenoAnalysis BSM_analysis(chain, HistoOutputFile, theDirectory, nDir);

}

using namespace std;
PhenoAnalysis::PhenoAnalysis(TChain& chain, TFile* theFile, TDirectory *cdDir[], int nDir)
{

   ifstream inFile;
   inFile.open ("config.in", ios::in); 
 
   if (!inFile)
     {
       cerr << "ERROR: Can't open input file: " << endl;
       exit (1);
     }
 
   string inputType = ""; 

   // This set of lines are used to open and read the "config.in" file. 
   /////////////////////////////////////////////////////////////////////// 
   TEnv *params = new TEnv ("config_file");
   params->ReadFile ("config.in", kEnvChange);
   
   float lead_jet_pt      = params->GetValue ("lead_jet_pt", 100.);
   float slead_jet_pt     = params->GetValue ("slead_jet_pt",100.);
   float jet_eta_max      = params->GetValue ("jet_eta_max", 5.0);
   float jj_delta_eta_min = params->GetValue ("jj_delta_eta_min", 3.5);  
   float jj_mass_min      = params->GetValue ("jj_mass_min", 100.);  
   float tau1_pt_min      = params->GetValue ("tau1_pt_min", 10.);
   float tau2_pt_min      = params->GetValue ("tau2_pt_min", 10.);
   float tau3_pt_min      = params->GetValue ("tau3_pt_min", 10.);
   float tau_eta_max      = params->GetValue ("tau_eta_max", 2.5);
    
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
       TLorentzVector Jet_subleading_vec(0., 0., 0., 0.);
       TLorentzVector Tau1_vec (0., 0., 0., 0.);
       TLorentzVector Tau2_vec (0., 0., 0., 0.);
       TLorentzVector Tau3_vec (0., 0., 0., 0.);
       vector<int> tau_index;
       bool is_b_jet = false;
       treeReader->ReadEntry(entry);
       METpointer = (MissingET*) branchMissingET->At(0);
       Double_t MET = METpointer->MET;
   
       if(branchJet->GetEntries() > 0)
	 {
           // For Jets
           double jet_highest_pt = 0.;
           double jet_Secondhighest_pt = 0.;
           int lead_ljet_index = 10000;
           int lead_sjet_index = 10000;

           if (branchJet->GetEntriesFast() < 2) continue;

           for (int j = 0; j < branchJet->GetEntriesFast(); j++)
             {
               Jet *jet = (Jet*) branchJet->At(j);
               if ((jet->PT > 20.0) && (jet->BTag == 1)){is_b_jet = true;}
               if (jet->TauTag == 1){
                 if ((abs(jet->Eta) < tau_eta_max)){
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
               if ((j == lead_ljet_index) || (jet->TauTag == 1)) continue;
               if (jet->PT > jet_Secondhighest_pt){
                  jet_Secondhighest_pt = jet->PT; 
                  lead_sjet_index = j;
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
               if (lead_sjet_index == j){
                 Jet_subleading_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
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
                  if (tau_index.size() == 3){
                    if(i == 0) Tau1_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
                    if(i == 1) Tau2_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
                    if(i == 2) Tau3_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
                  }
                                    
                }
            }
            // Apply central event selection criteria
            bool pass_lead_jet_cuts = false;
            bool pass_sublead_jet_cuts = false;
            cout << "LEADING JET "<<Jet_leading_vec.Pt()<<"  SUB LEADING JET "<<Jet_subleading_vec.Pt()<<endl;

            if ((Tau1_vec.Pt() > tau1_pt_min) || (Tau2_vec.Pt() > tau2_pt_min) || (Tau3_vec.Pt() > tau3_pt_min)){
               pass_cuts[0] = 1; 
            }
            if ((pass_cuts[0] == 1) && (!is_b_jet)){pass_cuts[1] = 1;}
            if ((pass_cuts[1] == 1) && (MET > 10.)){pass_cuts[2] = 1;}
            if ((pass_cuts[2] == 1) && (tau_index.size() == 1)){pass_cuts[3] = 1;}
            if ((pass_cuts[2] == 1) && (tau_index.size() == 2)){pass_cuts[4] = 1;}
            if ((pass_cuts[2] == 1) && (tau_index.size() == 3)){pass_cuts[5] = 1;}

            // Apply VBF selections
            if ((Jet_leading_vec.Pt() > lead_jet_pt) && (Jet_subleading_vec.Pt() > slead_jet_pt)){pass_cuts[6] = 1;} 
            if ((pass_cuts[6] == 1) && (abs(Jet_leading_vec.Eta()) < jet_eta_max) && (abs(Jet_subleading_vec.Eta()) < jet_eta_max)){pass_cuts[7] = 1;}
            if (pass_cuts[7] == 1){
               double eta_jj_product = Jet_leading_vec.Eta()*Jet_subleading_vec.Eta();
               if (eta_jj_product < 0.){pass_cuts[8] = 1;}
            }
            if (pass_cuts[8] == 1){
              double Delta_eta_jj = abs(Jet_leading_vec.Eta() - Jet_subleading_vec.Eta());
              if (Delta_eta_jj > jj_delta_eta_min){pass_cuts[9] = 1;}
            }
            if (pass_cuts[9] == 1){
              TLorentzVector jj_LV_sum = Jet_leading_vec + Jet_subleading_vec;
              double mjj = jj_LV_sum.M();
              if (mjj > jj_mass_min){pass_cuts[10] = 1;}
            }
            if ((pass_cuts[3] == 1) && (pass_cuts[10] == 1)){pass_cuts[11] = 1;}
            if ((pass_cuts[4] == 1) && (pass_cuts[10] == 1)){pass_cuts[12] = 1;}
            if ((pass_cuts[5] == 1) && (pass_cuts[10] == 1)){pass_cuts[13] = 1;}

	 }
       for (int i = 0; i < nDir; i++){
         if ( pass_cuts[i] == 1){
            _hmap_lead_jet_pT[i]->Fill(Jet_leading_vec.Pt());
            _hmap_lead_jet_eta[i]->Fill(Jet_leading_vec.Eta());
            _hmap_lead_jet_phi[i]->Fill(Jet_leading_vec.Phi());
            _hmap_sublead_jet_pT[i]->Fill(Jet_subleading_vec.Pt());
            _hmap_sublead_jet_eta[i]->Fill(Jet_subleading_vec.Eta());
            _hmap_sublead_jet_phi[i]->Fill(Jet_subleading_vec.Phi());
            double delta_eta_jj = abs(Jet_leading_vec.Eta() - Jet_subleading_vec.Eta());
            _hmap_delta_eta_jj[i]->Fill(delta_eta_jj);
            double mjj = (Jet_leading_vec + Jet_subleading_vec).M();
            _hmap_jet_mjj[i]->Fill(mjj);
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
            if(Tau3_vec.Pt() > 10.0){
              _hmap_tau3_pT[i]->Fill(Tau3_vec.Pt());
              _hmap_tau3_eta[i]->Fill(Tau3_vec.Eta());
              _hmap_tau3_phi[i]->Fill(Tau3_vec.Phi());
            }
             
         }
       }  
       Electron *elec1, *elec2;
     }

     theFile->cd();
     for (int d = 0; d < nDir; d++)
       {
         cdDir[d]->cd();
         _hmap_lead_jet_pT[d]->Write();
         _hmap_lead_jet_eta[d]->Write();
         _hmap_lead_jet_phi[d]->Write();
         _hmap_sublead_jet_pT[d]->Write();
         _hmap_sublead_jet_eta[d]->Write();
         _hmap_sublead_jet_phi[d]->Write();
         _hmap_delta_eta_jj[d]->Write();
         _hmap_jet_mjj[d]->Write();
         _hmap_tau1_pT[d]->Write();
         _hmap_tau1_eta[d]->Write();
         _hmap_tau1_phi[d]->Write();
         _hmap_tau2_pT[d]->Write();
         _hmap_tau2_eta[d]->Write();
         _hmap_tau2_phi[d]->Write();
         _hmap_tau3_pT[d]->Write();
         _hmap_tau3_eta[d]->Write();
         _hmap_tau3_phi[d]->Write();
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
       _hmap_sublead_jet_pT[i]  = new TH1F("jet_sublead_pT",   "j p_{T}", 1000, 0., 1000.);
       _hmap_sublead_jet_eta[i] = new TH1F("jet_sublead_eta",  "j #eta", 100, -5.0, 5.0);
       _hmap_sublead_jet_phi[i] = new TH1F("jet_sub_lead_phi", "j #phi", 72, -3.6, 3.6);
       _hmap_delta_eta_jj[i]    = new TH1F("delta_eta_jj",     "#Delta #eta_{jj}", 100, 0., 10.0);
       _hmap_jet_mjj[i]         = new TH1F("dijet_mjj",        "m_{jj}", 500, 0., 5000.0); 
       _hmap_tau1_pT[i]         = new TH1F("tau1_pT",          "p_{T}(#tau_{1})", 300, 0., 300.);
       _hmap_tau1_eta[i]        = new TH1F("tau1_eta",         "#eta(#tau_{1})", 50, -2.5, 2.5);
       _hmap_tau1_phi[i]        = new TH1F("tau1_phi",         "#phi(#tau_{1})", 72, -3.6, 3.6);
       _hmap_tau2_pT[i]         = new TH1F("tau2_pT",          "p_{T}(#tau_{2})", 300, 0., 300.);
       _hmap_tau2_eta[i]        = new TH1F("tau2_eta",         "#eta(#tau_{2})", 50, -2.5, 2.5);
       _hmap_tau2_phi[i]        = new TH1F("tau2_phi",         "#phi(#tau_{2})", 72, -3.6, 3.6);
       _hmap_tau3_pT[i]         = new TH1F("tau3_pT",          "p_{T}(#tau_{3})", 300, 0., 300.);
       _hmap_tau3_eta[i]        = new TH1F("tau3_eta",         "#eta(#tau_{3})", 50, -2.5, 2.5);
       _hmap_tau3_phi[i]        = new TH1F("tau3_phi",         "#phi(#tau_{3})", 72, -3.6, 3.6);
     }
}
