////////////////////////////////////////////////////////////////
//                                                            //
// Author: Andrés Flórez, Universidad de los Andes, Colombia  //
//                                                            //
////////////////////////////////////////////////////////////////



#include <iostream>
#include "ROOTFunctions.h"
#include "DelphesFunctions.h"
#include "PhenoAnalyzer.h"
#include <time.h>

int main(int argc, char *argv[]) {
  
  //TApplication app("App",&argc, argv);
  TChain chain("Delphes");
  chain.Add(argv[1]);
  TFile * HistoOutputFile = new TFile(argv[2], "RECREATE");
  int nDir = 23;
  TDirectory *theDirectory[nDir];
  theDirectory[0]  = HistoOutputFile->mkdir("Staus_No_cuts");
  theDirectory[1]  = HistoOutputFile->mkdir("Staus_with_one_Tau");
  theDirectory[2]  = HistoOutputFile->mkdir("Staus_with_two_Taus");
  theDirectory[3]  = HistoOutputFile->mkdir("Staus_After_Tau_pt_min_1Tau");
  theDirectory[4]  = HistoOutputFile->mkdir("Staus_After_Tau_pt_min_2Taus");
  theDirectory[5]  = HistoOutputFile->mkdir("Staus_After_Tau_Eta_min_1Tau");
  theDirectory[6]  = HistoOutputFile->mkdir("Staus_After_Tau_Eta_min_2Taus");
  theDirectory[7]  = HistoOutputFile->mkdir("Staus_After_b_jet_veto_1Tau");
  theDirectory[8]  = HistoOutputFile->mkdir("Staus_After_b_jet_veto_2Taus");
  theDirectory[9]  = HistoOutputFile->mkdir("Staus_After_jet_pt_cut_1Tau");
  theDirectory[10]  = HistoOutputFile->mkdir("Staus_After_jet_pt_cut_2Taus");
  theDirectory[11]  = HistoOutputFile->mkdir("Staus_After_jet_Eta_cut_1Taus");
  theDirectory[12]  = HistoOutputFile->mkdir("Staus_After_jet_Eta_cut_2Taus");
  theDirectory[13]  = HistoOutputFile->mkdir("Staus_Dphi_tau_jet_cut_1Tau");
  theDirectory[14]  = HistoOutputFile->mkdir("Staus_Dphi_tau_jet_cut_2Taus");
  theDirectory[15]  = HistoOutputFile->mkdir("Staus_Transverse_mass_cut_1Tau");
  theDirectory[16]  = HistoOutputFile->mkdir("Staus_Transverse_mass_cut_2Taus");
  theDirectory[17]  = HistoOutputFile->mkdir("Staus_After_MET_1Tau");
  theDirectory[18]  = HistoOutputFile->mkdir("Staus_After_MET_2Taus");
  theDirectory[19]  = HistoOutputFile->mkdir("Staus_Dphi_jet_met_1Tau");
  theDirectory[20]  = HistoOutputFile->mkdir("Staus_Dphi_jet_met_2Taus"); 
  theDirectory[21]  = HistoOutputFile->mkdir("Staus_Efective_mass_cut_1Tau");
  theDirectory[22]  = HistoOutputFile->mkdir("Staus_Efective_mass_cut_2Taus");
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

  //This set of lines are used to open and read the "config.in" file. 
  /////////////////////////////////////////////////////////////////////// 
  TEnv *params = new TEnv ("config_file");
  params->ReadFile ("config.in", kEnvChange);
  
  int BGQCD               = params->GetValue ("BGQCD", 0);
  double lead_jet_pt      = params->GetValue ("lead_jet_pt", 100.);
  double jet_eta_max      = params->GetValue ("jet_eta_max", 5.0);
  double tau1_pt_min      = params->GetValue ("tau1_pt_min", 10.);
  double tau2_pt_min      = params->GetValue ("tau2_pt_min", 10.);
  double tau_eta_max      = params->GetValue ("tau_eta_max", 2.1);
  double b_jet_pt_min     = params->GetValue ("b_jet_pt_min", 20.0);
  double DR_jet_elec_max  = params->GetValue ("DR_jet_elec_max", 0.3);
  double met_min          = params->GetValue ("met_min", 10.);
  double deltajetmetphi   = params->GetValue ("deltajetmetphi",1.5);
  double transversemassmin = params->GetValue ("transversemassmin",0.);
  double efectivemassmin =  params->GetValue ("efectivemassmin",0.);
  double Dphitaujetmax   =  params->GetValue ("Dphitaujetmax",0.); 
  crateHistoMasps(nDir);
  
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");  
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle"); 
  
  MissingET *METpointer; 
  TH1 *Numbertaus_D1     = new TH1F("Numbertaus_D1",  "Numbertaus_D1", 30, -5.5, 24.5); 		 
  TH1 *Numbertaus_D2     = new TH1F("Numbertaus_D2",  "Numbertaus_D2", 10, 0, 10);
  TH1 *Numbertaus_2      = new TH1F("Numbertaus_2",   "Numbertaus_D2", 30, -5.5, 24.5);
  TH1 *Mothers           = new TH1F("Mothers",        "Mothers", 510, -10, 500);
  
  int counter_elec_muon  = 0;

  
    for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  //for(Int_t entry = 0; entry < 100; ++entry)
    {
      int pass_cuts[nDir] = {0};
      
      TLorentzVector Jet_leading_vec(0., 0., 0., 0.);
      TLorentzVector electron (0., 0., 0., 0.);
      TLorentzVector Tau1_1cand_vec (0., 0., 0., 0.);
      TLorentzVector Tau1_2cand_vec(0., 0., 0., 0.);
      TLorentzVector Tau1_3cand_vec(0., 0., 0., 0.);
      TLorentzVector Tau1_4cand_vec(0., 0., 0., 0.);
      TLorentzVector Tau2_1cand_vec (0., 0., 0., 0.);
      TLorentzVector Tau2_2cand_vec (0., 0., 0., 0.);
      TLorentzVector Tau2_3cand_vec (0., 0., 0., 0.);
      TLorentzVector Tau2_4cand_vec (0., 0., 0., 0.);
      TLorentzVector Tau1cand_vec (0., 0., 0., 0.);
      TLorentzVector Tau2cand_vec(0., 0., 0., 0.);
      TLorentzVector Tau3cand_vec (0., 0., 0., 0.);
      TLorentzVector Tau4cand_vec (0., 0., 0., 0.);
      TLorentzVector Tau1Had_vec (0., 0., 0., 0.);
      TLorentzVector Tau2Had_vec (0., 0., 0., 0.);      
      TLorentzVector Tau3Had_vec (0., 0., 0., 0.);
      TLorentzVector Tau4Had_vec (0., 0., 0., 0.);
      TLorentzVector Tau1HadTLV (0., 0., 0., 0.);
      TLorentzVector Tau2HadTLV (0., 0., 0., 0.);    
      TLorentzVector Tau3HadTLV (0., 0., 0., 0.);
      TLorentzVector Tau4HadTLV (0., 0., 0., 0.);
      TLorentzVector TauLep1cand_vec (0., 0., 0., 0.);
      TLorentzVector TauLep2cand_vec (0., 0., 0., 0.);
      TLorentzVector TauLep3cand_vec (0., 0., 0., 0.);
      TLorentzVector TauLep4cand_vec (0., 0., 0., 0.);
      TLorentzVector NeuTau1 (0., 0., 0., 0.);
      TLorentzVector NeuTau2 (0., 0., 0., 0.);
      TLorentzVector NeuTau3 (0., 0., 0., 0.);
      TLorentzVector NeuTau4 (0., 0., 0., 0.);
      TLorentzVector NeuTau1cand_vec (0., 0., 0., 0.);
      TLorentzVector NeuTau2cand_vec (0., 0., 0., 0.);
      TLorentzVector NeuTau3cand_vec (0., 0., 0., 0.);
      TLorentzVector NeuTau4cand_vec (0., 0., 0., 0.);
      TLorentzVector NeuTau1Had_vec (0., 0., 0., 0.);
      TLorentzVector NeuTau2Had_vec (0., 0., 0., 0.);
      TLorentzVector NeuTau3Had_vec (0., 0., 0., 0.);
      TLorentzVector NeuTau4Had_vec (0., 0., 0., 0.);
      TLorentzVector NeuTau1HadTLV (0., 0., 0., 0.);
      TLorentzVector NeuTau2HadTLV (0., 0., 0., 0.);
      TLorentzVector NeuTau3HadTLV (0., 0., 0., 0.);
      TLorentzVector NeuTau4HadTLV (0., 0., 0., 0.);
      TLorentzVector NeuLep1cand_vec (0., 0., 0., 0.);
      TLorentzVector NeuLep2cand_vec (0., 0., 0., 0.);
      TLorentzVector NeuLep3cand_vec (0., 0., 0., 0.);
      TLorentzVector NeuLep4cand_vec (0., 0., 0., 0.);
      TLorentzVector Lepton1_vec (0., 0., 0., 0.);
      TLorentzVector Lepton2_vec (0., 0., 0., 0.);
      TLorentzVector Lepton3_vec (0., 0., 0., 0.);
      TLorentzVector Lepton4_vec (0., 0., 0., 0.);
      TLorentzVector Neu_Tau1_1cand_vec(0., 0., 0., 0.);
      TLorentzVector Neu_Tau1_2cand_vec(0., 0., 0., 0.);
      TLorentzVector Neu_Tau1_3cand_vec(0., 0., 0., 0.);
      TLorentzVector Neu_Tau1_4cand_vec(0., 0., 0., 0.);
      TLorentzVector Neu_Tau2_1cand_vec(0., 0., 0., 0.);
      TLorentzVector Neu_Tau2_2cand_vec(0., 0., 0., 0.);
      TLorentzVector Neu_Tau2_3cand_vec(0., 0., 0., 0.);
      TLorentzVector Neu_Tau2_4cand_vec(0., 0., 0., 0.);
      bool is_b_jet = false;
      treeReader->ReadEntry(entry);
      int tau_conter=0;
      METpointer = (MissingET*) branchMissingET->At(0);
      double MET = METpointer->MET;
      double MET_phi = METpointer->Phi;
      int njets_counter = 0;
      int ntau_counter = 0; 
      pass_cuts[0] = 1;
      int countertau=0;
      //bool passed_jet_tau = false;
      //bool fill_tau1 = false;
      //bool fill_tau2 = false;
      //bool fill_tau3 = false;
      if(branchJet->GetEntries() > 0)
	{
	  // Take first jet
	  Jet *firstjet = (Jet*) branchJet->At(0);
	  // Plot jet transverse momentum
	  // For Jets
	  double jet_min_pt = 20.;
	  
	  // We need at least 2 jets, 1 (2) tau jets and the ISR jet.
	  
	  if (branchJet->GetEntriesFast() > 0) 
	    {
	      TLorentzVector jet_i(0., 0., 0., 0.);
	      TLorentzVector elec_i(0., 0., 0., 0.);
              double N_muons = 0.;
	      double N_elec  = 0.; 
	      for (int muo = 0; muo < branchMuon->GetEntriesFast(); muo++){
		Muon *muon = (Muon*) branchMuon->At(muo);
		if ((muon->PT > 10.) && (abs(muon->Eta) < 2.5)){N_muons++;}            
	      }
	      
              int index_tau1=0;
              int index_tau2=0;
              int index_tau3=0;
	      
	      //////Jet-tau fake rate///////////////////
	      
	      bool filled_tau1_vec = false;
	      bool filled_tau2_vec = false;
	      bool filled_tau3_vec = false;
	      bool filled_tau4_vec = false;
	      
	      for (int j = 0; j < branchJet->GetEntriesFast(); j++) {
                bool passed_jet_tau = false;
                Jet *jet = (Jet*) branchJet->At(j);
                double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
                jet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
                passed_jet_tau = TauIDJet(jet_i);
                //cout << "passed_jet_tau "<<passed_jet_tau<<endl;
                if((passed_jet_tau == true) && (filled_tau1_vec == false)){
		  Tau1Had_vec.SetPtEtaPhiE(jet_i.Pt(), jet_i.Eta(), jet_i.Phi(), jet_i.E());
		  filled_tau1_vec = true;
		  continue;
                }
                if((passed_jet_tau == true) && (filled_tau2_vec == false) && (filled_tau1_vec == true)){
		  Tau2Had_vec.SetPtEtaPhiE(jet_i.Pt(), jet_i.Eta(), jet_i.Phi(), jet_i.E());
		  filled_tau2_vec = true;
		  continue;
                }
                if((passed_jet_tau == true) && (filled_tau3_vec == false) && (filled_tau1_vec == true) && (filled_tau2_vec == true)){
		  Tau3Had_vec.SetPtEtaPhiE(jet_i.Pt(), jet_i.Eta(), jet_i.Phi(), jet_i.E());
		  filled_tau3_vec = true;
		  continue;
                }
                if((passed_jet_tau == true) && (filled_tau4_vec == false) && (filled_tau1_vec == true) && (filled_tau2_vec == true) && 
                   (filled_tau3_vec == true)){
		  Tau4Had_vec.SetPtEtaPhiE(jet_i.Pt(), jet_i.Eta(), jet_i.Phi(), jet_i.E());
		  filled_tau4_vec = true;
                }
              }
	      
              ////////////////////////////////////
	      for (int j = 0; j < branchJet->GetEntriesFast(); j++)
		{
		  bool passed_jet_tau = false;
                  
                  bool is_jet_elec_overlap = false;
		  Jet *jet = (Jet*) branchJet->At(j);
		  double jet_energy = calculateE(jet->Eta, jet->PT, jet->Mass);
		  jet_i.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
                  
		  // Remove overlaps of all jets with electrons
		  for (int el = 0; el < branchElectron->GetEntriesFast(); el++){
		    Electron *elec = (Electron*) branchElectron->At(el);
		    double elec_energy = calculateE(elec->Eta, elec->PT, 0.000510998902);
		    elec_i.SetPtEtaPhiE(elec->PT, elec->Eta, elec->Phi, elec_energy);
		    double DR_jet_elec = jet_i.DeltaR(elec_i);
		    if (DR_jet_elec < DR_jet_elec_max){ 
		      is_jet_elec_overlap = true;
		      break;
		    }   
		    if ((elec->PT > 10.) && (abs(elec->Eta)) < 2.5){N_elec++;}
		  }         
		  
		  if ((jet->PT > b_jet_pt_min) && (jet->BTag == 1)){is_b_jet = true;}
		  if ((jet->PT > jet_min_pt) && (!is_jet_elec_overlap) && ( jet->TauTag == 0 ) && (jet->BTag == 0)){
		    njets_counter++;
		    jet_min_pt = jet->PT;
                    Jet_leading_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy); 
		    Jet_leading_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
		    if ((Tau1Had_vec.Pt() > 2.0) && (Jet_leading_vec.DeltaR(Tau1Had_vec) < 0.3)){
		      Jet_leading_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		    }
		    if ((Tau2Had_vec.Pt() > 2.0) && (Jet_leading_vec.DeltaR(Tau2Had_vec) < 0.3)){
		      Jet_leading_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		    }
		    if ((Tau3Had_vec.Pt() > 2.0) && (Jet_leading_vec.DeltaR(Tau3Had_vec) < 0.3)){
		      Jet_leading_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		    }
		    if ((Tau4Had_vec.Pt() > 2.0) && (Jet_leading_vec.DeltaR(Tau4Had_vec) < 0.3)){
		      Jet_leading_vec.SetPtEtaPhiE(0., 0., 0., 0.); 
		    }
                    
		  }
		  // }
		}
	      // check if there is at least one jet tagged as a tau
	      
	      bool found_first_tau_lep = false;
	      bool found_second_tau_lep = false;
	      bool found_third_tau_lep = false;
              bool first_lep = false;
              bool second_lep = false;
              bool third_lep = false;
	      bool found_first_neu = false;
	      bool found_second_neu = false;
	      bool found_third_neu = false;
	      bool found_first_neu_lep = false;
	      bool found_second_neu_lep = false;
	      bool found_third_neu_lep = false;
	      bool tau1 = false;
	      bool tau1_2 = false;
	      bool tau1_3 = false;
	      bool tau1_4 = false;
	      bool tau2 = false;
	      bool tau2_2 = false;
	      bool tau2_3 = false;
	      bool tau2_4 = false;	       
	      //cout << "------------------------------------------Inicio de Evento--------------------------------------"<<endl;
              // Loop over GenTau particles and emulate the identification and 
              // reconstruction of hadronic taus. 
	      
	      for(int tau_i = 0; tau_i < branchGenParticle->GetEntriesFast(); tau_i++){
		// Look for the particle PDGID
		GenParticle *tau_mother = (GenParticle*) branchGenParticle->At(tau_i);
		//cout <<tau_mother->PID<<", ";
		if(((abs(tau_mother->PID) == 23) || (abs(tau_mother->PID) == 24) || (abs(tau_mother->PID) == 15) ||  
		    (abs(tau_mother->PID) == 1000015) || (abs(tau_mother->PID) == 1000023)) && (tau_mother->Status == 2) ){
		  // Search for the daughters
		  GenParticle *daughter_1 = (GenParticle*) branchGenParticle->At(tau_mother->D1);
		  GenParticle *daughter_2 = (GenParticle*) branchGenParticle->At(tau_mother->D2);
		  //cout << "hija 1: "<< daughter_1->PID<<"  Hija 2: "<<daughter_2->PID<<endl;
		  // Find tau candidates that come from D1               
		  if((abs(daughter_1->PID) == 15)){
		    // Energy is not stored in the Delphes ntuple, so we calculate it.
		    double tau1cand_energy = calculateE(daughter_1->Eta, daughter_1->PT, daughter_1->Mass);
                    // Need to look at the tau daughters in order to identify 
                    // the associated neutrinos
		    GenParticle *daughther1_Tau1cand = (GenParticle*) branchGenParticle->At(daughter_1->D1);
		    GenParticle *daughther2_Tau1cand = (GenParticle*) branchGenParticle->At(daughter_1->D2);
                    double Neu_tau1_energy_d1 = 0.0;
                    double Neu_tau1_energy_d2 = 0.0;
		    
                    if (daughther1_Tau1cand->Status == 1){
		      Neu_tau1_energy_d1 = calculateE(daughther1_Tau1cand->Eta, daughther1_Tau1cand->PT, daughther1_Tau1cand->Mass);
		    }
		    if (daughther2_Tau1cand->Status == 1){
		      Neu_tau1_energy_d2 = calculateE(daughther2_Tau1cand->Eta, daughther2_Tau1cand->PT, daughther2_Tau1cand->Mass);
		    }
		    
                    // These booleans are used to keep track of the filled TLorentz vectors 
                    // that store the information of the tau_had candidates
		    
		    if((tau1 == false) && (tau1_2 == false) && (tau1_3 == false)){
		      
		      tau1 = true;
		      
                      // Fill the first TLorentz vector. This tau might be leptonic or hadronic 
                      // Later in the code we determine its nature.
		      Tau1_1cand_vec.SetPtEtaPhiE(daughter_1->PT, daughter_1->Eta, daughter_1->Phi, tau1cand_energy);
                      if((Tau1_1cand_vec.DeltaR(Tau2_1cand_vec) < 0.2) || (Tau1_1cand_vec.DeltaR(Tau2_2cand_vec) < 0.2) ||
                         (Tau1_1cand_vec.DeltaR(Tau2_3cand_vec) < 0.2) || (Tau1_1cand_vec.DeltaR(Tau2_4cand_vec) < 0.2)){
			Tau1_1cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
			tau1 = false;
                      }
		      
                      // Save the information for tau neutrinos. This is important to identidy later on 
                      // if the tau us hadronic or leptonic 
                      
		      if((abs(daughther1_Tau1cand->PID)==16) && (daughther1_Tau1cand->Status == 1)){
			Neu_Tau1_1cand_vec.SetPtEtaPhiE(daughther1_Tau1cand->PT, daughther1_Tau1cand->Eta, daughther1_Tau1cand->Phi, Neu_tau1_energy_d1);
			if((Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_1cand_vec) < 0.2) || (Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2) ||
			   (Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) || (Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2)){
                          Neu_Tau1_1cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
			}  
		      }
		      if((abs(daughther2_Tau1cand->PID)==16) && (daughther2_Tau1cand->Status == 1)){ 
                        Neu_Tau1_1cand_vec.SetPtEtaPhiE(daughther2_Tau1cand->PT, daughther2_Tau1cand->Eta, daughther2_Tau1cand->Phi, Neu_tau1_energy_d2);                       
                        if((Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_1cand_vec) < 0.2) || (Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2) ||
			   (Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) || (Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2)){
                          Neu_Tau1_1cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
			}  
		      }
		    }
                    // If the first TLorentz vector is been already filled, then 
                    // fill the next vector. The code down below is the same as above. 
		    if((tau1 == true) && (tau1_2 == false) && (tau1_3 == false)){
		      Tau1_2cand_vec.SetPtEtaPhiE(daughter_1->PT, daughter_1->Eta, daughter_1->Phi, tau1cand_energy);
                      tau1_2 = true;
                      if((Tau1_2cand_vec.DeltaR(Tau1_1cand_vec) < 0.2) || (Tau1_2cand_vec.DeltaR(Tau2_1cand_vec) < 0.2) || 
                         (Tau1_2cand_vec.DeltaR(Tau2_2cand_vec) < 0.2) || (Tau1_2cand_vec.DeltaR(Tau2_3cand_vec) < 0.2) || 
                         (Tau1_2cand_vec.DeltaR(Tau2_4cand_vec) < 0.2)){
			Tau1_2cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
			tau1_2 = false;
                      }
		      
		      if((abs(daughther1_Tau1cand->PID)==16) && (daughther1_Tau1cand->Status == 1)){
                        Neu_Tau1_2cand_vec.SetPtEtaPhiE(daughther1_Tau1cand->PT, daughther1_Tau1cand->Eta, daughther1_Tau1cand->Phi, Neu_tau1_energy_d1);
                        if((Neu_Tau1_2cand_vec.DeltaR(Neu_Tau1_1cand_vec) < 0.2) || (Neu_Tau1_2cand_vec.DeltaR(Neu_Tau2_1cand_vec) < 0.2) ||
			   (Neu_Tau1_2cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2) || (Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) ||
			   (Neu_Tau1_2cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2)){
			  Neu_Tau1_2cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
			}
		      }
		      if((abs(daughther2_Tau1cand->PID)==16) && (daughther2_Tau1cand->Status == 1)){
			Neu_Tau1_2cand_vec.SetPtEtaPhiE(daughther2_Tau1cand->PT, daughther2_Tau1cand->Eta, daughther2_Tau1cand->Phi, Neu_tau1_energy_d2);
		        if((Neu_Tau1_2cand_vec.DeltaR(Neu_Tau1_1cand_vec) < 0.2) || (Neu_Tau1_2cand_vec.DeltaR(Neu_Tau2_1cand_vec) < 0.2) ||
			   (Neu_Tau1_2cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2) || (Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) ||
			   (Neu_Tau1_2cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2)){
			  Neu_Tau1_2cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
			}
		      }
                    }
                    // Same for TLV3
                    if((tau1 == true) && (tau1_2 == true) && (tau1_3 == false)){
                      Tau1_3cand_vec.SetPtEtaPhiE(daughter_1->PT, daughter_1->Eta, daughter_1->Phi, tau1cand_energy);
                      tau1_3 = true;
		      if((Tau1_3cand_vec.DeltaR(Tau1_1cand_vec) < 0.2) || (Tau1_3cand_vec.DeltaR(Tau1_2cand_vec) < 0.2) || 
                         (Tau1_3cand_vec.DeltaR(Tau2_1cand_vec) < 0.2) || (Tau1_3cand_vec.DeltaR(Tau2_2cand_vec) < 0.2) || 
                         (Tau1_3cand_vec.DeltaR(Tau2_3cand_vec) < 0.2) || (Tau1_3cand_vec.DeltaR(Tau2_4cand_vec) < 0.2) ){
			Tau1_3cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
			tau1_3 = false;
		      }
		      
                      if((abs(daughther1_Tau1cand->PID)==16) && (daughther1_Tau1cand->Status == 1)){
			Neu_Tau1_3cand_vec.SetPtEtaPhiE(daughther1_Tau1cand->PT, daughther1_Tau1cand->Eta, daughther1_Tau1cand->Phi, Neu_tau1_energy_d1);
		        if((Neu_Tau1_3cand_vec.DeltaR(Neu_Tau1_1cand_vec) < 0.2) || (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau1_2cand_vec) < 0.2) || 
                           (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau2_1cand_vec) < 0.2) || (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2) || 
                           (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) || (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2) ){
			  Neu_Tau1_3cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
			}
                      }
		      if((abs(daughther2_Tau1cand->PID)==16) && (daughther2_Tau1cand->Status == 1)){
			Neu_Tau1_3cand_vec.SetPtEtaPhiE(daughther2_Tau1cand->PT, daughther2_Tau1cand->Eta, daughther2_Tau1cand->Phi, Neu_tau1_energy_d2);
                        if((Neu_Tau1_3cand_vec.DeltaR(Neu_Tau1_1cand_vec) < 0.2) || (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau1_2cand_vec) < 0.2) || 
                           (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau2_1cand_vec) < 0.2) || (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2) ||
                           (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) || (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2) ){ 
			  Neu_Tau1_3cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
			}
		      }
                    }
                    // Same for TLV4
		    if((tau1 == true) && (tau1_2 == true) && (tau1_3 == true)){ 
                      Tau1_4cand_vec.SetPtEtaPhiE(daughter_1->PT, daughter_1->Eta, daughter_1->Phi, tau1cand_energy);
		      if((Tau1_4cand_vec.DeltaR(Tau1_1cand_vec) < 0.2) || (Tau1_4cand_vec.DeltaR(Tau1_2cand_vec) < 0.2) || 
                         (Tau1_4cand_vec.DeltaR(Tau1_3cand_vec) < 0.2) || (Tau1_4cand_vec.DeltaR(Tau2_1cand_vec) < 0.2) || 
                         (Tau1_4cand_vec.DeltaR(Tau2_2cand_vec) < 0.2) || (Tau1_4cand_vec.DeltaR(Tau2_3cand_vec) < 0.2) ||
                         (Tau1_4cand_vec.DeltaR(Tau2_4cand_vec) < 0.2) ){
			Tau1_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);  
		      }
                      if((abs(daughther1_Tau1cand->PID)==16) && (daughther1_Tau1cand->Status == 1)){
			Neu_Tau1_4cand_vec.SetPtEtaPhiE(daughther1_Tau1cand->PT, daughther1_Tau1cand->Eta, daughther1_Tau1cand->Phi, Neu_tau1_energy_d1);
                        if((Neu_Tau1_4cand_vec.DeltaR(Neu_Tau1_1cand_vec) < 0.2) || (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau1_2cand_vec) < 0.2) ||
                           (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau1_3cand_vec) < 0.2) || (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau2_1cand_vec) < 0.2) || 
                           (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2) || (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) || 
                           (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2) ){
			  Neu_Tau1_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);  
                        }
		      }
		      if((abs(daughther2_Tau1cand->PID)==16) && (daughther2_Tau1cand->Status == 1)){
			Neu_Tau1_4cand_vec.SetPtEtaPhiE(daughther2_Tau1cand->PT, daughther2_Tau1cand->Eta, daughther2_Tau1cand->Phi, Neu_tau1_energy_d2);
		        if((Neu_Tau1_4cand_vec.DeltaR(Neu_Tau1_1cand_vec) < 0.2) || (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau1_2cand_vec) < 0.2) ||
                           (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau1_3cand_vec) < 0.2) || (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau2_1cand_vec) < 0.2) || 
                           (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2) || (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) || 
                           (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2) ){
			  Neu_Tau1_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);  
			}
		      }
		    }
		  } // Close: if((abs(daughter_1->PID) == 15))
                  // Now do the same as for the "daughter_1" case, but for "daughter_2".
                  // We will fill a different set of TLVs and later on we will look which vectors
                  // are filled and we will select at most 4 taus. If there are more than 4 taus, 
                  // we will selecto those with the highest pT.
		  if((abs(daughter_2->PID) == 15)){
		    double tau2cand_energy = calculateE(daughter_2->Eta, daughter_2->PT, daughter_2->Mass);
		    GenParticle *daughther1_Tau2cand = (GenParticle*) branchGenParticle->At(daughter_2->D1);
		    GenParticle *daughther2_Tau2cand = (GenParticle*) branchGenParticle->At(daughter_2->D2);
		    double Neu_tau2_energy_d1 = 0.0;
		    double Neu_tau2_energy_d2 = 0.0;
		    if (daughther1_Tau2cand->Status == 1){
		      Neu_tau2_energy_d1 = calculateE(daughther1_Tau2cand->Eta, daughther1_Tau2cand->PT, daughther1_Tau2cand->Mass);
		    }
		    if (daughther2_Tau2cand->Status == 1){
		      Neu_tau2_energy_d2 = calculateE(daughther2_Tau2cand->Eta, daughther2_Tau2cand->PT, daughther2_Tau2cand->Mass); 
		    }
		    if((tau2 == false) && (tau2_2 == false) && (tau2_3 == false)){
		      tau2 = true;
		      Tau2_1cand_vec.SetPtEtaPhiE(daughter_2->PT, daughter_2->Eta, daughter_2->Phi, tau2cand_energy);
		      if((Tau1_1cand_vec.DeltaR(Tau2_1cand_vec) < 0.2) || (Tau1_2cand_vec.DeltaR(Tau2_1cand_vec) < 0.2) || 
                         (Tau1_3cand_vec.DeltaR(Tau2_1cand_vec) < 0.2) || (Tau1_4cand_vec.DeltaR(Tau2_1cand_vec) < 0.2)){
			Tau2_1cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
			tau2 = false;
		      }
		      if((abs(daughther1_Tau2cand->PID)==16) && (daughther1_Tau2cand->Status == 1)){
			Neu_Tau2_1cand_vec.SetPtEtaPhiE(daughther1_Tau2cand->PT, daughther1_Tau2cand->Eta, daughther1_Tau2cand->Phi, Neu_tau2_energy_d1);
			if((Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_1cand_vec) < 0.2) || (Neu_Tau1_2cand_vec.DeltaR(Neu_Tau2_1cand_vec) < 0.2) ||
			   (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau2_1cand_vec) < 0.2) || (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau2_1cand_vec) < 0.2)){
			  Neu_Tau2_1cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
			}
		      }
		      
		      if((abs(daughther2_Tau2cand->PID)==16) && (daughther2_Tau2cand->Status == 1)){
			Neu_Tau2_1cand_vec.SetPtEtaPhiE(daughther2_Tau2cand->PT, daughther2_Tau2cand->Eta, daughther2_Tau2cand->Phi, Neu_tau2_energy_d2);
                        if((Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_1cand_vec) < 0.2) || (Neu_Tau1_2cand_vec.DeltaR(Neu_Tau2_1cand_vec) < 0.2) ||
			   (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau2_1cand_vec) < 0.2) || (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau2_1cand_vec) < 0.2)){
			  Neu_Tau2_1cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
			}
		      }   
                    }
		    
		    if((tau2 == true) && (tau2_2 == false) && (tau2_3 == false)){
		      Tau2_2cand_vec.SetPtEtaPhiE(daughter_2->PT, daughter_2->Eta, daughter_2->Phi, tau2cand_energy);
		      tau2_2 = true;
		      if((Tau2_2cand_vec.DeltaR(Tau2_1cand_vec) < 0.2) || (Tau2_2cand_vec.DeltaR(Tau1_1cand_vec) < 0.2) || 
                         (Tau2_2cand_vec.DeltaR(Tau1_2cand_vec) < 0.2) || (Tau2_2cand_vec.DeltaR(Tau1_3cand_vec) < 0.2) || 
                         (Tau2_2cand_vec.DeltaR(Tau1_4cand_vec) < 0.2)){
			Tau2_2cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
			tau2_2 = false;
		      }
		      
		      if((abs(daughther1_Tau2cand->PID)==16) && (daughther1_Tau2cand->Status == 1)){
			Neu_Tau2_2cand_vec.SetPtEtaPhiE(daughther1_Tau2cand->PT, daughther1_Tau2cand->Eta, daughther1_Tau2cand->Phi, Neu_tau2_energy_d1);
                        if((Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2) || (Neu_Tau1_2cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2) ||
                           (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2) || (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2) || 
                           (Neu_Tau2_1cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2)){
			  Neu_Tau2_2cand_vec.SetPtEtaPhiE(0., 0., 0., 0.); 		      
			}
		      }
		      if((abs(daughther2_Tau2cand->PID)==16) && (daughther2_Tau2cand->Status == 1)){
			Neu_Tau2_2cand_vec.SetPtEtaPhiE(daughther2_Tau2cand->PT, daughther2_Tau2cand->Eta, daughther2_Tau2cand->Phi, Neu_tau2_energy_d2);
		        if((Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2) || (Neu_Tau1_2cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2) ||
                           (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2) || (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2) || 
                           (Neu_Tau2_1cand_vec.DeltaR(Neu_Tau2_2cand_vec) < 0.2)){
			  Neu_Tau2_2cand_vec.SetPtEtaPhiE(0., 0., 0., 0.); 
			}
		      }
		    }
		    
		    if((tau2 == true) && (tau2_2 == true) && (tau2_3 == false)){
		      Tau2_3cand_vec.SetPtEtaPhiE(daughter_2->PT, daughter_2->Eta, daughter_2->Phi, tau2cand_energy);
		      tau2_3 = true;
		      if((Tau2_3cand_vec.DeltaR(Tau2_1cand_vec) < 0.2) || (Tau2_3cand_vec.DeltaR(Tau2_2cand_vec) < 0.2) || 
                         (Tau2_3cand_vec.DeltaR(Tau1_1cand_vec) < 0.2) || (Tau2_3cand_vec.DeltaR(Tau1_2cand_vec) < 0.2) || 
                         (Tau2_3cand_vec.DeltaR(Tau1_3cand_vec) < 0.2) || (Tau2_3cand_vec.DeltaR(Tau1_4cand_vec) < 0.2)){
			Tau2_3cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
			tau2_3 = false;
		      }
		      if((abs(daughther1_Tau2cand->PID)==16) && (daughther1_Tau2cand->Status == 1)){
			Neu_Tau2_3cand_vec.SetPtEtaPhiE(daughther1_Tau2cand->PT, daughther1_Tau2cand->Eta, daughther1_Tau2cand->Phi, Neu_tau2_energy_d1);
		        if((Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) || (Neu_Tau1_2cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) ||
                           (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) || (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) || 
                           (Neu_Tau2_1cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) || (Neu_Tau2_2cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2)){
			  Neu_Tau2_3cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
			}
                      }
		      if((abs(daughther2_Tau2cand->PID)==16) && (daughther2_Tau2cand->Status == 1)){
			Neu_Tau2_3cand_vec.SetPtEtaPhiE(daughther2_Tau2cand->PT, daughther2_Tau2cand->Eta, daughther2_Tau2cand->Phi, Neu_tau2_energy_d2);
		        if((Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) || (Neu_Tau1_2cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) ||
                           (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) || (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) || 
                           (Neu_Tau2_1cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2) || (Neu_Tau2_2cand_vec.DeltaR(Neu_Tau2_3cand_vec) < 0.2)){
			  Neu_Tau2_3cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
                        }
                      }
		    }
		    
		    if((tau2 == true) && (tau2_2 == true) && (tau2_3 == true)){
		      Tau2_4cand_vec.SetPtEtaPhiE(daughter_2->PT, daughter_2->Eta, daughter_2->Phi, tau2cand_energy);
		      if((Tau2_4cand_vec.DeltaR(Tau2_1cand_vec) < 0.2) || (Tau2_4cand_vec.DeltaR(Tau2_2cand_vec) < 0.2) || 
                         (Tau2_4cand_vec.DeltaR(Tau2_3cand_vec) < 0.2) || (Tau2_4cand_vec.DeltaR(Tau1_1cand_vec) < 0.2) || 
                         (Tau2_4cand_vec.DeltaR(Tau1_2cand_vec) < 0.2) || (Tau2_4cand_vec.DeltaR(Tau1_3cand_vec) < 0.2) ||
			 (Tau2_4cand_vec.DeltaR(Tau1_4cand_vec) < 0.2)){
			Tau2_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		      }
		      if((abs(daughther1_Tau2cand->PID)==16) && (daughther1_Tau2cand->Status == 1)){
			Neu_Tau2_4cand_vec.SetPtEtaPhiE(daughther1_Tau2cand->PT, daughther1_Tau2cand->Eta, daughther1_Tau2cand->Phi, Neu_tau2_energy_d1);
                        if((Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2) || (Neu_Tau1_2cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2) ||
                           (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2) || (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2) || 
                           (Neu_Tau2_1cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2) || (Neu_Tau2_2cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2) || 
                           (Neu_Tau2_3cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2)){
			  Neu_Tau2_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);	      
			}
		      }
		      if((abs(daughther2_Tau2cand->PID)==16) && (daughther2_Tau2cand->Status == 1)){
			Neu_Tau2_4cand_vec.SetPtEtaPhiE(daughther2_Tau2cand->PT, daughther2_Tau2cand->Eta, daughther2_Tau2cand->Phi, Neu_tau2_energy_d2);
			if((Neu_Tau1_1cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2) || (Neu_Tau1_2cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2) ||
			   (Neu_Tau1_3cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2) || (Neu_Tau1_4cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2) || 
			   (Neu_Tau2_1cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2) || (Neu_Tau2_2cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2) || 
			   (Neu_Tau2_3cand_vec.DeltaR(Neu_Tau2_4cand_vec) < 0.2)){
			  Neu_Tau2_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
			}
		      }
		    }
		  }  // close:  if((abs(daughter_2->PID) == 15))
		}   // close: if(tau_mother->PID = 23, 24, 1000015, 1000023)
		
		if (((abs(tau_mother->PID) == 11) || (abs(tau_mother->PID) == 13)) && (tau_mother->Status == 1)){
		  GenParticle *mother_elec = (GenParticle*) branchGenParticle->At(tau_mother->M1);  
		  if (abs(mother_elec->PID) == 15){
		    double lep_e = calculateE(tau_mother->Eta, tau_mother->PT, tau_mother->Mass);
		    double MotherElcand_energy = calculateE(mother_elec->Eta, mother_elec->PT, mother_elec->Mass);
		    if (found_first_tau_lep == false && found_second_tau_lep == false && found_third_tau_lep == false) { 
		      TauLep1cand_vec.SetPtEtaPhiE(mother_elec->PT, mother_elec->Eta, mother_elec->Phi, MotherElcand_energy);
		      Lepton1_vec.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, lep_e);
                      found_first_tau_lep == true;    
		    }
		    if(found_first_tau_lep == true && found_second_tau_lep == false && found_third_tau_lep == false) { 
		      TauLep2cand_vec.SetPtEtaPhiE(mother_elec->PT, mother_elec->Eta, mother_elec->Phi, MotherElcand_energy);
		      Lepton2_vec.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, lep_e);
                      found_second_tau_lep == true;
		    }
                    if(found_first_tau_lep == true && found_second_tau_lep == true && found_third_tau_lep == false) { 
		      TauLep3cand_vec.SetPtEtaPhiE(mother_elec->PT, mother_elec->Eta, mother_elec->Phi, MotherElcand_energy);
		      Lepton3_vec.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, lep_e);
                      found_third_tau_lep == true;
		    }
                    if(found_first_tau_lep == true && found_second_tau_lep == true && found_third_tau_lep == true) { 
		      TauLep4cand_vec.SetPtEtaPhiE(mother_elec->PT, mother_elec->Eta, mother_elec->Phi, MotherElcand_energy);
                      Lepton4_vec.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, lep_e);
		    }
		  }
		}
		if ((abs(tau_mother->PID) == 16) &&  (tau_mother->Status == 1)){
		  double TauNeu_energy = calculateE(tau_mother->Eta, tau_mother->PT, tau_mother->Mass);
		  if ((found_first_neu == false) && (found_second_neu == false) && (found_third_neu == false)){
		    NeuTau1.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, TauNeu_energy);
		    found_first_neu = true;
		  } 
		  if((found_first_neu == true) && (found_second_neu == false) && (found_third_neu == false)) {
		    NeuTau2.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, TauNeu_energy);
		    found_second_neu == true;
		  }
		  if((found_first_neu == true) && (found_second_neu == true) && (found_third_neu == false)) {
		    NeuTau3.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, TauNeu_energy);
		    found_third_neu == true;
		  }
                  if((found_first_neu == true) && (found_second_neu == true) && (found_third_neu == true)) {
		    NeuTau4.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, TauNeu_energy);
		  }
		}
                
		if(((abs(tau_mother->PID) == 12) || (abs(tau_mother->PID) == 14)) && (tau_mother->Status == 1)){
		  double LepNeu_energy = calculateE(tau_mother->Eta, tau_mother->PT, tau_mother->Mass);
		  if (found_first_neu_lep == false && (found_second_neu_lep == false) && (found_third_neu_lep == false)){
		    NeuLep1cand_vec.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, LepNeu_energy);
		    found_first_neu_lep = true;
		  } 
		  if((found_first_neu_lep == true) && (found_second_neu_lep == false) && (found_third_neu_lep == false)){
		    NeuLep2cand_vec.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, LepNeu_energy);
		    found_second_neu_lep = true; 
		  }
		  if((found_first_neu_lep == true) && (found_second_neu_lep == true) && (found_third_neu_lep == false)){
		    NeuLep3cand_vec.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, LepNeu_energy);
		    found_third_neu_lep = true; 
		  }
		  if((found_first_neu_lep == true) && (found_second_neu_lep == true) && (found_third_neu_lep == true)){
		    NeuLep4cand_vec.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, LepNeu_energy);
		  } 
		}
	      }
	      
	      //cout << "Momento neutrino tau1Cand: "<<Neu_Tau1_1cand_vec.Pt()<<"  Momento neutrino tau2Cand: "<<Neu_Tau2_1cand_vec.Pt()<<endl;
	      //cout <<endl;
              //cout <<"------------------------------------------Fin del evento--------------------------------"<<endl;
              //cout <<endl;
              /*cout <<"Tau1_1cand_vec.Pt "<<Tau1_1cand_vec.Pt()<<endl;
		cout <<"Tau1_2cand_vec.Pt "<<Tau1_2cand_vec.Pt()<<endl;
		cout <<"Tau1_3cand_vec.Pt "<<Tau1_3cand_vec.Pt()<<endl;
		cout <<"Tau1_4cand_vec.Pt "<<Tau1_4cand_vec.Pt()<<endl;
		cout <<"Tau2_1cand_vec.Pt "<<Tau2_1cand_vec.Pt()<<endl;
		cout <<"Tau2_2cand_vec.Pt "<<Tau2_2cand_vec.Pt()<<endl;
		cout <<"Tau2_3cand_vec.Pt "<<Tau2_3cand_vec.Pt()<<endl;
		cout <<"Tau2_4cand_vec.Pt "<<Tau2_4cand_vec.Pt()<<endl;
		cout <<"--------------Taus leptonicos----------------- "<<endl;
		cout << "TauLep1cand_vec "<<TauLep1cand_vec.Pt()<<endl;
		cout << "TauLep2cand_vec "<<TauLep2cand_vec.Pt()<<endl;
		cout << "TauLep3cand_vec "<<TauLep3cand_vec.Pt()<<endl;
		cout << "TauLep4cand_vec "<<TauLep4cand_vec.Pt()<<endl;
	      */
              bool tau1_1fill = false;
              bool tau1_2fill = false;
              bool tau1_3fill = false;
              bool tau1_4fill = false;
              bool tau2_1fill = false;
              bool tau2_2fill = false;
              bool tau2_3fill = false;
              bool tau2_4fill = false;
	      
              if(Tau1_1cand_vec.Pt() > 2.0 ){tau1_1fill = true;}
              if(Tau1_2cand_vec.Pt() > 2.0 ){tau1_2fill = true;}
              if(Tau1_3cand_vec.Pt() > 2.0 ){tau1_3fill = true;}
              if(Tau1_4cand_vec.Pt() > 2.0 ){tau1_4fill = true;}  
              if(Tau2_1cand_vec.Pt() > 2.0 ){tau2_1fill = true;}
              if(Tau2_2cand_vec.Pt() > 2.0 ){tau2_2fill = true;}
              if(Tau2_3cand_vec.Pt() > 2.0 ){tau2_3fill = true;}
              if(Tau2_4cand_vec.Pt() > 2.0 ){tau2_4fill = true;}
              
	      /////Si tenemos un tau en el evento se tiene 2 posibilidades ///////
	      
	      if((tau1_1fill == true) && (tau1_2fill == false) && (tau1_3fill == false) && (tau1_4fill == false) && (tau2_1fill == false) && 
		 (tau2_2fill == false) && (tau2_3fill == false) && (tau2_4fill == false)){
		Tau1cand_vec = Tau1_1cand_vec;
		NeuTau1cand_vec = Neu_Tau1_1cand_vec;
		Tau2cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Tau3cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Tau4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		NeuTau2cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		NeuTau3cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		NeuTau4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.); 
	      }
              
	      if((tau1_1fill == false) && (tau1_2fill == false) && (tau1_3fill == false) && (tau1_4fill == false) && (tau2_1fill == true) && 
		 (tau2_2fill == false) && (tau2_3fill == false) && (tau2_4fill == false)){
		Tau1cand_vec = Tau2_1cand_vec;
		NeuTau1cand_vec = Neu_Tau2_1cand_vec;  
		Tau2cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Tau3cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Tau4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		NeuTau2cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		NeuTau3cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		NeuTau4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);  
	      }
	      ////////////////////////////////////////////////////////////////////
	      //////Si tenemos dos taus en el evento se tiene tres posibilidades ///////////////////
	      if((tau1_1fill == true) && (tau1_2fill == true) && (tau1_3fill == false) && (tau1_4fill == false) && (tau2_1fill == false) && 
		 (tau2_2fill == false) && (tau2_3fill == false) && (tau2_4fill == false)){
		Tau1_3cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Tau1_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Neu_Tau1_3cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Neu_Tau1_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);   
		Mayor_menorPt(Tau1_1cand_vec, Tau1_2cand_vec, Tau1_3cand_vec, Tau1_4cand_vec, Neu_Tau1_1cand_vec, Neu_Tau1_2cand_vec, Neu_Tau1_3cand_vec, Neu_Tau1_4cand_vec,
                              &Tau1cand_vec, &Tau2cand_vec, &Tau3cand_vec, &Tau4cand_vec, &NeuTau1cand_vec, &NeuTau2cand_vec, &NeuTau3cand_vec, &NeuTau4cand_vec);
	      }
	      
	      if((tau1_1fill == true) && (tau1_2fill == false) && (tau1_3fill == false) && (tau1_4fill == false) && (tau2_1fill == true) && 
		 (tau2_2fill == false) && 
                 (tau2_3fill == false) && (tau2_4fill == false)){
		Tau1_3cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Tau1_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Neu_Tau1_3cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Neu_Tau1_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Mayor_menorPt(Tau1_1cand_vec, Tau2_1cand_vec, Tau1_3cand_vec, Tau1_4cand_vec, Neu_Tau1_1cand_vec, Neu_Tau2_1cand_vec, Neu_Tau1_3cand_vec, Neu_Tau1_4cand_vec, 
                              &Tau1cand_vec, &Tau2cand_vec, &Tau3cand_vec, &Tau4cand_vec, &NeuTau1cand_vec, &NeuTau2cand_vec, &NeuTau3cand_vec, &NeuTau4cand_vec);
	      }
	      if((tau1_1fill == false) && (tau1_2fill == false) && (tau1_3fill == false) && (tau1_4fill == false) && (tau2_1fill == true) && 
		 (tau2_2fill == true) && 
                 (tau2_3fill == false) && (tau2_4fill == false)){
		Tau1_3cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Tau1_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Neu_Tau1_3cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Neu_Tau1_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Mayor_menorPt(Tau2_1cand_vec, Tau2_2cand_vec, Tau1_3cand_vec, Tau1_4cand_vec, Neu_Tau2_1cand_vec, Neu_Tau2_2cand_vec, Neu_Tau1_3cand_vec, Neu_Tau1_4cand_vec, 
                              &Tau1cand_vec, &Tau2cand_vec, &Tau3cand_vec, &Tau4cand_vec, &NeuTau1cand_vec, &NeuTau2cand_vec, &NeuTau3cand_vec, &NeuTau4cand_vec);
	      }
	      ///////////////////////////////////////////////////////////////////
	      //// Si tenemos tres taus en el evento tenemos cuatro posibilidades ///////////////////////////
	      if((tau1_1fill == true) && (tau1_2fill == true) && (tau1_3fill == true) && (tau1_4fill == false) && (tau2_1fill == false) && (tau2_2fill == false) && 
                 (tau2_3fill == false) && (tau2_4fill == false)){
		Tau1_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Neu_Tau1_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.); 
		Mayor_menorPt(Tau1_1cand_vec, Tau1_2cand_vec, Tau1_3cand_vec, Tau1_4cand_vec, Neu_Tau1_1cand_vec, Neu_Tau1_2cand_vec, Neu_Tau1_3cand_vec, Neu_Tau1_4cand_vec, 
                              &Tau1cand_vec, &Tau2cand_vec, &Tau3cand_vec, &Tau4cand_vec, &NeuTau1cand_vec, &NeuTau2cand_vec, &NeuTau3cand_vec, &NeuTau4cand_vec);
	      }
	      
	      if((tau1_1fill == true) && (tau1_2fill == true) && (tau1_3fill == false) && (tau1_4fill == false) && (tau2_1fill == true) && (tau2_2fill == false) && 
                 (tau2_3fill == false) && (tau2_4fill == false)){
		Tau1_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Neu_Tau1_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);  
		Mayor_menorPt(Tau1_1cand_vec, Tau1_2cand_vec, Tau2_1cand_vec, Tau1_4cand_vec, Neu_Tau1_1cand_vec, Neu_Tau1_2cand_vec, Neu_Tau2_1cand_vec, Neu_Tau1_4cand_vec, 
                              &Tau1cand_vec, &Tau2cand_vec, &Tau3cand_vec, &Tau4cand_vec, &NeuTau1cand_vec, &NeuTau2cand_vec, &NeuTau3cand_vec, &NeuTau4cand_vec); 
	      }
	      if((tau1_1fill == true) && (tau1_2fill == false) && (tau1_3fill == false) && (tau1_4fill == false) && (tau2_1fill == true) && (tau2_2fill == true) && 
                 (tau2_3fill == false) && (tau2_4fill == false)){
		Tau1_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Neu_Tau1_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Mayor_menorPt(Tau1_1cand_vec, Tau2_1cand_vec, Tau2_2cand_vec, Tau1_4cand_vec, Neu_Tau1_1cand_vec, Neu_Tau2_1cand_vec, Neu_Tau2_2cand_vec, Neu_Tau1_4cand_vec, 
                              &Tau1cand_vec, &Tau2cand_vec, &Tau3cand_vec, &Tau4cand_vec, &NeuTau1cand_vec, &NeuTau2cand_vec, &NeuTau3cand_vec, &NeuTau4cand_vec);
	      }
	      if((tau1_1fill == false) && (tau1_2fill == false) && (tau1_3fill == false) && (tau1_4fill == false) && (tau2_1fill == true) && (tau2_2fill == true) && 
                 (tau2_3fill == true) && (tau2_4fill == false)){
		Tau1_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Neu_Tau1_4cand_vec.SetPtEtaPhiE(0., 0., 0., 0.);
		Mayor_menorPt(Tau2_1cand_vec, Tau2_2cand_vec, Tau2_3cand_vec, Tau1_4cand_vec, Neu_Tau2_1cand_vec, Neu_Tau2_2cand_vec, Neu_Tau2_3cand_vec, Neu_Tau1_4cand_vec, 
                              &Tau1cand_vec, &Tau2cand_vec, &Tau3cand_vec, &Tau4cand_vec, &NeuTau1cand_vec, &NeuTau2cand_vec, &NeuTau3cand_vec, &NeuTau4cand_vec);  
	      }
	      
	      ///////////////Si hay cuatro taus se tienen 5 posibilidades//////////////////////////////////////////////////
	      if((tau1_1fill == true) && (tau1_2fill == true) && (tau1_3fill == true) && (tau1_4fill == true) && (tau2_1fill == false) && (tau2_2fill == false) && 
                 (tau2_3fill == false) && (tau2_4fill == false)){
		Mayor_menorPt(Tau1_1cand_vec, Tau1_2cand_vec, Tau1_3cand_vec, Tau1_4cand_vec, Neu_Tau1_1cand_vec, Neu_Tau1_2cand_vec, Neu_Tau1_3cand_vec, Neu_Tau1_4cand_vec, 
                              &Tau1cand_vec, &Tau2cand_vec, &Tau3cand_vec, &Tau4cand_vec, &NeuTau1cand_vec, &NeuTau2cand_vec, &NeuTau3cand_vec, &NeuTau4cand_vec);
	      }
	      
	      if((tau1_1fill == true) && (tau1_2fill == true) && (tau1_3fill == true) && (tau1_4fill == false) && (tau2_1fill == true) && (tau2_2fill == false) &&
                 (tau2_3fill == false) && (tau2_4fill == false)){
		Mayor_menorPt(Tau1_1cand_vec, Tau1_2cand_vec, Tau1_3cand_vec, Tau2_1cand_vec, Neu_Tau1_1cand_vec, Neu_Tau1_2cand_vec, Neu_Tau1_3cand_vec, Neu_Tau2_1cand_vec, 
                              &Tau1cand_vec, &Tau2cand_vec, &Tau3cand_vec, &Tau4cand_vec, &NeuTau1cand_vec, &NeuTau2cand_vec, &NeuTau3cand_vec, &NeuTau4cand_vec);  
	      }
	      
	      if((tau1_1fill == true) && (tau1_2fill == true) && (tau1_3fill == false) && (tau1_4fill == false) && (tau2_1fill == true) && (tau2_2fill == true) &&
                 (tau2_3fill == false) && (tau2_4fill == false)){
                Mayor_menorPt(Tau1_1cand_vec, Tau1_2cand_vec, Tau2_1cand_vec, Tau2_2cand_vec, Neu_Tau1_1cand_vec, Neu_Tau1_2cand_vec, Neu_Tau2_1cand_vec, Neu_Tau2_2cand_vec, 
                              &Tau1cand_vec, &Tau2cand_vec, &Tau3cand_vec, &Tau4cand_vec, &NeuTau1cand_vec, &NeuTau2cand_vec, &NeuTau3cand_vec, &NeuTau4cand_vec);
	      }
	      
	      if((tau1_1fill == true) && (tau1_2fill == false) && (tau1_3fill == false) && (tau1_4fill == false) && (tau2_1fill == true) && (tau2_2fill == true) &&
                 (tau2_3fill == true) && (tau2_4fill == false)){
                Mayor_menorPt(Tau1_1cand_vec, Tau2_1cand_vec, Tau2_2cand_vec, Tau2_3cand_vec, Neu_Tau1_1cand_vec, Neu_Tau2_1cand_vec, Neu_Tau2_2cand_vec, Neu_Tau2_3cand_vec, 
                              &Tau1cand_vec, &Tau2cand_vec, &Tau3cand_vec, &Tau4cand_vec, &NeuTau1cand_vec, &NeuTau2cand_vec, &NeuTau3cand_vec, &NeuTau4cand_vec);
	      }
	      
	      if((tau1_1fill == false) && (tau1_2fill == false) && (tau1_3fill == false) && (tau1_4fill == false) && (tau2_1fill == true) && (tau2_2fill == true) &&
                 (tau2_3fill == true) && (tau2_4fill == true)){
                Mayor_menorPt(Tau2_1cand_vec, Tau2_2cand_vec, Tau2_3cand_vec, Tau2_4cand_vec, Neu_Tau2_1cand_vec, Neu_Tau2_2cand_vec, Neu_Tau2_3cand_vec, Neu_Tau2_4cand_vec, 
                              &Tau1cand_vec, &Tau2cand_vec, &Tau3cand_vec, &Tau4cand_vec, &NeuTau1cand_vec, &NeuTau2cand_vec, &NeuTau3cand_vec, &NeuTau4cand_vec); 
	      }
              /*cout << "TLorentz que se han llenado "<<endl;
		cout <<"Tau1cand_vec.Pt "<<Tau1cand_vec.Pt()<<"  NeuTau1cand_vec.Pt "<<NeuTau1cand_vec.Pt()<<endl;
		cout <<"Tau2cand_vec.Pt "<<Tau2cand_vec.Pt()<<"  NeuTau2cand_vec.Pt "<<NeuTau2cand_vec.Pt()<<endl;
		cout <<"Tau3cand_vec.Pt "<<Tau3cand_vec.Pt()<<"  NeuTau3cand_vec.Pt "<<NeuTau3cand_vec.Pt()<<endl;
		cout <<"Tau4cand_vec.Pt "<<Tau4cand_vec.Pt()<<"  NeuTau4cand_vec.Pt "<<NeuTau4cand_vec.Pt()<<endl;
              */
	      bool NeuTau1_isHad = false;
	      bool NeuTau2_isHad = false;
	      
	      TausHadronicos(Tau1cand_vec, Tau2cand_vec, Tau3cand_vec, Tau4cand_vec, TauLep1cand_vec, TauLep2cand_vec, TauLep3cand_vec, TauLep4cand_vec, NeuTau1cand_vec, 
			     NeuTau2cand_vec, NeuTau3cand_vec, NeuTau4cand_vec, &Tau1Had_vec, &Tau2Had_vec, &Tau3Had_vec, &Tau4Had_vec, &NeuTau1Had_vec,  &NeuTau2Had_vec, 
			     &NeuTau3Had_vec, &NeuTau4Had_vec);
	      
              /*cout << "--------------Taus Hadronicos sin restar momento neutrino-----------"<<endl;
		cout <<"Tau1Had_vec.Pt "<<Tau1Had_vec.Pt()<<" NeuTau1Had_vec.Pt() "<<NeuTau1Had_vec.Pt()<<endl;
		cout <<"Tau2Had_vec.Pt "<<Tau2Had_vec.Pt()<<" NeuTau2Had_vec.Pt() "<<NeuTau2Had_vec.Pt()<<endl;
		cout <<"Tau3Had_vec.Pt "<<Tau3Had_vec.Pt()<<" NeuTau3Had_vec.Pt() "<<NeuTau3Had_vec.Pt()<<endl;
		cout <<"Tau4Had_vec.Pt "<<Tau4Had_vec.Pt()<<" NeuTau4Had_vec.Pt() "<<NeuTau4Had_vec.Pt()<<endl; 
		
		cout << "--------------Taus Hadronicos-----------"<<endl;
		cout <<"Tau1Had_vec.Pt "<<Tau1Had_vec.Pt()<<" Tau1Had_vec.Eta() "<<Tau1Had_vec.Eta()<<endl;
		cout <<"Tau2Had_vec.Pt "<<Tau2Had_vec.Pt()<<" Tau2Had_vec.Eta() "<<Tau2Had_vec.Eta()<<endl;
		cout <<"Tau3Had_vec.Pt "<<Tau3Had_vec.Pt()<<" Tau3Had_vec.Eta() "<<Tau3Had_vec.Eta()<<endl;
		cout <<"Tau4Had_vec.Pt "<<Tau4Had_vec.Pt()<<" Tau4Had_vec.Eta() "<<Tau4Had_vec.Eta()<<endl; 
              */
	      
              if(NeuTau1Had_vec.Pt()>0){
		Tau1Had_vec.SetPtEtaPhiE(abs(Tau1Had_vec.Pt()-NeuTau1Had_vec.Pt()), Tau1Had_vec.Eta(), Tau1Had_vec.Phi(), Tau1Had_vec.E());
	      }
	      
              if(NeuTau2Had_vec.Pt()>0){
		Tau2Had_vec.SetPtEtaPhiE(abs(Tau2Had_vec.Pt()-NeuTau2Had_vec.Pt()), Tau2Had_vec.Eta(), Tau2Had_vec.Phi(), Tau2Had_vec.E());
	      }
	      
	      if(NeuTau3Had_vec.Pt()>0){
		Tau3Had_vec.SetPtEtaPhiE(abs(Tau3Had_vec.Pt()-NeuTau3Had_vec.Pt()), Tau3Had_vec.Eta(), Tau3Had_vec.Phi(), Tau3Had_vec.E());
	      } 
	      
	      if(NeuTau4Had_vec.Pt()>0){
		Tau4Had_vec.SetPtEtaPhiE(abs(Tau4Had_vec.Pt()-NeuTau4Had_vec.Pt()), Tau4Had_vec.Eta(), Tau4Had_vec.Phi(), Tau4Had_vec.E());
	      }
              
	/*	cout << "--------------Taus Hadronicos una vez se resta el momento-----------"<<endl;
		cout <<"Tau1Had_vec.Pt "<<Tau1Had_vec.Pt()<<" Tau1Had_vec.Eta() "<<Tau1Had_vec.Eta()<<endl;
		cout <<"Tau2Had_vec.Pt "<<Tau2Had_vec.Pt()<<" Tau2Had_vec.Eta() "<<Tau2Had_vec.Eta()<<endl;
		cout <<"Tau3Had_vec.Pt "<<Tau3Had_vec.Pt()<<" Tau3Had_vec.Eta() "<<Tau3Had_vec.Eta()<<endl;
		cout <<"Tau4Had_vec.Pt "<<Tau4Had_vec.Pt()<<" Tau4Had_vec.Eta() "<<Tau4Had_vec.Eta()<<endl;              
          */    
              //Se deben organizar nuevamente los TLV de los taus ya que al restar el momento de los neutrinos puede cambiar el orden.
              Mayor_menorPt(Tau1Had_vec, Tau2Had_vec, Tau3Had_vec, Tau4Had_vec, NeuTau1Had_vec, NeuTau2Had_vec, NeuTau3Had_vec, NeuTau4Had_vec, &Tau1HadTLV, &Tau2HadTLV, 
                            &Tau3HadTLV, &Tau4HadTLV, &NeuTau1HadTLV,  &NeuTau2HadTLV,  &NeuTau3HadTLV,  &NeuTau4HadTLV);
		  /*cout << "--------------Taus Hadronicos una vez se resta el momento-----------"<<endl;
		  cout <<"Tau1HadTLV.Pt "<<Tau1HadTLV.Pt()<<" Tau1HadTLV.Eta() "<<Tau1HadTLV.Eta()<<endl;
		  cout <<"Tau2HadTLV.Pt "<<Tau2HadTLV.Pt()<<" Tau2HadTLV.Eta() "<<Tau2HadTLV.Eta()<<endl;
		  cout <<"Tau3HadTLV.Pt "<<Tau3HadTLV.Pt()<<" Tau3HadTLV.Eta() "<<Tau3HadTLV.Eta()<<endl;
		  cout <<"Tau4HadTLV.Pt "<<Tau4HadTLV.Pt()<<" Tau4HadTLV.Eta() "<<Tau4HadTLV.Eta()<<endl;
	          */

	      // Needed to generate a unique seed for the random number function used in the TauID function.
	      srand48(entry);
	      bool passed_tau1Cand =  TauID(Tau1HadTLV);
	      if (passed_tau1Cand) {
		Tau1HadTLV = TauSmearing(Tau1HadTLV);
		ntau_counter++;
	      } else {
		Tau1HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
	      }
	      srand48(entry+numberOfEntries);
	      bool passed_tau2Cand =  TauID(Tau2HadTLV);
	      if (passed_tau2Cand) { 
		Tau2HadTLV = TauSmearing(Tau2HadTLV);
		ntau_counter++;
	      } else {
		Tau2HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
	      }
              
              //cout << "Tau1HadTLV "<< Tau1HadTLV.Pt() << "  Tau2HadTLV "<< Tau2HadTLV.Pt() <<endl;
	      
	      if(Tau1HadTLV.Pt() == 0. && Tau2HadTLV.Pt() > 0.){
		Tau1HadTLV = Tau2HadTLV;
		Tau2HadTLV.SetPtEtaPhiE(0.,0.,0.,0.);
              }
	      
              //cout << "Tau1HadTLV "<< Tau1HadTLV.Pt() << "  Tau2HadTLV "<< Tau2HadTLV.Pt() <<endl;
	      
              /////Numeros de taus en cada evento/////
	      //cout << "Numero de taus "<<ntau_counter++<<endl;
	      //////////////////////////////////////
	      bool pass_lead_jet_cuts = false;
	      // Events with no cuts
	      
	      //pass_cuts[1] = 1;
	      
	      
	      // for events with exactly 1 Tau
	      if ((pass_cuts[0] == 1) && (ntau_counter == 1)){
		pass_cuts[1] = 1; 
	      }
              ///Exactly 2 Taus
              if ((pass_cuts[0] == 1) && (ntau_counter == 2)){
                pass_cuts[2] = 1;
              }  
              // for events where the tau(s) passes a pT cut
              if ((pass_cuts[1] == 1) && (Tau1HadTLV.Pt() > tau1_pt_min)){
		pass_cuts[3] = 1;
	      } 
              if ((pass_cuts[2] == 1) && (Tau1HadTLV.Pt() > tau1_pt_min) && (Tau2HadTLV.Pt() > tau2_pt_min)){
                pass_cuts[4] = 1;
              }
              // for events where the tau(s) passes an eta cut
	      if((pass_cuts[3] == 1) && (abs(Tau1HadTLV.Eta()) < tau_eta_max)){
		pass_cuts[5] = 1;
	      }
              if((pass_cuts[4] == 1) && (abs(Tau1HadTLV.Eta()) < tau_eta_max) && (abs(Tau2HadTLV.Eta()) < tau_eta_max)){
                pass_cuts[6] = 1;
              }
	      // events passing also b-jet veto requirement
	      if ( (pass_cuts[5]==1) && (!is_b_jet)){pass_cuts[7] = 1;}
	      if ( (pass_cuts[6]==1) && (!is_b_jet)){pass_cuts[8] = 1;}
              // events passing jet pt cut
	      if ((pass_cuts[7]==1) && (Jet_leading_vec.Pt() > lead_jet_pt)){pass_cuts[9] = 1;}
	      if ((pass_cuts[8]==1) && (Jet_leading_vec.Pt() > lead_jet_pt)){pass_cuts[10] = 1;}
              // events passing jet eta cut
	      if ((pass_cuts[9]==1) && (abs(Jet_leading_vec.Eta()) < jet_eta_max) ){pass_cuts[11] =1;}
              if ((pass_cuts[10]==1) && (abs(Jet_leading_vec.Eta()) < jet_eta_max) ){pass_cuts[12] =1;}
 	      // events passing dphitaujet cut 
	      double tau1_jet_dphi = abs(normalizedDphi(Jet_leading_vec.Phi() - Tau1HadTLV.Phi()));
	      if ((pass_cuts[11] == 1) && (tau1_jet_dphi > Dphitaujetmax)){pass_cuts[13] = 1;}
	      if ((pass_cuts[12] == 1) && (tau1_jet_dphi > Dphitaujetmax)){pass_cuts[14] = 1;}
              // events passing mt cut
	      double transmass = TMath::Sqrt(TMath::Abs(2*Tau1HadTLV.Pt()*MET*(1-TMath::Cos(TMath::Abs(Tau1HadTLV.Phi() - MET_phi)))));
	      if ((pass_cuts[11] == 1) && (transmass > transversemassmin)){pass_cuts[15] = 1;}
              if ((pass_cuts[12] == 1) && (transmass > transversemassmin)){pass_cuts[16] = 1;}
              // events passing also MET cut
	      if ((pass_cuts[11] == 1) && (MET > met_min)){pass_cuts[17] = 1;}
	      if ((pass_cuts[12] == 1) && (MET > met_min)){pass_cuts[18] = 1;}
              // events passing dphijetmet cut
	      double jet_met_dphi = abs(normalizedDphi(Jet_leading_vec.Phi() - MET_phi));
	      if ((pass_cuts[11] == 1) && (jet_met_dphi > deltajetmetphi)){pass_cuts[19]=1;}
              if ((pass_cuts[12] == 1) && (jet_met_dphi > deltajetmetphi)){pass_cuts[20]=1;}
	      // events passing effective mass cut
	      double efective_mass = TMath::Sqrt(Jet_leading_vec.Pt()*Jet_leading_vec.Pt()+Tau1HadTLV.Pt()*Tau1HadTLV.Pt()+MET*MET);
	      if ((pass_cuts[11] == 1) && (efective_mass > efectivemassmin)){pass_cuts[21] = 1;}        
	      if ((pass_cuts[12] == 1) && (efective_mass > efectivemassmin)){pass_cuts[22] = 1;}
	      
	      // events passing previous cuts with exactly 1 tau
	    }
	} //if(branchJet->GetEntries() > 0)
      
      // save some important histograms
      // _hmap_lead_jet_pT[0]->Fill(Jet_leading_vec.Pt());  
      
      for (int i = 0; i < nDir; i++){
	_hmap_Nevents[i]->Fill(0.0);
	_hmap_n_jets[i]->Fill(njets_counter);
	_hmap_n_tau[i]->Fill(ntau_counter);
	
	if ( pass_cuts[i] == 1){
	  if (Jet_leading_vec.Pt() > 5.0) {
	    double jet_met_dphi = abs(normalizedDphi(Jet_leading_vec.Phi() - MET_phi));
	    _hmap_jet_met_Dphi[i]->Fill(abs(jet_met_dphi));
	    _hmap_jet_met_metDphi[i]->Fill(abs(jet_met_dphi),MET);
	    _hmap_Nevents[i]->Fill(1.0);
	  }
	  
	  _hmap_met[i]->Fill(MET);
	  if(Jet_leading_vec.Pt() > 5.0){
	    _hmap_lead_jet_pT[i]->Fill(Jet_leading_vec.Pt());
	    _hmap_lead_jet_eta[i]->Fill(Jet_leading_vec.Eta());
	    _hmap_lead_jet_phi[i]->Fill(Jet_leading_vec.Phi());
	  }
	  
	  if(Tau1HadTLV.Pt() > 5.0){
	    _hmap_tau1_pT[i]->Fill(Tau1HadTLV.Pt());
	    _hmap_tau1_eta[i]->Fill(Tau1HadTLV.Eta());
	    _hmap_tau1_phi[i]->Fill(Tau1HadTLV.Phi());
	    double tau1_met_dphi = abs(normalizedDphi(Tau1HadTLV.Phi() - MET_phi));
	    _hmap_tau1_met_Dphi[i]->Fill(tau1_met_dphi);
	    
	    double transmass = TMath::Sqrt(TMath::Abs(2*Tau1HadTLV.Pt()*MET*(1-TMath::Cos(TMath::Abs(Tau1HadTLV.Phi() - MET_phi)))));
	    
	    _hmap_transverse_mass[i]->Fill(transmass);
	    double tau1_jet_Dphi = TMath::Abs(Jet_leading_vec.Phi() - Tau1HadTLV.Phi());
	    _hmap_tau1_jet_Dphi[i]->Fill(tau1_jet_Dphi);      
	    double efective_mass = TMath::Sqrt(Jet_leading_vec.Pt()*Jet_leading_vec.Pt()+Tau1HadTLV.Pt()*Tau1HadTLV.Pt()+MET*MET);
	    _hmap_efective_mass[i]->Fill(efective_mass);       
	    
	  }
	  if(Tau2HadTLV.Pt() > 5.0){
	    _hmap_tau2_pT[i]->Fill(Tau2HadTLV.Pt());
	    _hmap_tau2_eta[i]->Fill(Tau2HadTLV.Eta());
	    _hmap_tau2_phi[i]->Fill(Tau2HadTLV.Phi());
	    double tau2_met_dphi = abs(normalizedDphi(Tau2HadTLV.Phi() - MET_phi));
	    _hmap_tau2_met_Dphi[i]->Fill(tau2_met_dphi);
	  }
	}
      } 
    }
  
  theFile->cd();
  Numbertaus_D1->Write();  
  Numbertaus_D2->Write();
  Numbertaus_2->Write();
  Mothers->Write();
  for (int d = 0; d < nDir; d++)
    {
      cdDir[d]->cd();
      
      _hmap_Nevents[d]->Write();
      _hmap_n_jets[d]->Write();
      _hmap_n_tau[d]->Write();
      _hmap_lead_jet_pT[d]->Write();
      _hmap_lead_jet_eta[d]->Write();
      _hmap_lead_jet_phi[d]->Write();
      _hmap_tau1_pT[d]->Write();
      _hmap_tau1_eta[d]->Write();
      _hmap_tau1_phi[d]->Write();
      _hmap_tau2_pT[d]->Write();
      _hmap_tau2_eta[d]->Write();
      _hmap_tau2_phi[d]->Write();
      _hmap_jet_met_Dphi[d]->Write();
      _hmap_met[d]->Write();
      _hmap_tau1_met_Dphi[d]->Write();
      _hmap_tau2_met_Dphi[d]->Write();
      _hmap_jet_met_metDphi[d]->Write();
      _hmap_tau1_jet_Dphi[d]->Write();
      _hmap_transverse_mass[d]->Write();
      _hmap_efective_mass[d]->Write();
      
    }
  
  //   cdDir[0]->cd();
  theFile->Close();
  
}



PhenoAnalysis::~PhenoAnalysis()
{
  // do anything here that needs to be done at desctruction time
}

void PhenoAnalysis::Mayor_menorPt(TLorentzVector Tau1Cand, TLorentzVector Tau2Cand, TLorentzVector Tau3Cand, TLorentzVector Tau4Cand, TLorentzVector NeuTau1Cand, TLorentzVector NeuTau2Cand, TLorentzVector NeuTau3Cand, TLorentzVector NeuTau4Cand, TLorentzVector *Tau5Cand, TLorentzVector *Tau6Cand, TLorentzVector *Tau7Cand, TLorentzVector *Tau8Cand, TLorentzVector *NeuTau5Cand,  TLorentzVector *NeuTau6Cand, TLorentzVector *NeuTau7Cand, TLorentzVector *NeuTau8Cand) {
  
  if((Tau1Cand.Pt() > Tau2Cand.Pt()) && (Tau1Cand.Pt() > Tau3Cand.Pt()) && (Tau1Cand.Pt() > Tau4Cand.Pt())){
    if((Tau2Cand.Pt() > Tau3Cand.Pt()) && (Tau2Cand.Pt() > Tau4Cand.Pt())){
      if(Tau3Cand.Pt() > Tau4Cand.Pt()){
	*Tau5Cand = Tau1Cand;
	*Tau6Cand = Tau2Cand;
	*Tau7Cand = Tau3Cand;
	*Tau8Cand = Tau4Cand;
	*NeuTau5Cand = NeuTau1Cand;
	*NeuTau6Cand = NeuTau2Cand;
	*NeuTau7Cand = NeuTau3Cand;
	*NeuTau8Cand = NeuTau4Cand;
      }
      else{
	*Tau5Cand = Tau1Cand;
	*Tau6Cand = Tau2Cand;
	*Tau7Cand = Tau4Cand;
	*Tau8Cand = Tau3Cand;
	*NeuTau5Cand = NeuTau1Cand;
	*NeuTau6Cand = NeuTau2Cand;
	*NeuTau7Cand = NeuTau4Cand;
	*NeuTau8Cand = NeuTau3Cand;
      }
    }
    if((Tau3Cand.Pt() > Tau2Cand.Pt()) && (Tau3Cand.Pt() > Tau4Cand.Pt())){
      if(Tau2Cand.Pt() > Tau4Cand.Pt()){
	*Tau5Cand = Tau1Cand;
	*Tau6Cand = Tau3Cand;
	*Tau7Cand = Tau2Cand;
	*Tau8Cand = Tau4Cand;
	*NeuTau5Cand = NeuTau1Cand;
	*NeuTau6Cand = NeuTau3Cand;
	*NeuTau7Cand = NeuTau2Cand;
	*NeuTau8Cand = NeuTau4Cand;
      }
      else{
	*Tau5Cand = Tau1Cand;
	*Tau6Cand = Tau3Cand;
	*Tau7Cand = Tau4Cand;
	*Tau8Cand = Tau2Cand;
	*NeuTau5Cand = NeuTau1Cand;
	*NeuTau6Cand = NeuTau3Cand;
	*NeuTau7Cand = NeuTau4Cand;
	*NeuTau8Cand = NeuTau2Cand; 
      }
    }
    if((Tau4Cand.Pt() > Tau2Cand.Pt()) && (Tau4Cand.Pt() > Tau3Cand.Pt())){
      if(Tau2Cand.Pt() > Tau3Cand.Pt()){
	*Tau5Cand = Tau1Cand;
	*Tau6Cand = Tau4Cand;
	*Tau7Cand = Tau2Cand;
	*Tau8Cand = Tau3Cand;
	*NeuTau5Cand = NeuTau1Cand;
	*NeuTau6Cand = NeuTau4Cand;
	*NeuTau7Cand = NeuTau2Cand;
	*NeuTau8Cand = NeuTau3Cand;
      }
      else{
	*Tau5Cand = Tau1Cand;
	*Tau6Cand = Tau4Cand;
	*Tau7Cand = Tau4Cand;
	*Tau8Cand = Tau2Cand;
	*NeuTau5Cand = NeuTau1Cand;
	*NeuTau6Cand = NeuTau4Cand;
	*NeuTau7Cand = NeuTau4Cand;
	*NeuTau8Cand = NeuTau2Cand;
      }
    }
  }
  ///////////////////////Si Tau2 tiene pt mayor////////////////////////
  if((Tau2Cand.Pt() > Tau1Cand.Pt()) && (Tau2Cand.Pt() > Tau3Cand.Pt()) && (Tau2Cand.Pt() > Tau4Cand.Pt())){
    if((Tau1Cand.Pt() > Tau3Cand.Pt()) && (Tau1Cand.Pt() > Tau4Cand.Pt())){
      if(Tau3Cand.Pt() > Tau4Cand.Pt()){
	*Tau5Cand = Tau2Cand;
	*Tau6Cand = Tau1Cand;
	*Tau7Cand = Tau3Cand;
	*Tau8Cand = Tau4Cand;
	*NeuTau5Cand = NeuTau2Cand;
	*NeuTau6Cand = NeuTau1Cand;
	*NeuTau7Cand = NeuTau3Cand;
	*NeuTau8Cand = NeuTau4Cand;
        
      }
      else{
	*Tau5Cand = Tau2Cand;
	*Tau6Cand = Tau1Cand;
	*Tau7Cand = Tau4Cand;
	*Tau8Cand = Tau3Cand;
	*NeuTau5Cand = NeuTau2Cand;
	*NeuTau6Cand = NeuTau1Cand;
	*NeuTau7Cand = NeuTau4Cand;
	*NeuTau8Cand = NeuTau3Cand; 
      }
    }
    if((Tau3Cand.Pt() > Tau1Cand.Pt()) && (Tau3Cand.Pt() > Tau4Cand.Pt())){
      if(Tau1Cand.Pt() > Tau4Cand.Pt()){
	*Tau5Cand = Tau2Cand;
	*Tau6Cand = Tau3Cand;
	*Tau7Cand = Tau1Cand;
	*Tau8Cand = Tau4Cand;
	*NeuTau5Cand = NeuTau2Cand;
	*NeuTau6Cand = NeuTau3Cand;
	*NeuTau7Cand = NeuTau1Cand;
	*NeuTau8Cand = NeuTau4Cand;
      }
      else{
	*Tau5Cand = Tau2Cand;
	*Tau6Cand = Tau3Cand;
	*Tau7Cand = Tau4Cand;
	*Tau8Cand = Tau1Cand;
	*NeuTau5Cand = NeuTau2Cand;
	*NeuTau6Cand = NeuTau3Cand;
	*NeuTau7Cand = NeuTau4Cand;
	*NeuTau8Cand = NeuTau1Cand;
      }
    }
    if((Tau4Cand.Pt() > Tau1Cand.Pt()) && (Tau4Cand.Pt() > Tau3Cand.Pt())){
      if(Tau1Cand.Pt() > Tau3Cand.Pt()){
	*Tau5Cand = Tau2Cand;
	*Tau6Cand = Tau4Cand;
	*Tau7Cand = Tau1Cand;
	*Tau8Cand = Tau3Cand;
	*NeuTau5Cand = NeuTau2Cand;
	*NeuTau6Cand = NeuTau4Cand;
	*NeuTau7Cand = NeuTau1Cand;
	*NeuTau8Cand = NeuTau3Cand;
      }
      else{
	*Tau5Cand = Tau2Cand;
	*Tau6Cand = Tau4Cand;
	*Tau7Cand = Tau3Cand;
	*Tau8Cand = Tau1Cand;
	*NeuTau5Cand = NeuTau2Cand;
	*NeuTau6Cand = NeuTau4Cand;
	*NeuTau7Cand = NeuTau3Cand;
	*NeuTau8Cand = NeuTau1Cand; 
      }
    }
  }
  ///////////////////////Si Tau3 tiene pt mayor////////////////////////
  if((Tau3Cand.Pt() > Tau1Cand.Pt()) && (Tau3Cand.Pt() > Tau2Cand.Pt()) && (Tau3Cand.Pt() > Tau4Cand.Pt())){
    if((Tau1Cand.Pt() > Tau2Cand.Pt()) && (Tau1Cand.Pt() > Tau4Cand.Pt())){
      if(Tau2Cand.Pt() > Tau4Cand.Pt()){
	*Tau5Cand = Tau3Cand;
	*Tau6Cand = Tau1Cand;
	*Tau7Cand = Tau2Cand;
	*Tau8Cand = Tau4Cand;
	*NeuTau5Cand = NeuTau3Cand;
	*NeuTau6Cand = NeuTau1Cand;
	*NeuTau7Cand = NeuTau2Cand;
	*NeuTau8Cand = NeuTau4Cand; 
      }
      else{
	*Tau5Cand = Tau3Cand;
	*Tau6Cand = Tau1Cand;
	*Tau7Cand = Tau4Cand;
	*Tau8Cand = Tau2Cand;
	*NeuTau5Cand = NeuTau3Cand;
	*NeuTau6Cand = NeuTau1Cand;
	*NeuTau7Cand = NeuTau4Cand;
	*NeuTau8Cand = NeuTau2Cand; 
      }
    }
    if((Tau2Cand.Pt() > Tau1Cand.Pt()) && (Tau2Cand.Pt() > Tau4Cand.Pt())){
      if(Tau1Cand.Pt() > Tau4Cand.Pt()){
	*Tau5Cand = Tau3Cand;
	*Tau6Cand = Tau2Cand;
	*Tau7Cand = Tau1Cand;
	*Tau8Cand = Tau4Cand;
	*NeuTau5Cand = NeuTau3Cand;
	*NeuTau6Cand = NeuTau2Cand;
	*NeuTau7Cand = NeuTau1Cand;
	*NeuTau8Cand = NeuTau4Cand; 
      }
      else{
	*Tau5Cand = Tau3Cand;
	*Tau6Cand = Tau2Cand;
	*Tau7Cand = Tau4Cand;
	*Tau8Cand = Tau1Cand;
	*NeuTau5Cand = NeuTau3Cand;
	*NeuTau6Cand = NeuTau2Cand;
	*NeuTau7Cand = NeuTau4Cand;
	*NeuTau8Cand = NeuTau1Cand; 
      }
    }
    if((Tau4Cand.Pt() > Tau1Cand.Pt()) && (Tau4Cand.Pt() > Tau2Cand.Pt())){
      if(Tau1Cand.Pt() > Tau2Cand.Pt()){
	*Tau5Cand = Tau3Cand;
	*Tau6Cand = Tau4Cand;
	*Tau7Cand = Tau1Cand;
	*Tau8Cand = Tau2Cand;
	*NeuTau5Cand = NeuTau3Cand;
	*NeuTau6Cand = NeuTau4Cand;
	*NeuTau7Cand = NeuTau1Cand;
	*NeuTau8Cand = NeuTau2Cand;
      }
      else{
	*Tau5Cand = Tau3Cand;
	*Tau6Cand = Tau4Cand;
	*Tau7Cand = Tau2Cand;
	*Tau8Cand = Tau1Cand;
	*NeuTau5Cand = NeuTau3Cand;
	*NeuTau6Cand = NeuTau4Cand;
	*NeuTau7Cand = NeuTau2Cand;
	*NeuTau8Cand = NeuTau1Cand;
      }
    }
  }
  ///////////////////////Si Tau4 tiene pt mayor////////////////////////
  if((Tau4Cand.Pt() > Tau1Cand.Pt()) && (Tau4Cand.Pt() > Tau2Cand.Pt()) && (Tau4Cand.Pt() > Tau3Cand.Pt())){
    if((Tau1Cand.Pt() > Tau2Cand.Pt()) && (Tau1Cand.Pt() > Tau3Cand.Pt())){
      if(Tau2Cand.Pt() > Tau3Cand.Pt()){
	*Tau5Cand = Tau4Cand;
	*Tau6Cand = Tau1Cand;
	*Tau7Cand = Tau2Cand;
	*Tau8Cand = Tau3Cand;
	*NeuTau5Cand = NeuTau4Cand;
	*NeuTau6Cand = NeuTau1Cand;
	*NeuTau7Cand = NeuTau2Cand;
	*NeuTau8Cand = NeuTau3Cand;
      }
      else{
	*Tau5Cand = Tau4Cand;
	*Tau6Cand = Tau1Cand;
	*Tau7Cand = Tau3Cand;
	*Tau8Cand = Tau2Cand;
	*NeuTau5Cand = NeuTau4Cand;
	*NeuTau6Cand = NeuTau1Cand;
	*NeuTau7Cand = NeuTau3Cand;
	*NeuTau8Cand = NeuTau2Cand;  
      }
    }
    if((Tau2Cand.Pt() > Tau1Cand.Pt()) && (Tau2Cand.Pt() > Tau3Cand.Pt())){
      if(Tau1Cand.Pt() > Tau3Cand.Pt()){
	*Tau5Cand = Tau4Cand;
	*Tau6Cand = Tau2Cand;
	*Tau7Cand = Tau1Cand;
	*Tau8Cand = Tau3Cand;
	*NeuTau5Cand = NeuTau4Cand;
	*NeuTau6Cand = NeuTau2Cand;
	*NeuTau7Cand = NeuTau1Cand;
	*NeuTau8Cand = NeuTau3Cand; 
      }
      else{
	*Tau5Cand = Tau4Cand;
	*Tau6Cand = Tau2Cand;
	*Tau7Cand = Tau3Cand;
	*Tau8Cand = Tau1Cand;
	*NeuTau5Cand = NeuTau4Cand;
	*NeuTau6Cand = NeuTau2Cand;
	*NeuTau7Cand = NeuTau3Cand;
	*NeuTau8Cand = NeuTau1Cand; 
      }
    }
    if((Tau3Cand.Pt() > Tau1Cand.Pt()) && (Tau3Cand.Pt() > Tau2Cand.Pt())){
      if(Tau1Cand.Pt() > Tau2Cand.Pt()){
	*Tau5Cand = Tau4Cand;
	*Tau6Cand = Tau3Cand;
	*Tau7Cand = Tau1Cand;
	*Tau8Cand = Tau2Cand;
	*NeuTau5Cand = NeuTau4Cand;
	*NeuTau6Cand = NeuTau3Cand;
	*NeuTau7Cand = NeuTau1Cand;
	*NeuTau8Cand = NeuTau2Cand; 
      }
      else{
	*Tau5Cand = Tau4Cand;
	*Tau6Cand = Tau3Cand;
	*Tau7Cand = Tau2Cand;
	*Tau8Cand = Tau1Cand;
	*NeuTau5Cand = NeuTau4Cand;
	*NeuTau6Cand = NeuTau3Cand;
	*NeuTau7Cand = NeuTau2Cand;
	*NeuTau8Cand = NeuTau1Cand; 
      }
    }
  }
}

void PhenoAnalysis::TausHadronicos(TLorentzVector Tau1cand, TLorentzVector Tau2cand, TLorentzVector Tau3cand, TLorentzVector Tau4cand, TLorentzVector TauLep1, TLorentzVector TauLep2, TLorentzVector TauLep3, TLorentzVector TauLep4, TLorentzVector NeuTau1vec, TLorentzVector NeuTau2vec, TLorentzVector NeuTau3vec, TLorentzVector NeuTau4vec, TLorentzVector *Tau1Had, TLorentzVector *Tau2Had, TLorentzVector *Tau3Had, TLorentzVector *Tau4Had, TLorentzVector *NeuTau1Had,  TLorentzVector *NeuTau2Had, TLorentzVector *NeuTau3Had, TLorentzVector *NeuTau4Had) {
  bool Tau1cand_isLep = false;
  bool Tau2cand_isLep = false;
  bool Tau3cand_isLep = false;
  bool Tau4cand_isLep = false;
  
  if(TauLep1.Pt() > 5.){
    if (Tau1cand.Eta() == TauLep1.Eta()){Tau1cand_isLep = true;}
    if (Tau2cand.Eta() == TauLep1.Eta()){Tau2cand_isLep = true;}
    if (Tau3cand.Eta() == TauLep1.Eta()){Tau3cand_isLep = true;}
    if (Tau4cand.Eta() == TauLep1.Eta()){Tau4cand_isLep = true;}
  }
  
  if (TauLep2.Pt() > 5.){
    if (Tau1cand.Eta() == TauLep2.Eta()){Tau1cand_isLep = true;}
    if (Tau2cand.Eta() == TauLep2.Eta()){Tau2cand_isLep = true;}
    if (Tau3cand.Eta() == TauLep2.Eta()){Tau3cand_isLep = true;}
    if (Tau4cand.Eta() == TauLep2.Eta()){Tau4cand_isLep = true;}
  }
  
  if (TauLep3.Pt() > 5.){
    if (Tau1cand.Eta() == TauLep3.Eta()){Tau1cand_isLep = true;}
    if (Tau2cand.Eta() == TauLep3.Eta()){Tau2cand_isLep = true;}
    if (Tau3cand.Eta() == TauLep3.Eta()){Tau3cand_isLep = true;}
    if (Tau4cand.Eta() == TauLep3.Eta()){Tau4cand_isLep = true;}
  }  	       
  
  if (TauLep4.Pt() > 5.){
    if (Tau1cand.Eta() == TauLep4.Eta()){Tau1cand_isLep = true;}
    if (Tau2cand.Eta() == TauLep4.Eta()){Tau2cand_isLep = true;}
    if (Tau3cand.Eta() == TauLep4.Eta()){Tau3cand_isLep = true;}
    if (Tau4cand.Eta() == TauLep4.Eta()){Tau4cand_isLep = true;}
  }
  
  TLorentzVector Tau_nulo(0., 0., 0., 0.);
  
  if((Tau1cand_isLep == false) && (Tau2cand_isLep == true) && (Tau3cand_isLep == true) && (Tau4cand_isLep == true)){
    if (Tau1cand.Pt() > Tau1Had->Pt()){*Tau1Had = Tau1cand; *NeuTau1Had = NeuTau1vec;}
  }
  if((Tau1cand_isLep == true) && (Tau2cand_isLep == false) && (Tau3cand_isLep == true) && (Tau4cand_isLep == true)){
    if (Tau2cand.Pt() > Tau1Had->Pt()){*Tau1Had = Tau2cand; *NeuTau1Had = NeuTau2vec;}
  }
  if((Tau1cand_isLep == true) && (Tau2cand_isLep == true) && (Tau3cand_isLep == false) && (Tau4cand_isLep == true)){
    if (Tau3cand.Pt() > Tau1Had->Pt()){*Tau1Had = Tau3cand; *NeuTau1Had = NeuTau3vec;} 
  }
  if((Tau1cand_isLep == true) && (Tau2cand_isLep == true) && (Tau3cand_isLep == true) && (Tau4cand_isLep == false)){
    if (Tau4cand.Pt() > Tau1Had->Pt()){*Tau1Had = Tau4cand; *NeuTau1Had = NeuTau4vec;}
  }   
  ////////////////////DOS TAUS HADRONICOS ///////////////////               
  if((Tau1cand_isLep == false) && (Tau2cand_isLep == false) && (Tau3cand_isLep == true) && (Tau4cand_isLep == true)){
    if (Tau1cand.Pt() > Tau1Had->Pt()){*Tau1Had = Tau1cand; *NeuTau1Had = NeuTau1vec;}
    if (Tau2cand.Pt() > Tau2Had->Pt()){*Tau2Had = Tau2cand; *NeuTau2Had = NeuTau2vec;}
  }
  if((Tau1cand_isLep == false) && (Tau2cand_isLep == true) && (Tau3cand_isLep == false) && (Tau4cand_isLep == true)){
    if (Tau1cand.Pt() > Tau1Had->Pt()){*Tau1Had = Tau1cand; *NeuTau1Had = NeuTau1vec;}
    if (Tau3cand.Pt() > Tau2Had->Pt()){*Tau2Had = Tau3cand; *NeuTau2Had = NeuTau3vec;}
  }  
  if((Tau1cand_isLep == false) && (Tau2cand_isLep == true) && (Tau3cand_isLep == true) && (Tau4cand_isLep == false)){
    if (Tau1cand.Pt() > Tau1Had->Pt()){*Tau1Had = Tau1cand; *NeuTau1Had = NeuTau1vec;}
    if (Tau4cand.Pt() > Tau2Had->Pt()){*Tau2Had = Tau4cand; *NeuTau2Had = NeuTau4vec;}
  }
  if((Tau1cand_isLep == true) && (Tau2cand_isLep == false) && (Tau3cand_isLep == false) && (Tau4cand_isLep == true)){
    if (Tau2cand.Pt() > Tau1Had->Pt()){*Tau1Had = Tau2cand; *NeuTau1Had = NeuTau2vec;}
    if (Tau3cand.Pt() > Tau2Had->Pt()){*Tau2Had = Tau3cand; *NeuTau2Had = NeuTau3vec;}
  }
  if((Tau1cand_isLep == true) && (Tau2cand_isLep == false) && (Tau3cand_isLep == true) && (Tau4cand_isLep == false)){
    if (Tau2cand.Pt() > Tau1Had->Pt()){*Tau1Had = Tau2cand; *NeuTau1Had = NeuTau2vec;}
    if (Tau4cand.Pt() > Tau2Had->Pt()){*Tau2Had = Tau4cand; *NeuTau2Had = NeuTau4vec;}
  }
  if((Tau1cand_isLep == true) && (Tau2cand_isLep == true) && (Tau3cand_isLep == false) && (Tau4cand_isLep == false)){
    if (Tau3cand.Pt() > Tau1Had->Pt()){*Tau1Had = Tau3cand; *NeuTau1Had = NeuTau3vec;}
    if (Tau4cand.Pt() > Tau2Had->Pt()){*Tau2Had = Tau4cand; *NeuTau2Had = NeuTau4vec;}
  }    
  ///////////////////TRES TAUS HADRONICOS///////////////////////////////
  if((Tau1cand_isLep == false) && (Tau2cand_isLep == false) && (Tau3cand_isLep == false) && (Tau4cand_isLep == true)){
    if (Tau1cand.Pt() > Tau1Had->Pt()){*Tau1Had = Tau1cand; *NeuTau1Had = NeuTau1vec;}
    if (Tau2cand.Pt() > Tau2Had->Pt()){*Tau2Had = Tau2cand; *NeuTau2Had = NeuTau2vec;}
    if (Tau3cand.Pt() > Tau3Had->Pt()){*Tau3Had = Tau3cand; *NeuTau3Had = NeuTau3vec;}
  }
  if((Tau1cand_isLep == false) && (Tau2cand_isLep == false) && (Tau3cand_isLep == true) && (Tau4cand_isLep == false)){
    if (Tau1cand.Pt() > Tau1Had->Pt()){*Tau1Had = Tau1cand; *NeuTau1Had = NeuTau1vec;}
    if (Tau2cand.Pt() > Tau2Had->Pt()){*Tau2Had = Tau2cand; *NeuTau2Had = NeuTau2vec;}
    if (Tau4cand.Pt() > Tau3Had->Pt()){*Tau3Had = Tau4cand; *NeuTau3Had = NeuTau4vec;}
  }
  if((Tau1cand_isLep == false) && (Tau2cand_isLep == true) && (Tau3cand_isLep == false) && (Tau4cand_isLep == false)){
    if (Tau1cand.Pt() > Tau1Had->Pt()){*Tau1Had = Tau1cand; *NeuTau1Had = NeuTau1vec;}
    if (Tau3cand.Pt() > Tau2Had->Pt()){*Tau2Had = Tau3cand; *NeuTau2Had = NeuTau3vec;}
    if (Tau4cand.Pt() > Tau3Had->Pt()){*Tau3Had = Tau4cand; *NeuTau3Had = NeuTau4vec;}
  }
  if((Tau1cand_isLep == true) && (Tau2cand_isLep == false) && (Tau3cand_isLep == false) && (Tau4cand_isLep == false)){
    if (Tau2cand.Pt() > Tau1Had->Pt()){*Tau1Had = Tau2cand; *NeuTau1Had = NeuTau2vec;}
    if (Tau3cand.Pt() > Tau2Had->Pt()){*Tau2Had = Tau3cand; *NeuTau2Had = NeuTau3vec;}
    if (Tau4cand.Pt() > Tau3Had->Pt()){*Tau3Had = Tau4cand; *NeuTau3Had = NeuTau4vec;}
    
  }
  if((Tau1cand_isLep == false) && (Tau2cand_isLep == false) && (Tau3cand_isLep == false) && (Tau4cand_isLep == false)){
    if (Tau1cand.Pt() > Tau1Had->Pt()){*Tau1Had = Tau1cand; *NeuTau1Had = NeuTau1vec;}
    if (Tau2cand.Pt() > Tau2Had->Pt()){*Tau2Had = Tau2cand; *NeuTau2Had = NeuTau2vec;}
    if (Tau3cand.Pt() > Tau3Had->Pt()){*Tau3Had = Tau3cand; *NeuTau3Had = NeuTau3vec;}
    if (Tau4cand.Pt() > Tau4Had->Pt()){*Tau4Had = Tau4cand; *NeuTau4Had = NeuTau4vec;}
  }
  //////////////////////////////////////////////////////////////////////
  
}


bool PhenoAnalysis::TauID(TLorentzVector TauCand) {
  
  bool passedTauID = false;
  if (TauCand.Pt() > 0.){
    double r = drand48();
    if (r <= 0.8){passedTauID = true;}
  }
  return passedTauID;
}

bool PhenoAnalysis::TauIDJet(TLorentzVector jet){
  bool passedTauIDJet = false;
  if(jet.Pt() > 0.){
  double x = drand48();
  if(x < 0.03){passedTauIDJet = true;}
  }
  return passedTauIDJet;
}

TLorentzVector PhenoAnalysis::TauSmearing(TLorentzVector TauCand) {
  
  if (TauCand.Pt() > 5.){
    double smeared_pt  = 0.03*TauCand.Pt() + TauCand.Pt();
    double smeared_eta = TauCand.Eta();
    double smeared_phi = TauCand.Phi();
    double smeared_e   = TauCand.E(); 
    TauCand.SetPtEtaPhiE(smeared_pt, smeared_eta, smeared_phi, smeared_e);
  }
  return TauCand;
}


double PhenoAnalysis::calculateE(double eta, double pt, double mass){
  
  double theta = TMath::ATan(TMath::Exp(-eta)); 
  double cos_theta = TMath::Cos(2*theta);
  double p= pt/cos_theta;
  double e = sqrt(pow(p, 2) + pow(mass, 2));
  
  return e;
  
}

double PhenoAnalysis::normalizedDphi(double phi){
  const double PI  = 3.141592653589793238463;
  double twoPI = 2.0*PI;
  if ( phi < -PI ){phi += twoPI;}
  if ( phi > PI ){phi -= twoPI;}
  return phi;
}
void PhenoAnalysis::crateHistoMasps (int directories)
{
  for (int i = 0; i < directories; i++)
    {
      _hmap_Nevents[i]       = new TH1F("Nevents", "Nevents", 3,0.,3); 
      _hmap_lead_jet_pT[i]   = new TH1F("jet_lead_pT",    "Pt leading jet", 100, 0., 1000.);
      _hmap_lead_jet_eta[i]  = new TH1F("jet_lead_eta",   "#eta jet", 50, -5.0, 5.0);
      _hmap_lead_jet_phi[i]  = new TH1F("jet_lead_phi",   "#phi jet", 70, -3.6, 3.6); 
      _hmap_n_jets[i]        = new TH1F("N_jets",         "N(jet)", 4, 0., 4);
      _hmap_n_tau[i]         = new TH1F("N_tau",          "N(tau)", 4, 0., 4);
      _hmap_tau1_pT[i]       = new TH1F("tau1_pT",        "p_{T}(#tau_{1})", 60, 0., 600.);
      _hmap_tau1_eta[i]      = new TH1F("tau1_eta",       "#eta(#tau_{1})", 50, -3.5, 3.5);
      _hmap_tau1_phi[i]      = new TH1F("tau1_phi",       "#phi(#tau_{1})", 70, -3.6, 3.6);
      _hmap_tau2_pT[i]       = new TH1F("tau2_pT",        "p_{T}(#tau_{2})", 75, 0., 300.);
      _hmap_tau2_eta[i]      = new TH1F("tau2_eta",       "#eta(#tau_{2})", 50, -3.5, 3.5);
      _hmap_tau2_phi[i]      = new TH1F("tau2_phi",       "#phi(#tau_{2})", 70, -3.6, 3.6);
      _hmap_jet_met_Dphi[i]  = new TH1F("jet_met_Dphi",   "#Delta #phi(jet, MET)", 32, 0, 3.2);
      _hmap_transverse_mass[i] = new TH1F("transverse_mass",   "Transverse Mass", 60, 0, 600);
      _hmap_met[i]           = new TH1F("met",            "Met", 120, 0., 1200.); 
      _hmap_tau1_met_Dphi[i]  = new TH1F("tau_met_Dphi",   "#Delta #phi(#tau_{1}, MET)", 32, 0, 3.2);
      _hmap_tau2_met_Dphi[i]  = new TH1F("tau_met_Dphi",   "#Delta #phi(#tau_{2}, MET)", 32, 0, 3.2);
      _hmap_tau1_jet_Dphi[i]  = new TH1F("tau1_jet_Dphi",   "#Delta #phi(#tau_{1}, jet)", 32, 0, 3.2);
      _hmap_efective_mass[i]  = new TH1F("efective_mass", "Efective_mass", 80, 0, 800); 
      _hmap_jet_met_metDphi[i] = new TH2F("jet_met_metDphi", "#Delta #phi(jet,MET), MET", 32, 0, 3.2, 120, 0, 1200);
    }
}
