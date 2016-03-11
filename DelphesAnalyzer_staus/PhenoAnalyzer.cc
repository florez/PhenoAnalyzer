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
  int nDir = 12;
  TDirectory *theDirectory[nDir+1];
  theDirectory[0]  = HistoOutputFile->mkdir("No_cuts");
  theDirectory[1]  = HistoOutputFile->mkdir("At_least_one_Tau");
  theDirectory[2]  = HistoOutputFile->mkdir("After_Tau_pt_min");
  theDirectory[3]  = HistoOutputFile->mkdir("After_Tau_Eta_min");
  theDirectory[4]  = HistoOutputFile->mkdir("After_b_jet_veto");
  theDirectory[5]  = HistoOutputFile->mkdir("After_jet_pt_cut");
  theDirectory[6]  = HistoOutputFile->mkdir("After_jet_Eta_cut");
  theDirectory[7]  = HistoOutputFile->mkdir("Dphi_tau_jet_cut");
  theDirectory[8]  = HistoOutputFile->mkdir("Transverse_mass_cut");
  theDirectory[9]  = HistoOutputFile->mkdir("After_MET");
  theDirectory[10]  = HistoOutputFile->mkdir("Dphi_jet_met");   
  theDirectory[11]  = HistoOutputFile->mkdir("Efective_mass_cut");
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
  TH1 *histJetPT         = new TH1F("jet_PT",         "j p_{T}", 3000, 0., 1000.);
  TH1 *histJetEta        = new TH1F("jet_Eta",        "j #eta", 100, -5.0, 5.0);
  TH1 *histJetPhi        = new TH1F("Jet_Phi",        "j #phi", 72, -3.6, 3.6);
  TH1 *Numbertaus_D1     = new TH1F("Numbertaus_D1",  "Numbertaus_D1", 30, -5.5, 24.5); 		 
  TH1 *Numbertaus_D2     = new TH1F("Numbertaus_D2",  "Numbertaus_D2", 10, 0, 10);
  TH1 *Numbertaus_2      = new TH1F("Numbertaus_2",   "Numbertaus_D2", 30, -5.5, 24.5);
  TH1 *Mothers           = new TH1F("Mothers",        "Mothers", 510, -10, 500);
  
  int counter_elec_muon  = 0;

  
  //for(Int_t entry = 0; entry < numberOfEntries; ++entry)
    for(Int_t entry = 0; entry < 100; ++entry)
    {
      int pass_cuts[nDir] = {0};
       cout << " "<<endl;
       cout << "=================  Event: " << entry <<"===================="<<endl;
       cout << "PID\t"<<"M1\t"<<"D1\t"<<"D2\t"<<endl;
       cout << " "<<endl;

      TLorentzVector Jet_leading_vec(0., 0., 0., 0.);
      TLorentzVector electron (0., 0., 0., 0.);
      TLorentzVector Tau1cand_vec (0., 0., 0., 0.);
      TLorentzVector Tau2cand_vec (0., 0., 0., 0.);
      TLorentzVector TauLep1cand_vec (0., 0., 0., 0.);
      TLorentzVector TauLep2cand_vec (0., 0., 0., 0.);
      TLorentzVector NeuTau1cand_vec (0., 0., 0., 0.);
      TLorentzVector NeuTau2cand_vec (0., 0., 0., 0.);
      TLorentzVector lepton_vec (0., 0., 0., 0.);

      bool is_b_jet = false;
      treeReader->ReadEntry(entry);
      int tau_conter=0;
      METpointer = (MissingET*) branchMissingET->At(0);
      double MET = METpointer->MET;
      double MET_phi = METpointer->Phi;
      int njets_counter = 0;
      int ntau_counter = 0; 
      pass_cuts[0]=1;
      int countertau=0;


      if(branchJet->GetEntries() > 0)
	{
	  // Take first jet
	  Jet *firstjet = (Jet*) branchJet->At(0);
	  // Plot jet transverse momentum
	  histJetPT->Fill(firstjet->PT);
	  histJetEta->Fill(firstjet->Eta);
	  histJetPhi->Fill(firstjet->Phi);
	  // For Jets
	  double jet_min_pt = 20.;
	  
          // We need at least 2 jets, 1 (2) tau jets and the ISR jet.
	  
	  if (branchJet->GetEntriesFast() > 0) 
	    {
	      TLorentzVector jet_i(0., 0., 0., 0.);
	      TLorentzVector elec_i(0., 0., 0., 0.);
	      
	      double N_muons = 0;
	      double N_elec  = 0; 
	      for (int muo = 0; muo < branchMuon->GetEntriesFast(); muo++){
		Muon *muon = (Muon*) branchMuon->At(muo);
		if ((muon->PT > 10.) && (abs(muon->Eta) < 2.5)){N_muons++;}            
	      }
	      
	      for (int j = 0; j < branchJet->GetEntriesFast(); j++)
		{
		  bool is_jet_elec_overlap = false;
		  Jet *jet = (Jet*) branchJet->At(j);
                  //cout << "BTag: "<<jet->BTag<<endl;
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
		  
		  //if ((!is_jet_elec_overlap) && (N_muons == 0) && (N_elec == 0)){
		  if ((jet->PT > b_jet_pt_min) && (jet->BTag == 1)){is_b_jet = true;}
		  // if the jet is tagged as a tau jet, then save the ID 
		  // to use it later
                  
		  else if (jet->PT > jet_min_pt){
		    njets_counter++;
		    jet_min_pt = jet->PT; 
		    Jet_leading_vec.SetPtEtaPhiE(jet->PT, jet->Eta, jet->Phi, jet_energy);
		  }
		  // }
		}
	      // check if there is at least one jet tagged as a tau

              bool found_first_lep = false;
              bool found_first_neu = false;

	      for(int tau_i = 0; tau_i < branchGenParticle->GetEntriesFast(); tau_i++){
		GenParticle *tau_mother = (GenParticle*) branchGenParticle->At(tau_i);
		if( (abs(tau_mother->PID) == 23) && (tau_mother->Status == 2) ){
		  GenParticle *daughter_1 = (GenParticle*) branchGenParticle->At(tau_mother->D1);
		  GenParticle *daughter_2 = (GenParticle*) branchGenParticle->At(tau_mother->D2);
		  GenParticle *mother_1 = (GenParticle*) branchGenParticle->At(tau_mother->M1);
		  // Find tau candidates                  
		  if ( (abs(daughter_1->PID) == 15) ) {
		    double tau1cand_energy = calculateE(daughter_1->Eta, daughter_1->PT, daughter_1->Mass);
		    Tau1cand_vec.SetPtEtaPhiE(daughter_1->PT, daughter_1->Eta, daughter_1->Phi, tau1cand_energy);
		    cout<<tau_mother->PID<<"\t"<<mother_1->PID<<"\t"<<daughter_1->PID<<"\t"<<daughter_2->PID<<endl; 
                    
		  }
		  if ( (abs(daughter_2->PID) == 15) ) {
		    double tau2cand_energy = calculateE(daughter_2->Eta, daughter_2->PT, daughter_2->Mass);
		    Tau2cand_vec.SetPtEtaPhiE(daughter_2->PT, daughter_2->Eta, daughter_2->Phi, tau2cand_energy);
		  }
		}
		if (((abs(tau_mother->PID) == 11) || (abs(tau_mother->PID) == 13)) && (tau_mother->Status == 1)){
                  double lepton_energy = calculateE(tau_mother->Eta, tau_mother->PT, tau_mother->Mass);
		  lepton_vec.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, lepton_energy);
                  GenParticle *mother_elec = (GenParticle*) branchGenParticle->At(tau_mother->M1);  
		  if (abs(mother_elec->PID) == 15){
		    double MotherElcand_energy = calculateE(mother_elec->Eta, mother_elec->PT, mother_elec->Mass);
		    if (found_first_lep == false ) { 
		      TauLep1cand_vec.SetPtEtaPhiE(mother_elec->PT, mother_elec->Eta, mother_elec->Phi, MotherElcand_energy); 
		      found_first_lep == true;    
		    } else { 
		      TauLep2cand_vec.SetPtEtaPhiE(mother_elec->PT, mother_elec->Eta, mother_elec->Phi, MotherElcand_energy);
		    }
		  }
		}

                if ((abs(tau_mother->PID) == 16) &&  (tau_mother->Status == 1)){
                   double TauNeu_energy = calculateE(tau_mother->Eta, tau_mother->PT, tau_mother->Mass);
                   if (found_first_neu == false){
                     NeuTau1cand_vec.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, TauNeu_energy);
                     found_first_neu = true;
                   } else {
                     NeuTau2cand_vec.SetPtEtaPhiE(tau_mother->PT, tau_mother->Eta, tau_mother->Phi, TauNeu_energy);
                   }
                }

	      }
	      bool Tau1cand_isLep = false;
              bool Tau2cand_isLep = false;
 
              if (TauLep1cand_vec.Pt() > 5.){
                if (Tau1cand_vec.Eta() == TauLep1cand_vec.Eta()){Tau1cand_isLep = true;}
                if (Tau2cand_vec.Eta() == TauLep1cand_vec.Eta()){Tau2cand_isLep = true;}
              }

              if (TauLep2cand_vec.Pt() > 5.){
                if (Tau1cand_vec.Eta() == TauLep2cand_vec.Eta()){Tau1cand_isLep = true;}
                if (Tau2cand_vec.Eta() == TauLep2cand_vec.Eta()){Tau2cand_isLep = true;}
              }
              if ( (Tau1cand_isLep == false) && (Tau2cand_isLep == true) && (Tau1cand_vec.Pt() > 0.)){
                cout << "Tau1Candi is Had"<<endl;
                cout <<"TauLep2cand_vec "<< TauLep2cand_vec.Pt() <<"  lepton Pt+ NeuTau1cand Pt  "<< lepton_vec.Pt()+ NeuTau1cand_vec.Pt()<<endl;
                cout <<"TauLep2cand_vec "<< TauLep2cand_vec.Pt() <<"  lepton Pt+ NeuTau1cand Pt  "<< lepton_vec.Pt()+ NeuTau1cand_vec.Pt()<<endl;
              }
              if ( (Tau1cand_isLep == true) && (Tau2cand_isLep == false) && (Tau2cand_vec.Pt() > 0.)){
                cout << "Tau2Candi is Had"<<endl;
                Tau1cand_vec.SetPtEtaPhiE(Tau2cand_vec.Pt(), Tau2cand_vec.Eta(), Tau2cand_vec.Phi(), Tau2cand_vec.E());
              }
              if ( (Tau1cand_isLep == true) && (Tau2cand_isLep == true)){
                cout << "No Tau Had"<<endl;
                Tau1cand_vec.SetPtEtaPhiE(0.,0.,0.,0.);
                Tau2cand_vec.SetPtEtaPhiE(0.,0.,0.,0.);
              }
              if ( (Tau1cand_isLep == false) && (Tau2cand_isLep == false) && (Tau1cand_vec.Pt() > 0.) && (Tau2cand_vec.Pt() > 0.)){
                cout << "Two Tau Had"<<endl;
              }

              // Needed to generate a unique seed for the random number function used in the TauID function.
              srand48(entry);
              bool passed_tau1Cand =  TauID(Tau1cand_vec);
              if (passed_tau1Cand) {
                Tau1cand_vec = TauSmearing(Tau1cand_vec);
                ntau_counter++;
              } else {
                Tau1cand_vec.SetPtEtaPhiE(0.,0.,0.,0.);
              }
              srand48(entry+numberOfEntries);
              bool passed_tau2Cand =  TauID(Tau2cand_vec);
              if (passed_tau2Cand) { 
                Tau2cand_vec = TauSmearing(Tau2cand_vec);
                ntau_counter++;
              } else {
                Tau2cand_vec.SetPtEtaPhiE(0.,0.,0.,0.);
              }           
    
	      bool pass_lead_jet_cuts = false;
	      // Events with no cuts
	      
	      pass_cuts[1] = 1;
	      
	      // for events where the tau(s) passes a pT cut with exactly 1 Tau
	      if ((pass_cuts[1] == 1) && (Tau1cand_vec.Pt() > tau1_pt_min) && (ntau_counter==1)){
		pass_cuts[2] = 1; 
	      }
	      // for events where the tau(s) passes an eta cut
	      if ((pass_cuts[2] == 1) && (abs(Tau1cand_vec.Eta()) < tau_eta_max)){
		pass_cuts[3] = 1;
	      }
	      // events passing also b-jet veto requirement
	      if ( (pass_cuts[3]==1) && (!is_b_jet)){pass_cuts[4] = 1;}
	      // events passing jet pt cut
	      if ((Jet_leading_vec.Pt() > lead_jet_pt) && (pass_cuts[4] == 1)){pass_cuts[5] =1;}
	      // events passing jet eta cut
	      if ((abs(Jet_leading_vec.Eta()) < jet_eta_max) && (pass_cuts[5] == 1)){pass_cuts[6] =1;}
	      // events passing dphitaujet cut 
	      double tau1_jet_Dphi = TMath::Abs(Jet_leading_vec.Phi() - Tau1cand_vec.Phi());
	      if ((pass_cuts[6] == 1) && (tau1_jet_Dphi < Dphitaujetmax)){pass_cuts[7] = 1;}
	      // events passing mt cut
	      double transmass = TMath::Sqrt(TMath::Abs(2*Tau1cand_vec.Pt()*MET*(1-TMath::Cos(TMath::Abs(Tau1cand_vec.Phi() - MET_phi)))));
	      if ((pass_cuts[7] == 1) && (transmass > transversemassmin)){pass_cuts[8] = 1;}
	      
	      // events passing also MET cut
	      if ((pass_cuts[8] == 1) && (MET > met_min)){pass_cuts[9] = 1;}
	      // events passing dphijetmet cut
	      double jet_met_dphi = TMath::Abs(Jet_leading_vec.Phi() - MET_phi);
	      if ((pass_cuts[9] == 1) && (jet_met_dphi > deltajetmetphi)){pass_cuts[10]=1;}
	      // events passing effective mass cut
	      
	      double efective_mass = TMath::Sqrt(Jet_leading_vec.Pt()*Jet_leading_vec.Pt()+Tau1cand_vec.Pt()*Tau1cand_vec.Pt()+MET*MET);
	      if ((pass_cuts[10] == 1) && (efective_mass > efectivemassmin)){pass_cuts[11] = 1;}        
	          
	          
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
          if (Jet_leading_vec.Pt() > 1) {
            double jet_met_dphi = TMath::Abs(Jet_leading_vec.Phi() - MET_phi);
            _hmap_jet_met_Dphi[i]->Fill(abs(jet_met_dphi));
            _hmap_jet_met_metDphi[i]->Fill(abs(jet_met_dphi),MET);
	    _hmap_Nevents[i]->Fill(1.0);
          }
	  
          _hmap_met[i]->Fill(MET);
	  if(Jet_leading_vec.Pt() > 1){
	    _hmap_lead_jet_pT[i]->Fill(Jet_leading_vec.Pt());
	    _hmap_lead_jet_eta[i]->Fill(Jet_leading_vec.Eta());
	    _hmap_lead_jet_phi[i]->Fill(Jet_leading_vec.Phi());
	  }
	  if(Tau1cand_vec.Pt() > 2){
	    _hmap_tau1_pT[i]->Fill(Tau1cand_vec.Pt());
	    _hmap_tau1_eta[i]->Fill(Tau1cand_vec.Eta());
	    _hmap_tau1_phi[i]->Fill(Tau1cand_vec.Phi());
            double tau1_met_dphi = TMath::Abs(Tau1cand_vec.Phi() - MET_phi);
            _hmap_tau1_met_Dphi[i]->Fill(tau1_met_dphi);
	    
	    double transmass = TMath::Sqrt(TMath::Abs(2*Tau1cand_vec.Pt()*MET*(1-TMath::Cos(TMath::Abs(Tau1cand_vec.Phi() - MET_phi)))));
	    
 	    _hmap_transverse_mass[i]->Fill(transmass);
            double tau1_jet_Dphi = TMath::Abs(Jet_leading_vec.Phi() - Tau1cand_vec.Phi());
            _hmap_tau1_jet_Dphi[i]->Fill(tau1_jet_Dphi);      
	    double efective_mass = TMath::Sqrt(Jet_leading_vec.Pt()*Jet_leading_vec.Pt()+Tau1cand_vec.Pt()*Tau1cand_vec.Pt()+MET*MET);
            _hmap_efective_mass[i]->Fill(efective_mass);       
	    
	  }
	  if(Tau2cand_vec.Pt() > 2){
	    _hmap_tau2_pT[i]->Fill(Tau2cand_vec.Pt());
	    _hmap_tau2_eta[i]->Fill(Tau2cand_vec.Eta());
	    _hmap_tau2_phi[i]->Fill(Tau2cand_vec.Phi());
            double tau2_met_dphi = TMath::Abs(Tau2cand_vec.Phi() - MET_phi);
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

bool PhenoAnalysis::TauID(TLorentzVector TauCand) {

  bool passedTauID = false;
  if (TauCand.Pt() > 0.){
    double r = drand48();
    if (r <= 0.8){passedTauID = true;}
    TauCand.SetPtEtaPhiE(TauCand.Pt(), TauCand.Eta(), TauCand.Phi(), TauCand.E());
  }
  return passedTauID;

}

TLorentzVector PhenoAnalysis::TauSmearing(TLorentzVector TauCand) {

  if (TauCand.Pt() > 0.){
    double smeared_pt  = 0.03*TauCand.Pt() + TauCand.Pt();
    double smeared_eta = 0.03*TauCand.Eta() + TauCand.Eta();
    double smeared_phi = 0.03*TauCand.Phi() + TauCand.Phi();
    double smeared_e   = 0.03*TauCand.E() + TauCand.E(); 
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
