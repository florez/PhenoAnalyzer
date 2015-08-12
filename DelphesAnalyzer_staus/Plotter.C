{
  TFile* f1 = new TFile ("param_sps1a_LSP000_Stau025_Chargino100_run_1_output.root");
  TH1F* jet_pt_1_cut5 = (TH1F*)f1->Get("After_cut_5/jet_pT");
  TH1F* jet_eta_led_1_cut5 = (TH1F*)f1->Get("After_cut_5/jet_eta");
  TH1F* jet_phi_1_cut5 = (TH1F*)f1->Get("After_cut_5/jet_phi");

  TFile* f2 = new TFile ("param_sps1a_LSP000_Stau150_Chargino300_run_1_output.root");
  TH1F* jet_pt_2_cut5 = (TH1F*)f2->Get("After_cut_5/jet_pT");
  TH1F* jet_eta_led_2_cut5 = (TH1F*)f2->Get("After_cut_5/jet_eta");
  TH1F* jet_phi_2_cut5 = (TH1F*)f2->Get("After_cut_5/jet_phi");
  
  TFile* f3 = new TFile ("param_sps1a_LSP250_Stau295_Chargino300_run_1_output.root");
  TH1F* jet_pt_3_cut5 = (TH1F*)f3->Get("After_cut_5/jet_pT");
  TH1F* jet_eta_led_3_cut5 = (TH1F*)f3->Get("After_cut_5/jet_eta");
  TH1F* jet_phi_3_cut5 = (TH1F*)f3->Get("After_cut_5/jet_phi");

  TCanvas *c_jet = new TCanvas("c_jet", "c_jet");
  jet_pt_1_cut5->SetLineColor(kRed);
  jet_pt_2_cut5->SetLineColor(kBlue);
  jet_pt_3_cut5->SetLineColor(kBlack);
  jet_pt_1_cut5->Draw();
  jet_pt_2_cut5->Draw("same");
  jet_pt_3_cut5->Draw("same");

    TLegend* leg = new TLegend(0.6, 0.6, 0.9, 0.9);
  leg->AddEntry(jet_pt_1_cut5, "CUT 5 m(#tilde{#tau}) = 25 GeV, m(#tilde{#chi}^{0}_{1}) = 0 GeV, m(#tilde{#chi}^{#pm}) = 100 GeV");
  leg->AddEntry(jet_pt_2_cut5, "CUT 5 m(#tilde{#tau}) = 50 GeV, m(#tilde{#chi}^{0}_{1}) = 0 GeV, m(#tilde{#chi}^{#pm}) = 100 GeV");
  leg->AddEntry(jet_pt_3_cut5, "CUT 5 m(#tilde{#tau}) = 75 GeV, m(#tilde{#chi}^{0}_{1}) = 0 Gev,  m(#tilde{#chi}^{#pm}) = 100 GeV");
  leg->Draw();

   TCanvas *c_jet_eta_led = new TCanvas("c_jet_eta_led", "c_jet_eta_led");
  jet_eta_led_1_cut5->SetLineColor(kRed);
  jet_eta_led_2_cut5->SetLineColor(kBlack);
  jet_eta_led_3_cut5->SetLineColor(kBlue);
  jet_eta_led_1_cut5->Draw();
  jet_eta_led_2_cut5->Draw("same");
  jet_eta_led_3_cut5->Draw("same");

   leg->Draw();

   TCanvas *c_jet_phi = new TCanvas("c_jet_phi", "c_jet_phi");
  jet_phi_1_cut5->SetLineColor(kRed);
  jet_phi_2_cut5->SetLineColor(kBlack);
  jet_phi_3_cut5->SetLineColor(kBlue);
  jet_phi_1_cut5->Draw();
  jet_phi_2_cut5->Draw("same");
  jet_phi_3_cut5->Draw("same");

  leg->Draw();
  
}
