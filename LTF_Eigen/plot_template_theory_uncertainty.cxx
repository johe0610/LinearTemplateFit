void plot_template_theory_uncertainty() {

  vector<string> variables = {
    "dr_bl"    ,
    "dr_bw"    ,
    "eta_lep1" ,
    "mT_lep1"  ,
    "m_bl"     ,
    "m_bw"     ,
    "m_minimax",
    "m_wbbl"   ,
    "mblmin"   ,
    "mtblmin"  ,
    "pT_Wbbl"  ,
    "pT_Wbblnu",
    "pT_bjet1" ,
    "pT_bjet2" ,
    "pT_bl"    ,
    "pT_blnu"  ,
    "pT_bw"    ,
    "pT_lep1"  ,
    "pT_met"   ,
    "pT_owj1"  ,
    "pT_owj2"  ,
    "pT_whad"  ,
    "r_bj"     ,
    "r_bq"     ,
    "y_whad"
  };

  TCanvas c1("canvas", "Histograms with Ratio", 800, 800);
  c1.Print("theory_uncertainty.pdf[");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TFile *f = TFile::Open("files/template_theory_uncertainty.root");

  bool doDrawAll = false;
  
  for (auto var: variables) {
    TLegend* legend0 = new TLegend(0.68, 0.78, 0.80, 0.84);
    TLegend* legend1 = new TLegend(0.68, 0.60, 0.80, 0.78);
    TLegend* legend2 = new TLegend(0.81, 0.60, 0.93, 0.84);

    c1.Divide(1, 2, 0, 0);
    c1.cd(1)->SetBottomMargin(0);
    c1.cd(1)->SetTopMargin(0.15);
    c1.cd(1)->SetRightMargin(0.05);

    c1.cd(2)->SetTopMargin(0);
    c1.cd(2)->SetBottomMargin(0.15);
    c1.cd(2)->SetRightMargin(0.05);

    c1.cd(1)->SetLogy();
    TH1 *h = (TH1*)f->Get(var.c_str());
    TH1* h_err_up   = (TH1*)h->Clone(Form("%s_err_up",var.c_str()));
    TH1* h_err_down = (TH1*)h->Clone(Form("%s_err_down",var.c_str()));
    h_err_up->Reset("C");
    h_err_down->Reset("C");

    h->SetLineColor(kWhite);
    h->GetYaxis()->SetTitle("d#sigma/dx [pb]");
    h->Draw("hist");

    c1.cd(2);

    TH1 *h_isr_aS_up = (TH1*)f->Get(Form("%s_templ_unc_isr_aS_up", var.c_str()));
    TH1 *h_isr_aS_down = (TH1*)f->Get(Form("%s_templ_unc_isr_aS_down", var.c_str()));
    h_isr_aS_up->SetLineColor(632);
    h_isr_aS_down->SetLineColor(632);
    h_isr_aS_up->Divide(h);
    h_isr_aS_down->Divide(h);
    h_isr_aS_up->GetYaxis()->SetRangeUser(-.53, .53);
    h_isr_aS_up->GetYaxis()->SetTitle("Relative uncertainty");
    h_isr_aS_up->GetXaxis()->SetTitle(var.c_str());
    h_isr_aS_up->Draw("hist");
    h_isr_aS_down->Draw("hist same");
    legend1->AddEntry(h_isr_aS_up, "ISR #alpha_{s}", "L");
    
    // Scale variations, keep largest uncertainty
    vector<string> scale_vars = {"muR_05_muF_05", "muR_05_muF_10", "muR_10_muF_05", "muR_10_muF_20", "muR_20_muF_10", "muR_20_muF_20"};
    TH1 *h_scale_up = (TH1*)h->Clone("h_scale_up");
    TH1 *h_scale_down = (TH1*)h->Clone("h_scale_down");
    h_scale_up->Reset("C");
    h_scale_down->Reset("C");
    for ( int i = 0; i<scale_vars.size(); i++ ) {
      TH1* hist = (TH1*)f->Get(Form("%s_templ_unc_%s", var.c_str(), scale_vars[i].c_str()));
      for ( int j = 1; j <= hist->GetNbinsX(); j++ ) {
	h_scale_up->SetBinContent(j, max(h_scale_up->GetBinContent(j), hist->GetBinContent(j)));
	h_scale_down->SetBinContent(j, min(h_scale_down->GetBinContent(j), hist->GetBinContent(j)));
      }
    }
    h_scale_up->SetLineColor(634);
    h_scale_down->SetLineColor(634);
    h_scale_up->Divide(h);
    h_scale_down->Divide(h);
    h_scale_up->Draw("hist same");
    h_scale_down->Draw("hist same");
    legend1->AddEntry(h_scale_up, "Scale Var.", "L");
    
    // For aS, take difference between errors and symmetrize
    TH1 *h_aS_up = (TH1*)f->Get(Form("%s_templ_unc_aS_up", var.c_str()));
    TH1 *h_aS_down = (TH1*)f->Get(Form("%s_templ_unc_aS_down", var.c_str()));
    h_aS_up->SetLineColor(416);
    h_aS_up->Add(h_aS_down, -1.);
    h_aS_up->Divide(h);
    h_aS_up->Draw("hist same");
    h_aS_up->Scale(-1);
    h_aS_up->Draw("hist same");
    legend1->AddEntry(h_aS_up, "#alpha_{S} var.", "L");
    
    TH1 *h_fsr_scale_up = (TH1*)f->Get(Form("%s_templ_unc_fsr_scale_up", var.c_str()));
    TH1 *h_fsr_scale_down = (TH1*)f->Get(Form("%s_templ_unc_fsr_scale_down", var.c_str()));
    h_fsr_scale_up->SetLineColor(800);
    h_fsr_scale_down->SetLineColor(800);
    h_fsr_scale_up->Divide(h);
    h_fsr_scale_down->Divide(h);
    h_fsr_scale_up->Draw("hist same");
    h_fsr_scale_down->Draw("hist same");
    legend1->AddEntry(h_fsr_scale_up, "FSR scale", "L");

    // Underlying event from A14 Var1 variation
    TH1 *h_ue_up = (TH1*)f->Get(Form("%s_templ_unc_A14Var1Up", var.c_str()));
    TH1 *h_ue_down = (TH1*)f->Get(Form("%s_templ_unc_A14Var1Down", var.c_str()));
    h_ue_up->SetLineColor(616);
    h_ue_down->SetLineColor(616);
    h_ue_up->Divide(h);
    h_ue_down->Divide(h);
    h_ue_up->Draw("hist same");
    h_ue_down->Draw("hist same");
    legend1->AddEntry(h_ue_up, "UE", "L");

    // Color reconnection. Take the larger from the CR1 and CR2 variations
    // Should this be symmetrized?
    TH1 *h_cr_up = (TH1*)h->Clone("h_cr_up");
    h_cr_up->Reset("C");
    {
      TH1 *h_cr1= (TH1*)f->Get(Form("%s_templ_unc_CR1", var.c_str()));
      TH1 *h_cr2 = (TH1*)f->Get(Form("%s_templ_unc_CR2", var.c_str()));
      for ( int i = 1; i<=h_cr1->GetNbinsX(); i++ ) {
	h_cr_up->SetBinContent(i, max(abs(h_cr1->GetBinContent(i)), abs(h_cr2->GetBinContent(i))));
      }
    }
    h_cr_up->SetLineColor(432);
    h_cr_up->Divide(h);
    TH1 *h_cr_down = (TH1*)h_cr_up->Clone("h_cr_down");
    h_cr_down->Scale(-1);
    h_cr_up->Draw("hist same");
    h_cr_down->Draw("hist same");
    legend1->AddEntry(h_cr_up, "CR", "L");

    // Sum quadratically over all PDF contributions
    TH1 *h_PDF_up = (TH1*)f->Get(Form("%s_templ_unc_PDF_90900", var.c_str()))->Clone("h_PDF_up"); // histo has bin entries zero
    for ( int i = 1; i<= 30; i++ ) {
      string index = Form("0%d", i);
      if ( i > 9 ) index = Form("%d", i);
      TH1 *h_PDF_var_up = (TH1*)f->Get(Form("%s_templ_unc_PDF_909%s", var.c_str(), index.c_str()));
      for ( int j = 1; j <= h_PDF_up->GetNbinsX(); j++ ) {
	h_PDF_up->SetBinContent(j, std::sqrt(std::pow(h_PDF_up->GetBinContent(j), 2) + std::pow(h_PDF_var_up->GetBinContent(j), 2)));
      }
    }
    h_PDF_up->SetLineColor(400);
    h_PDF_up->Divide(h);
    TH1 *h_PDF_down = (TH1*)h_PDF_up->Clone("h_PDF_down");
    h_PDF_down->Scale(-1.);
    h_PDF_up->Draw("hist same");
    h_PDF_down->Draw("hist same");
    legend1->AddEntry(h_PDF_up, "PDF", "L");

    vector<Color_t> colors = {//kCyan-9, kGray+1, kYellow-6, kOrange+7, kOrange-7, kViolet+1, kGray, kRed+1, kOrange-3, kBlue+7};
      820, 840, 860,
      880, 900, 600,
      602, 801, 920
    };

    vector<TString> sys_list = {
      "hdamp", "herwig",
      "pthard", "PS",
      "recoilToTop", "bquark_frag",
      "lineshape", "DS",
      "NNLO_tt"
    };
    
    for ( int i = 0; i<sys_list.size(); i++ ) {
      TH1* hist = (TH1*)f->Get(Form("%s_templ_unc_%s", var.c_str(), sys_list[i].Data()));
      hist->SetLineColor(colors[i]);
      hist->Divide(h);
      TH1* hist_sym = (TH1*)hist->Clone(Form("%s_templ_unc_sym_%s", var.c_str(), sys_list[i].Data()));
      hist->Draw("hist same");
      hist_sym->Scale(-1);
      hist_sym->Draw("hist same");

      //cout<<sys_list[i]<<endl;
      //hist->Print("All");
      
      for ( int j = 1; j <= h->GetNbinsX(); j++ ) {
	h_err_up->SetBinContent(j, sqrt(pow(h_err_up->GetBinContent(j),2)+pow(hist->GetBinContent(j),2)));
	h_err_down->SetBinContent(j, sqrt(pow(h_err_down->GetBinContent(j),2)+pow(hist->GetBinContent(j),2)));
      }
      legend2->AddEntry(hist, sys_list[i], "L");
    }

    
    for ( int i = 1; i <= h->GetNbinsX(); i++ ) {
      
      // aS error symmetrized
      h_err_up->SetBinContent(i, sqrt(pow(h_err_up->GetBinContent(i),2)+pow(h_aS_up->GetBinContent(i),2)));
      h_err_down->SetBinContent(i, sqrt(pow(h_err_down->GetBinContent(i),2)+pow(h_aS_up->GetBinContent(i),2)));
      
      if ( h_isr_aS_up->GetBinContent(i) > 0 )   h_err_up->SetBinContent(i, sqrt(pow(h_err_up->GetBinContent(i),2)+pow(h_isr_aS_up->GetBinContent(i),2)));
      else  h_err_down->SetBinContent(i, sqrt(pow(h_err_down->GetBinContent(i),2)+pow(h_isr_aS_up->GetBinContent(i),2)));
      if ( h_isr_aS_down->GetBinContent(i) > 0 ) h_err_up->SetBinContent(i, sqrt(pow(h_err_up->GetBinContent(i),2)+pow(h_isr_aS_down->GetBinContent(i),2)));
      else  h_err_down->SetBinContent(i, sqrt(pow(h_err_down->GetBinContent(i),2)+pow(h_isr_aS_down->GetBinContent(i),2)));
      
      if ( h_scale_up->GetBinContent(i) > 0 )   h_err_up->SetBinContent(i, sqrt(pow(h_err_up->GetBinContent(i),2)+pow(h_scale_up->GetBinContent(i),2)));
      else  h_err_down->SetBinContent(i, sqrt(pow(h_err_down->GetBinContent(i),2)+pow(h_scale_up->GetBinContent(i),2)));
      if ( h_scale_down->GetBinContent(i) > 0 ) h_err_up->SetBinContent(i, sqrt(pow(h_err_up->GetBinContent(i),2)+pow(h_scale_down->GetBinContent(i),2)));
      else  h_err_down->SetBinContent(i, sqrt(pow(h_err_down->GetBinContent(i),2)+pow(h_scale_down->GetBinContent(i),2)));

      if ( h_fsr_scale_up->GetBinContent(i) > 0 )   h_err_up->SetBinContent(i, sqrt(pow(h_err_up->GetBinContent(i),2)+pow(h_fsr_scale_up->GetBinContent(i),2)));
      else  h_err_down->SetBinContent(i, sqrt(pow(h_err_down->GetBinContent(i),2)+pow(h_fsr_scale_up->GetBinContent(i),2)));
      if ( h_fsr_scale_down->GetBinContent(i) > 0 ) h_err_up->SetBinContent(i, sqrt(pow(h_err_up->GetBinContent(i),2)+pow(h_fsr_scale_down->GetBinContent(i),2)));
      else  h_err_down->SetBinContent(i, sqrt(pow(h_err_down->GetBinContent(i),2)+pow(h_fsr_scale_down->GetBinContent(i),2)));
      
      if ( h_PDF_up->GetBinContent(i) > 0 )   h_err_up->SetBinContent(i, sqrt(pow(h_err_up->GetBinContent(i),2)+pow(h_PDF_up->GetBinContent(i),2)));
      else  h_err_down->SetBinContent(i, sqrt(pow(h_err_down->GetBinContent(i),2)+pow(h_PDF_up->GetBinContent(i),2)));
      if ( h_PDF_down->GetBinContent(i) > 0 ) h_err_up->SetBinContent(i, sqrt(pow(h_err_up->GetBinContent(i),2)+pow(h_PDF_down->GetBinContent(i),2)));
      else  h_err_down->SetBinContent(i, sqrt(pow(h_err_down->GetBinContent(i),2)+pow(h_PDF_down->GetBinContent(i),2)));

      if ( h_cr_up->GetBinContent(i) > 0 )   h_err_up->SetBinContent(i, sqrt(pow(h_err_up->GetBinContent(i),2)+pow(h_cr_up->GetBinContent(i),2)));
      else  h_err_down->SetBinContent(i, sqrt(pow(h_err_down->GetBinContent(i),2)+pow(h_cr_up->GetBinContent(i),2)));
      if ( h_cr_down->GetBinContent(i) > 0 ) h_err_up->SetBinContent(i, sqrt(pow(h_err_up->GetBinContent(i),2)+pow(h_cr_down->GetBinContent(i),2)));
      else  h_err_down->SetBinContent(i, sqrt(pow(h_err_down->GetBinContent(i),2)+pow(h_cr_down->GetBinContent(i),2)));

      if ( h_ue_up->GetBinContent(i) > 0 )   h_err_up->SetBinContent(i, sqrt(pow(h_err_up->GetBinContent(i),2)+pow(h_ue_up->GetBinContent(i),2)));
      else  h_err_down->SetBinContent(i, sqrt(pow(h_err_down->GetBinContent(i),2)+pow(h_ue_up->GetBinContent(i),2)));
      if ( h_ue_down->GetBinContent(i) > 0 ) h_err_up->SetBinContent(i, sqrt(pow(h_err_up->GetBinContent(i),2)+pow(h_ue_down->GetBinContent(i),2)));
      else  h_err_down->SetBinContent(i, sqrt(pow(h_err_down->GetBinContent(i),2)+pow(h_ue_down->GetBinContent(i),2)));
    }

    h_err_up->SetLineColor(kBlack);
    h_err_down->SetLineColor(kBlack);
    h_err_up->Draw("hist same");
    h_err_down->Scale(-1.0);
    h_err_down->Draw("hist same");
    
    c1.cd(1);
    TGraphAsymmErrors* graph = new TGraphAsymmErrors(h);

    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        graph->SetPoint(i-1, h->GetBinCenter(i), h->GetBinContent(i));
        graph->SetPointError(i-1, 0.0, 0.0, -1.0 * h_err_down->GetBinContent(i)*h->GetBinContent(i), h_err_up->GetBinContent(i)*h->GetBinContent(i));
	// Absolute value has to be provided for SetPointError, therefore -1 * yerr_down
    }
    graph->SetLineColor(kBlack);
    graph->SetMarkerColor(kBlack);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.5);
    graph->Draw("P same");
    
    legend0->AddEntry(graph, "Pseudo-data", "PE");
    legend0->AddEntry(h_err_up, "Total unc.", "L");

    legend0->SetMargin(0.12);
    legend0->SetBorderSize(0);
    legend0->SetFillColor(kWhite);
    legend0->SetTextSize(0.03);
    legend0->Draw("same");
    legend1->SetMargin(0.12);
    legend1->SetBorderSize(0);
    legend1->SetFillColor(kWhite);
    legend1->SetTextSize(0.03);
    legend1->Draw("same");
    legend2->SetMargin(0.12);
    legend2->SetBorderSize(0);
    legend2->SetFillColor(kWhite);
    legend2->SetTextSize(0.03);
    legend2->Draw("same");

    //c1.Update();
    c1.Print("theory_uncertainty.pdf");
    c1.Clear();
    c1.Update();

  }
  c1.Print("theory_uncertainty.pdf]");
  f->Close();
}
 
