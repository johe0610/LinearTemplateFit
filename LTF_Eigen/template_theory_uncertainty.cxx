// template_theory_uncertainty.C
void template_theory_uncertainty() {
  
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
  
  TFile *f_default = TFile::Open("files/template_nominal_with_weights.root");
  TFile *fout = TFile::Open("files/template_theory_uncertainty.root", "RECREATE");

    // For later: Add a plotting function that compares the nominal sample, including all uncertainties to the templates
    
    // Compare nominal PS to template PS
    {
      TFile *f_templ_25 = TFile::Open("files/template_172_25.root");
      TFile *f_templ_75 = TFile::Open("files/template_172_75.root");
      
      for (auto var: variables) {
	TH1 *h = (TH1*)f_default->Get(var.c_str());
	TH1* h_clone = (TH1*)h->Clone(var.c_str());
	TH1 *h_templ_25 = (TH1*)f_templ_25->Get(var.c_str());
	TH1 *h_templ_75 = (TH1*)f_templ_75->Get(var.c_str());
	
	if (!h || !h_templ_75 || !h_templ_25) continue;
	
        TH1 *h_diff_templ = (TH1*)h_templ_75->Clone(Form("%s_templ_unc_PS", var.c_str()));
	h_diff_templ->Add(h_templ_25);
	h_diff_templ->Scale(0.5); // Build average from template with m = 172.25 and m = 175.75
        h_diff_templ->Add(h, -1.0);  // h1 - h2
        //h_diff_templ->Divide(h);     // (h1 - h2) / h2
	fout->cd();
	//cout<<"Absolute error"<<endl;
	//h_diff_templ->Print("All");
	//cout<<"Cross section"<<endl;
	//h_clone->Print("All");
	//cout<<"Relative error"<<endl;
	//for ( int i = 1; i<=h_clone->GetNbinsX(); i++){
	//  cout<<h_diff_templ->GetBinContent(i) / h_clone->GetBinContent(i)<<endl;
	//}
	h_clone->Write();
	h_diff_templ->Write();
      }
      f_templ_25->Close();
      f_templ_75->Close();
    }

    // Renormalization and factorization scale variation, 7-point scale variation
    {
      for (auto var: variables) {
	vector<string> scale_vars = {"muR_05_muF_05", "muR_05_muF_10", "muR_10_muF_05", "muR_10_muF_20", "muR_20_muF_10", "muR_20_muF_20"};
	TH1 *h = (TH1*)f_default->Get(var.c_str());
	if (!h) continue;
	//TH1 *h_scale_up = (TH1*)h->Clone(Form("%s_templ_unc_scale_up", var.c_str()));
        //TH1 *h_scale_down = (TH1*)h->Clone(Form("%s_templ_unc_scale_down", var.c_str()));
	for (auto scale: scale_vars) {
	  TH1 *h_scale = (TH1*)f_default->Get(Form("%s[_%s_]", var.c_str(), scale.c_str()))->Clone(Form("%s_templ_unc_%s", var.c_str(), scale.c_str()));
	  if (!h_scale) continue;
	  h_scale->Add(h, -1.0);
	  //h_scale->Divide(h);
	  fout->cd();
	  h_scale->Write();

//	  
//	  for ( int i = 1; i <= h->GetNbinsX(); i++ ) {
//	    h_scale_up->SetBinContent(i, std::max(h_scale->GetBinContent(i), h_scale_up->GetBinContent(i)));
//	    h_scale_down->SetBinContent(i,std::min(h_scale->GetBinContent(i), h_scale_down->GetBinContent(i)));
//	    //if ( h_scale->GetBinContent(i) > h_scale_up->GetBinContent(i) ) h_scale_up->SetBinContent(i, h_scale->GetBinContent(i));
//	    //if ( h_scale->GetBinContent(i) < h_scale_down->GetBinContent(i) ) h_scale_down->SetBinContent(i, h_scale->GetBinContent(i));
//	  }
//	}
//        h_scale_up->Add(h, -1.0);
//        h_scale_up->Divide(h);
//	h_scale_down->Add(h, -1.0);
//        h_scale_down->Divide(h);
//	fout->cd();
//        h_scale_up->Write();
//        h_scale_down->Write();
	}
      }
    }
    // alpha_S uncertainty
    {
      for (auto var: variables) {
	TH1 *h = (TH1*)f_default->Get(var.c_str());
        if (!h) continue;

	TH1 *h_aS_up = (TH1*)f_default->Get(Form("%s%s", var.c_str(), "[_PDF_265000_]"))->Clone(Form("%s_templ_unc_aS_up", var.c_str()));;
	TH1 *h_aS_down = (TH1*)f_default->Get(Form("%s%s", var.c_str(), "[_PDF_266000_]"))->Clone(Form("%s_templ_unc_aS_down", var.c_str()));

	h_aS_up->Add(h, -1.0);
        //h_aS_up->Divide(h);
	h_aS_down->Add(h, -1.0);
        //h_aS_down->Divide(h);
	fout->cd();
        h_aS_up->Write();
        h_aS_down->Write();
      }
    }

    // PDF uncertainties
    // Write out all systematic variations
    // Average over uncertainties after the fit
    // Choose between NNPDF3.0 and PDF4LHC
    {
      for (auto var: variables) {

        TH1 *h = (TH1*)f_default->Get(Form("%s[_PDF_90900_]",var.c_str()));
        if (!h) continue;
	TIter nextKey(f_default->GetListOfKeys());
	TKey* key;
	while ((key = (TKey*)nextKey())) {
	  std::string keyName = key->GetName();
	  if ( keyName.find(var.c_str()) == std::string::npos ) {
            continue;
          }
	  if ( keyName.find("blnu") != std::string::npos && var.find("blnu") == std::string::npos) {
            continue;
          }
	  if ( keyName.find("_fine") != std::string::npos ) {
            continue;
          }
	  if ( keyName.find("[_PDF_909") == std::string::npos ) {
            continue;
	  }
	  string name=keyName.substr(keyName.length() - 8).substr(0,6);
	  TH1D* h_PDF = (TH1D*) f_default->Get(keyName.c_str())->Clone(Form("%s_templ_unc_PDF%s", var.c_str(), name.c_str()));
	  h_PDF->Add(h, -1.0);
	  //h_PDF->Divide(h);
	  fout->cd();
	  h_PDF->Write();
	  h_PDF->Delete();
	}
      }
    }

    // ISR alpha_s variation is evaluated with Var3c variation
    {
      for (auto var: variables) {
        TH1 *h = (TH1*)f_default->Get(var.c_str());
	if (!h) continue;

        TH1 *h_isr_aS_up = (TH1*)f_default->Get(Form("%s%s", var.c_str(), "[Var3cUp]"))->Clone(Form("%s_templ_unc_isr_aS_up", var.c_str()));;
        TH1 *h_isr_aS_down = (TH1*)f_default->Get(Form("%s%s", var.c_str(), "[Var3cDown]"))->Clone(Form("%s_templ_unc_isr_aS_down", var.c_str()));

        h_isr_aS_up->Add(h, -1.0);
        //h_isr_aS_up->Divide(h);
        h_isr_aS_down->Add(h, -1.0);
        //h_isr_aS_down->Divide(h);
        fout->cd();
        h_isr_aS_up->Write();
        h_isr_aS_down->Write();
      }
    }

    // FSR scale variation, in case of large fluctuations, use dedicated samples instead
    {
      for (auto var: variables) {
        TH1 *h = (TH1*)f_default->Get(var.c_str());
        if (!h) continue;

        TH1 *h_fsr_scale_up = (TH1*)f_default->Get(Form("%s%s", var.c_str(), "[isrmuRfac10_fsrmuRfac05]"))->Clone(Form("%s_templ_unc_fsr_scale_up", var.c_str()));;
        TH1 *h_fsr_scale_down = (TH1*)f_default->Get(Form("%s%s", var.c_str(), "[isrmuRfac10_fsrmuRfac20]"))->Clone(Form("%s_templ_unc_fsr_scale_down", var.c_str()));

        h_fsr_scale_up->Add(h, -1.0);
        //h_fsr_scale_up->Divide(h);
        h_fsr_scale_down->Add(h, -1.0);
        //h_fsr_scale_down->Divide(h);
        fout->cd();
        h_fsr_scale_up->Write();
        h_fsr_scale_down->Write();
      }
    }


    {
      for (auto var: variables) {
        TH1 *h = (TH1*)f_default->Get(var.c_str());
        if (!h) continue;

	TFile *f_herwig = TFile::Open("files/template_uncertainty_herwig.root");
	TH1* h_herwig = (TH1*)f_herwig->Get(var.c_str())->Clone(Form("%s_templ_unc_herwig", var.c_str()));
	h_herwig->Add(h, -1.0);
	fout->cd();
	h_herwig->Write();
	f_herwig->Close();
	delete f_herwig;

	TFile *f_hdamp = TFile::Open("files/template_uncertainty_hdamp.root");
	TH1* h_hdamp = (TH1*)f_hdamp->Get(var.c_str())->Clone(Form("%s_templ_unc_hdamp", var.c_str()));
	h_hdamp->Add(h, -1.0);
        fout->cd();
        h_hdamp->Write();
	f_hdamp->Close();
        delete f_hdamp;

	TFile *f_pthard = TFile::Open("files/template_uncertainty_pthard.root");
        TH1* h_pthard = (TH1*)f_pthard->Get(var.c_str())->Clone(Form("%s_templ_unc_pthard", var.c_str()));
        h_pthard->Add(h, -1.0);
        fout->cd();
        h_pthard->Write();
	f_pthard->Close();
        delete f_pthard;

	TFile *f_recoilToTop = TFile::Open("files/template_uncertainty_recoilToTop.root");
	TH1* h_recoilToTop = (TH1*)f_recoilToTop->Get(var.c_str())->Clone(Form("%s_templ_unc_recoilToTop", var.c_str()));
        h_recoilToTop->Add(h, -1.0);
        fout->cd();
        h_recoilToTop->Write();
        f_recoilToTop->Close();
        delete f_recoilToTop;
	
	TFile *f_A14Var1Up = TFile::Open("files/template_uncertainty_A14Var1Up.root");
        TH1* h_A14Var1Up = (TH1*)f_A14Var1Up->Get(var.c_str())->Clone(Form("%s_templ_unc_A14Var1Up", var.c_str()));
        h_A14Var1Up->Add(h, -1.0);
        fout->cd();
        h_A14Var1Up->Write();
	f_A14Var1Up->Close();
        delete f_A14Var1Up;
	
	TFile *f_A14Var1Down = TFile::Open("files/template_uncertainty_A14Var1Down.root");
        TH1* h_A14Var1Down = (TH1*)f_A14Var1Down->Get(var.c_str())->Clone(Form("%s_templ_unc_A14Var1Down", var.c_str()));
        h_A14Var1Down->Add(h, -1.0);
        fout->cd();
        h_A14Var1Down->Write();
        f_A14Var1Down->Close();
        delete f_A14Var1Down;

	TFile *f_bquark_frag = TFile::Open("files/template_uncertainty_bquark_frag.root");
        TH1* h_bquark_frag = (TH1*)f_bquark_frag->Get(var.c_str())->Clone(Form("%s_templ_unc_bquark_frag", var.c_str()));
        h_bquark_frag->Add(h, -1.0);
        fout->cd();
        h_bquark_frag->Write();
	f_bquark_frag->Close();
        delete f_bquark_frag;

	TFile *f_lineshape = TFile::Open("files/template_uncertainty_lineshape.root");
        TH1* h_lineshape = (TH1*)f_lineshape->Get(var.c_str())->Clone(Form("%s_templ_unc_lineshape", var.c_str()));
        h_lineshape->Add(h, -1.0);
        fout->cd();
        h_lineshape->Write();
	f_lineshape->Close();
        delete f_lineshape;

	TFile *f_DS = TFile::Open("files/template_uncertainty_DS.root");
        TH1* h_DS = (TH1*)f_DS->Get(var.c_str())->Clone(Form("%s_templ_unc_DS", var.c_str()));
        h_DS->Add(h, -1.0);
        fout->cd();
        h_DS->Write();
	f_DS->Close();
        delete f_DS;

	TFile *f_CR0 = TFile::Open("files/template_uncertainty_CR0.root");
	TH1* h_CR0 = (TH1*)f_CR0->Get(var.c_str());
	TFile *f_CR1 = TFile::Open("files/template_uncertainty_CR1.root");
        TH1* h_CR1 = (TH1*)f_CR1->Get(var.c_str())->Clone(Form("%s_templ_unc_CR1", var.c_str()));
        h_CR1->Add(h_CR0, -1.0);
        fout->cd();
        h_CR1->Write();
	f_CR1->Close();
        delete f_CR1;

	TFile *f_CR2 = TFile::Open("files/template_uncertainty_CR2.root");
	TH1* h_CR2 = (TH1*)f_CR2->Get(var.c_str())->Clone(Form("%s_templ_unc_CR2", var.c_str()));
	h_CR2->Add(h_CR0, -1.0);
	fout->cd();
        h_CR2->Write();
	f_CR2->Close();
        delete f_CR2;
	f_CR0->Close();
        delete f_CR0;

	//cout<<"Checking sign of CR variation for obs "<<var<<endl;
	//for( int i = 1; i<=h_CR2->GetNbinsX(); i++) {
	//  if ( h_CR1->GetBinContent(i) * h_CR2->GetBinContent(i) < 0 ) {
	//    cout<<"Obs "<<var<<" bin "<<i<<" CR1 "<<h_CR1->GetBinContent(i)/h->GetBinContent(i)<<" CR2 "<<h_CR2->GetBinContent(i)/h->GetBinContent(i)<<endl;
	//  }
	//}
	//cout<<endl;
	
	TFile *f_NNLO = TFile::Open("files/template_uncertainty_nnlo_reweighting.root");
        TH1* h_NNLO = (TH1*)f_DS->Get(var.c_str())->Clone(Form("%s_templ_unc_NNLO_tt", var.c_str()));
        h_NNLO->Add(h, -1.0);
        fout->cd();
        h_NNLO->Write();
        f_NNLO->Close();
        delete f_NNLO;

      }
    }
    
    fout->Close();
    f_default->Close();


}
