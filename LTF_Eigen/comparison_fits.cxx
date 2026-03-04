// plotMany.C
// Usage in ROOT: .x plotMany.C

#include <vector>
#include <string>
#include <TMath.h>

void comparison_fits() {

    // -----------------------------
    // User configuration
    // -----------------------------
    std::vector<std::string> filenames = {
      "summary/WbWb_Slurm_Template_162_5.root",
      "summary/WbWb_Slurm_Template_165.root",
      "summary/WbWb_Slurm_Template_167_5.root",
      "summary/WbWb_Slurm_Template_170.root",
      "summary/WbWb_Slurm_Template_172_5.root",
      "summary/WbWb_Slurm_Template_175.root",
      "summary/WbWb_Slurm_Template_177_5.root",
      "summary/WbWb_Slurm_Template_180.root",
      "summary/WbWb_Slurm_Template_182_5.root"
    };

    std::map<string, double> mass_points = {
      {"summary/WbWb_Slurm_Template_162_5.root", 162.5},
      {"summary/WbWb_Slurm_Template_165.root", 165},
      {"summary/WbWb_Slurm_Template_167_5.root", 167.5},
      {"summary/WbWb_Slurm_Template_170.root", 170},
      {"summary/WbWb_Slurm_Template_172_5.root", 172.5},
      {"summary/WbWb_Slurm_Template_175.root", 175},
      {"summary/WbWb_Slurm_Template_177_5.root", 177.5},
      {"summary/WbWb_Slurm_Template_180.root", 180},
      {"summary/WbWb_Slurm_Template_182_5.root", 182.5}
    };

    vector<Color_t> colors = {kCyan-9, kGray+1, kYellow-6, kOrange+7, kOrange-7, kViolet+1, kGray, kRed+1, kOrange-3, kBlue+7, kGreen, kRed};
    //std::vector<std::string> functions = {"pol0", "pol1", "pol2", "pol3", "pol4", "pol5", "pol6", "pol7", "pol8", "pol9",
    //"gaus", "expo", "[0]*TMath::Gaus(x, [1], [2], true) + [3]*TMath::Landau(x, [4], [5], true)"};
    std::vector<std::string> functions = {//"pol6", "pol7", "pol8", "pol9",
					  //"gaus",
      "cheb7", "cheb8", "cheb9"};
//      "[0]*TMath::Gaus(x, [1], [2], true) + [3]*TMath::Gaus(x, [4], [5], true)",
//					  "[0]*TMath::Gaus(x, [1], [2], true) + [3]*TMath::Gaus(x, [4], [5], true) + [6]*TMath::Gaus(x, [7], [8], true)",
//					  "[0]*TMath::Landau(x, [1], [2], true)",
//					  "[0]*TMath::Gaus(x, [1], [2], true) + [3]*TMath::Landau(x, [4], [5], true)"};

    std::string histName = "m_bl_fine";   // same histogram name in all files

    // -----------------------------
    // ROOT objects
    // -----------------------------
    
    TCanvas* c1 = new TCanvas("c1", "Overlayed Histograms", 800, 600);
    c1->Print("test.ps[");
    c1->SetLogy(false);


    TLegend* leg = new TLegend(0.7, 0.7, 0.88, 0.88);
    leg->SetBorderSize(0);

    std::vector<TFile*> files;
    std::vector<TGraphErrors*> graphs;
    std::map<string, std::map<double, vector<double>>> fitParams;
    std::map<string, std::map<double, vector<double>>> fitParamErrors;
    
    int icolor = 0;
    int iParamMax = -1;
    double xmin = 50;
    double xmax = 300;
    // -----------------------------
    // Loop over files
    // -----------------------------
    for (const auto& fname : filenames) {

        TFile* f = TFile::Open(fname.c_str(), "READ");
        if (!f || f->IsZombie()) {
            std::cerr << "Error opening file: " << fname << std::endl;
            continue;
        }

        TH1* h = dynamic_cast<TH1*>(f->Get(histName.c_str()));
        if (!h) {
            std::cerr << "Histogram " << histName
                      << " not found in file " << fname << std::endl;
            f->Close();
            continue;
        }

	TLegend* l = new TLegend(0.7, 0.7, 0.88, 0.88);
	l->SetBorderSize(0);
        // Clone so histogram survives file closure
        TH1* hclone = dynamic_cast<TH1*>(h->Clone());
        hclone->SetDirectory(nullptr);
	TGraphErrors* g = new TGraphErrors();
	bool scaleX = false;
        for (int i = 1; i <= hclone->GetNbinsX(); i++) {
	  double xnew;
	  if (scaleX ) xnew = (2*hclone->GetBinCenter(i) - (hclone->GetXaxis()->GetXmax() + hclone->GetXaxis()->GetXmin())) /
			 (hclone->GetXaxis()->GetXmax() + hclone->GetXaxis()->GetXmin());
          else xnew = hclone->GetBinCenter(i);
	    g->SetPoint(i-1, xnew, hclone->GetBinContent(i));
          g->SetPointError(i-1, 0, hclone->GetBinError(i));
	  cout<<"Set point "<<xnew<<"\t"<<hclone->GetBinError(i)<<endl;
        }
	g->GetYaxis()->SetRangeUser(0.4*hclone->GetMinimum(), 3*hclone->GetMaximum());
	g->SetLineColor(kBlack);
        g->SetLineWidth(2);
	g->SetMarkerColor(kBlack);
	g->SetMarkerStyle(20);
	g->SetMarkerSize(0.5);
	g->SetStats(0);
	g->Draw("PAE0");
	g->SetTitle(Form("Template m_{t}=%.1f",mass_points[fname]));
	l->AddEntry(g, Form("Template m_{t}=%.1f",mass_points[fname]), "P");
	icolor = 0;
	for (const auto& function: functions) {
	  TF1* f = new TF1("f", function.c_str(), xmin, xmax);
	  if ( function == "[0]*TMath::Gaus(x, [1], [2], true) + [3]*TMath::Landau(x, [4], [5], true)" ) {
	    f->SetParameters(1,    // Gauss amplitude
			     -0.8, // Gauss mean
			     0.2,  // Gauss sigma
			     1,    // Landau amplitude
			     -0.6, // Landau MPV
			     0.4   // Landau width
			     );
	  }
	  else if ( function == "[0]*TMath::Landau(x, [1], [2], true)" ) {
            f->SetParameters(1, -0.7, 0.2);
          }
	  else if ( function == "[0]*TMath::Gaus(x, [1], [2], true) + [3]*TMath::Gaus(x, [4], [5], true)" ) {
	    f->SetParameters(1, -0.7, 0.2, 1, -0.7, 0.2);
	  }
	  else if ( function == "[0]*TMath::Gaus(x, [1], [2], true) + [3]*TMath::Gaus(x, [4], [5], true) + [6]*TMath::Gaus(x, [7], [8], true)" ) {
	    f->SetParameters(1, -0.9, 0.2, 1, -0.7, 0.2, 0.1, 0, 0.5);
	  }
	  else if ( function == "cheb7" ) {
	    f->SetParameters(0.018, -0.021, -0.0017, 0.0196, -0.0018, 0.011, -0.002);
	  }
	  else if ( function == "cheb8" ) {
            f->SetParameters(0.013, -0.039, -0.0125, 0.004, -0.026, -0.004, -0.002, -0.001);
          }
	  else if ( function == "cheb9" ) {
            f->SetParameters(0.034, -0.01, 0.025, 0.03, 0.006, 0.012, 0.02, -0.0045, 0.011);
          }
	  TFitResultPtr fitRes = g->Fit(f, "SQ0RW", "", xmin, xmax);
	  if ( f->GetNpar() > iParamMax ) iParamMax = f->GetNpar();
	  for(int i=0; i<f->GetNpar(); i++) {
	    fitParams[function.c_str()][mass_points[fname]].push_back(f->GetParameter(i));
	    fitParamErrors[function.c_str()][mass_points[fname]].push_back(f->GetParError(i));
	  }
	  f->SetLineColor(colors[icolor]);
	  f->DrawClone("same");
	  std::ostringstream oss;
	  if ( function == "gaus") oss << "Gaussian";
	  else if ( function == "[0]*TMath::Landau(x, [1], [2], true)" ) oss <<"Landau";
	  else if ( function == "[0]*TMath::Gaus(x, [1], [2], true) + [3]*TMath::Landau(x, [4], [5], true)" ) oss << "Gaus+Landau";
	  else if ( function == "[0]*TMath::Gaus(x, [1], [2], true) + [3]*TMath::Gaus(x, [4], [5], true)" ) oss << "2 Gaussians";
	  else if ( function == "[0]*TMath::Gaus(x, [1], [2], true) + [3]*TMath::Gaus(x, [4], [5], true) + [6]*TMath::Gaus(x, [7], [8], true)" ) oss << "3 Gaussians";
	  else oss << function;
	  oss << "  (#chi^{2}/ndf = " << fitRes->Chi2() << " / " << fitRes->Ndf() << ")";
	  l->AddEntry(f, oss.str().c_str()  , "l");
	  icolor++;
	}
	g->Draw("same PE0");
	l->DrawClone();
	c1->Print("test.ps");
	c1->Clear();
        graphs.push_back(g);
        files.push_back(f);
    }

    // -----------------------------
    // Draw histograms
    // -----------------------------
    icolor=0;
    for (auto* g : graphs) {
      g->SetLineColor(colors[icolor]);
      g->SetMarkerColor(colors[icolor]);
      if (icolor == 0) {
	g->SetMarkerStyle(1);
	g->Draw("AP E0");
      } else {
	g->Draw("E0 SAME");
      }
      leg->AddEntry(g, g->GetTitle(), "l");

      icolor++;
    }

    leg->Draw();
    c1->Print("test.ps");

    c1->Clear();
    c1->SetLogy(false);
    c1->Divide(iParamMax, functions.size());
    // Draw fit parameters
    int irow=0;
    for (const auto& function: functions) {
      for (int i = 0; i<fitParamErrors[function.c_str()][mass_points["summary/WbWb_Slurm_Template_162_5.root"]].size(); i++){
	TGraphErrors* g = new TGraphErrors();
	int j = 0;
	for (const auto& fname : filenames) {
	  g->SetPoint(j, mass_points[fname], fitParams[function.c_str()][mass_points[fname]][i]);
	  g->SetPointError(j, 0, fitParamErrors[function.c_str()][mass_points[fname]][i]);
	  j++;
	}
	c1->cd(irow*iParamMax+i+1);
	g->SetMarkerColor(kBlack);
	g->SetMarkerSize(.3);
	g->SetMarkerStyle(20);
	g->Draw("APE");

	TFitResultPtr resPol1 = g->Fit("pol1","SQ0W");
	//TMatrixDSym covPol1 = resPol1->GetCovarianceMatrix();
	TF1* pol1 = (TF1*)g->GetFunction("pol1")->Clone("pol1");
	pol1->SetLineColor(kBlue+3);
	pol1->SetLineStyle(1);
	pol1->SetLineWidth(1);
	pol1->Draw("same");

	TFitResultPtr resPol2 = g->Fit("pol2","SQ0W");
        //TMatrixDSym covPol1 = resPol1->GetCovarianceMatrix();
        TF1* pol2 = (TF1*)g->GetFunction("pol2")->Clone("pol2");
        pol2->SetLineColor(kBlue);
        pol2->SetLineStyle(1);
        pol2->SetLineWidth(1);
	pol2->Draw("same");

	g->Draw("PEsame");

	TLegend* l = new TLegend(0.15, 0.65, 0.48, 0.87);
	l->AddEntry(g, Form("param%d",i), "l");
	//l->AddEntry(pol1,Form("pol1: %.2f / %d",resPol1->Chi2(),resPol1->Ndf()),"l");
	//l->AddEntry(pol2,Form("pol2: %.2f / %d",resPol2->Chi2(),resPol2->Ndf()),"l");
	l->AddEntry(pol1,"pol1","l");
	l->AddEntry(pol2,"pol2","l");

	l->SetBorderSize(0);
	l->SetTextSize(0.05);
	l->DrawClone();
      }
      irow++;
    }
    c1->Update();

    c1->Print("test.ps");
    c1->Print("test.ps]");
    // -----------------------------
    // Cleanup (optional)
    // -----------------------------
    for (auto* f : files) {
        if (f) f->Close();
    }
}
