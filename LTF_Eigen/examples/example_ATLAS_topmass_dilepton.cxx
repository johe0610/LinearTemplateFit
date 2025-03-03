/*
  Copyright (c) 2021, D. Britzger, Max-Planck-Institute for Physics, Munich, Germany

  Permission is hereby granted, free of charge, to any person obtaining
  a copy of this software and associated documentation files (the
  "Software"), to deal in the Software without restriction, including
  without limitation the rights to use, copy, modify, merge, publish,
  distribute, sublicense, and/or sell copies of the Software, and to
  permit persons to whom the Software is furnished to do so, subject to
  the following conditions:

  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.

  The Software is provided "as is", without warranty of any kind,
  express or implied, including but not limited to the warranties of
  merchantability, fitness for a particular purpose and
  noninfringement. In no event shall the authors or copyright holders be
  liable for any claim, damages or other liability, whether in an action
  of contract, tort or otherwise, arising from, out of or in connection
  with the Software or the use or other dealings in the Software.
*/
// -------------------------------------------------------------------- //

#include <iostream>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TKey.h>
#include <LTF/LTF_ROOTTools.h>
#include <LTF/LTF.h>
#include <TH2D.h>

void PrintAsciiTable(const map<double,TH1D*>&, TH1D* data);

int fitMultipleObservables(const char* ps_name, const vector<TString> fit_vars, const vector<TString> fit_vars_rivet);

vector<vector<double > > TH2D_to_vecvec(TH2D* hist2D) {
   const int n = hist2D->GetNbinsX();
   vector<vector<double > > vecvec(n);
   for ( size_t i = 0 ; i<n ; i++ ) {
      vecvec[i].resize(n);
      for ( size_t j = 0 ; j<n ; j++ ) {
         vecvec[i][j] = hist2D->GetBinContent(i+1,j+1);
      }
   }
   return vecvec;
}

int example_ATLAS_topmass_dilepton() {


  //if (fitMultipleObservables("plots_dilepton/fit_minimax.ps", {"22"}, {"minimaxmbl"}) > 0) return 1;
  //if (fitMultipleObservables("plots_dilepton/fit_mT_bb4l.ps", {"5"}, {"mT_bb4l"}) > 0) return 1;
  if (fitMultipleObservables("plots_dilepton/fit_m_bbll.ps", {"3"}, {"m_bbll"}) > 0) return 1;
  
  return 0;
}



//! ------------------------------------------------------------------------ //
//! main function
#ifndef __CLING__

int example_ATLAS_topmass_dilepton();

int main(int ,const char **) {

   gROOT->SetBatch();

   return example_ATLAS_topmass_dilepton();
}
#endif

int fitMultipleObservables(const char* ps_name, const vector<TString> file_number, const vector<TString> fit_vars_rivet) {
   using namespace std;

#ifdef __CLING__
   TH1D::AddDirectory(false);
#endif
   TH1::SetDefaultSumw2(true);

   map<double,TH1D*> templates;
   
   const int     iRebin       = 1;
   const int     iRebinData   = 1;
   const TString datafile = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output_dilepton/HEPData-1739985874-v1-Table_"+file_number[0]+".root";
   //const TString datafile = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output_dilepton/S3beta_Cluster_WbWb_dilepton_170.root";
   //const TString covariancefile = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_170_1258_matrices.root";

   int bins_number = 0;
   //for ( auto& tmp: fit_vars_rivet ) {
   for ( auto& tmp: file_number ) { //johannes loop over files, not over observables (only true for data and systematics
     TString name = "Table "+tmp+"/Hist1D_y1";
     TH1F* tmp_data = TFile::Open(datafile)->Get<TH1F>(name);
     //TH1D* tmp_data = TFile::Open(datafile)->Get<TH1D>(tmp);

     //if ( !tmp_data ) { cerr<<"Could not find data histogram " << name <<endl; exit(1);}
     //else cout<<"Found data histogram "<<name<<endl;
     cout<<"Adding "<<tmp_data->GetNbinsX()-1<<" bins for variable "<<tmp<<endl;
     bins_number += tmp_data->GetNbinsX()-1;
     tmp_data->Clear();
   }
   
   TH1D* combined_data = new TH1D("combined_data", "combined_data", bins_number, 0, bins_number);
   int bin_offset = 0;
   for ( auto& tmp: file_number ) {
   TString name = "Table "+tmp+"/Hist1D_y1";
   TH1F* tmp_data = TFile::Open(datafile)->Get<TH1F>(name);
   //for ( auto& tmp: fit_vars_rivet ) {
   //TH1D* tmp_data = TFile::Open(datafile)->Get<TH1D>(tmp);
     for ( int i = 1; i <= tmp_data->GetNbinsX()-1; i++ ) {
       combined_data->SetBinContent(i+bin_offset, tmp_data->GetBinContent(i));
     }
     bin_offset += tmp_data->GetNbinsX()-1;
   }
   combined_data -> Rebin(iRebinData);
   //combined_data->Scale(scale);
   combined_data->SetLineColor(kBlack);
   combined_data->SetMarkerSize(1.8);
   cout<<"This is the new combined data hist"<<endl;
   combined_data->Print("All");
   
   {
     TH1D* combined_template_160 = new TH1D("combined_template_160", "combined_template_160", bins_number, 0, bins_number);
     TH1D* combined_template_165 = new TH1D("combined_template_165", "combined_template_165", bins_number, 0, bins_number);
     TH1D* combined_template_170 = new TH1D("combined_template_170", "combined_template_170", bins_number, 0, bins_number);
     TH1D* combined_template_175 = new TH1D("combined_template_175", "combined_template_175", bins_number, 0, bins_number);
     TH1D* combined_template_180 = new TH1D("combined_template_180", "combined_template_180", bins_number, 0, bins_number);
     int bin_offset = 0;
     //vector<double> lumi={803.612954972981, 904.407021431961, 1145.08363345207, 1176.18904529949, 1354.78090409675};
     vector<double> lumi={815.805084365513, 919.340378420208, 1170.27887095863, 1190.7914999182, 1370.88882772107};

     for ( auto& tmp: fit_vars_rivet ) {
       TH1D* h_tmp_160 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output_dilepton/S3beta_Cluster_WbWb_dilepton_160.root")->Get<TH1D>("/RAW"+tmp);
       TH1D* h_tmp_165 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output_dilepton/S3beta_Cluster_WbWb_dilepton_165.root")->Get<TH1D>("/RAW"+tmp);
       TH1D* h_tmp_170 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output_dilepton/S3beta_Cluster_WbWb_dilepton_170.root")->Get<TH1D>("/RAW"+tmp);
       TH1D* h_tmp_175 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output_dilepton/S3beta_Cluster_WbWb_dilepton_175.root")->Get<TH1D>("/RAW"+tmp);
       TH1D* h_tmp_180 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output_dilepton/S3beta_Cluster_WbWb_dilepton_180.root")->Get<TH1D>("/RAW"+tmp);
       for ( int i = 1; i<= h_tmp_160->GetNbinsX(); i++ ) {
	 combined_template_160->SetBinContent(i+bin_offset, h_tmp_160->GetBinContent(i) / h_tmp_160->GetXaxis()->GetBinWidth(i) / lumi[0] / 1000);
         combined_template_165->SetBinContent(i+bin_offset, h_tmp_165->GetBinContent(i) / h_tmp_160->GetXaxis()->GetBinWidth(i) / lumi[1] / 1000);
         combined_template_170->SetBinContent(i+bin_offset, h_tmp_170->GetBinContent(i) / h_tmp_160->GetXaxis()->GetBinWidth(i) / lumi[2] / 1000);
         combined_template_175->SetBinContent(i+bin_offset, h_tmp_175->GetBinContent(i) / h_tmp_160->GetXaxis()->GetBinWidth(i) / lumi[3] / 1000);
         combined_template_180->SetBinContent(i+bin_offset, h_tmp_180->GetBinContent(i) / h_tmp_160->GetXaxis()->GetBinWidth(i) / lumi[4] / 1000);

	 combined_template_160->SetBinError(i+bin_offset, h_tmp_160->GetBinError(i) / h_tmp_160->GetXaxis()->GetBinWidth(i));
         combined_template_165->SetBinError(i+bin_offset, h_tmp_165->GetBinError(i) / h_tmp_160->GetXaxis()->GetBinWidth(i));
         combined_template_170->SetBinError(i+bin_offset, h_tmp_170->GetBinError(i) / h_tmp_160->GetXaxis()->GetBinWidth(i));
         combined_template_175->SetBinError(i+bin_offset, h_tmp_175->GetBinError(i) / h_tmp_160->GetXaxis()->GetBinWidth(i));
         combined_template_180->SetBinError(i+bin_offset, h_tmp_180->GetBinError(i) / h_tmp_160->GetXaxis()->GetBinWidth(i));
       }
       bin_offset += h_tmp_160->GetNbinsX();
     }
     templates[160] = combined_template_160;
     templates[165] = combined_template_165;
     templates[170] = combined_template_170;
     templates[175] = combined_template_175;
     templates[180] = combined_template_180;
     }

   for ( auto [MM,hist] : templates ) {
      hist->Rebin(iRebin);
   }

   // ------------------------------------------------ //
   // ---  Do linear template fit
   // ------------------------------------------------ //
   // --- instantiate LTF object
   LTF ltf;
   ltf.SetGamma(vector<double>{1});
   ltf.UseNuisanceParameters(true);
   ltf.UseLogNormalUncertainties(false);

   // --- initialize templates
   for ( auto [MM,hist] : templates ) {
      ltf.AddTemplate(MM,  hist->GetNbinsX(),  hist->GetArray()+1 ); // set template
      ltf.AddTemplateErrorSquared("statY", MM , hist->GetNbinsX(), hist->GetSumw2()->GetArray()+1, 0.); // set template error dY
   }

   // --- initialize data
   ltf.SetData( combined_data->GetNbinsX(), combined_data->GetArray()+1);
   
   std::unique_ptr<TFile> file(TFile::Open(datafile));
   
   if (!file || file->IsOpen() == kFALSE) {
      std::cerr << "Error: Couldn't open the file!" << std::endl;
      return 1;
   }

   /*
   TH2D* combined_covariance = new TH2D("combined_cov", "combined_cov", bins_number, 0, bins_number, bins_number, 0, bins_number);
   bin_offset = 0;
   for ( auto& fit_var: fit_vars ) {
     TString histnameCovStat("unfolding_"+fit_var+"_NOSYS"); // "unfolding_mbl_selected_NOSYS"  ->  unfolding_covariance_matrix_ptl1_covariance_STAT_DATA
      histnameCovStat.ReplaceAll("unfolding_","unfolding_covariance_matrix_");
      histnameCovStat.ReplaceAll("_NOSYS","_covariance_STAT_DATA");
      TH2D* cov_stat_dat = file->Get<TH2D>(histnameCovStat);
      if ( !cov_stat_dat ) { cerr<<"Could not find covariance matrix " << histnameCovStat <<endl; exit(1);}
      else cout<<"Found covarinace matrix "<<histnameCovStat<<endl;
      TH1D* hist = file->Get<TH1D>("unfolding_error_"+fit_var+"_direct_envelope_STAT_DATA__1up");
      for ( int i = 1; i < cov_stat_dat->GetNbinsX(); i++ ) {
	for ( int j = 1; j < cov_stat_dat->GetNbinsY(); j++ ) {
	  double bin_width_x = hist->GetXaxis()->GetBinWidth(i);
	  double bin_width_y = hist->GetXaxis()->GetBinWidth(j);
	  combined_covariance->SetBinContent(i+bin_offset,j+bin_offset, cov_stat_dat->GetBinContent(i,j)/(bin_width_x*bin_width_y));
	  //cout<<"Filling histo with "<<cov_stat_dat->GetBinContent(i,j)<<" / ("<<bin_width_x<<" * "<<bin_width_y<<")"<<endl;
	}
      }
      bin_offset += cov_stat_dat->GetNbinsX() - 1;
      //cov_stat_dat->Print("All");
   }
   //combined_covariance->Print("All");
   vector<vector<double > > vecCov2_1 = TH2D_to_vecvec(combined_covariance);
   for ( auto& tmp_vec: vecCov2_1 ){
     for ( auto& tmp: tmp_vec ) cout<<tmp<<"\t";
     cout<<endl;
   }
   if ( fit_vars.size() > 1 ) {
     bin_offset = 0;
     for ( int v1 = 0; v1 < fit_vars_rivet.size(); v1++ ) {
       for ( int v2 = v1+1; v2 < fit_vars_rivet.size(); v2++ ) {
	 TString histname = fit_vars_rivet[v1]+"_"+fit_vars_rivet[v2];
	 TH2D* cov = TFile::Open(covariancefile)->Get<TH2D>(histname);
	 if ( !cov ) { cerr<<"Could not find covariance matrix " << histname <<endl; exit(1);}
	 else cout<<"Found covarinace matrix "<<histname<<endl;
	 TH1D* h_err_v1 = file->Get<TH1D>("unfolding_error_"+fit_vars[v1]+"_direct_envelope_STAT_DATA__1up");
	 TH1D* h_err_v2 = file->Get<TH1D>("unfolding_error_"+fit_vars[v2]+"_direct_envelope_STAT_DATA__1up"); // get relative error
	 //TH1D* hist_data = ... //get cross section in single bin // might have to divide by bin width later
	 bin_offset += h_err_v1->GetNbinsX()-1;
	 // Projections of correlation matrix
	 TH1D* projection_var1 = cov->ProjectionX("pro_v1",0,-1,"e");
	 TH1D* projection_var2 = cov->ProjectionY("pro_v2",0,-1,"e");
	 //projection_var1->Print("All");
	 //projection_var2->Print("All");
	 
	 // Cross section file
	 TH1D* h_data_var1 = file->Get<TH1D>("unfolding_"+fit_vars[v1]+"_NOSYS");
	 TH1D* h_data_var2 = file->Get<TH1D>("unfolding_"+fit_vars[v2]+"_NOSYS");
	 if ( !h_data_var1 ) { cerr<<"Could not find data for " << fit_vars[v1] <<endl; exit(1);}
         else cout<<"Found data for "<<fit_vars[v1]<<endl;
	 if ( !h_data_var2 ) { cerr<<"Could not find data for " << fit_vars[v2] <<endl; exit(1);}
         else cout<<"Found data for "<<fit_vars[v2]<<endl;
	 //cout<<"Var1 content ";
	 //for ( int i = 1; i <= h_data_var1->GetNbinsX(); i++ ) cout<<projection_var1->GetBinContent(i) / h_data_var1->GetBinContent(i)<<"\t";
	 //cout<<endl<<"Var1 error   ";
	 //for ( int i = 1; i <= h_data_var1->GetNbinsX(); i++ ) cout<<projection_var1->GetBinError(i) / h_data_var1->GetBinError(i)<<"\t";
	 //cout<<endl<<"Var2 content ";
	 //for ( int i = 1; i <= h_data_var2->GetNbinsX(); i++ ) cout<<projection_var2->GetBinContent(i) / h_data_var2->GetBinContent(i)<<"\t";
	 //cout<<endl<<"Var2 error   ";
	 //for ( int i = 1; i <= h_data_var2->GetNbinsX(); i++ ) cout<<projection_var2->GetBinError(i) / h_data_var2->GetBinError(i)<<"\t";
	 //cout<<endl;

	 for ( int i = 1; i < h_err_v1->GetNbinsX(); i++ ) {
	   for ( int j = 1; j < h_err_v2->GetNbinsX(); j++ ) {
	     // Calculate the covariance in the data from the one in Sherpa
	     // 1. Get the Sherpa covariance: std::pow(cov->GetBinError(i,j),2)
	     // 2. Divide by the std deviation from Sherpa: sigma_template1, sigma_template2
	     // 3. Multiply with the std deviation from the data: sigma_data1, sigma_data2
	     double sigma_data1 = h_data_var1->GetBinContent(i) * h_err_v1->GetBinContent(i);
	     double sigma_data2 = h_data_var2->GetBinContent(j) * h_err_v2->GetBinContent(j);
	     double sigma_template1 = projection_var1->GetBinError(i);
	     double sigma_template2 = projection_var2->GetBinError(j);
	     cout<<"Correlation "<<std::pow(cov->GetBinError(i,j),2)/ sigma_template1 / sigma_template2<<endl;
	     combined_covariance->SetBinContent(i, j+bin_offset, std::pow(cov->GetBinError(i,j),2)*(sigma_data1*sigma_data2) / (sigma_template1*sigma_template2));
	     combined_covariance->SetBinContent(j+bin_offset, i, std::pow(cov->GetBinError(i,j),2)*(sigma_data1*sigma_data2) / (sigma_template1*sigma_template2));
	   }
	 }
       }
     }
   }
   vector<vector<double > > vecCov2 = TH2D_to_vecvec(combined_covariance);
   for ( auto& tmp_vec: vecCov2 ){
     for ( auto& tmp: tmp_vec ) cout<<tmp<<"\t";
     cout<<endl;
   }
   ltf.AddErrorRelative("STAT_DATA", vecCov2, LTF::Uncertainty::Constrained);
*/


   vector<string> uncertainties = {"stat", "systematics", "total uncertainty", "MC stat", //will not be used
				   "Electron energy resolution", "Electron energy scale AF2", "Electron energy scale", "Muon ID momentum resolution",
				   "Muon MS momentum resolution", "Muon sagitta residual bias", "Muon sagitta rho", "Muon momentum scale",
				   "Electron ID efficiency", "Electron isolation efficiency", "Electron reconstruction efficiency",
				   "Electron trigger efficiency", "Muon ID efficiency (stat.)", "Muon ID efficiency (syst.)", "Muon isolation efficiency (stat.)",
				   "Muon isolation efficiency (syst.)", "Muon trigger efficiency (stat.)", "Muon trigger efficiency (syst.)",
				   "Muon TTVA efficiency (stat.)", "Muon TTVA efficiency (syst.)", // lepton uncertainties
				   "Jet vertex tagger efficiency","b-jet energy scale", "Jet energy scale (detector NP1)", "Jet energy scale (detector NP2)",
				   "Jet energy scale (mixed NP1)", "Jet energy scale (mixed NP2)", "Jet energy scale (mixed NP3)", "Jet energy scale (modelling NP1)",
				   "Jet energy scale (modelling NP2)", "Jet energy scale (modelling NP3)", "Jet energy scale (modelling NP4)", "Jet energy scale (statistical NP1)",
				   "Jet energy scale (statistical NP2)", "Jet energy scale (statistical NP3)", "Jet energy scale (statistical NP4)",
				   "Jet energy scale (statistical NP5)", "Jet energy scale (statistical NP6)", "Jet energy scale $\eta$ intercalib (modelling)",
				   "Jet energy scale $\eta$ intercalib (non-closure 2018 data)", "Jet energy scale $\eta$ intercalib (non-closure high $E$)",
				   "Jet energy scale $\eta$ intercalib (non-closure neg $\eta$)", "Jet energy scale $\eta$ intercalib (non-closure pos $\eta$)",
				   "Jet energy scale $\eta$ intercalib (stat.)","Jet energy scale (flavour composition)",
				   "Jet energy scale (flavour response)", // JES uncertainties
				   "Jet energy resolution (data vs MC)", "Jet energy resolution (NP1)", "Jet energy resolution (NP10)", "Jet energy resolution (NP11)",
				   "Jet energy resolution (NP12restTerm)", "Jet energy resolution (NP2)", "Jet energy resolution (NP3)", "Jet energy resolution (NP4)", "Jet energy resolution (NP5)",
				   "Jet energy resolution (NP6)", "Jet energy resolution (NP7)", "Jet energy resolution (NP8)", "Jet energy resolution (NP9)", // JER uncertainties
				   "Jet energy scale (pileup offset)", "Jet energy scale (pileup NPV)", "Jet energy scale (pileup pT term)", "Jet energy scale (pileup rho topology",
				   "Jet energy scale (punch-through)", "Jet energy scale (single particle high pT)", "MET soft term resolution (parallel)",
				   "MET soft term resolution (perpendicular)", "MET soft term scale", // pile-up and MET
				   "b-jet efficiency (eigenvar 0)",
				   "b-jet efficiency (eigenvar 1)", "b-jet efficiency (eigenvar 10)", "b-jet efficiency (eigenvar 11)", "b-jet efficiency (eigenvar 12)",
				   "b-jet efficiency (eigenvar 13)", "b-jet efficiency (eigenvar 14)", "b-jet efficiency (eigenvar 15)", "b-jet efficiency (eigenvar 16)",
				   "b-jet efficiency (eigenvar 17)", "b-jet efficiency (eigenvar 18)", "b-jet efficiency (eigenvar 19)", "b-jet efficiency (eigenvar 2)",
				   "b-jet efficiency (eigenvar 20)", "b-jet efficiency (eigenvar 21)", "b-jet efficiency (eigenvar 22)", "b-jet efficiency (eigenvar 23)",
				   "b-jet efficiency (eigenvar 24)", "b-jet efficiency (eigenvar 25)", "b-jet efficiency (eigenvar 26)", "b-jet efficiency (eigenvar 27)",
				   "b-jet efficiency (eigenvar 28)", "b-jet efficiency (eigenvar 29)", "b-jet efficiency (eigenvar 3)", "b-jet efficiency (eigenvar 30)",
				   "b-jet efficiency (eigenvar 31)", "b-jet efficiency (eigenvar 32)", "b-jet efficiency (eigenvar 33)", "b-jet efficiency (eigenvar 34)",
				   "b-jet efficiency (eigenvar 35)", "b-jet efficiency (eigenvar 36)", "b-jet efficiency (eigenvar 37)", "b-jet efficiency (eigenvar 38)",
				   "b-jet efficiency (eigenvar 39)", "b-jet efficiency (eigenvar 4)", "b-jet efficiency (eigenvar 40)", "b-jet efficiency (eigenvar 41)",
				   "b-jet efficiency (eigenvar 42)", "b-jet efficiency (eigenvar 43)", "b-jet efficiency (eigenvar 44)", "b-jet efficiency (eigenvar 5)",
				   "b-jet efficiency (eigenvar 6)", "b-jet efficiency (eigenvar 7)", "b-jet efficiency (eigenvar 8)", "b-jet efficiency (eigenvar 9)",
				   "c-jet efficiency (eigenvar 0)", "c-jet efficiency (eigenvar 1)", "c-jet efficiency (eigenvar 10)", "c-jet efficiency (eigenvar 11)",
				   "c-jet efficiency (eigenvar 12)", "c-jet efficiency (eigenvar 13)", "c-jet efficiency (eigenvar 14)", "c-jet efficiency (eigenvar 15)",
				   "c-jet efficiency (eigenvar 16)", "c-jet efficiency (eigenvar 17)", "c-jet efficiency (eigenvar 18)", "c-jet efficiency (eigenvar 19)",
				   "c-jet efficiency (eigenvar 2)", "c-jet efficiency (eigenvar 3)", "c-jet efficiency (eigenvar 4)", "c-jet efficiency (eigenvar 5)",
				   "c-jet efficiency (eigenvar 6)", "c-jet efficiency (eigenvar 7)", "c-jet efficiency (eigenvar 8)", "c-jet efficiency (eigenvar 9)",
				   "light-jet efficiency (eigenvar 0)", "light-jet efficiency (eigenvar 1)", "light-jet efficiency (eigenvar 10",
				   "light-jet efficiency (eigenvar 11", "light-jet efficiency (eigenvar 12", "light-jet efficiency (eigenvar 13", "light-jet efficiency (eigenvar 14",
				   "light-jet efficiency (eigenvar 15", "light-jet efficiency (eigenvar 16", "light-jet efficiency (eigenvar 17", "light-jet efficiency (eigenvar 18",
				   "light-jet efficiency (eigenvar 19", "light-jet efficiency (eigenvar 2)", "light-jet efficiency (eigenvar 3)", "light-jet efficiency (eigenvar 4)",
				   "light-jet efficiency (eigenvar 5)", "light-jet efficiency (eigenvar 6)", "light-jet efficiency (eigenvar 7)", "light-jet efficiency (eigenvar 8)",
				   "light-jet efficiency (eigenvar 9)", //b-tagging uncertainties
				   "pileup reweighting", "$ttV$ normalisation", "Diboson normalisation", "Fakes normalisation", "Z+jets",
				   "ttbar normalisation", "tW normalisation", "tW DS vs DR", "Luminosity", // background uncertainties
				   "h_{damp}", "FSR #mu_{R}", "Scale #mu_{R}", "Scale #mu_{F}", "ISR #alpha_{S} Var3c", "Parton shower",
				   "Matching", "Recoil to top", "PDF", "top mass" //modeling uncertainties
   };
   

   
   // Systematical uncertainties
   for ( int i = 5; i <= 175; i++ ) {
     vector<double> combined_error;
     for ( auto& tmp: file_number ) { //johannes loop over files, not variables
       TH1F* hist = file->Get<TH1F>("Table "+tmp+"/Hist1D_y1_e"+std::to_string(i)+"plus");
       cout<<uncertainties[i-1]<<"\t";
       for (int j=1; j< hist->GetNbinsX(); j++) {
	 cout<<hist->GetBinContent(j)<<"\t";
	 combined_error.push_back(hist->GetBinContent(j));
       }
       cout<<endl;
     }
     double corr = 1.0;
     if ( combined_error.size() > 0 ) ltf.AddError(uncertainties[i-1], combined_error, corr, LTF::Uncertainty::Constrained);
     combined_error.clear();
   }
   
   // Stat uncertainties
   {
     vector<double> combined_error;
     for ( auto& tmp: file_number ) {
       TH1F* hist = file->Get<TH1F>("Table "+tmp+"/Hist1D_y1_e1"); // data stat
       for (int i=1; i< hist->GetNbinsX(); i++) {
         combined_error.push_back(hist->GetBinContent(i));
       }
     }
     double corr = 0.0;
     if ( combined_error.size() > 0 ) ltf.AddError("STAT_DATA", combined_error, corr, LTF::Uncertainty::Constrained);
     combined_error.clear();
   }
   {
     vector<double> combined_error;
     for ( auto& tmp: file_number ) {
       TH1F* hist = file->Get<TH1F>("Table "+tmp+"/Hist1D_y1_e4plus"); // MC stat
       for (int i=1; i< hist->GetNbinsX(); i++) {
         combined_error.push_back(hist->GetBinContent(i));
       }
     }
     double corr = 0.0;
     if ( combined_error.size() > 0 ) ltf.AddError("STAT_MC", combined_error, corr, LTF::Uncertainty::Constrained);
     combined_error.clear();
   }

   
//   // External uncertainties -> e2 (sys unc) and e3 (total unc)
//   for ( auto& uncertainty: external_uncertainties ) {
//     vector<double> combined_error;
//     for ( auto& fit_variable: fit_vars ) {
//       TH1D* hist = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_"+uncertainty+"__1up");
//       for (int i=1; i< hist->GetNbinsX(); i++) {
//         combined_error.push_back(hist->GetBinContent(i));
//       }
//     }
//     double corr = 1;
//     if ( combined_error.size() > 0 ) ltf.AddErrorRelative(uncertainty, combined_error, corr, LTF::Uncertainty::External);
//     combined_error.clear();
//   }

   PrintAsciiTable(templates,combined_data);

   LTF::LiTeFit fit = ltf.DoLiTeFit();
   fit.PrintFull();

   vector<double> bins{};
   for(int i = 1; i <= combined_data->GetNbinsX(); i++) {
      bins.push_back(combined_data->GetBinLowEdge(i));
   }
   bins.push_back(combined_data->GetXaxis()->GetBinUpEdge(combined_data->GetNbinsX()));
   string label = "";
   for (TString tmp: fit_vars_rivet) label += tmp+"\t"; 

   LTF_ROOTTools::plotLiTeFit(fit, bins, ps_name, "1/#sigma d#sigma/dx", label,"m_{t} [GeV]"); // Johannes take care of plotting later

   return 0;
}


//! ------------------------------------------------------------------------ //
//! --- write templates, and data, to ascii file
void PrintAsciiTable(const map<double,TH1D*>& templates, TH1D* data){
   cout<<endl;
   printf(" %11s %11s",Form("Data"),Form("Stat"));
   for ( auto [mean,hist] : templates ) 
      printf(" %11s %11s",Form("Tmpl_%5.2f",mean),Form("Stat_%5.2f",mean));
   cout<<endl;
   for ( int i=1; i<=data->GetNbinsX() ;i++ ) {
      printf(" %11.4f %11.4f",data->GetBinContent(i),data->GetBinError(i));
   for ( auto [mean,hist] : templates ) 
         printf(" %11.4f %11.4f",hist->GetBinContent(i),hist->GetBinError(i));
      cout<<endl;
   }
   cout<<endl;
}

