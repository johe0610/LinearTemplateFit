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

int fitMultipleObservables(const char* ps_name, const vector<TString> fit_vars, const vector<TString> fit_vars_short);

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


  if (fitMultipleObservables("plots/fit_mbl.ps", {"mbl_selected"},    {"m_bl"}) > 0) return 1;

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

int fitMultipleObservables(const char* ps_name, const vector<TString> fit_vars, const vector<TString> fit_vars_short) {
   using namespace std;

#ifdef __CLING__
   TH1D::AddDirectory(false);
#endif
   TH1::SetDefaultSumw2(true);

   map<double,TH1D*> templates;
   
   const int     iRebin       = 1;
   const int     iRebinData   = 1;
   const TString datafile = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/unfolding_SR_Whad_Final_l_Whad_particle_TUnfoldStandalone_OptionA_data_nonClosureAlternative.root";
   //const TString covariancefile = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_170_1258_matrices.root";

   int bins_number = 0;
   for ( auto& tmp: fit_vars ) { //johannes loop over files, not over observables (only true for data and systematics
     TString name = "Table 23/Hist1D_"+tmp;
     TH1D* tmp_data = TFile::Open(datafile)->Get<TH1D>(name);
     if ( !tmp_data ) { cerr<<"Could not find data histogram " << name <<endl; exit(1);}
     else cout<<"Found data histogram "<<name<<endl;
     cout<<"Adding "<<tmp_data->GetNbinsX()-1<<" bins for variable "<<tmp<<endl;
     bins_number += tmp_data->GetNbinsX()-1;
     tmp_data->Clear();
   }
   
   TH1D* combined_data = new TH1D("combined_data", "combined_data", bins_number, 0, bins_number);
   int bin_offset = 0;
   for ( auto& tmp: fit_vars ) {
     TString name = "Table 23/Hist1D_"+tmp;
     TH1D* tmp_data = TFile::Open(datafile)->Get<TH1D>(name);
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
     vector<double> scaledBy = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
     for ( auto& tmp: fit_vars_short ) {
       TH1D* h_tmp_160 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_160_1256.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_165 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_165_1246.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_170 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_170_1248.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_175 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_175_1250.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_180 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_180_1252.root")->Get<TH1D>(tmp);
       for ( int i = 1; i<= h_tmp_160->GetNbinsX(); i++ ) {
	 combined_template_160->SetBinContent(i+bin_offset, h_tmp_160->GetBinContent(i)*scaledBy[1] / h_tmp_160->GetXaxis()->GetBinWidth(i) * 1000); // check if cross sec needs to be scaled
	 combined_template_165->SetBinContent(i+bin_offset, h_tmp_165->GetBinContent(i)*scaledBy[2] / h_tmp_160->GetXaxis()->GetBinWidth(i) * 1000);
	 combined_template_170->SetBinContent(i+bin_offset, h_tmp_170->GetBinContent(i)*scaledBy[3] / h_tmp_160->GetXaxis()->GetBinWidth(i) * 1000);
	 combined_template_175->SetBinContent(i+bin_offset, h_tmp_175->GetBinContent(i)*scaledBy[4] / h_tmp_160->GetXaxis()->GetBinWidth(i) * 1000);
	 combined_template_180->SetBinContent(i+bin_offset, h_tmp_180->GetBinContent(i)*scaledBy[5] / h_tmp_160->GetXaxis()->GetBinWidth(i) * 1000);

         combined_template_160->SetBinError(i+bin_offset, h_tmp_160->GetBinError(i)*scaledBy[1] / h_tmp_160->GetXaxis()->GetBinWidth(i) * 1000);
         combined_template_165->SetBinError(i+bin_offset, h_tmp_165->GetBinError(i)*scaledBy[2] / h_tmp_160->GetXaxis()->GetBinWidth(i) * 1000);
         combined_template_170->SetBinError(i+bin_offset, h_tmp_170->GetBinError(i)*scaledBy[3] / h_tmp_160->GetXaxis()->GetBinWidth(i) * 1000);
         combined_template_175->SetBinError(i+bin_offset, h_tmp_175->GetBinError(i)*scaledBy[4] / h_tmp_160->GetXaxis()->GetBinWidth(i) * 1000);
         combined_template_180->SetBinError(i+bin_offset, h_tmp_180->GetBinError(i)*scaledBy[5] / h_tmp_160->GetXaxis()->GetBinWidth(i) * 1000);
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
     for ( int v1 = 0; v1 < fit_vars_short.size(); v1++ ) {
       for ( int v2 = v1+1; v2 < fit_vars_short.size(); v2++ ) {
	 TString histname = fit_vars_short[v1]+"_"+fit_vars_short[v2];
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

   
   // Systematical uncertainties
   for ( int i = 5; i <= 175; i++ ) {
     vector<double> combined_error;
     for ( auto& fit_variable: fit_vars ) { //johannes loop over files, not variables
       TH1D* hist = file->Get<TH1D>("Table 23/Hist1D_"+fit_variable+"_e"+std::to_string(i)+"plus");
       for (int j=1; j< hist->GetNbinsX(); j++) {
	 combined_error.push_back(hist->GetBinContent(j));
       }
     }
     double corr = 1.0;
     if ( combined_error.size() > 0 ) ltf.AddErrorRelative("e"+i, combined_error, corr, LTF::Uncertainty::Constrained);
     combined_error.clear();
   }
   
   // Stat uncertainties
   {
     vector<double> combined_error;
     for ( auto& fit_variable: fit_vars ) {
       TH1D* hist = file->Get<TH1D>("Table 23/Hist1D_"+fit_variable+"_e1"); // data stat
       for (int i=1; i< hist->GetNbinsX(); i++) {
         combined_error.push_back(hist->GetBinContent(i));
       }
     }
     double corr = 0.0;
     if ( combined_error.size() > 0 ) ltf.AddErrorRelative("e1", combined_error, corr, LTF::Uncertainty::Constrained);
     combined_error.clear();
   }
   {
     vector<double> combined_error;
     for ( auto& fit_variable: fit_vars ) {
       TH1D* hist = file->Get<TH1D>("Table 23/Hist1D_"+fit_variable+"_e4plus"); // MC stat
       for (int i=1; i< hist->GetNbinsX(); i++) {
         combined_error.push_back(hist->GetBinContent(i));
       }
     }
     double corr = 0.0;
     if ( combined_error.size() > 0 ) ltf.AddErrorRelative("e4", combined_error, corr, LTF::Uncertainty::Constrained);
     combined_error.clear();
   }

   
//   // External uncertainties
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
   for (TString tmp: fit_vars_short) label += tmp+"\t"; 
   //LTF_ROOTTools::plotLiTeFit(fit, bins, ps_name, "1/#sigma d#sigma/dx", label,"m_{t} [GeV]"); // Johannes take care of plotting later

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

