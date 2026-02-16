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



//Johannes: To do
//  Update templates, run the fit again for all observables, also combinations
//  Possibly exclude templates or bins from the fit
//  Run the LINEAR template fit and check if observables where edm is small are well described by chisq parabola






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


int example_ATLAS_topmass() {

  if (fitMultipleObservables("plots/fit_mbl.ps", {"mbl_selected"},    {"m_bl"}) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_mbw.ps", {"mbwhad_selected"},    {"m_bw"}) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptl1.ps", {"ptl1"},    {"pT_lep1"}) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_drbl.ps", {"dRbl_selected"},    {"dr_bl"}) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_drbw.ps", {"dRbwhad_selected"},    {"dr_bw"}) > 0) return 1;
//  //if (fitMultipleObservables("plots/fit_etal1.ps", {"etal1"},    {"eta_lep1"}) > 0) return 1;
//  //if (fitMultipleObservables("plots/fit_mtlepmet.ps", {"mtlepmet"},    {"mT_lep1met"}) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_minimax.ps", {"minimax_whadbbl"},    {"m_minimax"}) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_mwbbl.ps", {"mwhadbbl"},    {"m_wbbl"}) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptb1.ps", {"ptb1"},    {"pT_bjet1"}) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptb2.ps", {"ptb2"},    {"pT_bjet2"}) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptbl.ps", {"ptbl_selected"},    {"pT_bl"}) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptbw.ps", {"ptbwhad_selected"},    {"pT_bw"}) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptj1.ps", {"ptj1"},    {"pT_jet1"}) > 0) return 1;
//  //if (fitMultipleObservables("plots/fit_ptj2.ps", {"ptj2"},    {"pT_jet2"}) > 0) return 1;
//  //if (fitMultipleObservables("plots/fit_ptmet.ps", {"met"},    {"pT_met"}) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptowj1.ps", {"ptOWj1"},    {"pT_owj1"}) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptowj2.ps", {"ptOWj2"},    {"pT_owj2"}) > 0) return 1;
//  //ptwhadbbl
//  if (fitMultipleObservables("plots/fit_ywhad.ps", {"rapiditywhad"},    {"y_whad"}) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptwhad.ps", {"ptwhad"},    {"pT_whad"}) > 0) return 1;
//
//  if (fitMultipleObservables("plots/fit_mbl_mbw.ps", {"mbl_selected", "mbwhad_selected"}, {"m_bl", "m_bw"}) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_mbl_ptl1.ps", {"mbl_selected", "ptl1"}, {"m_bl", "pT_lep1"}) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_mbw_ptl1.ps", {"mbwhad_selected", "ptl1"}, {"m_bw", "pT_lep1"}) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_mbl_mbw_ptl1.ps", {"mbl_selected", "mbwhad_selected", "ptl1"}, {"m_bl", "m_bw", "pT_lep1"}) > 0) return 1;
//  //  if (fitMultipleObservables("plots/fit_mbl_mbw_ptl1_ptw.ps", {"mbl_selected", "mbwhad_selected", "ptl1", "ptwhad"}, {"m_bl", "m_bw", "pT_lep1", "pT_whad"}) > 0) return 1;
  
  return 0;
}



//! ------------------------------------------------------------------------ //
//! main function
#ifndef __CLING__

int example_ATLAS_topmass();

int main(int ,const char **) {

   gROOT->SetBatch();
   return example_ATLAS_topmass();
}
#endif

int fitMultipleObservables(const char* ps_name, const vector<TString> fit_vars, const vector<TString> fit_vars_short) {
   using namespace std;
  
#ifdef __CLING__
   TH1D::AddDirectory(false);
#endif
   TH1::SetDefaultSumw2(true);

   map<double,TH1D*> templates;
   bool doPseudo = true;
   bool doJesStressTest = false;

   const int     iRebin       = 1;
   const int     iRebinData   = 1;
   const int     iRemoveBins  = 0;
   const TString datafile = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/unfolding_SR_Whad_Final_l_Whad_particle_TUnfoldStandalone_OptionA_data_nonClosureAlternative.root";
   TString pseudodatafile;
   if ( doJesStressTest ) {
     pseudodatafile = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/summary/WbWb_Slurm_Template_171_JES_095.root";
   } else {
     pseudodatafile = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/summary/WbWb_Slurm_Template_171.root";
   }
   //const TString pseudodatafile = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/summary/WbWb_Slurm_Template_171_JES_095.root";
   //const TString covariancefile = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_170_1258_matrices.root";
   //const TString covariancefile = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/summary/WbWb_Slurm_Template_171.root";
   const TString covariancefile = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/summary/matrices.root";
   int bins_number = 0;
   std::vector<int> matrix_dimension = {0}; // Needed to put covariance matrix together
   for ( auto& tmp: fit_vars ) {
     TString name = "unfolding_"+tmp+"_NOSYS";
     TH1D* tmp_data = TFile::Open(datafile)->Get<TH1D>(name);
     if ( !tmp_data ) { cerr<<"Could not find data histogram " << name <<" in file "<<datafile<<endl; exit(1);}
     bins_number += tmp_data->GetNbinsX()-1 - iRemoveBins;
     matrix_dimension.push_back(tmp_data->GetNbinsX()-1 - iRemoveBins);
     tmp_data->Clear();
   }
   
   TH1D* combined_data = new TH1D("combined_data", "combined_data", bins_number, 0, bins_number);
   int bin_offset = 0;
   if ( doPseudo ) {
     for ( auto& tmp: fit_vars_short ) {
       TH1D* tmp_data = TFile::Open(pseudodatafile)->Get<TH1D>(tmp);
       for ( int i = 1; i <= tmp_data->GetNbinsX() - iRemoveBins; i++ ) {
	 double bin_width = tmp_data->GetXaxis()->GetBinWidth(i);
	 combined_data->SetBinContent(i+bin_offset, tmp_data->GetBinContent(i)*bin_width);
	 combined_data->SetBinError(i+bin_offset, tmp_data->GetBinError(i)*bin_width);
       }
       bin_offset += tmp_data->GetNbinsX() - iRemoveBins;
     }
   }
   else {
     for ( auto& tmp: fit_vars ) {
       TString name = "unfolding_"+tmp+"_NOSYS";
       TH1D* tmp_data = TFile::Open(datafile)->Get<TH1D>(name);
       // Loop only to NbinsX-1, because last bin is overflow bin
       for ( int i = 1; i <= tmp_data->GetNbinsX()-1 - iRemoveBins; i++ ) {
	 combined_data->SetBinContent(i+bin_offset, tmp_data->GetBinContent(i)/tmp_data->GetXaxis()->GetBinWidth(i));
	 cout<<"For bin "<<i<<" add content "<<tmp_data->GetBinContent(i)<<" (data)"<<endl;
       }
       bin_offset += tmp_data->GetNbinsX()-1 - iRemoveBins;
     }
   }
   combined_data->Rebin(iRebinData);
   combined_data->SetLineColor(kBlack);
   combined_data->SetMarkerSize(1.8);
   /*   
   TH1D* data = TFile::Open(pseudodatafile)->Get<TH1D>(histname); // pseudo data, use Sherpa 3 with m_t = 170 GeV for now
   double binning[9] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0};
   //TH1D* data_tmp      = TFile::Open(pseudodatafile)->Get<TH1D>(histname); // pseudo data, use Sherpa 3 with m_t = 170 GeV for now
   //TH1D* data      = new TH1D("data", "data", 8, binning);
   //for ( int i =1; i <= data->GetNbinsX(); i++) {
   //   data->SetBinContent(i, data_tmp->GetBinContent(i));
   //   data->SetBinError(i, data_tmp->GetBinError(i));
   //}
   data->SetLineColor(kBlack);
   data->SetMarkerSize(1.8);

   if ( var_name_index == 4 || var_name_index == 5 ) {
      
      TH1D* template_1_tmp = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_155_1258.root")->Get<TH1D>(histname);
      TH1D* template_2_tmp = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_160_1256.root")->Get<TH1D>(histname);
      TH1D* template_3_tmp = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_165_1246.root")->Get<TH1D>(histname);
      TH1D* template_4_tmp = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_170_1248.root")->Get<TH1D>(histname);
      TH1D* template_5_tmp = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_175_1250.root")->Get<TH1D>(histname);
      TH1D* template_6_tmp = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_180_1252.root")->Get<TH1D>(histname);
      TH1D* template_7_tmp = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/output/Ana_S3beta_Cluster_H_mtop_185_1254.root")->Get<TH1D>(histname);
      
      TH1D* template_1   = new TH1D("template_155", "template_155", 8, binning);
      TH1D* template_2   = new TH1D("template_160", "template_160", 8, binning);
      TH1D* template_3   = new TH1D("template_165", "template_165", 8, binning);
      TH1D* template_4   = new TH1D("template_170", "template_170", 8, binning);
      TH1D* template_5   = new TH1D("template_175", "template_175", 8, binning);
      TH1D* template_6   = new TH1D("template_180", "template_180", 8, binning);
      TH1D* template_7   = new TH1D("template_185", "template_185", 8, binning);
      for ( int i =1; i <= data->GetNbinsX(); i++) {
         template_1->SetBinContent(i, template_1_tmp->GetBinContent(i));
         template_1->SetBinError(  i, template_1_tmp->GetBinError(i));
         template_2->SetBinContent(i, template_2_tmp->GetBinContent(i));
         template_2->SetBinError(  i, template_2_tmp->GetBinError(i));
         template_3->SetBinContent(i, template_3_tmp->GetBinContent(i));
         template_3->SetBinError(  i, template_3_tmp->GetBinError(i));
         template_4->SetBinContent(i, template_4_tmp->GetBinContent(i));
         template_4->SetBinError(  i, template_4_tmp->GetBinError(i));
         template_5->SetBinContent(i, template_5_tmp->GetBinContent(i));
         template_5->SetBinError(  i, template_5_tmp->GetBinError(i));
         template_6->SetBinContent(i, template_6_tmp->GetBinContent(i));
         template_6->SetBinError(  i, template_6_tmp->GetBinError(i));
         template_7->SetBinContent(i, template_7_tmp->GetBinContent(i));
         template_7->SetBinError(  i, template_7_tmp->GetBinError(i));
      }
      
      templates[155] = template_1;
      templates[160] = template_2;
      templates[165] = template_3;
      templates[170] = template_4;
      templates[175] = template_5;
      templates[180] = template_6;
      templates[185] = template_7;
   }
   else */
   {
     TH1D* combined_template_150   = new TH1D("combined_template_150", "combined_template_150", bins_number, 0, bins_number);
     TH1D* combined_template_162_5 = new TH1D("combined_template_162_5", "combined_template_162_5", bins_number, 0, bins_number);
     TH1D* combined_template_165   = new TH1D("combined_template_165", "combined_template_165", bins_number, 0, bins_number);
     TH1D* combined_template_167_5 = new TH1D("combined_template_167_5", "combined_template_167_5", bins_number, 0, bins_number);
     TH1D* combined_template_170   = new TH1D("combined_template_170", "combined_template_170", bins_number, 0, bins_number);
     TH1D* combined_template_172_5 = new TH1D("combined_template_172_5", "combined_template_172_5", bins_number, 0, bins_number);
     TH1D* combined_template_175   = new TH1D("combined_template_175", "combined_template_175", bins_number, 0, bins_number);
     TH1D* combined_template_177_5 = new TH1D("combined_template_177_5", "combined_template_177_5", bins_number, 0, bins_number);
     TH1D* combined_template_180   = new TH1D("combined_template_180", "combined_template_180", bins_number, 0, bins_number);
     TH1D* combined_template_182_5 = new TH1D("combined_template_182_5", "combined_template_182_5", bins_number, 0, bins_number);
     TH1D* combined_template_200   = new TH1D("combined_template_200", "combined_template_200", bins_number, 0, bins_number);

     int bin_offset = 0;
     for ( auto& tmp: fit_vars_short ) {
       TH1D* h_tmp_150   = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/summary/WbWb_Slurm_Template_150.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_162_5 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/summary/WbWb_Slurm_Template_162_5.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_165   = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/summary/WbWb_Slurm_Template_165.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_167_5 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/summary/WbWb_Slurm_Template_167_5.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_170   = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/summary/WbWb_Slurm_Template_170.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_172_5 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/summary/WbWb_Slurm_Template_172_5.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_175   = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/summary/WbWb_Slurm_Template_175.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_177_5 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/summary/WbWb_Slurm_Template_177_5.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_180   = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/summary/WbWb_Slurm_Template_180.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_182_5 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/summary/WbWb_Slurm_Template_182_5_new.root")->Get<TH1D>(tmp);
       TH1D* h_tmp_200   = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/summary/WbWb_Slurm_Template_200.root")->Get<TH1D>(tmp);

       double kFactor =  1.0; // Define k-factor when comparing templates to data
       for ( int i = 1; i<= h_tmp_165->GetNbinsX() - iRemoveBins; i++ ) {
	 double bin_width_x = h_tmp_165->GetXaxis()->GetBinWidth(i);
         combined_template_150->SetBinContent(  i+bin_offset, h_tmp_150->GetBinContent(i)  * kFactor*bin_width_x);
         combined_template_162_5->SetBinContent(i+bin_offset, h_tmp_162_5->GetBinContent(i)* kFactor*bin_width_x);
	 combined_template_165->SetBinContent(  i+bin_offset, h_tmp_165->GetBinContent(i)  * kFactor*bin_width_x);
         combined_template_167_5->SetBinContent(i+bin_offset, h_tmp_167_5->GetBinContent(i)* kFactor*bin_width_x);
	 combined_template_170->SetBinContent(  i+bin_offset, h_tmp_170->GetBinContent(i)  * kFactor*bin_width_x);
	 combined_template_172_5->SetBinContent(i+bin_offset, h_tmp_172_5->GetBinContent(i)* kFactor*bin_width_x);
	 combined_template_175->SetBinContent(  i+bin_offset, h_tmp_175->GetBinContent(i)  * kFactor*bin_width_x);
         combined_template_177_5->SetBinContent(i+bin_offset, h_tmp_177_5->GetBinContent(i)* kFactor*bin_width_x);
         combined_template_180->SetBinContent(  i+bin_offset, h_tmp_180->GetBinContent(i)  * kFactor*bin_width_x);
         combined_template_182_5->SetBinContent(i+bin_offset, h_tmp_182_5->GetBinContent(i)* kFactor*bin_width_x);
         combined_template_200->SetBinContent(  i+bin_offset, h_tmp_200->GetBinContent(i)  * kFactor*bin_width_x);

         combined_template_150->SetBinError(   i+bin_offset, h_tmp_150->GetBinError(i)   * kFactor * bin_width_x);
	 combined_template_162_5->SetBinError( i+bin_offset, h_tmp_162_5->GetBinError(i) * kFactor * bin_width_x);
	 combined_template_165->SetBinError(   i+bin_offset, h_tmp_165->GetBinError(i)   * kFactor * bin_width_x);
         combined_template_167_5->SetBinError( i+bin_offset, h_tmp_167_5->GetBinError(i) * kFactor * bin_width_x);
	 combined_template_170->SetBinError(   i+bin_offset, h_tmp_170->GetBinError(i)   * kFactor * bin_width_x);
	 combined_template_172_5->SetBinError( i+bin_offset, h_tmp_172_5->GetBinError(i) * kFactor * bin_width_x);
         combined_template_175->SetBinError(   i+bin_offset, h_tmp_175->GetBinError(i)   * kFactor * bin_width_x);
	 combined_template_177_5->SetBinError( i+bin_offset, h_tmp_177_5->GetBinError(i) * kFactor * bin_width_x);
         combined_template_180->SetBinError(   i+bin_offset, h_tmp_180->GetBinError(i)   * kFactor * bin_width_x);
         combined_template_182_5->SetBinError( i+bin_offset, h_tmp_182_5->GetBinError(i) * kFactor * bin_width_x);
	 combined_template_200->SetBinError(   i+bin_offset, h_tmp_200->GetBinError(i)   * kFactor * bin_width_x);
       }
       bin_offset += h_tmp_165->GetNbinsX() - iRemoveBins;
     }
     //templates[150] = combined_template_150;
     templates[162.5] = combined_template_162_5;
     templates[165] = combined_template_165;
     templates[167.5] = combined_template_167_5;
     templates[170] = combined_template_170;
     templates[172.5] = combined_template_172_5;
     templates[175] = combined_template_175;
     templates[177.5] = combined_template_177_5;
     templates[180] = combined_template_180;
     templates[182.5] = combined_template_182_5;
     //templates[200] = combined_template_200;

   }

   for ( auto [MM,hist] : templates ) {
      hist->Rebin(iRebin);
   }

   
   // Add predictions to see if chisquare is as expected
   /*
   const TString aMCatNLO_ttbar = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/ttbar_enhanced_aMcAtNlo_fast.Theory.root";
   const TString aMCatNLO_single = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/singletop_inclusive_Wt_DR_aMcAtNlo_fast.Theory.root";
   const TString pythia_single_DR = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/singletop_inclusive_Wt_DR.Theory.root";
   const TString MiNNLO_ttbar_alt = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/ttbar_enhanced_MiNNLO_alt.Theory.root";
   const TString MiNNLO_ttbar = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/examples/data/ttbar_enhanced_MiNNLO.Theory.root";
   {
     constexpr double LuminosityMC16a = 3244.54 + 33402.2;
     constexpr double LuminosityMC16d = 44630.6;
     constexpr double LuminosityMC16e = 58791.6;
     constexpr double LuminosityFull = LuminosityMC16a + LuminosityMC16d + LuminosityMC16e;
     constexpr double LuminosityInverse = 1.0 / LuminosityFull;

     TH1D* combined_template_190 = new TH1D("combined_template_190", "combined_template_190", bins_number, 0, bins_number);
     TH1D* h_tmp_aMCNLO_ttbar  = TFile::Open(aMCatNLO_ttbar)->Get<TH1D>( "THEORY_NOSYS_signal/hist/l_Whad_particle_"+fit_vars[0]);
     TH1D* h_tmp_aMCNLO_single = TFile::Open(aMCatNLO_single)->Get<TH1D>("THEORY_NOSYS_signal/hist/l_Whad_particle_"+fit_vars[0]);
     h_tmp_aMCNLO_ttbar->Add(h_tmp_aMCNLO_single);
     
     TH1D* combined_template_195 = new TH1D("combined_template_195", "combined_template_195", bins_number, 0, bins_number);
     TH1D* h_tmp_MiNNLO_ttbar  = TFile::Open(MiNNLO_ttbar)->Get<TH1D>( "THEORY_NOSYS_signal/hist/l_Whad_particle_"+fit_vars[0]);
     TH1D* h_tmp_pythia_DR_single = TFile::Open(pythia_single_DR)->Get<TH1D>("THEORY_NOSYS_signal/hist/l_Whad_particle_"+fit_vars[0]);
     h_tmp_MiNNLO_ttbar->Add(h_tmp_pythia_DR_single);

     TH1D* combined_template_200 = new TH1D("combined_template_200", "combined_template_200", bins_number, 0, bins_number);
     TH1D* h_tmp_MiNNLO_ttbar_alt  = TFile::Open(MiNNLO_ttbar_alt)->Get<TH1D>( "THEORY_NOSYS_signal/hist/l_Whad_particle_"+fit_vars[0]);
     h_tmp_MiNNLO_ttbar_alt->Add(h_tmp_pythia_DR_single);
     
     for ( int i = 1; i <= h_tmp_aMCNLO_ttbar->GetNbinsX(); i++ ) {
       combined_template_190->SetBinContent(i, LuminosityInverse*h_tmp_aMCNLO_ttbar->GetBinContent(i) / h_tmp_aMCNLO_ttbar->GetXaxis()->GetBinWidth(i));
       combined_template_195->SetBinContent(i, LuminosityInverse*h_tmp_MiNNLO_ttbar->GetBinContent(i) / h_tmp_MiNNLO_ttbar->GetXaxis()->GetBinWidth(i));
       combined_template_200->SetBinContent(i, LuminosityInverse*h_tmp_MiNNLO_ttbar_alt->GetBinContent(i) / h_tmp_MiNNLO_ttbar_alt->GetXaxis()->GetBinWidth(i));

     }
     templates[190] = combined_template_190;
     templates[195] = combined_template_195;
     templates[200] = combined_template_200;

   }
   */
   
   // ------------------------------------------------ //
   // --- List of uncertainties
   // ------------------------------------------------ //
   vector<string> uncertainties = {
     "EG_RESOLUTION_ALL", "EG_SCALE_ALL", "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR", "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR",
     "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR", "EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR",
     "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR", "MUON_SAGITTA_DATASTAT", "MUON_SAGITTA_RESBIAS",
     "MUON_EFF_BADMUON_SYS", "MUON_EFF_ISO_STAT", "MUON_EFF_ISO_SYS", "MUON_EFF_RECO_STAT",
     "MUON_EFF_RECO_SYS", "MUON_EFF_TTVA_STAT", "MUON_EFF_TTVA_SYS", "MUON_EFF_TrigStatUncertainty",
     "MUON_EFF_TrigSystUncertainty", // Lepton uncertainties
     "JET_EffectiveNP_Detector1", "JET_EffectiveNP_Detector2", "JET_EffectiveNP_Mixed1", "JET_EffectiveNP_Mixed2",
     "JET_EffectiveNP_Mixed3", "JET_EffectiveNP_Modelling1", "JET_EffectiveNP_Modelling2",
     "JET_EffectiveNP_Modelling3", "JET_EffectiveNP_Modelling4", "JET_EffectiveNP_Statistical1",
     "JET_EffectiveNP_Statistical2", "JET_EffectiveNP_Statistical3", "JET_EffectiveNP_Statistical4",
     "JET_EffectiveNP_Statistical5", "JET_EffectiveNP_Statistical6", "JET_EtaIntercalibration_Modelling",
     "JET_EtaIntercalibration_NonClosure_2018data", "JET_EtaIntercalibration_NonClosure_highE",
     "JET_EtaIntercalibration_NonClosure_negEta", "JET_EtaIntercalibration_NonClosure_posEta",
     "JET_EtaIntercalibration_TotalStat", "JET_Flavor_Composition_prop", "JET_Flavor_Response_prop",
     "JET_Pileup_OffsetMu", "JET_Pileup_OffsetNPV", "JET_Pileup_PtTerm", "JET_Pileup_RhoTopology",
     "JET_PunchThrough_MC16", // JES uncertainty
     "JET_JER_DataVsMC_MC16_PseudoData", "JET_JER_EffectiveNP_1_PseudoData", "JET_JER_EffectiveNP_2_PseudoData",
     "JET_JER_EffectiveNP_3_PseudoData", "JET_JER_EffectiveNP_4_PseudoData", "JET_JER_EffectiveNP_5_PseudoData",
     "JET_JER_EffectiveNP_6_PseudoData", "JET_JER_EffectiveNP_7_PseudoData", "JET_JER_EffectiveNP_8_PseudoData",
     "JET_JER_EffectiveNP_9_PseudoData", "JET_JER_EffectiveNP_10_PseudoData", "JET_JER_EffectiveNP_11_PseudoData",
     "JET_JER_EffectiveNP_12restTerm_PseudoData", // JER uncertainty
     "JET_JvtEfficiency", "PRW_DATASF", "MET_SoftTrk_ResoPara", "MET_SoftTrk_ResoPerp", "MET_SoftTrk_Scale",
     "WMASS_VAR_signal", //MET+JVT+PileUp
     "FT_EFF_Eigen_B_0", "FT_EFF_Eigen_B_1", "FT_EFF_Eigen_B_2", "FT_EFF_Eigen_B_3", "FT_EFF_Eigen_B_4",
     "FT_EFF_Eigen_B_5", "FT_EFF_Eigen_B_6", "FT_EFF_Eigen_B_7", "FT_EFF_Eigen_B_8", "FT_EFF_Eigen_C_0",
     "FT_EFF_Eigen_C_1", "FT_EFF_Eigen_C_2", "FT_EFF_Eigen_C_3", "FT_EFF_Eigen_Light_0", "FT_EFF_Eigen_Light_1",
     "FT_EFF_Eigen_Light_2", "FT_EFF_Eigen_Light_3", "FT_EFF_extrapolation", "FT_EFF_extrapolation_from_charm", // b-tagging uncertainties
     "THEORY_CROSS_SECTION_signal", "THEORY_SHOWERING_HERWIG7_signal", "THEORY_SCALE_FACTORISATION_signal",
     "THEORY_SCALE_RENORMALISATION_signal", "THEORY_ISR_signal", "THEORY_FSR_signal", "THEORY_HDAMP_signal",
     "THEORY_PTHARD_signal", "THEORY_TOPRECOIL_signal", "THEORY_TOP_MASS_signal",
     "THEORY_PDF4LHC_VARIATION_signal", "THEORY_DR_DS_signal", // modelling uncertainties
     "THEORY_CROSS_SECTION_Wjets", "THEORY_SCALE_COMBINED_Wjets", "THEORY_PDF4LHC_VARIATION_Wjets", "THEORY_EWK_Wjets",
     "THEORY_SCALE_COMBINED_multiboson_noW", "THEORY_PDF4LHC_VARIATION_multiboson_noW", "THEORY_EWK_multiboson_noW",
     "THEORY_CROSS_SECTION_other_top_noWt", "FAKES_Electron", "FAKES_Muon", // bkgd uncertainties
     "LUMINOSITY" // lumi uncertainty
   };
   if ( doJesStressTest ) uncertainties.push_back("myJESUncertainty");
   
   
   vector<string> statistical_uncertainties = {"STAT_MC"};
   
   vector<string> external_uncertainties = {"FULL_SYS_SUM", "FULL_SYS_SUM_DETECTOR", "FULL_SYS_SUM_THEORY",
					    "TOTAL_SYSONLY", "TOTAL", "TOTAL_NO_DR_DS", "FULL_SYS_TOYS"};

   // ------------------------------------------------ //
   // ---  Do linear template fit
   // ------------------------------------------------ //
   // --- instantiate LTF object
   LTF ltf;
   ltf.SetGamma(vector<double>{1});
   ltf.UseNuisanceParameters(true);
   ltf.UseLogNormalUncertainties(false);// Johannes: Effect on fit is small for true or false

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

   if ( !doPseudo ) {
     vector<double> combined_error;
     for (int i=1; i<= combined_data->GetNbinsX(); i++ ) {
       combined_error.push_back(combined_data->GetBinError(i)/combined_data->GetBinContent(i));
     }
     ltf.AddErrorRelative("pseudoDataStat", combined_error, 0.0, LTF::Uncertainty::Constrained);
   }
   else {
   // Start building covariance matrix
   // Get data covariance matrix after unfolding from TUnfold
   // Covariance defined as rho_ij*sigma_i*sigma_j
   // sigma is the absolute error

   // In order to plot the correlation matrix in an interactive root session, just do
   // auto htext = (TH2D*)_file0->Get("combined_corr")
   // gStyle->SetPaintTextFormat("4.2f");
   // htext->SetMarkerSize(1.0);
   // htext->Draw("TEXT colz")
   TH2D* combined_covariance = new TH2D("combined_cov", "combined_cov", bins_number, 0, bins_number, bins_number, 0, bins_number);
   TH2D* combined_correlation = new TH2D("combined_corr", "combined_corr", bins_number, 0, bins_number, bins_number, 0, bins_number);
   bin_offset = 0;
   for ( auto& fit_var: fit_vars ) {
      TString histnameCovStat("unfolding_covariance_matrix_"+fit_var+"_covariance_STAT_DATA");

      TH2D* cov_stat_dat = file->Get<TH2D>(histnameCovStat);
      if ( !cov_stat_dat ) { cerr<<"Could not find covariance matrix " << histnameCovStat <<endl; exit(1);}
      else cout<<"Found covarinace matrix "<<histnameCovStat<<endl;
      TH1D* h_err  = file->Get<TH1D>("unfolding_error_"+fit_var+"_direct_envelope_STAT_DATA__1up");
      TH1D* h_data = file->Get<TH1D>("unfolding_"+fit_var+"_NOSYS");
      for ( int i = 1; i < cov_stat_dat->GetNbinsX() - iRemoveBins; i++ ) {
	double bin_width_x = h_err->GetXaxis()->GetBinWidth(i);
	double sigma_data_x = sqrt(cov_stat_dat->GetBinContent(i,i)); //h_data->GetBinContent(i) * h_err->GetBinContent(i);

	for ( int j = 1; j < cov_stat_dat->GetNbinsY() - iRemoveBins; j++ ) {
	  double bin_width_y = h_err->GetXaxis()->GetBinWidth(j);
	  double sigma_data_y = sqrt(cov_stat_dat->GetBinContent(j,j)); //h_data->GetBinContent(j) * h_err->GetBinContent(j);
	  combined_covariance->SetBinContent(i+bin_offset,j+bin_offset, cov_stat_dat->GetBinContent(i,j));// Divide by bin width, if templates are not divided by bin width!
	  combined_correlation->SetBinContent(i+bin_offset,j+bin_offset, cov_stat_dat->GetBinContent(i,j)/(sigma_data_x*sigma_data_y));
	}
      }
      bin_offset += cov_stat_dat->GetNbinsX() - 1 - iRemoveBins;
   }

   // Fill off-diagonal block matrices of the covariance matrix
   if ( fit_vars.size() > 1 ) {
     for ( int v1 = 0; v1 < fit_vars_short.size(); v1++ ) {
       TH1D* h_err_v1 = file->Get<TH1D>("unfolding_error_"+fit_vars[v1]+"_direct_envelope_STAT_DATA__1up");
       TH1D* h_data_var1 = file->Get<TH1D>("unfolding_"+fit_vars[v1]+"_NOSYS");
       if ( !h_data_var1 ) { cerr<<"Could not find data for " << fit_vars[v1] <<endl; exit(1);}
       for ( int v2 = v1+1; v2 < fit_vars_short.size(); v2++ ) {
	 TString histname = fit_vars_short[v1]+"_"+fit_vars_short[v2];
	 TH2D* cov = TFile::Open(covariancefile)->Get<TH2D>(histname);
	 if ( !cov ) { cerr<<"Could not find covariance matrix " << histname <<endl; exit(1);}
	 else cout<<"Found covarinace matrix "<<histname<<endl;
	 TH1D* h_err_v2 = file->Get<TH1D>("unfolding_error_"+fit_vars[v2]+"_direct_envelope_STAT_DATA__1up"); // get relative error
	 //TH1D* hist_data = ... //get cross section in single bin // might have to divide by bin width later
	 // Projections of correlation matrix
	 TH1D* projection_var1 = cov->ProjectionX("pro_v1",0,-1,"e");
	 TH1D* projection_var2 = cov->ProjectionY("pro_v2",0,-1,"e");
	 // Cross section file
	 TH1D* h_data_var2 = file->Get<TH1D>("unfolding_"+fit_vars[v2]+"_NOSYS");
	 if ( !h_data_var2 ) { cerr<<"Could not find data for " << fit_vars[v2] <<endl; exit(1);}
	 //cout<<"Var1 content ";
	 //for ( int i = 1; i <= h_data_var1->GetNbinsX(); i++ ) cout<<projection_var1->GetBinContent(i) / h_data_var1->GetBinContent(i)<<"\t";
	 //cout<<endl<<"Var1 error   ";
	 //for ( int i = 1; i <= h_data_var1->GetNbinsX(); i++ ) cout<<projection_var1->GetBinError(i) / h_data_var1->GetBinError(i)<<"\t";
	 //cout<<endl<<"Var2 content ";
	 //for ( int i = 1; i <= h_data_var2->GetNbinsX(); i++ ) cout<<projection_var2->GetBinContent(i) / h_data_var2->GetBinContent(i)<<"\t";
	 //cout<<endl<<"Var2 error   ";
	 //for ( int i = 1; i <= h_data_var2->GetNbinsX(); i++ ) cout<<projection_var2->GetBinError(i) / h_data_var2->GetBinError(i)<<"\t";
	 //cout<<endl;
	 int bin_offset_x = 0;
	 int bin_offset_y = 0;
	 for ( int v1_tmp = 0; v1_tmp <= v1; v1_tmp++ ) bin_offset_x += matrix_dimension[v1_tmp];
	 for ( int v2_tmp = 0; v2_tmp <= v2; v2_tmp++ ) bin_offset_y += matrix_dimension[v2_tmp];

	 TString histnameCovStatX("unfolding_covariance_matrix_"+fit_vars[v1]+"_covariance_STAT_DATA");
	 TH2D* cov_stat_dat_x = file->Get<TH2D>(histnameCovStatX);
	 TString histnameCovStatY("unfolding_covariance_matrix_"+fit_vars[v2]+"_covariance_STAT_DATA");
         TH2D* cov_stat_dat_y = file->Get<TH2D>(histnameCovStatY);
	 for ( int i = 1; i < h_err_v1->GetNbinsX() - iRemoveBins; i++ ) {
	   for ( int j = 1; j < h_err_v2->GetNbinsX() - iRemoveBins; j++ ) {
	     // Calculate the covariance in the data from the one in Sherpa
	     // 1. Get the Sherpa covariance: std::pow(cov->GetBinError(i,j),2)
	     // 2. Divide by the std deviation from Sherpa: sigma_template1, sigma_template2
	     // 3. Multiply with the std deviation from the data: sigma_data1, sigma_data2
	     double sigma_data1 = sqrt(cov_stat_dat_x->GetBinContent(i,i)); //h_data_var1->GetBinContent(i) * h_err_v1->GetBinContent(i);
	     double sigma_data2 = sqrt(cov_stat_dat_y->GetBinContent(j,j)); //h_data_var2->GetBinContent(j) * h_err_v2->GetBinContent(j);
	     double sigma_template1 = projection_var1->GetBinError(i);
	     double sigma_template2 = projection_var2->GetBinError(j);
	     //cout<<"Correlation "<<std::pow(cov->GetBinError(i,j),2)/ sigma_template1 / sigma_template2<<endl;
	     combined_covariance->SetBinContent(i+bin_offset_x, j+bin_offset_y, std::pow(cov->GetBinError(i,j),2)*(sigma_data1*sigma_data2) / (sigma_template1*sigma_template2));
	     combined_covariance->SetBinContent(j+bin_offset_y, i+bin_offset_x, std::pow(cov->GetBinError(i,j),2)*(sigma_data1*sigma_data2) / (sigma_template1*sigma_template2));
	     combined_correlation->SetBinContent(i+bin_offset_x, j+bin_offset_y, std::pow(cov->GetBinError(i,j),2) / (sigma_template1*sigma_template2));
             combined_correlation->SetBinContent(j+bin_offset_y, i+bin_offset_x, std::pow(cov->GetBinError(i,j),2) / (sigma_template1*sigma_template2));
	   }
	 }
       }
     }
   }
   const bool writeCovarianceMatrix = false;
   if ( writeCovarianceMatrix ) {
     TString name_outfile = "correlation";
     for ( auto tmp: fit_vars_short ) name_outfile += tmp;
     TFile *outfile = new TFile(name_outfile+".root", "RECREATE");
     combined_covariance->Write();
     combined_correlation->Write();
     outfile->Close();
     outfile->Delete();
   }

   cout<<"Johannes print final covariance matrix"<<endl;
   vector<vector<double > > vecCov2 = TH2D_to_vecvec(combined_covariance);
   for ( auto& tmp_vec: vecCov2 ){
     for ( auto& tmp: tmp_vec ) cout<<tmp<<"\t";
     cout<<endl;
   }
   ltf.AddError("STAT_DATA", vecCov2, LTF::Uncertainty::Constrained);

   cout<<"Johannes print final correltaion matrix"<<endl;
   vector<vector<double > > vecCor2 = TH2D_to_vecvec(combined_correlation);
   for ( auto& tmp_vec: vecCor2 ){
     for ( auto& tmp: tmp_vec ) cout<<tmp<<"\t";
     cout<<endl;
   }
  
   TH1D* total_error = new TH1D("total_error", "total_error", bins_number, 0, bins_number);

   // Systematical uncertainties
   for ( auto& uncertainty: uncertainties ) {
     vector<double> combined_error;

     if ( doJesStressTest && uncertainty == "myJESUncertainty" ) {
       cout<<"Now adding new uncertainty"<<endl;
       const TString datafile_nominal = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/summary/WbWb_Slurm_Template_171.root";
       for ( auto& tmp: fit_vars_short ) {
	 TH1D* tmp_data = TFile::Open(datafile_nominal)->Get<TH1D>(tmp);
	 cout<<"Retrieved the file"<<endl;
	 tmp_data->Print("All");
	 for ( int i = 1; i <= tmp_data->GetNbinsX() - iRemoveBins; i++ ) {
	   double error = (tmp_data->GetBinContent(i)-combined_data->GetBinContent(i))/combined_data->GetBinContent(i);
	   cout<<"Adding err "<<error<<" = "<<tmp_data->GetBinContent(i)<<" / "<<combined_data->GetBinContent(i)<<endl;
	   combined_error.push_back(error);
	   total_error->SetBinContent(i, std::sqrt(std::pow(total_error->GetBinContent(i),2)+std::pow(error,2)));
	 }
       }
       cout<<"Added new uncertainty"<<endl;
     }
     else {
       for ( auto& fit_variable: fit_vars ) {
	 TH1D* hist_up   = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_"+uncertainty+"__1up");
	 //TH1D* hist_down = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_"+uncertainty+"__1down");
	 for (int i=1; i< hist_up->GetNbinsX() - iRemoveBins; i++) {
	   //combined_error.push_back(std::max(abs(hist_up->GetBinContent(i)), abs(hist_down->GetBinContent(i))));
	   //total_error->SetBinContent(i, std::sqrt(std::pow(total_error->GetBinContent(i),2)+std::pow(std::max(abs(hist_up->GetBinContent(i)),abs(hist_down->GetBinContent(i))),2)));
	   combined_error.push_back(hist_up->GetBinContent(i));
	   total_error->SetBinContent(i, std::sqrt(std::pow(total_error->GetBinContent(i),2)+std::pow(hist_up->GetBinContent(i),2))); // Use up variation
	   //combined_error.push_back(hist_down->GetBinContent(i));
	   //total_error->SetBinContent(i, std::sqrt(std::pow(total_error->GetBinContent(i),2)+std::pow(hist_down->GetBinContent(i),2))); // Use down variation
	 }
       }
     }
     double corr = 1.0;
     if ( combined_error.size() > 0 ) ltf.AddErrorRelative(uncertainty, combined_error, corr, LTF::Uncertainty::Constrained);
     combined_error.clear();
   }
   total_error->Print("All");
   
   // Johannes added this unconstrained error to correct for offset in normalization 
   {
     vector<double> combined_error;
     for ( auto& fit_variable: fit_vars ) {
       TH1D* hist_up   = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_"+uncertainties[0]+"__1up");
       for (int i=1; i< hist_up->GetNbinsX() - iRemoveBins; i++) {
	 combined_error.push_back(0.1);
	 //total_error->SetBinContent(i, std::sqrt(std::pow(total_error->GetBinContent(i),2)+std::pow(0.1,2)));
       }
     }
     //if ( combined_error.size() > 0 ) ltf.AddErrorRelative("UnconstrainedError", combined_error, 1.0, LTF::Uncertainty::Unconstrained);
   }
   // Stat uncertainties
   for ( auto& uncertainty: statistical_uncertainties ) { // Only stat. unc. is MC stat. 
     vector<double> combined_error;
     for ( auto& fit_variable: fit_vars ) {
       TH1D* hist_up = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_"+uncertainty+"__1up");
       //TH1D* hist_down = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_"+uncertainty+"__1down");
       for (int i=1; i< hist_up->GetNbinsX() - iRemoveBins; i++) {
	 combined_error.push_back(hist_up->GetBinContent(i));
         total_error->SetBinContent(i, std::sqrt(std::pow(total_error->GetBinContent(i),2)+std::pow(hist_up->GetBinContent(i),2)));
	 cout<<uncertainty<<" "<<hist_up->GetBinContent(i)<<endl;
       }
     }
     double corr = 0.0; // Stat. unc. are uncorrelated
     if ( combined_error.size() > 0 ) ltf.AddErrorRelative(uncertainty, combined_error, corr, LTF::Uncertainty::Constrained);
     combined_error.clear();
   }
   total_error->Print("All");
   
   // When using pseudo-data, add the stat. uncertainties of the pseudo data as an additional uncertainty
   if ( doPseudo ) {
     vector<double> combined_error;
     for (int i=1; i<= combined_data->GetNbinsX(); i++ ) {
       combined_error.push_back(combined_data->GetBinError(i)/combined_data->GetBinContent(i));
       total_error->SetBinContent(i, std::sqrt(std::pow(total_error->GetBinContent(i),2)+std::pow(combined_data->GetBinError(i)/combined_data->GetBinContent(i),2)));

     }
     ltf.AddErrorRelative("pseudoDataStat", combined_error, 0.0, LTF::Uncertainty::Constrained);
   }
   total_error->Print("All");
   
   // External uncertainties
   for ( auto& uncertainty: external_uncertainties ) {
     vector<double> combined_error;
     for ( auto& fit_variable: fit_vars ) {
       TH1D* hist = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_"+uncertainty+"__1up");
       for (int i=1; i< hist->GetNbinsX() - iRemoveBins; i++) {
         combined_error.push_back(hist->GetBinContent(i));
       }
     }
     double corr = 1;
     if ( combined_error.size() > 0 ) ltf.AddErrorRelative(uncertainty, combined_error, corr, LTF::Uncertainty::External);
     combined_error.clear();
   }
   }
   PrintAsciiTable(templates,combined_data);
   
   LTF::LiTeFit fit = ltf.DoLiTeFit();
   //LTF::LiTeFit fit = ltf.DoQuadraticTemplateFit(3);
   fit.PrintFull();

   //fit.DoIterativeFitNewton(6,0.6,2,1);
   //fit.DoIterativeFitTaylor();
   //fit.PrintFull();
   vector<double> bins{};
   for(int i = 1; i <= combined_data->GetNbinsX(); i++) {
      bins.push_back(combined_data->GetBinLowEdge(i));
   }
   bins.push_back(combined_data->GetXaxis()->GetBinUpEdge(combined_data->GetNbinsX()));
   string label = "";
   for (TString tmp: fit_vars_short) label += tmp+"\t"; 
   LTF_ROOTTools::plotLiTeFit(fit, bins, ps_name, "d#sigma/dx [pb]", label,"m_{t} [GeV]");

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
      printf(" %11.6f %11.6f",data->GetBinContent(i),data->GetBinError(i));
   for ( auto [mean,hist] : templates ) 
     printf(" %11.6f %11.6f",hist->GetBinContent(i),hist->GetBinError(i));// johannes this was 11.4
      cout<<endl;
   }
   cout<<endl;
}

