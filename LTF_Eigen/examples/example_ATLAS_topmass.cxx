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
#include <TMatrixD.h>

void PrintAsciiTable(const map<double,TH1D*>&, TH1D* data);

int fitMultipleObservables(const char* ps_name, const vector<TString> fit_vars, const vector<TString> fit_vars_short, int excl_mbl, int excl_mbw, int excl_yw, int excl_mwbbl, int excl_extra);

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

  TFile *file = new TFile("fit_quality.root", "RECREATE");
  TH1D *h_chisq_prob_linear = new TH1D("h_chisq_prob_linear", "h_chisq_prob_linear", 10, 0.0, 1.0);
  TH1D *h_chisq_prob_quadratic = new TH1D("h_chisq_prob_quadratic", "h_chisq_prob_quadratic", 10, 0.0, 1.0);
  TH1D *h_chisq_prob_ratio = new TH1D("h_chisq_prob_ratio", "h_chisq_prob_ratio", 15, 0.0, 1.5);
  TH1D *h_chisq_prob_ratio_rel = new TH1D("h_chisq_prob_ratio_rel", "h_chisq_prob_ratio_rel", 25, -1.05, 1.45);
  TH1D *h_chisq_ratio = new TH1D("h_chisq_ratio", "h_chisq_ratio", 15, 0.0, 1.5);
  TH1D *h_cheb2_sign = new TH1D("h_cheb2_sign", "h_cheb2_sign", 15, 0.0, 3);
  TH1D *h_cheb3_sign = new TH1D("h_cheb3_sign", "h_cheb3_sign", 15, 0.0, 3);
  TH1D *h_cheb_sign_ratio = new TH1D("h_cheb_sign_ratio", "h_cheb_sign_ratio", 15, 0.0, 1.5);
  TH2D *h_cheb2_sign_cheb3_sign = new TH2D("h_cheb2_sign_cheb3_sign", "h_cheb2_sign_cheb3_sign", 10, 0.0, 3.0, 10, 0.0, 3);
  TH2D *h_chisq_ratio_cheb2_sign = new TH2D("h_chisq_ratio_cheb2_sign", "h_chisq_ratio_cheb2_sign", 10, 0.0, 1.5, 10, 0.0, 3);
  
  h_chisq_prob_linear->Write();
  h_chisq_prob_quadratic->Write();
  h_chisq_prob_ratio->Write();
  h_chisq_prob_ratio_rel->Write();
  h_chisq_ratio->Write();
  h_cheb2_sign->Write();
  h_cheb3_sign->Write();
  h_cheb_sign_ratio->Write();
  h_cheb2_sign_cheb3_sign->Write();
  h_chisq_ratio_cheb2_sign->Write();
  file->Close();

//  if (fitMultipleObservables("plots/fit_all.ps", {"mbl_selected", "mbwhad_selected", "dRbl_selected", "dRbwhad_selected", "ptl1", "ptwhad", "rbq", "rbj"}, {"m_bl", "m_bw", "dr_bl", "dr_bw", "pT_lep1", "pT_whad", "r_bq", "r_bj"}, 0, 0) > 0) return 1;
//  
//
//  if (fitMultipleObservables("plots/fit_drbl.ps",     {"dRbl_selected"},    {"dr_bl"}    , 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_drbw.ps",     {"dRbwhad_selected"}, {"dr_bw"}    , 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_etal1.ps",    {"etal1"},            {"eta_lep1"} , 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_mtlep.ps",    {"mtlep"},            {"mT_lep1"}  , 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_mbl.ps",      {"mbl_selected"},     {"m_bl"}     , 0, 0, 0, 0, 0) > 0) return 1; // 17 bins
//  if (fitMultipleObservables("plots/fit_mbw.ps",      {"mbwhad_selected"},  {"m_bw"}     , 0, 0, 0, 0, 0) > 0) return 1; // 11 bins
//  if (fitMultipleObservables("plots/fit_minimax.ps",  {"minimax_whadbbl"},  {"m_minimax"}, 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_mwbbl.ps",    {"mwhadbbl"},         {"m_wbbl"}   , 0, 0, 0, 0, 0) > 0) return 1; // 10 bins
//  if (fitMultipleObservables("plots/fit_mblmin.ps",   {"mblmin"},           {"mblmin"}   , 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_mtblmin.ps",  {"mtblmin"},          {"mtblmin"}  , 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptwbbl.ps",   {"ptwhadbbl"},        {"pT_Wbbl"}  , 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptwbblnu.ps", {"ptwhadbblnu"},      {"pT_Wbblnu"}, 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptb1.ps",     {"ptb1"},             {"pT_bjet1"} , 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptb2.ps",     {"ptb2"},             {"pT_bjet2"} , 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptbl.ps",     {"ptbl_selected"},    {"pT_bl"}    , 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptblnu.ps",   {"ptblnu_selected"},  {"pT_blnu"}  , 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptbw.ps",     {"ptbwhad_selected"}, {"pT_bw"}    , 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptl1.ps",     {"ptl1"},             {"pT_lep1"}  , 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptmet.ps",    {"met"},              {"pT_met"}   , 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptowj1.ps",   {"ptOWj1"},           {"pT_owj1"}  , 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptowj2.ps",   {"ptOWj2"},           {"pT_owj2"}  , 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_ptwhad.ps",   {"ptwhad"},           {"pT_whad"}  , 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_rbj.ps",      {"rbj"},              {"r_bj"}     , 0, 0, 0, 0, 0) > 0) return 1; // 7 bins
//  if (fitMultipleObservables("plots/fit_rbq.ps",      {"rbq"},              {"r_bq"}     , 0, 0, 0, 0, 0) > 0) return 1; // 7 bins

  if (fitMultipleObservables("plots_data/fit_dataratio_ywhad.ps",            {"rapiditywhad"},     {"y_whad"}   , 0, 0, 0, 0, 0) > 0) return 1; // 24 bins
  //if (fitMultipleObservables("plots_data/fit_dataratio_mbl_1_11.ps",         {"mbl_selected"},     {"m_bl"}     , 6, 0, 0, 0, 0) > 0) return 1;
  //if (fitMultipleObservables("plots_data/fit_dataratio_mbw_1_8.ps",          {"mbwhad_selected"},  {"m_bw"}     , 0, 3, 0, 0, 0) > 0) return 1;
  //if (fitMultipleObservables("plots_data/fit_dataratio_mwbbl_1_8.ps",        {"mwhadbbl"},         {"m_wbbl"}   , 0, 0, 0, 2, 0) > 0) return 1;
  //if (fitMultipleObservables("plots_data/fit_dataratio_mbl_mbw.ps",          {"mbl_selected", "mbwhad_selected"}, {"m_bl", "m_bw"}, 6, 3, 0, 2, 0) > 0) return 1;
  //if (fitMultipleObservables("plots_data/fit_dataratio_mbl_mbw_yw.ps",       {"mbl_selected", "mbwhad_selected", "rapiditywhad"}, {"m_bl", "m_bw", "y_whad"}, 6, 3, 0, 2, 0) > 0) return 1;
  if (fitMultipleObservables("plots_data/fit_dataratio_mbl_mbw_mwbbl_yw.ps", {"mbl_selected", "mbwhad_selected", "mwhadbbl", "rapiditywhad"}, {"m_bl", "m_bw", "m_wbbl", "y_whad"}, 6, 3, 0, 2, 0) > 0) return 1;

//  if (fitMultipleObservables("plots/fit_ptl1_asimov.ps",               {"ptl1"}, {"pT_lep1"}, 0, 0, 0, 0, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_mbl_mbw_ptl_asimov.ps",        {"mbl_selected", "mbwhad_selected", "ptl1"},{"m_bl", "m_bw", "pT_lep1"}, 6, 3, 0, 2, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_mbl_mbw_yw_ptl_asimov.ps",     {"mbl_selected", "mbwhad_selected", "ptl1", "rapiditywhad"},{"m_bl", "m_bw", "pT_lep1", "y_whad"}, 6, 3, 0, 2, 0) > 0) return 1;
//  if (fitMultipleObservables("plots/fit_mbl_mbw_yw_mwbbl_ptl_asimov.ps",  {"mbl_selected", "mbwhad_selected", "ptl1", "mwhadbbl", "rapiditywhad"},{"m_bl", "m_bw", "pT_lep1", "m_wbbl", "y_whad"}, 6, 3, 0, 2, 0) > 0) return 1;

  
//  for ( int i = 4; i < 8; i ++) {
//    for ( int j = 2; j < 6; j++) {
//      for ( int k = 0; k < 4; k++) {
//	for ( int l = 0; l < 2; l++) {
//	  TString name_mbl_mbw_mwbbl_ywhad = Form("plots/fit_mbl_%d_mbw_%d_mwbbl_%d_ywhad_%d_.ps", 17-i, 11-j, 10-k, 24-2*l);
//	  cout<<"Should create histo named "<<name_mbl_mbw_mwbbl_ywhad<<endl;
//	  if (fitMultipleObservables(name_mbl_mbw_mwbbl_ywhad, {"mbl_selected", "mbwhad_selected", "mwhadbbl", "rapiditywhad"},
//				     {"m_bl", "m_bw", "m_wbbl", "y_whad"}, i, j, k, l) > 0) return 1;
	
//	//TString name_mbl_mbw = Form("plots/fit_mbl_%d_mbw_%d.ps", i, j);
//	//cout<<"Should create histo named "<<name_mbl_mbw<<endl;
//	//if (fitMultipleObservables(name_mbl_mbw, {"mbl_selected", "mbwhad_selected"},
//	//			   {"m_bl", "m_bw"}, i, j, 0, 0) > 0) return 1;
//	//TString name_mbl_mbw_drbl = Form("plots/fit_mbl_%d_mbw_%d_drbl_%d.ps", i, j, l);
//	//cout<<"Should create histo named "<<name_mbl_mbw_drbl<<endl;
//	//if (fitMultipleObservables(name_mbl_mbw_drbl, {"mbl_selected", "mbwhad_selected", "dRbl_selected"},
//	//			   {"m_bl", "m_bw", "dr_bl"}, i, j, l, 0) > 0) return 1;
//	//TString name_mbl_mbw_drbw = Form("plots/fit_mbl_%d_mbw_%d_drbw_%d.ps", i, j, l);
//        //cout<<"Should create histo named "<<name_mbl_mbw_drbw<<endl;
//	//if (fitMultipleObservables(name_mbl_mbw_drbw, {"mbl_selected", "mbwhad_selected", "dRbwhad_selected"},
//        //                           {"m_bl", "m_bw", "dr_bw"}, i, j, 0, l) > 0) return 1;
//      
//	
//	//TString name_mbl_mbw_ptwhad = Form("plots/fit_mbl_%d_mbw_%d_ptwhad_%d_.ps", i, j, l);
//        //cout<<"Should create histo named "<<name_mbl_mbw_ptwhad<<endl;
//	//if (fitMultipleObservables(name_mbl_mbw_ptwhad, {"mbl_selected", "mbwhad_selected", "ptwhad"},
//        //                           {"m_bl", "m_bw", "pT_whad"}, i, j, 0, l) > 0) return 1;
//	TString name_mbl_mbw_mwbbl = Form("plots/fit_mbl_%d_mbw_%d_mwbbl_%d.ps", 17-i, 11-j, 10-l);
//        cout<<"Should create histo named "<<name_mbl_mbw_mwbbl<<endl;
//	if (fitMultipleObservables(name_mbl_mbw_mwbbl, {"mbl_selected", "mbwhad_selected", "mwhadbbl"},
//                                   {"m_bl", "m_bw", "m_wbbl"}, i, j, 0, l) > 0) return 1;
//	//TString name_mbl_mbw_etal = Form("plots/fit_mbl_%d_mbw_%d_etal_%d.ps", i, j, l);
//        //cout<<"Should create histo named "<<name_mbl_mbw_etal<<endl;
//	//if (fitMultipleObservables(name_mbl_mbw_etal, {"mbl_selected", "mbwhad_selected", "etal1"},
//        //                           {"m_bl", "m_bw", "eta_lep1"}, i, j, 0, l) > 0) return 1;
//
//	TString name_mbl_mbw_rbq = Form("plots/fit_mbl_%d_mbw_%d_rbq_%d_.ps", 17-i, 11-j, 7-l);
//        cout<<"Should create histo named "<<name_mbl_mbw_rbq<<endl;
//	if (fitMultipleObservables(name_mbl_mbw_rbq, {"mbl_selected", "mbwhad_selected", "rbq"},
//				   {"m_bl", "m_bw", "r_bq"}, i, j, 0, l) > 0) return 1;
//        TString name_mbl_mbw_rbj = Form("plots/fit_mbl_%d_mbw_%d_rbj_%d_.ps", 17-i, 11-j, 7-l);
//        cout<<"Should create histo named "<<name_mbl_mbw_rbj<<endl;
//	if (fitMultipleObservables(name_mbl_mbw_rbj, {"mbl_selected", "mbwhad_selected", "rbj"},
//				   {"m_bl", "m_bw", "r_bj"}, i, j, 0, l) > 0) return 1;
//        TString name_mbl_mbw_ywhad = Form("plots/fit_mbl_%d_mbw_%d_ywhad_%d_.ps", 17-i, 11-j, 24-2*l);
//        cout<<"Should create histo named "<<name_mbl_mbw_ywhad<<endl;
//	if (fitMultipleObservables(name_mbl_mbw_ywhad, {"mbl_selected", "mbwhad_selected", "rapiditywhad"},
//				   {"m_bl", "m_bw", "y_whad"}, i, j, 0, l) > 0) return 1;
//
//				   
//				   //if (fitMultipleObservables("plots/fit_mbl.ps", {"mbl_selected", "mbwhad_selected"},
//	//{"m_bl", "m_bw"}, i, j) > 0) return 1;
//	}
//      }
//    }
//  }

  //LTF_ROOTTools::plotFitComparison();
  
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

int fitMultipleObservables(const char* ps_name, const vector<TString> fit_vars, const vector<TString> fit_vars_short, int excl_mbl, int excl_mbw, int excl_yw, int excl_mwbbl, int excl_extra) {
   using namespace std;
  
#ifdef __CLING__
   TH1D::AddDirectory(false);
#endif
   TH1::SetDefaultSumw2(true);

   //std::map<double, std::unique_ptr<TH1D>> templates;
   map<double,TH1D*> templates;
   bool doPseudo = false;
   bool divideByBinWidth = false;
   bool doJesStressTest = false;

   const int     iRebin       = 1;
   const int     iRebinData   = 1;

   map<TString, int> rmbins_low_fit_vars = {
     {"dRbl_selected",     0},
     {"dRbwhad_selected",  0},
     {"etal1",             0},
     {"mtlep",             0},
     {"mbl_selected",      0},
     {"mbwhad_selected",   0},
     {"minimax_whadbbl",   0},
     {"mwhadbbl",          0},
     {"mblmin",            0},
     {"mtblmin",           0},
     {"ptwhadbbl",         0},
     {"ptwhadbblnu",       0},
     {"ptb1",              0},
     {"ptb2",              0},
     {"ptbl_selected",     0},
     {"ptblnu_selected",   0},
     {"ptbwhad_selected",  0},
     {"ptl1",              0},
     {"met",               0},
     {"ptOWj1",            0},
     {"ptOWj2",            0},
     {"ptwhad",            0},
     {"rbj",               0},
     {"rbq",               0},
     {"rapiditywhad",      excl_yw},
   };

   map<TString, int> rmbins_up_fit_vars = {
     {"dRbl_selected",     excl_extra},
     {"dRbwhad_selected",  excl_extra},
     {"etal1",             excl_extra},
     {"mtlep",             excl_extra},
     {"mbl_selected",      excl_mbl},
     {"mbwhad_selected",   excl_mbw},
     {"minimax_whadbbl",   excl_extra},
     {"mwhadbbl",          excl_mwbbl},
     {"mblmin",            excl_extra},
     {"mtblmin",           excl_extra},
     {"ptwhadbbl",         excl_extra},
     {"ptwhadbblnu",       excl_extra},
     {"ptb1",              excl_extra},
     {"ptb2",              excl_extra},
     {"ptbl_selected",     excl_extra},
     {"ptblnu_selected",   excl_extra},
     {"ptbwhad_selected",  excl_extra},
     {"ptl1",              excl_extra},
     {"met",               excl_extra},
     {"ptOWj1",            excl_extra},
     {"ptOWj2",            excl_extra},
     {"ptwhad",            excl_extra},
     {"rbj",               excl_extra},
     {"rbq",               excl_extra},
     {"rapiditywhad",      excl_yw},
   };
   
   map<TString, int> rmbins_low_fit_vars_short = {
     {"dr_bl"    ,  0},
     {"dr_bw"    ,  0},
     {"eta_lep1" ,  0},
     {"mT_lep1"  ,  0},
     {"m_bl"     ,  0},
     {"m_bw"     ,  0},
     {"m_minimax",  0},
     {"m_wbbl"   ,  0},
     {"mblmin"   ,  0},
     {"mtblmin"  ,  0},
     {"pT_Wbbl"  ,  0},
     {"pT_Wbblnu",  0},
     {"pT_bjet1" ,  0},
     {"pT_bjet2" ,  0},
     {"pT_bl"    ,  0},
     {"pT_blnu"  ,  0},
     {"pT_bw"    ,  0},
     {"pT_lep1"  ,  0},
     {"pT_met"   ,  0},
     {"pT_owj1"  ,  0},
     {"pT_owj2"  ,  0},
     {"pT_whad"  ,  0},
     {"r_bj"     ,  0},
     {"r_bq"     ,  0},
     {"y_whad"   ,  excl_yw},
   };

   map<TString, int> rmbins_up_fit_vars_short = {
     {"dr_bl"    ,  excl_extra},
     {"dr_bw"    ,  excl_extra},
     {"eta_lep1" ,  excl_extra},
     {"mT_lep1"  ,  excl_extra},
     {"m_bl"     ,  excl_mbl},
     {"m_bw"     ,  excl_mbw},
     {"m_minimax",  excl_extra},
     {"m_wbbl"   ,  excl_mwbbl},
     {"mblmin"   ,  excl_extra},
     {"mtblmin"  ,  excl_extra},
     {"pT_Wbbl"  ,  excl_extra},
     {"pT_Wbblnu",  excl_extra},
     {"pT_bjet1" ,  excl_extra},
     {"pT_bjet2" ,  excl_extra},
     {"pT_bl"    ,  excl_extra},
     {"pT_blnu"  ,  excl_extra},
     {"pT_bw"    ,  excl_extra},
     {"pT_lep1"  ,  excl_extra},
     {"pT_met"   ,  excl_extra},
     {"pT_owj1"  ,  excl_extra},
     {"pT_owj2"  ,  excl_extra},
     {"pT_whad"  ,  excl_extra},
     {"r_bj"     ,  excl_extra},
     {"r_bq"     ,  excl_extra},
     {"y_whad"   ,  excl_yw},
   };

   const TString datafile =     "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/files/unfolding_SR_Whad_Final_l_Whad_particle_TUnfoldStandalone_OptionA_data_nonClosureAlternative_paper.root";
   const TString datafile_alt = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/files/unfolding_SR_Whad_Final_l_Whad_particle_TUnfoldStandalone_OptionA_data_nonClosureAlternative_extra.root";
   const TString template_uncertainty_file = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/files/template_theory_uncertainty.root";
   
   const TString covariancefile = "/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/files/data.Nominal.root";
   int bins_number = 0;
   std::vector<int> matrix_dimension = {0}; // Needed to put covariance matrix together
   for ( auto& tmp: fit_vars ) {
     TString name = "unfolding_"+tmp+"_NOSYS";
     int nbins = -1;

     TFile* file1 = TFile::Open(datafile);
     if (!file1 || file1->IsZombie()) {cout<<"Could not open file"<<datafile<<endl; exit(1);}

     TH1D* tmp_data = file1->Get<TH1D>(name);
     if ( tmp_data ) {
       nbins = tmp_data->GetNbinsX();
     }
     file1->Close();
     delete file1;

     if ( nbins < 0 ) {
       TFile* file2 = TFile::Open(datafile_alt);
       if (!file2 || file2->IsZombie()) {cout<<"Could not open file"<<datafile_alt<<endl; exit(1);}
       TH1D* tmp_data2 = file2->Get<TH1D>(name);
       if (tmp_data2) {
	 nbins = tmp_data2->GetNbinsX();
       }
       file2->Close();
       delete file2;
     }
     if (nbins < 0 ) {
       cerr<<"Trying to get the number of bins for "<<name<<endl;
       cerr<<"Could not find the histogram! Exiting ..."<<endl;
       exit(1);
     }
     else {
       bins_number += nbins -1 - rmbins_low_fit_vars[tmp] - rmbins_up_fit_vars[tmp];
       matrix_dimension.push_back(nbins - 1 - rmbins_low_fit_vars[tmp] - rmbins_up_fit_vars[tmp]);
     }
   }
   
   TH1D* combined_data = new TH1D("combined_data", "combined_data", bins_number, 0, bins_number);
   int bin_offset = 0;
   if ( doPseudo ) {
     for ( auto& tmp: fit_vars_short ) {
       bin_offset -= rmbins_low_fit_vars_short[tmp];
       //TFile* f = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/files/template_172_75.root");
       TFile* f = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/files/template_nominal_with_weights.root");
       TH1D* h = f->Get<TH1D>(tmp);
       h->SetDirectory(nullptr);
       delete f;

       //TFile* f_templ_unc = TFile::Open(template_uncertainty_file); // johannes directly modify the input data
       //TH1D* h_templ_unc = f_templ_unc->Get<TH1D>(Form("%s_templ_unc_aS_down",tmp.Data()));
       int extrabin = 0;
       if ( tmp == "m_bw" ) extrabin =1;
       for ( int i = 1 + rmbins_low_fit_vars_short[tmp]; i <= h->GetNbinsX() - rmbins_up_fit_vars_short[tmp] - extrabin; i++ ) {
	 if ( divideByBinWidth ) {
	   combined_data->SetBinContent(i+bin_offset, h->GetBinContent(i));
           combined_data->SetBinError(i+bin_offset, h->GetBinError(i));
	 }
	 else {
	   double bin_width = h->GetXaxis()->GetBinWidth(i);
	   combined_data->SetBinContent(i+bin_offset, h->GetBinContent(i)*bin_width);
	   //combined_data->SetBinContent(i+bin_offset, (h->GetBinContent(i)+h_templ_unc->GetBinContent(i))*bin_width);
	   combined_data->SetBinError(i+bin_offset, h->GetBinError(i)*bin_width);
	 }
       }
       bin_offset += h->GetNbinsX() - rmbins_up_fit_vars_short[tmp] - extrabin;
     }
   }
   else {
     for ( auto& tmp: fit_vars ) {
       TString name = "unfolding_"+tmp+"_NOSYS";
       TH1D* tmp_data = TFile::Open(datafile)->Get<TH1D>(name);
       if (!tmp_data) {
	 tmp_data = TFile::Open(datafile_alt)->Get<TH1D>(name);
       }
       if (!tmp_data) {
	 cerr<<"Could not find data histogram "<<name<<endl;
	 cerr<<"Exiting..."<<endl;
	 exit(1);
       }
       // Loop only to NbinsX-1, because last bin is overflow bin
       bin_offset -= rmbins_low_fit_vars[tmp];
       for ( int i = 1 + rmbins_low_fit_vars[tmp]; i <= tmp_data->GetNbinsX()-1 - rmbins_up_fit_vars[tmp]; i++ ) {
	 if ( divideByBinWidth ) combined_data->SetBinContent(i+bin_offset, tmp_data->GetBinContent(i) / tmp_data->GetXaxis()->GetBinWidth(i)); 
	 else combined_data->SetBinContent(i+bin_offset, tmp_data->GetBinContent(i)); ///tmp_data->GetXaxis()->GetBinWidth(i));
       }
       bin_offset += tmp_data->GetNbinsX()-1 - rmbins_low_fit_vars[tmp] - rmbins_up_fit_vars[tmp];
     }
   }
   combined_data->Rebin(iRebinData);
   combined_data->SetLineColor(kBlack);
   combined_data->SetMarkerSize(1.8);

   {
     TH1D* combined_template_169    = new TH1D("combined_template_169",    "combined_template_169",    bins_number, 0, bins_number);
     TH1D* combined_template_171    = new TH1D("combined_template_171",    "combined_template_171",    bins_number, 0, bins_number);
     TH1D* combined_template_172_25 = new TH1D("combined_template_172_25", "combined_template_172_25", bins_number, 0, bins_number);
     TH1D* combined_template_172_75 = new TH1D("combined_template_172_75", "combined_template_172_75", bins_number, 0, bins_number);
     TH1D* combined_template_174    = new TH1D("combined_template_174",    "combined_template_174",    bins_number, 0, bins_number);
     TH1D* combined_template_176    = new TH1D("combined_template_176",    "combined_template_176",    bins_number, 0, bins_number);
     int bin_offset = 0;
     for ( auto& tmp: fit_vars_short ) {

       TFile* f_169 =   TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/files/template_169.root");
       TH1D* h_tmp_169 = f_169->Get<TH1D>(tmp);
       h_tmp_169->SetDirectory(nullptr);
       delete f_169;
       
       TFile* f_171 =	TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/files/template_171.root");
       TH1D* h_tmp_171 = f_171->Get<TH1D>(tmp);
       h_tmp_171->SetDirectory(nullptr);
       delete f_171;

       TFile* f_172_25 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/files/template_172_25.root");
       TH1D* h_tmp_172_25 = f_172_25->Get<TH1D>(tmp);
       h_tmp_172_25->SetDirectory(nullptr);
       delete f_172_25;

       TFile* f_172_75 = TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/files/template_172_75.root");
       TH1D* h_tmp_172_75 = f_172_75->Get<TH1D>(tmp);
       h_tmp_172_75->SetDirectory(nullptr);
       delete f_172_75;
       
       TFile* f_174 =   TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/files/template_174.root");
       TH1D* h_tmp_174 = f_174->Get<TH1D>(tmp);
       h_tmp_174->SetDirectory(nullptr);
       delete f_174;
              
       TFile* f_176 =   TFile::Open("/home/iwsatlas1/jhessler/LTF/LinearTemplateFit/LTF_Eigen/files/template_176.root");
       TH1D* h_tmp_176 = f_176->Get<TH1D>(tmp);
       h_tmp_176->SetDirectory(nullptr);
       delete f_176;

       double kFactor = 904.0/802.0; // Define k-factor when comparing templates to data
       bin_offset -= rmbins_low_fit_vars_short[tmp];
       int extrabin = 0;
       if ( tmp == "m_bw" ) extrabin = 1;
       for ( int i = 1 + rmbins_low_fit_vars_short[tmp]; i<= h_tmp_169->GetNbinsX() - rmbins_up_fit_vars_short[tmp] - extrabin; i++ ) {
	 if ( divideByBinWidth ) {
	   combined_template_169->SetBinContent(   i+bin_offset, h_tmp_169->GetBinContent(i)   * kFactor);
           combined_template_171->SetBinContent(   i+bin_offset, h_tmp_171->GetBinContent(i)   * kFactor);
           combined_template_172_25->SetBinContent(i+bin_offset, h_tmp_172_25->GetBinContent(i)* kFactor);
           combined_template_172_75->SetBinContent(i+bin_offset, h_tmp_172_75->GetBinContent(i)* kFactor);
           combined_template_174->SetBinContent(   i+bin_offset, h_tmp_174->GetBinContent(i)   * kFactor);
           combined_template_176->SetBinContent(   i+bin_offset, h_tmp_176->GetBinContent(i)   * kFactor);

           combined_template_169->SetBinError(    i+bin_offset, h_tmp_169->GetBinError(i)    * kFactor);
           combined_template_171->SetBinError(    i+bin_offset, h_tmp_171->GetBinError(i)    * kFactor);
           combined_template_172_25->SetBinError( i+bin_offset, h_tmp_172_25->GetBinError(i) * kFactor);
           combined_template_172_75->SetBinError( i+bin_offset, h_tmp_172_75->GetBinError(i) * kFactor);
           combined_template_174->SetBinError(    i+bin_offset, h_tmp_174->GetBinError(i)    * kFactor);
           combined_template_176->SetBinError(    i+bin_offset, h_tmp_176->GetBinError(i)    * kFactor);
	 }
	 else {
	   double bin_width_x = h_tmp_169->GetXaxis()->GetBinWidth(i);
	   combined_template_169->SetBinContent(   i+bin_offset, h_tmp_169->GetBinContent(i)   * kFactor*bin_width_x);
	   combined_template_171->SetBinContent(   i+bin_offset, h_tmp_171->GetBinContent(i)   * kFactor*bin_width_x);
	   combined_template_172_25->SetBinContent(i+bin_offset, h_tmp_172_25->GetBinContent(i)* kFactor*bin_width_x);
	   combined_template_172_75->SetBinContent(i+bin_offset, h_tmp_172_75->GetBinContent(i)* kFactor*bin_width_x);
	   combined_template_174->SetBinContent(   i+bin_offset, h_tmp_174->GetBinContent(i)   * kFactor*bin_width_x);
	   combined_template_176->SetBinContent(   i+bin_offset, h_tmp_176->GetBinContent(i)   * kFactor*bin_width_x);
	   
	   combined_template_169->SetBinError(    i+bin_offset, h_tmp_169->GetBinError(i)    * kFactor * bin_width_x);
	   combined_template_171->SetBinError(    i+bin_offset, h_tmp_171->GetBinError(i)    * kFactor * bin_width_x);
	   combined_template_172_25->SetBinError( i+bin_offset, h_tmp_172_25->GetBinError(i) * kFactor * bin_width_x);
	   combined_template_172_75->SetBinError( i+bin_offset, h_tmp_172_75->GetBinError(i) * kFactor * bin_width_x);
	   combined_template_174->SetBinError(    i+bin_offset, h_tmp_174->GetBinError(i)    * kFactor * bin_width_x);
	   combined_template_176->SetBinError(    i+bin_offset, h_tmp_176->GetBinError(i)    * kFactor * bin_width_x);
	 }
       }
       bin_offset += h_tmp_169->GetNbinsX() - rmbins_up_fit_vars_short[tmp] - extrabin;
     }
     templates[169]    = combined_template_169;
     templates[171]    = combined_template_171;
     templates[172.25] = combined_template_172_25;
     templates[172.75] = combined_template_172_75;
     templates[174]    = combined_template_174;
     templates[176]    = combined_template_176;
   }

   for ( auto [MM,hist] : templates ) {
      hist->Rebin(iRebin);
   }
   
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
     //"WMASS_VAR_signal", //MET+JVT+PileUp
     "FT_EFF_Eigen_B_0", "FT_EFF_Eigen_B_1", "FT_EFF_Eigen_B_2", "FT_EFF_Eigen_B_3", "FT_EFF_Eigen_B_4",
     "FT_EFF_Eigen_B_5", "FT_EFF_Eigen_B_6", "FT_EFF_Eigen_B_7", "FT_EFF_Eigen_B_8", "FT_EFF_Eigen_C_0",
     "FT_EFF_Eigen_C_1", "FT_EFF_Eigen_C_2", "FT_EFF_Eigen_C_3", "FT_EFF_Eigen_Light_0", "FT_EFF_Eigen_Light_1",
     "FT_EFF_Eigen_Light_2", "FT_EFF_Eigen_Light_3", "FT_EFF_extrapolation", "FT_EFF_extrapolation_from_charm", // b-tagging uncertainties
     "THEORY_CROSS_SECTION_signal", "THEORY_SHOWERING_HERWIG7_signal", "THEORY_SCALE_FACTORISATION_signal",
     "THEORY_SCALE_RENORMALISATION_signal", "THEORY_ISR_signal", "THEORY_FSR_signal", "THEORY_HDAMP_signal",
     "THEORY_PTHARD_signal", "THEORY_TOPRECOIL_signal", "THEORY_TOP_MASS_signal",
     "THEORY_PDF4LHC_VARIATION_signal", "THEORY_DR_DS_signal", // modelling uncertainties
     "THEORY_CROSS_SECTION_Wjets", "THEORY_SCALE_COMBINED_Wjets", "THEORY_PDF4LHC_VARIATION_Wjets", "THEORY_EWK_Wjets",
     "THEORY_SCALE_COMBINED_multiboson_noW", "THEORY_PDF4LHC_VARIATION_multiboson_noW", //"THEORY_EWK_multiboson_noW",
     "THEORY_CROSS_SECTION_other_top_noWt", "FAKES_Electron", "FAKES_Muon", // bkgd uncertainties
     "LUMINOSITY" // lumi uncertainty
   };
   if ( doJesStressTest ) uncertainties.push_back("myJESUncertainty");
   
   
   vector<string> statistical_uncertainties = {"STAT_MC"};
   
   vector<string> external_uncertainties = {"FULL_SYS_SUM", "FULL_SYS_SUM_DETECTOR", "FULL_SYS_SUM_THEORY",
					    "TOTAL_SYSONLY", "TOTAL", "TOTAL_NO_DR_DS", "FULL_SYS_TOYS"};

   vector<string> template_uncertainties_full = {"templ_unc_PS"};
   
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

   // --- Add template errors
   {
     vector<double> template_values = {169, 171, 172.25, 172.75, 174, 176};
     vector<string> template_uncertainties = {"PS", "aS_up", "aS_down", "isr_aS_up", "isr_aS_down", "fsr_scale_up", "fsr_scale_down", "herwig", "hdamp",
					      "pthard", "recoilToTop", "A14Var1Up", "A14Var1Down", "CR1", "CR2", "bquark_frag", "lineshape", "DS", "NNLO_tt",
					      "muR_05_muF_05", "muR_05_muF_10", "muR_10_muF_05", "muR_10_muF_20", "muR_20_muF_10", "muR_20_muF_20",
					      "PDF_90901", "PDF_90902", "PDF_90903", "PDF_90904", "PDF_90905", "PDF_90906", "PDF_90907", "PDF_90908",
					      "PDF_90909", "PDF_90910", "PDF_90911", "PDF_90912", "PDF_90913", "PDF_90914", "PDF_90915", "PDF_90916",
					      "PDF_90917", "PDF_90918", "PDF_90919", "PDF_90920", "PDF_90921", "PDF_90922", "PDF_90923", "PDF_90924",
					      "PDF_90925", "PDF_90926", "PDF_90927", "PDF_90928", "PDF_90929", "PDF_90930"
     };

     TFile* f_templ_unc = TFile::Open(template_uncertainty_file);
     cout<<"Adding template errors"<<endl;
     for (auto name: template_uncertainties) {
       vector<double> uncertainty;
       for ( auto& tmp: fit_vars_short ) {
	 //cout<<"Get histo "<<Form("%s_templ_unc_nom_%s",tmp.Data(), name.c_str())<<" from file "<<f_templ_unc->GetName()<<endl;
         TH1D* h = f_templ_unc->Get<TH1D>(Form("%s_templ_unc_%s",tmp.Data(), name.c_str()));
	 int extrabin = 0;
	 if ( tmp == "m_bw" ) extrabin = 1;
         for ( int i = 1 + rmbins_low_fit_vars_short[tmp]; i <= h->GetNbinsX() - rmbins_up_fit_vars_short[tmp] - extrabin; i++ ) {
	   if ( h->GetXaxis()->GetBinUpEdge(i) == 1000 ) continue;
           uncertainty.push_back(h->GetBinContent(i) * h->GetBinWidth(i));
         }
       }
       for (auto topmass: template_values) { // same error will be passed for each template
	 ltf.AddTemplateError(name, topmass, uncertainty, 1.0 );
       }
     }
     f_templ_unc->Close();
     delete f_templ_unc;
   }

   // --- initialize data
   ltf.SetData( combined_data->GetNbinsX(), combined_data->GetArray()+1);
   
   // Start building covariance matrix
   // Get data covariance matrix after unfolding from TUnfold
   // Covariance defined as rho_ij*sigma_i*sigma_j
   // sigma is the absolute error

   // In order to plot the correlation matrix in an interactive root session, just do
   // auto htext = (TH2D*)_file0->Get("combined_corr")
   // gStyle->SetPaintTextFormat("4.2f");
   // htext->SetMarkerSize(1.0);
   // htext->Draw("TEXT45 colz")
   TH2D* combined_covariance = new TH2D("combined_cov", "combined_cov", bins_number, 0, bins_number, bins_number, 0, bins_number);
   TH2D* combined_correlation = new TH2D("combined_corr", "combined_corr", bins_number, 0, bins_number, bins_number, 0, bins_number);
   bin_offset = 0;
   for ( auto& fit_var: fit_vars ) {

     TFile* file = TFile::Open(datafile);
     if (!file || file->IsOpen() == kFALSE) {
       std::cerr << "Error: Couldn't open the file "<<datafile<< std::endl;
       return 1;
     }
     TString histnameCovStat("unfolding_covariance_matrix_"+fit_var+"_covariance_STAT_DATA");
     TH2D* cov_stat_dat = file->Get<TH2D>(histnameCovStat);
     if ( !cov_stat_dat ) {
        file->Close();
	file =TFile::Open(datafile_alt);
	cov_stat_dat = file->Get<TH2D>(histnameCovStat);
     }
     if ( !cov_stat_dat ) {
	cerr<<"Could not find covariance matrix " << histnameCovStat <<endl; exit(1);
     }
     
      TH1D* h_err  = file->Get<TH1D>("unfolding_error_"+fit_var+"_direct_envelope_STAT_DATA__1up");
      TH1D* h_data = file->Get<TH1D>("unfolding_"+fit_var+"_NOSYS");
      bin_offset -= rmbins_low_fit_vars[fit_var];
      for ( int i = 1 + rmbins_low_fit_vars[fit_var]; i < cov_stat_dat->GetNbinsX() - rmbins_up_fit_vars[fit_var]; i++ ) {
	double bin_width_x = h_err->GetXaxis()->GetBinWidth(i);
	double sigma_data_x = sqrt(cov_stat_dat->GetBinContent(i,i)); //h_data->GetBinContent(i) * h_err->GetBinContent(i);

	for ( int j = 1 + rmbins_low_fit_vars[fit_var]; j < cov_stat_dat->GetNbinsY() - rmbins_up_fit_vars[fit_var]; j++ ) {
	  double bin_width_y = h_err->GetXaxis()->GetBinWidth(j);
	  double sigma_data_y = sqrt(cov_stat_dat->GetBinContent(j,j)); //h_data->GetBinContent(j) * h_err->GetBinContent(j);
	  if ( divideByBinWidth ) {
	    combined_covariance->SetBinContent(i+bin_offset,j+bin_offset, cov_stat_dat->GetBinContent(i,j) / (bin_width_x*bin_width_y));
	  }
	  else {
	    combined_covariance->SetBinContent(i+bin_offset,j+bin_offset, cov_stat_dat->GetBinContent(i,j));
	  }
	  combined_correlation->SetBinContent(i+bin_offset,j+bin_offset, cov_stat_dat->GetBinContent(i,j)/(sigma_data_x*sigma_data_y));
	}
      }
      bin_offset += cov_stat_dat->GetNbinsX() - 1 - rmbins_up_fit_vars[fit_var];
      file->Close();
      delete file;
   }
   
   // Fill off-diagonal block matrices of the covariance matrix
   if ( fit_vars.size() > 1 ) {
     for ( int v1 = 0; v1 < fit_vars_short.size(); v1++ ) {
       
       TFile* file = TFile::Open(datafile);
       TFile* file_alt = TFile::Open(datafile_alt);
       file->cd();
       if (!file || file->IsOpen() == kFALSE) {
	 std::cerr << "Error: Couldn't open the file "<<datafile<< std::endl;
	 return 1;
       }
       TH1D* h_err_v1 = file->Get<TH1D>("unfolding_error_"+fit_vars[v1]+"_direct_envelope_STAT_DATA__1up");
       TH1D* h_data_var1 = file->Get<TH1D>("unfolding_"+fit_vars[v1]+"_NOSYS");
       TString histnameCovStatX("unfolding_covariance_matrix_"+fit_vars[v1]+"_covariance_STAT_DATA");
       TH2D* cov_stat_dat_x = file->Get<TH2D>(histnameCovStatX);
       if ( !h_err_v1 || !h_data_var1 ) {
	 file_alt->cd();
	 h_err_v1 = file_alt->Get<TH1D>("unfolding_error_"+fit_vars[v1]+"_direct_envelope_STAT_DATA__1up");
	 h_data_var1 = file_alt->Get<TH1D>("unfolding_"+fit_vars[v1]+"_NOSYS");
	 cov_stat_dat_x = file_alt->Get<TH2D>(histnameCovStatX);
       }
       if ( !h_err_v1 || !h_data_var1 ) {
	 cerr<<"Could not find statistical uncertainty for variable " << fit_vars[v1] <<endl; exit(1);
       }
       
       //TH1D* h_err_v1 = file->Get<TH1D>("unfolding_error_"+fit_vars[v1]+"_direct_envelope_STAT_DATA__1up");
       //TH1D* h_data_var1 = file->Get<TH1D>("unfolding_"+fit_vars[v1]+"_NOSYS");

       for ( int v2 = v1+1; v2 < fit_vars_short.size(); v2++ ) {
	 file->cd();
	 TH1D* h_err_v2 = file->Get<TH1D>("unfolding_error_"+fit_vars[v2]+"_direct_envelope_STAT_DATA__1up"); // get relative error
	 TH1D* h_data_var2 = file->Get<TH1D>("unfolding_"+fit_vars[v2]+"_NOSYS");
	 if ( !h_err_v2 || !h_data_var2 ) {
	   file_alt->cd();
	   h_err_v2 = file_alt->Get<TH1D>("unfolding_error_"+fit_vars[v2]+"_direct_envelope_STAT_DATA__1up");
	   h_data_var2 = file_alt->Get<TH1D>("unfolding_"+fit_vars[v2]+"_NOSYS");
	 }
	 
         if ( !h_data_var2 ) {
	   cerr<<"Could not find data for " << fit_vars[v2] <<endl;
	   exit(1);
	 }
	 TFile* f_cov = TFile::Open(covariancefile);
	 TString histname = "NOSYS/hist/l_Whad_particle_"+fit_vars[v1]+"_"+fit_vars[v2];
	 TH2D* cov_orig = f_cov->Get<TH2D>(histname);
	 if ( !cov_orig ) { cerr<<"Could not find correlation matrix " << histname <<endl; exit(1);}
	 
	 std::vector<double> binEdgesX;
	 int nx = h_err_v1->GetXaxis()->GetNbins();
	 binEdgesX.reserve(nx + 1);
	 for (int i = 1; i <= nx; ++i) binEdgesX.push_back(h_err_v1->GetXaxis()->GetBinLowEdge(i));
	 binEdgesX.push_back(h_err_v1->GetXaxis()->GetBinUpEdge(nx));

	 std::vector<double> binEdgesY;
         int ny = h_err_v2->GetXaxis()->GetNbins();
         binEdgesY.reserve(ny + 1);
         for (int i = 1; i <= ny; ++i) binEdgesY.push_back(h_err_v2->GetXaxis()->GetBinLowEdge(i));
         binEdgesY.push_back(h_err_v2->GetXaxis()->GetBinUpEdge(ny));
	 
	 TH2D* cov = LTF_ROOTTools::rebin(cov_orig, binEdgesX, binEdgesY, "cov");
	 
	 // Projections of correlation matrix
	 TH1D* projection_var1 = cov->ProjectionX("pro_v1",0,-1,"e");
	 TH1D* projection_var2 = cov->ProjectionY("pro_v2",0,-1,"e");

	 int bin_offset_x = 0;
	 int bin_offset_y = 0;
	 for ( int v1_tmp = 0; v1_tmp <= v1; v1_tmp++ ) bin_offset_x += matrix_dimension[v1_tmp];
	 for ( int v2_tmp = 0; v2_tmp <= v2; v2_tmp++ ) bin_offset_y += matrix_dimension[v2_tmp];

	 TString histnameCovStatY("unfolding_covariance_matrix_"+fit_vars[v2]+"_covariance_STAT_DATA");
	 file->cd();
	 TH2D* cov_stat_dat_y = file->Get<TH2D>(histnameCovStatY);
	 if ( !cov_stat_dat_y ) {
	   file_alt->cd();
	   cov_stat_dat_y = file_alt->Get<TH2D>(histnameCovStatY);
	 }
	 if ( !cov_stat_dat_y )	{
	   cerr<<"Could not find the matrix "<<histnameCovStatY<<endl;
	   exit(1);
	 }
	 bin_offset_x -= rmbins_low_fit_vars[fit_vars[v1]];
	 bin_offset_y -= rmbins_low_fit_vars[fit_vars[v2]]; //johannes check if this is correct
	 for ( int i = 1 + rmbins_low_fit_vars[fit_vars[v1]]; i < h_err_v1->GetNbinsX() - rmbins_up_fit_vars[fit_vars[v1]]; i++ ) {
	   for ( int j = 1 + rmbins_low_fit_vars[fit_vars[v2]]; j < h_err_v2->GetNbinsX() - rmbins_up_fit_vars[fit_vars[v2]]; j++ ) {
	     // Calculate the covariance
	     // 1. Get the covariance matrix measured on detector level
	     // 2. Divide by the std deviation to get the correlation: sigma_detector1, sigma_detector2
	     // 3. Multiply with the std deviation from the unfolded data: sigma_data1, sigma_data2
	     double sigma_data1 = sqrt(cov_stat_dat_x->GetBinContent(i,i)); //h_data_var1->GetBinContent(i) * h_err_v1->GetBinContent(i);
	     double sigma_data2 = sqrt(cov_stat_dat_y->GetBinContent(j,j)); //h_data_var2->GetBinContent(j) * h_err_v2->GetBinContent(j);
	     double sigma_detector1 = projection_var1->GetBinError(i);
	     double sigma_detector2 = projection_var2->GetBinError(j);
	     if ( abs(sigma_detector1) > 0 && abs(sigma_detector2) > 0 ) {
	       combined_covariance->SetBinContent(i+bin_offset_x, j+bin_offset_y, std::pow(cov->GetBinError(i,j),2)*(sigma_data1*sigma_data2) / (sigma_detector1*sigma_detector2));
	       combined_covariance->SetBinContent(j+bin_offset_y, i+bin_offset_x, std::pow(cov->GetBinError(i,j),2)*(sigma_data1*sigma_data2) / (sigma_detector1*sigma_detector2));
	       combined_correlation->SetBinContent(i+bin_offset_x, j+bin_offset_y, std::pow(cov->GetBinError(i,j),2) / (sigma_detector1*sigma_detector2));
	       combined_correlation->SetBinContent(j+bin_offset_y, i+bin_offset_x, std::pow(cov->GetBinError(i,j),2) / (sigma_detector1*sigma_detector2));
	     }
	     else {
	       cout<<"ERROR in example_ATLAS_topmass! Migration martix needs more statistics! Set values to 0.0 for now"<<endl;
	       combined_covariance->SetBinContent(i+bin_offset_x, j+bin_offset_y, 0);
               combined_covariance->SetBinContent(j+bin_offset_y, i+bin_offset_x, 0);
               combined_correlation->SetBinContent(i+bin_offset_x, j+bin_offset_y, 0);
               combined_correlation->SetBinContent(j+bin_offset_y, i+bin_offset_x, 0);
	     }
	   }
	   //f_cov->Close(); //johannes
	   //delete f_cov;
	 }
       }
       delete file;
       delete file_alt;
     }
   }
   const bool writeCovarianceMatrix = true;
   if ( writeCovarianceMatrix ) {
     TString name_outfile = "matrices/correlation";
     for ( auto tmp: fit_vars_short ) name_outfile += tmp;
     TFile *outfile = new TFile(name_outfile+".root", "RECREATE");
     combined_covariance->Write();
     combined_correlation->Write();
     delete outfile;
   }

   vector<vector<double > > vecCov2 = TH2D_to_vecvec(combined_covariance);
   ltf.AddError("STAT_DATA", vecCov2, LTF::Uncertainty::Constrained);

   int rows = vecCov2.size();
   int cols = vecCov2[0].size();
   TMatrixD m(rows, cols);
   for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
        m(i, j) = vecCov2[i][j];
    }
   }
   
   TH1D* total_error = new TH1D("total_error", "total_error", bins_number, 0, bins_number);
   //johannes begin
   // Systematical uncertainties
   for ( auto& uncertainty: uncertainties ) {
     vector<double> combined_error, combined_error_up, combined_error_down;

     for ( auto& fit_variable: fit_vars ) {
       //cout<<"Adding systematical uncertaintiy "<<uncertainty<<endl;
       TFile* file = TFile::Open(datafile);
       if (!file || file->IsOpen() == kFALSE) {
	 std::cerr << "Error: Couldn't open the file "<<datafile<< std::endl;
	 return 1;
       }
       TH1D* hist_nom  = file->Get<TH1D>("unfolding_"+fit_variable+"_NOSYS");
       TH1D* hist_up   = file->Get<TH1D>("unfolding_"+fit_variable+"_direct_"+uncertainty+"__1up");
       
       if ( !hist_up ) {
	 file->Close();
	 file =TFile::Open(datafile_alt);
	 hist_nom  = file->Get<TH1D>("unfolding_"+fit_variable+"_NOSYS");
	 hist_up = file->Get<TH1D>("unfolding_"+fit_variable+"_direct_"+uncertainty+"__1up");
       }
       if ( !hist_up ) {
	 cerr<<"Could not find histogram for systematic uncertainty " << uncertainty <<endl; exit(1);
       }
       TH1D* hist_down = file->Get<TH1D>("unfolding_"+fit_variable+"_direct_"+uncertainty+"__1down");
       for (int i=1 + rmbins_low_fit_vars[fit_variable]; i< hist_up->GetNbinsX() - rmbins_up_fit_vars[fit_variable]; i++) {
	 
	 double err_up = (hist_up->GetBinContent(i) - hist_nom->GetBinContent(i)) / hist_nom->GetBinContent(i);
	 double err_down = (hist_down->GetBinContent(i)- hist_nom->GetBinContent(i)) / hist_nom->GetBinContent(i);

	 combined_error.push_back(err_up);
	 combined_error_up.push_back(err_up);
         combined_error_down.push_back(err_down);
	 total_error->SetBinContent(i, std::sqrt(std::pow(total_error->GetBinContent(i),2)+std::pow(err_up,2))); // Use up variation
	 //total_error->SetBinContent(i, std::sqrt(std::pow(total_error->GetBinContent(i),2)+std::pow(hist_down->GetBinContent(i),2))); // Use down variation
       }
       delete file;
     }
     //}
     double corr = 1.0;
     if ( combined_error.size() > 0 ) {
       ltf.AddErrorRelative(uncertainty, combined_error, corr, LTF::Uncertainty::Constrained);
       ltf.AddErrorRelative(Form("%s_up",uncertainty.c_str()), combined_error_up, corr, LTF::Uncertainty::External);
       ltf.AddErrorRelative(Form("%s_down",uncertainty.c_str()), combined_error_down, corr, LTF::Uncertainty::External);
     }
     combined_error.clear();
     combined_error_up.clear();
     combined_error_down.clear();
   }
//johannes end   
   // Stat uncertainties
   for ( auto& uncertainty: statistical_uncertainties ) { // Only stat. unc. is MC stat. 
     vector<double> combined_error;
     for ( auto& fit_variable: fit_vars ) {
       TFile* file = TFile::Open(datafile);
       if (!file || file->IsOpen() == kFALSE) {
	 std::cerr << "Error: Couldn't open the file "<<datafile<< std::endl;
	 return 1;
       }
       TH1D* hist_up = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_"+uncertainty+"__1up");

       if ( !hist_up ) {
	 file->Close();
	 file =TFile::Open(datafile_alt);
	 hist_up = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_"+uncertainty+"__1up");
       }
       if ( !hist_up ) {
	 cerr<<"Could not find histogram for systematic uncertainty " << uncertainty <<endl; exit(1);
       }
       //TH1D* hist_down = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_"+uncertainty+"__1down");
       for (int i=1 + rmbins_low_fit_vars[fit_variable]; i< hist_up->GetNbinsX() - rmbins_up_fit_vars[fit_variable]; i++) {
	 combined_error.push_back(hist_up->GetBinContent(i));
         total_error->SetBinContent(i, std::sqrt(std::pow(total_error->GetBinContent(i),2)+std::pow(hist_up->GetBinContent(i),2)));
       }
       delete file;
     }
     double corr = 0.0; // Stat. unc. are uncorrelated
     if ( combined_error.size() > 0 ) ltf.AddErrorRelative(uncertainty, combined_error, corr, LTF::Uncertainty::Constrained);
     combined_error.clear();
   }

   // Unconstrained uncertainty to account for normalization
   {
     vector<double> combined_error;
     for ( auto& fit_variable: fit_vars ) {
       TFile* file = TFile::Open(datafile);
       if (!file || file->IsOpen() == kFALSE) {
         std::cerr << "Error: Couldn't open the file "<<datafile<< std::endl;
         return 1;
       }
       // Get some random histo to get the correct binning
       TH1D* hist_up = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_STAT_MC__1up");
       if ( !hist_up ) {
         file->Close();
         file =TFile::Open(datafile_alt);
         hist_up = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_STAT_MC__1up");
       }
       if ( !hist_up ) {
         cerr<<"Could not find histogram for systematic uncertainty Normalization" << endl; exit(1);
       }
 
       for (int i=1 + rmbins_low_fit_vars[fit_variable]; i< hist_up->GetNbinsX() - rmbins_up_fit_vars[fit_variable]; i++) {
         combined_error.push_back(0.001); // 0.1% of normalization error
         //total_error->SetBinContent(i, std::sqrt(std::pow(total_error->GetBinContent(i),2)+std::pow(hist_up->GetBinContent(i),2)));
       }
       delete file;
     }
     double corr = 1.0;
     if ( combined_error.size() > 0 ) ltf.AddErrorRelative("Normalization", combined_error, corr, LTF::Uncertainty::Unconstrained);
     combined_error.clear();
   }
   
   // External uncertainties
   for ( auto& uncertainty: external_uncertainties ) {
     vector<double> combined_error;
     for ( auto& fit_variable: fit_vars ) {
       TFile* file = TFile::Open(datafile);
       if (!file || file->IsOpen() == kFALSE) {
         std::cerr << "Error: Couldn't open the file "<<datafile<< std::endl;
         return 1;
       }
       TH1D* hist_up = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_"+uncertainty+"__1up");

       if ( !hist_up ) {
         file->Close();
         file =TFile::Open(datafile_alt);
         hist_up = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_"+uncertainty+"__1up");
       }
       if ( !hist_up ) {
         cerr<<"Could not find histogram for systematic uncertainty " << uncertainty <<endl; exit(1);
       }
       cout<<"Number of bins "<<hist_up->GetNbinsX()<<endl;
       //TH1D* hist = file->Get<TH1D>("unfolding_error_"+fit_variable+"_direct_envelope_"+uncertainty+"__1up");
       for (int i=1 + rmbins_low_fit_vars[fit_variable]; i< hist_up->GetNbinsX() - rmbins_up_fit_vars[fit_variable]; i++) {
         combined_error.push_back(hist_up->GetBinContent(i));
       }
       delete file;
     }
     cout<<"Combined error size "<<combined_error.size()<<endl;
     double corr = 1;
     if ( combined_error.size() > 0 ) ltf.AddErrorRelative(uncertainty, combined_error, corr, LTF::Uncertainty::External);
     combined_error.clear();
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

   combined_data->Clear();
   for (auto& [key, hist] : templates) {
     delete hist;  // This also deregisters it from gDirectory
   }
   templates.clear();
   total_error->Delete();
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
