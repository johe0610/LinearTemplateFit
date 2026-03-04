//#include "include/LTF/LTF.h"
#include "LTF/LTF.h"
#include "LTF/LTF_ROOTTools.h"

#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TMatrixDfwd.h>
#include <TMatrixD.h>
#include <TMatrixTSym.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSystem.h>
#include <TFile.h>
#include <string>
#include <fstream>
#include <TRandom3.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TF1.h>
#include <TF2.h>
#include <TLine.h>
#include <TRandom3.h>
#include <TGraph2DErrors.h>
#include <TGraph2D.h>
#include "TColor.h"
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <set>

using namespace std;

// __________________________________________________________________________________ //
//!
//! read_input_table2()
//!
//! read input data table from file 
//! 
std::map < std::string, std::vector<double> > LTF_ROOTTools::read_input_table2(std::string filename, int ncol ) 
{
   std::map < std::string, std::vector<double> > ret;
   std::vector<std::string> cols;
   // open file for reading                                           
   std::ifstream istrm(filename.c_str(), std::ios::binary);
   if (!istrm.is_open()) {
      std::cout << "failed to open " << filename << std::endl;
      exit(1);
   }
   for ( int c=0 ; c<ncol ; c++ ) {
      std::string colname;
      istrm >> colname;
      ret[colname] = vector<double>();
      cols.push_back(colname);
   }
   while ( istrm.good()) {
      double value;
      for ( int c=0 ; c<ncol ; c++ ) {
         istrm >> value;
         if ( !istrm.good() ) break;
         std::string colname = cols[c];
         ret[colname].push_back(value);
      }
   }
   cout<<"Info.  [read_input_table]  Read "<<ret.size()<<" rows."<<endl;
   return ret;
}



// __________________________________________________________________________________ //
//! 
//!  MakeTGraph
//!
//!  make a TGraph for plotting purposes
//!
TGraphErrors* LTF_ROOTTools::MakeTGraph(const Eigen::VectorXd& xvalues, int ibin, const Eigen::MatrixXd& Y,
                         const std::map<std::string,Eigen::MatrixXd >& VSysY ) 
{

   TGraphErrors* graph = new TGraphErrors();
   for ( int i = 0 ; i<xvalues.size() ; i++ ) {
      double yvalue = Y(ibin,i);
      graph->SetPoint(i,xvalues(i),yvalue);
      double ey2 = 0;
      for ( auto [name,Vy] : VSysY ) {
	ey2 += pow(Vy(ibin,i),2); // johannes change this back!!
	//ey2 += pow(Vy(ibin,1),2);
      }
      graph->SetPointError(i,0,sqrt(ey2));
   }
   return graph;
}


// __________________________________________________________________________________ //
//! 
//!  MakeTGraph2D
//!
//!  make a TGraph2D for plotting purposes
//!
TGraph2DErrors* LTF_ROOTTools::MakeTGraph2D(const Eigen::VectorXd& xvalues, const Eigen::VectorXd& yvalues,
                                            int ibin, const Eigen::MatrixXd& Y,
                                            const std::map<std::string,Eigen::MatrixXd >& VSysY ) 
{

   TGraph2DErrors* graph = new TGraph2DErrors();
   for ( int i = 0 ; i<xvalues.size() ; i++ ) {
      double zvalue = Y(ibin,i);
      graph->SetPoint(i,xvalues(i),yvalues(i),zvalue);
      double ey2 = 0;
      for ( auto [name,Vy] : VSysY ) {
         ey2 += pow(Vy(ibin,i),2);
      }
      graph->SetPointError(i,0,0,sqrt(ey2));
   }
   return graph;
}


// __________________________________________________________________________________ //
//!
//! MakeHistogram
//!
//! make a histogram and fill it with random events according to a gauss
//! distribution around M
//!
TH1D* LTF_ROOTTools::MakeHistogram(int nEvents, int seed, double mean, double sigma, vector<double> bins ) 
{
   TRandom3 rn(seed);
   
   TH1D* hist = new TH1D("hist","hist",bins.size()-1, &bins[0]);
   hist->Sumw2();
   for ( int i=0 ; i<nEvents; i++ ) {
      hist->Fill( rn.Gaus(mean,sigma), 1 );
   }
   return hist;
}

double linTransform(double x)
{
  return (x-172.5)/10;
}


// __________________________________________________________________________________ //
//!
//!  MakeHistogram
//!
//!  make a histogram from an Eigen::Vector for plotting purposes
//!
TH1D* LTF_ROOTTools::MakeHistogram(const Eigen::VectorXd& values, vector<double> bins, const std::vector<std::pair<std::string,Eigen::MatrixXd > >& V ) 
{
   TH1D* hist = bins.empty() ?
      new TH1D("hist","hist",values.size(),0,values.size() ) :
      new TH1D("hist","hist",bins.size()-1, &bins[0]);
   if ( bins.size() && int(values.size()+1) != int(bins.size()) ) {cout<<"ERROR! binning and number of entries does not fit!"<<endl;exit(1);}
   for ( size_t i = 0 ;i<bins.size()-1; i++ ) {
      hist->SetBinContent(i+1, values(i));
      if ( V.size() ) {
         double e2 = 0;
         for ( auto apair : V ) {
            e2 += apair.second(i,i);
         }
         hist->SetBinError(i+1, sqrt(e2) );
      }
   }
   return hist;
}

// __________________________________________________________________________________ //
//!
//!  MakeHistogram
//!
//!  make histogram from an Eigen::Vector including stat. uncertainties of the template for plotting purposes
//!
TH1D* LTF_ROOTTools::MakeHistogram(const Eigen::VectorXd& values, const Eigen::VectorXd& errors, vector<double> bins )
{
   TH1D* hist = bins.empty() ?
      new TH1D("hist","hist",values.size(),0,values.size() ) :
      new TH1D("hist","hist",bins.size()-1, &bins[0]);
   if ( bins.size() && int(values.size()+1) != int(bins.size()) ) {cout<<"ERROR! binning and number of entries does not fit!"<<endl;exit(1);}
   if ( bins.size() && errors.size() && int(errors.size()+1) != int(bins.size()) ) {cout<<"ERROR! binning and number of errors does not fit!"<<endl;exit(1);}
   for ( size_t i = 0 ;i<bins.size()-1; i++ ) {
      hist->SetBinContent(i+1, values(i));
      if ( errors.size() > 0 ) hist->SetBinError(i+1, errors(i));
   }
   return hist;
}


// __________________________________________________________________________________ //
//!
//!  GetSigCheb2
//!
double LTF_ROOTTools::GetSigCheb2(const vector<double>& xvals, const vector<double>& yvals, vector<double> yerr) {
   //! Calculate the significance of the 2nd parameter
   //! of a fitted Chebychev polynomial of 2nd order
   //! to a set of data points with/without uncertainties
   //!
   //! Input
   //!    xvals:  x-values of the graph
   //!    yvals:  y-values of the graph
   //!    yerr:   uncertainties of the values (optional)
   //!            All values are considered to uncorrelated
   //!
   //! return
   //!    significance of 2nd Chebyshev parameter 

   const bool ErrorIsCorrelated = false ;
   
   // ---------------------------------------------- //
   // convert monomial to Chebyshev base
   static const TMatrixD MtoC(3,3, &vector<double>{
      1. , 0., 0.5, 
      0. , 1., 0. , 
      0. , 0., 0.5, 
      }[0]);
   static const TMatrixD MtoCT(TMatrixD::kTransposed,MtoC);

   // ---------------------------------------------- //
   // check input
   const int npar = int(xvals.size());
   if ( int(yvals.size()) != npar ) { cout<<"ERROR in GetSigCheb2()! Number of x and y values is not identical. Exiting..."<<endl; exit(1); }
   if ( int(yerr.size()) && ( int(yerr.size()) != npar ) ) { cout<<"ERROR in GetSigCheb2()! Number of y-errors and y-values is not identical. Exiting..."<<endl; exit(1); }

   
   // ---------------------------------------------- //
   // construct fit matrix
   TMatrixD M2(npar,3);
   for ( int i = 0 ; i<npar ; i ++ ) {
      M2(i,0) = 1.;
      M2(i,1) = xvals[i];
      M2(i,2) = xvals[i]*xvals[i];
   }
   
   TMatrixD M2T(TMatrixD::kTransposed,M2);


   TMatrixD yval(npar, 1, &yvals[0]);
		 
   // ---------------------------------------------- //
   // construct covariance matrix
   TMatrixD V(npar,npar);
   for ( int i = 0 ; i<npar ; i++ ) {
      if ( yerr.size() ) {
	 if ( !ErrorIsCorrelated )
	    V(i,i) = yerr[i]*yerr[i];
	 else {
	    cout<<"Error in GetSigCheb2()! Fit is undefined, when the errors are fully correlated."<<endl; exit(1);
	    for ( int j = 0 ; i<npar ; i++ )
	       V(j,i) = yerr[i]*yerr[j];
	 }
      }
      else
	 V(i,i) = 1.;
   }

   TMatrixD Vinv(V); Vinv.Invert();

   // ---------------------------------------------- //
   // fit matrix
   TMatrixD M2x = (M2T*Vinv*M2).Invert() * M2T*Vinv;
   TMatrixD M2xT (TMatrixD::kTransposed,M2x);

   // ---------------------------------------------- //
   // fit results
   TMatrixD VV(M2x * V * M2xT );
   TMatrixD chat(M2x * yval);
   //cout<<"Before transform p2 = "<<chat(2,0)<<" +/- "<<sqrt(VV(2,2))<<" ratio "<<chat(2,0)/sqrt(VV(2,2))<<endl;
   // ---------------------------------------------- //
   // fit results in Chebyshev base
   TMatrixD VVC(MtoC * VV * MtoCT ); // or is it (MtoCT * VV * MtoC )  ??
   TMatrixD cHatC(MtoC * chat);
         
   // ---------------------------------------------- //
   // significance of 2nd Chebyshev parameter
   double c2    = cHatC(2,0);
   double c2err = sqrt(VVC(2,2));
   //cout<<"After transform: p2 = "<<c2<<" / "<<c2err<<" ratio "<<c2/c2err<<endl;
   // ---------------------------------------------- //
   // return significance
   return c2/c2err;
}

// __________________________________________________________________________________ //
//!
//! GetSigCheb3
//!
double LTF_ROOTTools::GetSigCheb3(const vector<double>& xvals, const vector<double>& yvals, vector<double> yerr = vector<double>()) {
   //! Calculate the significance of the 3rd parameter
   //! of a fitted Chebychev polynomial of 2nd order
   //! to a set of data points with/without uncertainties
   //!
   //! Input
   //!    xvals:  x-values of the graph
   //!    yvals:  y-values of the graph
   //!    yerr:   uncertainties of the values (optional)
   //!            All values are considered to uncorrelated
   //!
   //! return
   //!    significance of 2nd Chebyshev parameter


   const bool ErrorIsCorrelated = false ;

   // ---------------------------------------------- //
   // check input
   const int npar = int(xvals.size());
   if ( int(yvals.size()) != npar ) { cout<<"ERROR in GetSigCheb3()! Number of x and y values is not identical. Exiting..."<<endl; exit(1); }
   if ( int(yerr.size()) && ( int(yerr.size()) != npar ) ) { cout<<"ERROR in GetSigCheb3()! Number of y-errors and y-values is not identical. Exiting..."<<endl; exit(1); }
   if ( int(yerr.size()) < 3 ) { cout<<"ERROR in GetSigCheb3()! Too few input values for a third order fit.. Exiting..."<<endl; exit(1); }

   TMatrixD yval(npar, 1, &yvals[0]);
   
   // ---------------------------------------------- //
   // construct fit matrix
   TMatrixD CM3(npar,4);
   for ( int i = 0 ; i<npar ; i ++ ) {
      CM3(i,0) = 1.; // T0
      CM3(i,1) = xvals[i]; // T1
      CM3(i,2) = 2*xvals[i]*xvals[i] - 1; //T2
      CM3(i,3) = 4*xvals[i]*xvals[i]*xvals[i] - 3*xvals[i]; //T3
   }
   TMatrixD CM3T(TMatrixD::kTransposed,CM3);


   // ---------------------------------------------- //
   // construct covariance matrix
   TMatrixD V(npar,npar);
   for ( int i = 0 ; i<npar ; i++ ) {
      V(i,i) = yerr[i]*yerr[i];
   }
   TMatrixD Vinv(V); Vinv.Invert();

   // ---------------------------------------------- //
   // fit in Cheb space
   TMatrixD CM3x = (CM3T*Vinv*CM3).Invert() * CM3T*Vinv;
   TMatrixD CM3xT (TMatrixD::kTransposed,CM3x);  
   
   TMatrixD C3VVC( CM3x * V * CM3xT );
   TMatrixD C3chat( CM3x * yval);
   
   // ---------------------------------------------- //
   // significance of 2nd Chebyshev parameter
   double c3    = C3chat(3,0);
   double c3err = sqrt(C3VVC(3,3));
   return c3/c3err;
}



// __________________________________________________________________________________ // 
//!
//!  makeErrorPlot
//!
//!  Add one line 
//!
double LTF_ROOTTools::makeErrorPlot(TCanvas& c, const char* ps_name, const char* title, const LTF::LiTeFit& fit, const vector<string> &uncertainties) {
   bool useNuisanceParameter = true;
   int nPar = 1; //M.cols()-1;
   double sum_error = 0; // this needs to be a vector in the case of more than 1 parameter

   if ( !useNuisanceParameter ) {
     for ( int i = 0 ; i<nPar ; i++ ) {
         TH1D* h  = new TH1D(title, title, uncertainties.size()+1, 0, uncertainties.size()+1);
	 
         for ( const string &source: uncertainties ) {
            double error = 0;
            //if (source.find("stat.")!= std::string::npos ) error = std::sqrt(fabs(fit.Vsource.find(source)->second(1,1)));
            //else error = std::sqrt(fabs(fit.Vsource.find("unfolding_error_"+variable+"_direct_envelope_"+source+"__1up")->second(i,i)));
            h->Fill(source.c_str(), error);
            sum_error += pow(error,2);
         }
         h->SetBinContent(h->GetNbinsX(), std::sqrt(sum_error));
         h->GetXaxis()->SetBinLabel(h->GetNbinsX(), "Total unc.");
         h->SetBarWidth(0.85);
         h->GetYaxis()->SetTitle("Uncertainty [GeV]");
         h->GetXaxis()->SetTickLength(0);
         h->Draw("hbar");
      }
   }
   else {
      for ( int i = 0 ; i<nPar ; i++ ) {
	 c.SetLeftMargin(0.6);
	 

	 TH1D* h  = new TH1D(title, "", uncertainties.size()+1, 0, uncertainties.size()+1);
         for ( const string &source: uncertainties ) {
           double error = 0;
	   if (source.find("STAT_DATA")!= std::string::npos ) error =std::sqrt(fit.Vsource.find(source)->second(0,0));
	   else if (source.find("STAT_MC")!= std::string::npos ) error =std::sqrt(fit.Vsource.find(source)->second(0,0));
	   else if (source.find("pseudoDataStat")!= std::string::npos ) error =std::sqrt(fit.Vsource.find(source)->second(0,0));
	   else    error = fabs(fit.DeltaSys.find(source)->second(i));
	   
           h->Fill(source.c_str(), error);
           sum_error += pow(error,2);
         }
	 if ( uncertainties.size() > 30 ) {
	   h->GetXaxis()->SetLabelSize(6);
           h->GetXaxis()->SetTitleSize(6);
	   h->GetYaxis()->SetLabelSize(6);
           h->GetYaxis()->SetTitleSize(6);
	   h->GetYaxis()->SetTitleOffset(1.2);
	   h->GetYaxis()->SetRangeUser(0.0, 1.5);
         }
         else {
	   h->GetYaxis()->SetRangeUser(0.0, 1.5);
	   h->GetYaxis()->SetTitleOffset(1.2);
	 }
         h->SetBinContent(h->GetNbinsX(), std::sqrt(sum_error));
         h->GetXaxis()->SetBinLabel(h->GetNbinsX(), "Total unc.");
         h->SetBarOffset(0.1);
         h->SetBarWidth(0.8);
         h->GetYaxis()->SetTitle("Uncertainty [GeV]");
         h->GetXaxis()->SetTickLength(0);
         h->Draw("hbar");

	 c.Print(ps_name);
         c.Clear();
	 //c.SetLeftMargin(0.6);
	 int nbins = uncertainties.size()+2;
	 double x[nbins], y[nbins], xerr[nbins], xerr2[nbins], yerr[nbins];
	 for (int i = 0; i <= nbins; i++) {
	   y[i] = h->GetBinCenter(i);
	   x[i] = 0;
	   xerr[i] = 1.0;
	   xerr2[i] = 2.0;
	   yerr[i] = 1.0;
	 }
	 TGraphErrors *band = new TGraphErrors(nbins, x, y, xerr, yerr);
	 band->SetFillColor(kGreen);
	 TGraphErrors *band2 = new TGraphErrors(nbins, x, y, xerr2, yerr);
         band2->SetFillColor(kYellow);

	 TH1D* h1  = new TH1D("NP", "", uncertainties.size(), 0, uncertainties.size());
         TGraphErrors* g = new TGraphErrors(uncertainties.size());
         for ( long unsigned int j = 0; j < uncertainties.size(); j++ ) {
	    h1->Fill(uncertainties[j].c_str(), fit.map_nuisance.find(uncertainties[j])->second.first);
            h1->SetBinError(j,fit.map_nuisance.find(uncertainties[j])->second.second);
            g->SetPoint(j, fit.map_nuisance.find(uncertainties[j])->second.first, j+0.5);
	    g->SetPointError(j, fit.map_nuisance.find(uncertainties[j])->second.second, 0);
	    //g->SetPoint(j, j+0.5, fit.map_nuisance.find(uncertainties[j])->second.first);
            //g->SetPointError(j, 0, fit.map_nuisance.find(uncertainties[j])->second.second);
	    //g->SetPointError(j, fit.map_nuisance.find(uncertainties[j])->second.second / fit.DeltaSys.find(uncertainties[j])->second(i), 0);
         }
         //h1->SetBinContent(h1->GetNbinsX(), 0.);
         //h1->GetXaxis()->SetBinLabel(h1->GetNbinsX(), "");
         gStyle->SetHistMinimumZero();
	 
         h1->SetBarOffset(0.95);
         h1->SetBarWidth(0);
         h1->SetLineColor(10);
         h1->SetFillColor(10);
	 h1->SetMarkerColor(10);
         h1->SetLineColorAlpha(10,0);
	 g->SetMarkerStyle(20);
         g->SetMarkerColor(kBlue);
	 g->SetMarkerSize(0.8);
	 if ( uncertainties.size() > 30 ) {
           h1->GetXaxis()->SetLabelSize(6);
           h1->GetXaxis()->SetTitleSize(6);
           h1->GetYaxis()->SetLabelSize(6);
           h1->GetYaxis()->SetTitleSize(6);
           h1->GetYaxis()->SetTitleOffset(1.2);
           g->SetMarkerSize(0.5);
	 }
	 
	 h1->GetYaxis()->SetTitle("Nuisance parameter");
         h1->GetXaxis()->SetTickLength(0);
	 h1->GetYaxis()->SetRangeUser(-3,3);

         h1->Draw("hbar e");
	 band2->Draw("same 2");
	 band->Draw("same 2");
	 g->Draw("same PE");
	 //band->Draw("same 3");//("3 SAME");
	 //g->Draw("same PE");

         TLine* l = new TLine(0, 0, 0, h1->GetNbinsX());
         l->SetLineStyle(3);
         l->Draw("same");
	 gPad->RedrawAxis();
	 //gPad->RedrawAxis("G");
	 TLine l1;
	 l1.DrawLine(-3, h1->GetXaxis()->GetXmax (), 3, h1->GetXaxis()->GetXmax ());
	 
	 c.Print(ps_name);
         c.Clear();
         h->Delete();
         h1->Delete();
         g->Delete();
	 band->Delete();
	 band2->Delete();
	 l->Delete();
	 l1.Delete();
      }
   }
   return std::sqrt(sum_error);
}

void LTF_ROOTTools::makeErrorPlotDilepton(TCanvas& c1, const char* ps_name, const LTF::LiTeFit& fit) {

  c1.Clear();
  c1.SetLogy(0);
  c1.SetLeftMargin(0.2);

  vector<string> lepton_uncertainties = {"Electron energy resolution",
					 "Electron energy scale AF2",
					 "Electron energy scale",
					 "Muon ID momentum resolution",
					 "Muon MS momentum resolution",
					 "Muon sagitta residual bias",
					 "Muon sagitta rho",
					 "Muon momentum scale",
					 "Electron ID efficiency",
					 "Electron isolation efficiency",
					 "Electron reconstruction efficiency",
					 "Electron trigger efficiency",
					 "Muon ID efficiency (stat.)",
					 "Muon ID efficiency (syst.)",
					 "Muon isolation efficiency (stat.)",
					 "Muon isolation efficiency (syst.)",
					 "Muon trigger efficiency (stat.)",
					 "Muon trigger efficiency (syst.)",
					 "Muon TTVA efficiency (stat.)",
					 "Muon TTVA efficiency (syst.)"};
  
  vector<string> jes_uncertainties = {"Jet vertex tagger efficiency",
                                      "b-jet energy scale",
				      "Jet energy scale (detector NP1)",
				      "Jet energy scale (detector NP2)",
				      "Jet energy scale (mixed NP1)",
				      "Jet energy scale (mixed NP2)",
				      "Jet energy scale (mixed NP3)",
				      "Jet energy scale (modelling NP1)",
				      "Jet energy scale (modelling NP2)",
				      "Jet energy scale (modelling NP3)",
				      "Jet energy scale (modelling NP4)",
				      "Jet energy scale (statistical NP1)",
				      "Jet energy scale (statistical NP2)",
				      "Jet energy scale (statistical NP3)",
				      "Jet energy scale (statistical NP4)",
				      "Jet energy scale (statistical NP5)",
				      "Jet energy scale (statistical NP6)",
				      "Jet energy scale $\eta$ intercalib (modelling)",
				      "Jet energy scale $\eta$ intercalib (non-closure 2018 data)",
				      "Jet energy scale $\eta$ intercalib (non-closure high $E$)",
				      "Jet energy scale $\eta$ intercalib (non-closure neg $\eta$)",
				      "Jet energy scale $\eta$ intercalib (non-closure pos $\eta$)",
				      "Jet energy scale $\eta$ intercalib (stat.)",
				      "Jet energy scale (flavour composition)",
				      "Jet energy scale (flavour response)"};
  
  vector<string> jer_uncertainties = {"Jet energy resolution (data vs MC)",
				      "Jet energy resolution (NP1)",
				      "Jet energy resolution (NP10)",
				      "Jet energy resolution (NP11)",
				      "Jet energy resolution (NP12restTerm)",
				      "Jet energy resolution (NP2)",
				      "Jet energy resolution (NP3)",
				      "Jet energy resolution (NP4)",
				      "Jet energy resolution (NP5)",
				      "Jet energy resolution (NP6)",
				      "Jet energy resolution (NP7)",
				      "Jet energy resolution (NP8)",
				      "Jet energy resolution (NP9)"};
  vector<string> jvt_pileup_met_uncertainties = { "Jet energy scale (pileup offset)",
						  "Jet energy scale (pileup NPV)",
						  "Jet energy scale (pileup pT term)",
						  "Jet energy scale (pileup rho topology",
						  "Jet energy scale (punch-through)",
						  "Jet energy scale (single particle high pT)",
						  "MET soft term resolution (parallel)",
						  "MET soft term resolution (perpendicular)",
						  "MET soft term scale"};

 
  vector<string> b_tagging_uncertainties = {"b-jet efficiency (eigenvar 0)",
					    "b-jet efficiency (eigenvar 1)",
					    "b-jet efficiency (eigenvar 10)",
					    "b-jet efficiency (eigenvar 11)",
					    "b-jet efficiency (eigenvar 12)",
					    "b-jet efficiency (eigenvar 13)",
					    "b-jet efficiency (eigenvar 14)",
					    "b-jet efficiency (eigenvar 15)",
					    "b-jet efficiency (eigenvar 16)",
					    "b-jet efficiency (eigenvar 17)",
					    "b-jet efficiency (eigenvar 18)",
					    "b-jet efficiency (eigenvar 19)",
					    "b-jet efficiency (eigenvar 2)",
					    "b-jet efficiency (eigenvar 20)",
					    "b-jet efficiency (eigenvar 21)",
					    "b-jet efficiency (eigenvar 22)",
					    "b-jet efficiency (eigenvar 23)",
					    "b-jet efficiency (eigenvar 24)",
					    "b-jet efficiency (eigenvar 25)",
					    "b-jet efficiency (eigenvar 26)",
					    "b-jet efficiency (eigenvar 27)",
					    "b-jet efficiency (eigenvar 28)",
					    "b-jet efficiency (eigenvar 29)",
					    "b-jet efficiency (eigenvar 3)",
					    "b-jet efficiency (eigenvar 30)",
					    "b-jet efficiency (eigenvar 31)",
					    "b-jet efficiency (eigenvar 32)",
					    "b-jet efficiency (eigenvar 33)",
					    "b-jet efficiency (eigenvar 34)",
					    "b-jet efficiency (eigenvar 35)",
					    "b-jet efficiency (eigenvar 36)",
					    "b-jet efficiency (eigenvar 37)",
					    "b-jet efficiency (eigenvar 38)",
					    "b-jet efficiency (eigenvar 39)",
					    "b-jet efficiency (eigenvar 4)",
					    "b-jet efficiency (eigenvar 40)",
					    "b-jet efficiency (eigenvar 41)",
					    "b-jet efficiency (eigenvar 42)",
					    "b-jet efficiency (eigenvar 43)",
					    "b-jet efficiency (eigenvar 44)",
					    "b-jet efficiency (eigenvar 5)",
					    "b-jet efficiency (eigenvar 6)",
					    "b-jet efficiency (eigenvar 7)",
					    "b-jet efficiency (eigenvar 8)",
					    "b-jet efficiency (eigenvar 9)",
					    "c-jet efficiency (eigenvar 0)",
					    "c-jet efficiency (eigenvar 1)",
					    "c-jet efficiency (eigenvar 10)",
					    "c-jet efficiency (eigenvar 11)",
					    "c-jet efficiency (eigenvar 12)",
					    "c-jet efficiency (eigenvar 13)",
					    "c-jet efficiency (eigenvar 14)",
					    "c-jet efficiency (eigenvar 15)",
					    "c-jet efficiency (eigenvar 16)",
					    "c-jet efficiency (eigenvar 17)",
					    "c-jet efficiency (eigenvar 18)",
					    "c-jet efficiency (eigenvar 19)",
					    "c-jet efficiency (eigenvar 2)",
					    "c-jet efficiency (eigenvar 3)",
					    "c-jet efficiency (eigenvar 4)",
					    "c-jet efficiency (eigenvar 5)",
					    "c-jet efficiency (eigenvar 6)",
					    "c-jet efficiency (eigenvar 7)",
					    "c-jet efficiency (eigenvar 8)",
					    "c-jet efficiency (eigenvar 9)",
					    "light-jet efficiency (eigenvar 0)",
					    "light-jet efficiency (eigenvar 1)",
					    "light-jet efficiency (eigenvar 10",
					    "light-jet efficiency (eigenvar 11",
					    "light-jet efficiency (eigenvar 12",
					    "light-jet efficiency (eigenvar 13",
					    "light-jet efficiency (eigenvar 14",
					    "light-jet efficiency (eigenvar 15",
					    "light-jet efficiency (eigenvar 16",
					    "light-jet efficiency (eigenvar 17",
					    "light-jet efficiency (eigenvar 18",
					    "light-jet efficiency (eigenvar 19",
					    "light-jet efficiency (eigenvar 2)",
					    "light-jet efficiency (eigenvar 3)",
					    "light-jet efficiency (eigenvar 4)",
					    "light-jet efficiency (eigenvar 5)",
					    "light-jet efficiency (eigenvar 6)",
					    "light-jet efficiency (eigenvar 7)",
					    "light-jet efficiency (eigenvar 8)",
					    "light-jet efficiency (eigenvar 9)"};

  vector<string> bkgd_uncertainties = {"pileup reweighting",
				       "$ttV$ normalisation",
				       "Diboson normalisation",
				       "Fakes normalisation",
				       "Z+jets",
				       "ttbar normalisation",
				       "tW normalisation",
				       "tW DS vs DR",
				       "Luminosity"};
  
  vector<string> modelling_uncertainties = {"h_{damp}",
					    "FSR #mu_{R}",
					    "Scale #mu_{R}",
					    "Scale #mu_{F}",
					    "ISR #alpha_{S} Var3c",
					    "Parton shower",
					    "Matching",
					    "Recoil to top",
					    "PDF",
					    "top mass"};
  vector<string> all_uncertainties;
  all_uncertainties.insert(all_uncertainties.end(), lepton_uncertainties.begin(), lepton_uncertainties.end());
  all_uncertainties.insert(all_uncertainties.end(), jes_uncertainties.begin(), jes_uncertainties.end());

  std::map<string, double> error_summary;
  error_summary.insert({"Lepton",           makeErrorPlot(c1, ps_name, "Lepton uncertainties", fit, lepton_uncertainties)});
  error_summary.insert({"JES",              makeErrorPlot(c1, ps_name, "JES uncertainties", fit, jes_uncertainties)});
  error_summary.insert({"JER",              makeErrorPlot(c1, ps_name, "JER uncertainties", fit, jer_uncertainties)});
  error_summary.insert({"JVT+PileUp+MET",   makeErrorPlot(c1, ps_name, "JVT+PileUp+MET uncertainties", fit, jvt_pileup_met_uncertainties)});
  error_summary.insert({"b-tagging",        makeErrorPlot(c1, ps_name, "b-tagging uncertainties", fit, b_tagging_uncertainties)});
  error_summary.insert({"theory modelling", makeErrorPlot(c1, ps_name, "modelling uncertainties", fit, modelling_uncertainties)});
  error_summary.insert({"bkgd modelling",   makeErrorPlot(c1, ps_name, "background uncertainties", fit, bkgd_uncertainties)});
  makeErrorPlot(c1, ps_name, "all uncertainties", fit, all_uncertainties); 
  //error_summary.insert({"Stat.+Lumi",       makeErrorPlot(c1, ps_name, "statistical uncertainties", fit, other_uncertainties)});

  TH1D* h  = new TH1D("Full error breakdown", "", error_summary.size()+1, 0, error_summary.size()+1);
  double sum_error_sq = 0;
  for( auto& tmp_err: error_summary ) {
    h->Fill(tmp_err.first.c_str(), tmp_err.second);
    sum_error_sq += pow(tmp_err.second,2);
  }

  h->SetBinContent(h->GetNbinsX(), std::sqrt(sum_error_sq));
  h->GetXaxis()->SetBinLabel(h->GetNbinsX(), "Total unc.");
  h->SetBarWidth(0.85);
  //h->GetYaxis()->SetLabelSize(0.03);
  h->GetYaxis()->SetTitle("Uncertainty [GeV]");
  //h->GetXaxis()->SetLabelSize(0.02);
  h->GetXaxis()->SetTickLength(0);
  h->GetXaxis()->LabelsOption("<");
  h->Draw("hbar");
  c1.Print(ps_name);
  c1.Clear();
}

void LTF_ROOTTools::makeErrorPlotSingle(TCanvas& c1, const char* ps_name, const LTF::LiTeFit& fit) {

  c1.Clear();
  c1.SetLogy(0);
  c1.SetLeftMargin(0.2);

  bool doJesStressTest = false;
  
  vector<string> lepton_uncertainties = {"EG_RESOLUTION_ALL",
                                         "EG_SCALE_ALL",
                                         "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR",
                                         "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR",
                                         "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR",
                                         "EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR",
                                         "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR",
                                         "MUON_SAGITTA_DATASTAT",
                                         "MUON_SAGITTA_RESBIAS",
                                         "MUON_EFF_BADMUON_SYS",
                                         "MUON_EFF_ISO_STAT",
                                         "MUON_EFF_ISO_SYS",
                                         "MUON_EFF_RECO_STAT",
                                         "MUON_EFF_RECO_SYS",
                                         "MUON_EFF_TTVA_STAT",
                                         "MUON_EFF_TTVA_SYS",
                                         "MUON_EFF_TrigStatUncertainty",
                                         "MUON_EFF_TrigSystUncertainty"};

  vector<string> jes_uncertainties = {"JET_EffectiveNP_Detector1",
                                      "JET_EffectiveNP_Detector2",
                                      "JET_EffectiveNP_Mixed1",
                                      "JET_EffectiveNP_Mixed2",
                                      "JET_EffectiveNP_Mixed3",
                                      "JET_EffectiveNP_Modelling1",
                                      "JET_EffectiveNP_Modelling2",
                                      "JET_EffectiveNP_Modelling3",
                                      "JET_EffectiveNP_Modelling4",
                                      "JET_EffectiveNP_Statistical1",
                                      "JET_EffectiveNP_Statistical2",
                                      "JET_EffectiveNP_Statistical3",
                                      "JET_EffectiveNP_Statistical4",
                                      "JET_EffectiveNP_Statistical5",
                                      "JET_EffectiveNP_Statistical6",
                                      "JET_EtaIntercalibration_Modelling",
                                      "JET_EtaIntercalibration_NonClosure_2018data",
                                      "JET_EtaIntercalibration_NonClosure_highE",
                                      "JET_EtaIntercalibration_NonClosure_negEta",
                                      "JET_EtaIntercalibration_NonClosure_posEta",
                                      "JET_EtaIntercalibration_TotalStat",
                                      "JET_Flavor_Composition_prop",
                                      "JET_Flavor_Response_prop",
                                      "JET_Pileup_OffsetMu",
                                      "JET_Pileup_OffsetNPV",
                                      "JET_Pileup_PtTerm",
                                      "JET_Pileup_RhoTopology",
                                      "JET_PunchThrough_MC16"};
  if (doJesStressTest) jes_uncertainties.push_back("myJESUncertainty");
  
  vector<string> jer_uncertainties = {"JET_JER_DataVsMC_MC16_PseudoData",
                                      "JET_JER_EffectiveNP_1_PseudoData",
                                      "JET_JER_EffectiveNP_2_PseudoData",
                                      "JET_JER_EffectiveNP_3_PseudoData",
                                      "JET_JER_EffectiveNP_4_PseudoData",
                                      "JET_JER_EffectiveNP_5_PseudoData",
                                      "JET_JER_EffectiveNP_6_PseudoData",
                                      "JET_JER_EffectiveNP_7_PseudoData",
                                      "JET_JER_EffectiveNP_8_PseudoData",
                                      "JET_JER_EffectiveNP_9_PseudoData",
                                      "JET_JER_EffectiveNP_10_PseudoData",
                                      "JET_JER_EffectiveNP_11_PseudoData",
                                      "JET_JER_EffectiveNP_12restTerm_PseudoData"};

  vector<string> jvt_pileup_met_uncertainties = {"JET_JvtEfficiency",
                                                 "PRW_DATASF",
                                                 "MET_SoftTrk_ResoPara",
                                                 "MET_SoftTrk_ResoPerp",
                                                 "MET_SoftTrk_Scale",
                                                 "WMASS_VAR_signal"};

  vector<string> b_tagging_uncertainties = {"FT_EFF_Eigen_B_0",
                                            "FT_EFF_Eigen_B_1",
                                            "FT_EFF_Eigen_B_2",
                                            "FT_EFF_Eigen_B_3",
                                            "FT_EFF_Eigen_B_4",
                                            "FT_EFF_Eigen_B_5",
                                            "FT_EFF_Eigen_B_6",
                                            "FT_EFF_Eigen_B_7",
                                            "FT_EFF_Eigen_B_8",
                                            "FT_EFF_Eigen_C_0",
                                            "FT_EFF_Eigen_C_1",
                                            "FT_EFF_Eigen_C_2",
                                            "FT_EFF_Eigen_C_3",
                                            "FT_EFF_Eigen_Light_0",
                                            "FT_EFF_Eigen_Light_1",
                                            "FT_EFF_Eigen_Light_2",
                                            "FT_EFF_Eigen_Light_3",
                                            "FT_EFF_extrapolation",
                                            "FT_EFF_extrapolation_from_charm"};

  vector<string> modelling_uncertainties = {"THEORY_CROSS_SECTION_signal",
                                            "THEORY_SHOWERING_HERWIG7_signal",
                                            "THEORY_SCALE_FACTORISATION_signal",
                                            "THEORY_SCALE_RENORMALISATION_signal",
                                            "THEORY_ISR_signal",
                                            "THEORY_FSR_signal",
                                            "THEORY_HDAMP_signal",
                                            "THEORY_PTHARD_signal",
                                            "THEORY_TOPRECOIL_signal",
                                            "THEORY_TOP_MASS_signal",
                                            "THEORY_PDF4LHC_VARIATION_signal",
                                            "THEORY_DR_DS_signal"};

  vector<string> bkgd_uncertainties = {"THEORY_CROSS_SECTION_Wjets",
                                       "THEORY_SCALE_COMBINED_Wjets",
                                       "THEORY_PDF4LHC_VARIATION_Wjets",
                                       "THEORY_EWK_Wjets",
                                       "THEORY_SCALE_COMBINED_multiboson_noW",
                                       "THEORY_PDF4LHC_VARIATION_multiboson_noW",
                                       "THEORY_EWK_multiboson_noW",
                                       "THEORY_CROSS_SECTION_other_top_noWt",
                                       "FAKES_Electron",
                                       "FAKES_Muon"};

  vector<string> other_uncertainties = {"LUMINOSITY",
					"pseudoDataStat",
					"STAT_MC",
					"STAT_DATA"};

  vector<string> all_uncertainties;
  all_uncertainties.insert(all_uncertainties.end(), lepton_uncertainties.begin(), lepton_uncertainties.end());
  all_uncertainties.insert(all_uncertainties.end(), jes_uncertainties.begin(), jes_uncertainties.end());
  all_uncertainties.insert(all_uncertainties.end(), jer_uncertainties.begin(), jer_uncertainties.end());
  all_uncertainties.insert(all_uncertainties.end(), jvt_pileup_met_uncertainties.begin(), jvt_pileup_met_uncertainties.end());
  all_uncertainties.insert(all_uncertainties.end(), b_tagging_uncertainties.begin(), b_tagging_uncertainties.end());
  all_uncertainties.insert(all_uncertainties.end(), modelling_uncertainties.begin(), modelling_uncertainties.end());
  all_uncertainties.insert(all_uncertainties.end(), bkgd_uncertainties.begin(), bkgd_uncertainties.end());
  all_uncertainties.insert(all_uncertainties.end(), other_uncertainties.begin(), other_uncertainties.end());
  
  std::map<string, double> error_summary;
  error_summary.insert({"Lepton",           makeErrorPlot(c1, ps_name, "Lepton uncertainties", fit, lepton_uncertainties)});
  error_summary.insert({"JES",              makeErrorPlot(c1, ps_name, "JES uncertainties", fit, jes_uncertainties)});
  error_summary.insert({"JER",              makeErrorPlot(c1, ps_name, "JER uncertainties", fit, jer_uncertainties)});
  error_summary.insert({"JVT+PileUp+MET",   makeErrorPlot(c1, ps_name, "JVT+PileUp+MET uncertainties", fit, jvt_pileup_met_uncertainties)});
  error_summary.insert({"b-tagging",        makeErrorPlot(c1, ps_name, "b-tagging uncertainties", fit, b_tagging_uncertainties)});
  error_summary.insert({"theory modelling", makeErrorPlot(c1, ps_name, "modelling uncertainties", fit, modelling_uncertainties)});
  error_summary.insert({"bkgd modelling",   makeErrorPlot(c1, ps_name, "background uncertainties", fit, bkgd_uncertainties)});
  error_summary.insert({"Stat.+Lumi",       makeErrorPlot(c1, ps_name, "statistical uncertainties", fit, other_uncertainties)});
  makeErrorPlot(c1, ps_name, "all uncertainties", fit, all_uncertainties);
  
  TH1D* h  = new TH1D("Full error breakdown", "", error_summary.size()+1, 0, error_summary.size()+1);
  double sum_error_sq = 0;
  for( auto& tmp_err: error_summary ) {
    h->Fill(tmp_err.first.c_str(), tmp_err.second);
    sum_error_sq += pow(tmp_err.second,2);
  }

  h->SetBinContent(h->GetNbinsX(), std::sqrt(sum_error_sq));
  h->GetXaxis()->SetBinLabel(h->GetNbinsX(), "Total unc.");
  h->SetBarWidth(0.85);
  h->GetYaxis()->SetTitle("Uncertainty [GeV]");
  h->GetXaxis()->SetTickLength(0);
  h->GetXaxis()->LabelsOption("<");
  h->Draw("hbar");
  c1.Print(ps_name);
  c1.Clear();
}




// __________________________________________________________________________________ //
//!
//!
//!  Plot a LiTeFit object using ROOT
//!
//!  The binning needs to be provided to the plotting function,
//!  since this is not included in LTF::LiTeFit
//! 
void LTF_ROOTTools::plotLiTeFit(LTF::LiTeFit& fit, const vector<double>& bins, 
				const char*   ps_name,
				const string& yaxistitle,
				const string& xaxistitle,
				const string& referencename)
{
   gStyle->SetOptStat(0);
   gStyle->SetLabelFont(43, "XYZ");
   gStyle->SetTitleFont(43, "XYZ");
   gStyle->SetLegendFont(43);
   gStyle->SetLabelSize(14, "XYZ");
   gStyle->SetTitleSize(14, "XYZ");
   gStyle->SetTitleFontSize(14);
   gStyle->SetLegendTextSize(14);
   gStyle->SetTextFont(43);
   gStyle->SetTextSize(14);
   gROOT->ForceStyle();

   gSystem->mkdir("plots");
   auto& M = fit.M;
   
   // sanity check
   if ( M.cols() != 2 ) {cout<<"Error! only 1-dim plotting implemented."<<endl;exit(1);}
   Eigen::VectorXd reference_values = M.col(1);
   //TH1D* hist = new TH1D("hist","hist",bins.size()-1, &bins[0]);
   
   //const vector<string> var_name_latex = {"m_{bl}", "m_{bW}", "m_{Wbbl}", "m_{Wb,bl}^{minimax}", "#Delta R(b,l)", "#Delta R(b,W)", "p_{T}(l)",  "p_{T}(b_{1})", "m(W_{had})", "y(W_{had})"};

   map<double,TH1D*> templates;
   for ( int iref = 0 ; iref<reference_values.size() ; iref++ ) {
     templates[iref] = MakeHistogram(fit.Y.col(iref), fit.SysY.at("statY").col(iref), bins); // templates including stat. uncertainties
   }

   TH1D* data    = MakeHistogram(fit.Dt,bins,fit.Vs);
   TH1D* TheoFit = MakeHistogram(fit.TheoFit,bins);
   if (fit.GetLogNormal()) for (int i =1; i<=TheoFit->GetNbinsX();i++) TheoFit->SetBinContent(i, std::log(TheoFit->GetBinContent(i)));
   
   TCanvas c1("c1","LTF plots",800,800);
   c1.SetRightMargin(0.05);
   c1.SetLeftMargin(0.15);
   c1.SetTopMargin(0.08);
   // c1.SetRightMargin(0.02);

   //const char* ps_name = fit.GetLogNormal() ?
   //   "plots/LTFlog_plots.ps" :
   //   "plots/LTF_plots.ps";
   c1.Print( (string(ps_name)+"[").c_str() );
   c1.Divide(1, 2, 0, 0);
   c1.cd(1)->SetBottomMargin(0);
   c1.cd(2)->SetTopMargin(0);
   if (!fit.GetLogNormal()) c1.cd(1)->SetLogy();
   // ---------------------------------------------- //
   // main plot
   // ---------------------------------------------- //

   data->SetMarkerStyle(20);
   data->SetMarkerSize(1.4);
   data->SetLineColor(kBlack);
   vector<Color_t> colors = {kCyan-9, kGray+1, kYellow-6, kOrange+7, kOrange-7, kViolet+1, kGray, kRed+1, kOrange-3, kBlue+7};
   for ( int iref = 0 ; iref<reference_values.size() ; iref++ ) {
      c1.cd(1);
      templates[iref]->SetLineWidth(2);
      if ( iref == 0 ) {
         templates[0]->SetTitle((";"+xaxistitle+";"+yaxistitle).c_str());
         templates[0]->SetLineColor(kRed+1);
         //if ( templates[0]->GetMaximum()>0 )templates[0]->SetMinimum(0);
         if ( !fit.GetLogNormal() ) {
	   templates[0]->SetMinimum(templates[0]->GetMinimum()*0.5);
	   templates[0]->SetMaximum(templates[0]->GetMaximum()*100);
	 }
	 else {
	   templates[0]->SetMinimum(templates[0]->GetMinimum()*1.2); // if cross section in one bin is smaller than 1
	   templates[0]->SetMaximum(templates[0]->GetMaximum()*3.);
	 }
         templates[0]->SetLineWidth(3);
         templates[0]->DrawClone("hist E");
      }
      else if ( iref==reference_values.size()-1) {
	 //if ( templates[0]->GetMaximum()>0 ) templates[iref]->SetFillColorAlpha(kBlue,0.15);
         templates[iref]->SetLineColor(kBlue+2);
         templates[iref]->SetLineWidth(3);
         templates[iref]->Draw("histsame E");
      }
      else {
	 templates[iref]->SetLineColor(colors[iref-1]);
         templates[iref]->SetLineWidth(2);
         templates[iref]->Draw("histsame E");
      }
      c1.cd(2);
      TH1D* tmp = (TH1D*)templates[iref]->Clone("tmp");
      for ( int i = 1; i <= data->GetNbinsX(); i++ ) {
	tmp->SetBinError(i, templates[iref]->GetBinError(i) / templates[iref]->GetBinContent(i) );
        tmp->SetBinContent(i, templates[iref]->GetBinContent(i) / TheoFit->GetBinContent(i) );
	//tmp->SetBinContent(i, templates[iref]->GetBinContent(i) / data->GetBinContent(i) );
      }
      if ( iref == 0 ) {
	tmp->GetYaxis()->SetRangeUser(0.45, 1.65);
	tmp->GetYaxis()->CenterTitle();
	tmp->GetYaxis()->SetTitle("Ratio to best model");
      }
      tmp->Draw("hist same E");
   }   
   TLine *line1 = new TLine(templates[0]->GetXaxis()->GetXmin(), 1.0, templates[0]->GetXaxis()->GetXmax(), 1.0);
   line1->SetLineColor(kBlack);
   line1->SetLineStyle(2);
   line1->SetLineWidth(2);
   line1->Draw("same");
   TH1D* data_clone = (TH1D*)data->Clone("data_clone");
   for ( int i = 1; i <= data->GetNbinsX(); i++ ) {
     data_clone->SetBinError(i, data_clone->GetBinError(i) / data_clone->GetBinContent(i) );
     data_clone->SetBinContent(i, data_clone->GetBinContent(i) / TheoFit->GetBinContent(i) );
   }
   data_clone->Draw("e0same");

   c1.cd(1);
   templates[0]->Draw("histsame");

   data->SetMarkerStyle(20);
   data->SetMarkerSize(1.4);
   data->SetLineColor(kBlack);
   data->Draw("e0same");

   TheoFit->SetLineColor(923);
   TheoFit->SetLineWidth(4);
   TheoFit->SetLineStyle(2);
   TheoFit->Draw("histsame");

   TLegend legend(0.18,0.74,0.94,0.97,"","NDC");
   legend.SetNColumns(3);
   legend.SetFillStyle(0);
   legend.SetBorderSize(0);
   legend.AddEntry(data,"Data","E0P");
   for ( int iref = 0 ; iref<reference_values.size() ; iref++ ) {
      legend.AddEntry(templates[iref],Form("Template m_{t}=%6.2f",reference_values[iref]),"FL");
   }
   legend.AddEntry(TheoFit,"Estimated best model","L");
   legend.Draw();

   TLatex text;
   text.SetNDC();
   text.SetTextAlign(13);
   text.DrawLatex(0.75,0.6,Form("m_{fit} =  %.2f +/- %.2f GeV",fit.ahat(0),fit.ahat_errorFit(0)));

   c1.Print(ps_name);
   c1.Clear();
   c1.SetRightMargin(0.02);
   c1.SetTopMargin(0.02);
   c1.SetLeftMargin(0.16);
   c1.SetBottomMargin(0.16);
   //c1.Print("plots/LTF_plot.pdf");


   // ---------------------------------------------- //
   // print relative size of all errors
   // ---------------------------------------------- //
   makeErrorPlotSingle(c1, ps_name, fit);
   //makeErrorPlotDilepton(c1, ps_name, fit);

   // ---------------------------------------------- //
   // print linear-functions in every bin
   // ---------------------------------------------- //
   c1.Clear();
   c1.SetRightMargin(0.02);
   c1.SetTopMargin(0.02);
   c1.SetLeftMargin(0.16);
   c1.SetBottomMargin(0.16);
   c1.SetLogy(false);

   gPad->SetTicky(1);

   //gStyle->SetLabelSize(0.05,"XYZ");
   //gStyle->SetTitleSize(0.05,"XYZ");
   gStyle->SetTitleOffset(1.4,"X");
   gStyle->SetTitleOffset(2.3,"Y");
   gStyle->SetMarkerSize(1.5);

   // Save ChiSq probabilities in separate file
   TFile* file = new TFile("fit_quality.root", "UPDATE");
   if (file->IsZombie()) {
        return;
    }

   TH1D* h_prob_linear = file->Get<TH1D>("h_chisq_prob_linear");
   TH1D* h_prob_quadratic = file->Get<TH1D>("h_chisq_prob_quadratic");
   TH1D* h_prob_ratio = file->Get<TH1D>("h_chisq_prob_ratio");
   TH1D* h_prob_ratio_rel = file->Get<TH1D>("h_chisq_prob_ratio_rel");
   TH1D* h_chisq_ratio = file->Get<TH1D>("h_chisq_ratio");
   TH1D* h_cheb2_sign= file->Get<TH1D>("h_cheb2_sign");
   TH1D* h_cheb3_sign= file->Get<TH1D>("h_cheb3_sign");
   TH1D* h_cheb_sign_ratio= file->Get<TH1D>("h_cheb_sign_ratio");
   TH2D* h_cheb2_sign_cheb3_sign= file->Get<TH2D>("h_cheb2_sign_cheb3_sign");
   TH2D* h_chisq_ratio_cheb2_sign= file->Get<TH2D>("h_chisq_ratio_cheb2_sign");

   bool doLinTransform = true;
   for ( int ibin = 0 ; ibin<fit.Dt.size() ; ibin++ ) {
      TGraphErrors* gdata = new TGraphErrors();
      if ( doLinTransform ) gdata->SetPoint(0,linTransform(fit.ahat(0)),fit.Dt(ibin));
      else gdata->SetPoint(0,linTransform(fit.ahat(0)),fit.Dt(ibin));
      gdata->SetPointError(0,0,data->GetBinError(ibin+1));
      gdata->SetMarkerStyle(20);

      TGraphErrors* graph;
      if (doLinTransform) {
	Eigen::VectorXd xvalues(fit.M.col(1).size());
	for( int i = 0; i < fit.M.col(1).size(); i++ ) {
	  xvalues(i) = linTransform(fit.M.col(1)(i));
	}
	graph = MakeTGraph(xvalues,ibin,fit.Y,fit.SysY);
      }
      else {
	graph = MakeTGraph(fit.M.col(1),ibin,fit.Y,fit.SysY);
      }
      
      TFitResultPtr resPol1 = graph->Fit("pol1","SQ0"); // weighted fit takes uncertainties into account, to neglect weights (set all weights to 1) use option "W"
      TMatrixDSym covPol1 = resPol1->GetCovarianceMatrix();
      TF1* pol1 = (TF1*)graph->GetFunction("pol1")->Clone("pol1");
      pol1->SetLineColor(kBlue+3);
      pol1->SetLineStyle(1);
      pol1->SetLineWidth(2);

      // Matrix for conversion to Chebyshev (orthogonal) polynomials
      //double matrixval[3][3] = {{1, 0, 0.5},{0, 1, 0}, {0, 0, 0.5}};
      //TMatrixD convToCheb(3,3,&matrixval[0][0]);
      //double matrixvalT[3][3] = {{1, 0, 0},{0, 1, 0}, {0.5, 0, 0.5}};
      //TMatrixD convToChebTrans(3,3,&matrixvalT[0][0]);
      //
      TFitResultPtr resPol2 = graph->Fit("pol2","SQ0");
      //TMatrixDSym covPol2 = resPol2->GetCovarianceMatrix();
      //cout<<"Covariancematrix"<<endl;
      //covPol2.Print();
     
      //TMatrixDSymEigen eig(covPol2);
      //TMatrixD V = eig.GetEigenVectors();
      //TMatrixD D = V.Invert()*covPol2*V;
      //TVectorD eigenvalues =eig.GetEigenValues();
      //cout<<"Matrix of Eigenvectors V"<<endl;
      //V.Print();
      //cout<<"Eigenvalues"<<endl;
      //eigenvalues.Print();
      //cout<<"Diagonal matrix with eigen values D"<<endl;
      //D.Print();
      //cout<<"V^-1*V should be unity"<<endl;
      //TMatrixD temp1 = V.Invert()*V;
      //temp1.Print();
      //cout<<"V*D*V^-1 should give A"<<endl;
      //TMatrixD temp2 = V * D * V.Invert();
      //temp2.Print();
      //TMatrixD temp3 = covPol2*V-D*V;
      //cout<<"Norm AV-VD "<<temp3.NormInf()<<" (should be zero)"<<endl;
      
      //eig.GetEigenValues().Print();
      //
      //cout<<"err0 "<<resPol2->ParError(0)<<" err1 "<<resPol2->ParError(1)<<" err2 "<<resPol2->ParError(2)<<endl;
      //cout<<"err02 "<<pow(resPol2->ParError(0),2)<<" err12 "<<pow(resPol2->ParError(1),2)<<" err22 "<<pow(resPol2->ParError(2),2)<<endl;

      //convToCheb.Print();
      //covPol2.Print();
      //convToChebTrans.Print();
      //
      //TMatrixD tmp =covPol2 * convToChebTrans;
      //tmp.Print();
      //
      //TMatrixD covCheb = convToCheb * covPol2 * convToChebTrans;
      //covCheb.Print();

      TF1* pol2 = (TF1*)graph->GetFunction("pol2")->Clone("pol2");
      pol2->SetLineColor(kOrange-3);
      pol2->SetLineStyle(1);
      pol2->SetLineWidth(2);
      

      //TF1* f2= new TF1("f2","[0]+[1]*pow(x,2)", -FLT_MIN,FLT_MAX );
      //graph->Fit(f2,"QW");
//      TF1* f1= new TF1("f1","[0]+[1]*pow(x,[2])", -FLT_MIN,FLT_MAX );
//      bool UseLTFOutput = true;
//      if ( UseLTFOutput ) {
//         Eigen::VectorXd ltfpol1param =(fit.Y*fit.Mc().transpose()).row(ibin);
//         //cout<<"M+*Y:"<<endl<< ltfpol1param <<endl;
//         pol1->SetParameter(0,ltfpol1param(0));
//         pol1->SetParameter(1,ltfpol1param(1));
//
//	 //graph->GetFunction("pol1")->SetParameter(0,ltfpol1param(0));
//         //graph->GetFunction("pol1")->SetParameter(1,ltfpol1param(1));
//	 //f1->SetParameter(0,ltfpol1param(0));
//         //f1->SetParameter(1,ltfpol1param(1));
//         //f1->SetParameter(2,fit.Gamma[0]);
//      }
//      if ( fit.Gamma[0]!=1 )   cout<<"Warning! Plotting with gamma factor !=1 not correctly implmeneted!"<<endl;

      graph->SetMarkerStyle(47);
      graph->SetMarkerColor(kRed+3);
      graph->SetLineColor(kRed+3);
      if ( fit.GetLogNormal() ) 
         graph->SetTitle((";"+referencename+";log("+yaxistitle+")").c_str()); // log(value/unit)
      else
         graph->SetTitle((";"+referencename+";"+yaxistitle).c_str());

      TF1* f1log = NULL;
      if ( fit.GetLogNormal() ) {
         //f1log = new TF1("pol1log","[0]+[1]*exp(x)", -FLT_MIN,FLT_MAX );
         f1log = new TF1("pol1log","log([0]+[1]*pow(x,1))", -FLT_MIN,FLT_MAX );

         TGraph* gexp = new TGraph();
         for ( int i = 0 ; i<graph->GetN() ; i++ ) 
            gexp->SetPoint(i,graph->GetX()[i], exp(graph->GetY()[i]));
         gexp->Fit("pol1","QW");

         Eigen::VectorXd ltfpol1param =(fit.Y*fit.Mc().transpose()).row(ibin);
         if ( graph->GetFunction("pol1") )
            graph->GetFunction("pol1")->SetParameter(0,ltfpol1param(0));

         f1log->SetParameter(0, gexp->GetFunction("pol1")->GetParameter(0));
         f1log->SetParameter(1, gexp->GetFunction("pol1")->GetParameter(1));
         if ( fit.Gamma[0]!=1 )   cout<<"Warning! Plotting with gamma factor !=1 not correctly implmeneted!"<<endl;

         // //f1 = new TF1("pol1","[0]+[1]*x", -FLT_MIN,FLT_MAX );
         // graph->Fit(f1log,"QW"); // "W": Ignore all point errors when fitting a TGraphErrors 
         f1log->SetLineColor(kBlue+1);
         f1log->SetLineStyle(7);
         f1log->SetLineWidth(2);
      }

      if ( gdata->GetY()[0] > 0 )
         graph->SetMinimum(0);

      graph->SetMaximum( max(gdata->GetY()[0],max(graph->GetY()[0],graph->GetY()[graph->GetN()-1]))*1.2);
      graph->GetYaxis()->SetRangeUser(min(gdata->GetY()[0],max(graph->GetY()[0],graph->GetY()[graph->GetN()-1]))*0.6,
				      max(gdata->GetY()[0],max(graph->GetY()[0],graph->GetY()[graph->GetN()-1]))*1.2);
      graph->Draw("APE0");
      if ( fit.GetLogNormal() )
         f1log->Draw("Lsame");
      else {
	//f1->Draw("Lsame");
	 pol1->Draw("Lsame");
	 pol2->Draw("Lsame");
      }
      //graph->GetFunction("pol1")->Draw("Lsame");
      if ( reference_values.size()+1<= 6 ) 
	graph->GetHistogram()->GetXaxis()->SetNdivisions(graph->GetN()+1);
      else
	graph->GetHistogram()->GetXaxis()->SetNdivisions(int(graph->GetN())+1+200);
      //f2->Draw("same");
      graph->Draw("PE0 same");
      gdata->Draw("PE0 same");

      double cheb2_sign = 0;
      double cheb3_sign = 0;
      {
        vector<double> xvals, yvals, yerr;
        for( int i = 0; i < fit.M.col(1).size(); i++ ) {
          xvals.push_back(linTransform(fit.M.col(1)(i)));
          yvals.push_back(fit.Y(ibin,i));
          double err = 0;
          for ( auto [name,Vy] : fit.SysY ) {
            err += pow(Vy(ibin,i),2);
	  }
          yerr.push_back(sqrt(err));
	}
        cheb2_sign = GetSigCheb2(xvals, yvals, yerr);
        cheb3_sign = GetSigCheb3(xvals, yvals, yerr);
      }

      
      //if ( ibin==0 ) {
      //double xmin = fit.GetLogNormal() ? 0.36 : 0.45;
      TLegend legend(0.55,0.68,0.96,0.94,"","NDC");
      legend.SetFillStyle(0);
      legend.SetBorderSize(0);
      legend.AddEntry(data,"Data","E0P");
      legend.AddEntry(graph,"Templates","PE0");
      if ( fit.GetLogNormal() ) {
	legend.AddEntry(graph->GetFunction("pol1"),"Linear log(model)","L");
	legend.AddEntry(f1log,"#scale[0.9]{Linearized model #scale[0.7]{(unused)}}","L");
      }
      else {
	legend.AddEntry(pol1,Form("Linearized model ( #chi^{2} / ndf = %.3f / %d)", pol1->GetChisquare(),graph->GetN()-2),"L");
	legend.AddEntry(pol2,Form("Quadratic model, ( #chi^{2} / ndf = %.3f / %d)", pol2->GetChisquare(),graph->GetN()-3),"L");
	legend.AddEntry(pol2,Form("p0 = %.3f +/- %.3f)", pol2->GetParameter(0), pol2->GetParError(0)),"");
	legend.AddEntry(pol2,Form("p1 = %.3f +/- %.3f)", pol2->GetParameter(1), pol2->GetParError(1)),"");
        legend.AddEntry(pol2,Form("p2 = %.3f +/- %.3f)", pol2->GetParameter(2), pol2->GetParError(2)),"");

	legend.AddEntry(pol2,Form(" (#chi^{2}_{quad}*ndf_{lin}) / (#chi^{2}_{lin}*ndf_{quad})  = %.3f )",
				  (pol2->GetChisquare()*(graph->GetN()-2))/(pol1->GetChisquare()*(graph->GetN()-3))),"");
	legend.AddEntry(pol2,Form("Cheb. 2: p2err/p2 = %.3f", cheb2_sign),"");
        legend.AddEntry(pol2,Form("Cheb. 3: p3err/p3 = %.3f", cheb3_sign),"");

	//legend.AddEntry(graph->GetFunction("pol1"),"Linearized model","L");
	//legend.AddEntry(pol1,"Weighted fit #scale[0.7]{(unused)}","L");
	//legend.AddEntry(pol2,"Weighted fit #scale[0.7]{(unused)}","L");
      }
      legend.DrawClone();
      
      TLatex text;
      text.SetNDC();
      text.SetTextAlign(11);
      //text.DrawLatex(0.20,0.93,Form("%3.1f_{ }<_{ }|y|_{ }<_{ }%3.1f",input_table["ylow"][ibin],input_table["yhigh"][ibin]));
      TString infotext = "_{ }<_{ }" + xaxistitle + "_{ }<_{ }";
      infotext.Prepend(Form("%3.0f", bins[ibin]));
      infotext.Append(Form("%3.0f", bins[ibin+1]));
      infotext.Append("_{}");
      text.DrawLatex(0.20,0.93, infotext);

      c1.Print(ps_name);


      h_prob_linear->Fill(resPol1->Prob());
      h_prob_quadratic->Fill(resPol2->Prob());
      h_prob_ratio->Fill(resPol1->Prob()/resPol2->Prob());
      h_prob_ratio_rel->Fill((resPol2->Prob()-resPol1->Prob())/resPol2->Prob());
      h_chisq_ratio->Fill((pol2->GetChisquare()*(graph->GetN()-2))/(pol1->GetChisquare()*(graph->GetN()-3)));
      if ( abs(cheb3_sign / cheb2_sign) > 1.5 ) h_cheb_sign_ratio->Fill(1.45);
      else h_cheb_sign_ratio->Fill(abs(cheb3_sign / cheb2_sign));
      if ( abs(cheb2_sign) > 2.9 ) cheb2_sign = 2.9;
      if ( abs(cheb3_sign) > 2.9 ) cheb3_sign = 2.9;
      h_cheb2_sign->Fill(cheb2_sign);
      h_cheb3_sign->Fill(cheb3_sign);
      h_cheb2_sign_cheb3_sign->Fill(cheb2_sign, cheb3_sign);
      h_chisq_ratio_cheb2_sign->Fill((pol2->GetChisquare()*(graph->GetN()-2))/(pol1->GetChisquare()*(graph->GetN()-3)), cheb2_sign);
     
    }
    h_prob_linear->Write("", TObject::kOverwrite);
    h_prob_quadratic->Write("", TObject::kOverwrite);
    h_prob_ratio->Write("", TObject::kOverwrite);
    h_prob_ratio_rel->Write("", TObject::kOverwrite);
    h_chisq_ratio->Write("", TObject::kOverwrite);
    h_cheb2_sign->Write("", TObject::kOverwrite);
    h_cheb3_sign->Write("", TObject::kOverwrite);
    h_cheb_sign_ratio->Write("", TObject::kOverwrite);
    h_cheb2_sign_cheb3_sign->Write("", TObject::kOverwrite);
    h_chisq_ratio_cheb2_sign->Write("", TObject::kOverwrite);
    file->Close();
    // ---------------------------------------------- //
    //   chisq plot
    // ---------------------------------------------- //
    c1.SetLogy(false);
    TGraphErrors* gChi2 = new TGraphErrors();
    TGraphErrors* gChi2_fit = new TGraphErrors();
    int ndf = (fit.Dt.rows()-(fit.M.cols()-1));
    for ( int itmpl = 0 ; itmpl<fit.chisq_y.size() ; itmpl++ )  {
      gChi2->SetPoint(itmpl, reference_values[itmpl], fit.chisq_y(itmpl)/ndf);
      gChi2->SetPointError(itmpl, 0, fit.chisq_y_error(itmpl)/ndf);
      gChi2_fit->SetPoint(itmpl, reference_values[itmpl], fit.chisq_fit(itmpl)/ndf);
      gChi2_fit->SetPointError(itmpl, 0, fit.chisq_y_error(itmpl)/ndf);
    }
    gChi2->Fit("pol2","Q");
    gChi2->SetMarkerStyle(20);
    gChi2->SetMarkerSize(1);
    gChi2_fit->SetMarkerStyle(20);
    gChi2_fit->SetMarkerSize(1);
    gChi2_fit->SetMarkerColor(kCyan-7);
    gChi2_fit->SetLineColor(kCyan-7);

    TGraphErrors* gChi2LTF = new TGraphErrors();
    gChi2LTF->SetPoint(0, fit.ahat(0), fit.chisq/ndf);
    gChi2LTF->SetPointError(0, 0, fit.chisq_error/ndf);
    gChi2LTF->SetMarkerStyle(29);
    gChi2LTF->SetMarkerSize(2.5);
    gChi2LTF->SetMarkerColor(kViolet+2);
    gChi2LTF->SetLineColor(kViolet+2);

    TGraph* gChi2chk = new TGraph();
    gChi2chk->SetPoint(0, fit.achk(0), fit.achk_chisq/ndf);
    gChi2chk->SetMarkerSize(2.2);
    gChi2chk->SetMarkerStyle(24);
    gChi2chk->SetMarkerColor(kRed);
    
    //gChi2LTF->Print("all");
    
    //gChi2->SetTitle(";#alpha_{0} [unit];#chi^{2}/ndf");
    gChi2->SetTitle((";"+referencename+";#chi^{2}/ndf").c_str());
    gChi2->SetMinimum(-1.0);
    gChi2->SetMaximum(5.0); //johannes set these values again once the fit is stable 

    gChi2->Draw("ape");
    if ( reference_values.size()+1<= 8 ) 
       gChi2->GetHistogram()->SetNdivisions(reference_values.size()+1+300,"X");
    else
       gChi2->GetHistogram()->SetNdivisions(int(reference_values.size()/2)+1+200,"X");
    
    TLine line;
    line.SetLineColor(920);
    line.SetLineStyle(3);
    line.DrawLine( 
	 gChi2->GetHistogram()->GetXaxis()->GetXmin(),(fit.chisq+1.)/ndf,
	 gChi2->GetHistogram()->GetXaxis()->GetXmax(),(fit.chisq+1)/ndf);

    //line.DrawLine( 
    //   fit.achk_chisq/ndf+1,1.,
    //   fit.achk_chisq/ndf+1,1.);
    gChi2chk->Draw("Psame");
    gChi2LTF->Draw("PEsame");
    gChi2_fit->Draw("PEsame");
    {
       TLegend legend(0.4,0.80,0.8,0.97,"","NDC");
       //legend.SetNColumns(3);
       legend.SetFillStyle(0);
       legend.SetBorderSize(0);
       legend.AddEntry(gChi2LTF,"#hat#chi^{2} of the Linear Template Fit","P");
       legend.AddEntry(gChi2,   "#chi^{2}_{#font[12]{j}} of the individual templates","P");
       legend.AddEntry(gChi2->GetFunction("pol2"),"Parabola","L");
       legend.AddEntry(gChi2chk,"Minimum of #chi^{2} parabola #scale[0.8]{(#check#chi^{2})}","P"); //  (#check#chi^{2})
       legend.DrawClone();
    }
    {
      TF1 *myfit = gChi2->GetFunction("pol2");
      TLatex latex;
      latex.SetNDC();
      text.SetTextAlign(11);
      latex.DrawLatex(0.47, 0.750, "Quadratic fit p_{0}+p_{1}*x+p_{2}*x^{2}");
      latex.DrawLatex(0.47, 0.725, Form("p_{0} = %.3f +/- %.3f", myfit->GetParameter(0), myfit->GetParError(0)));
      latex.DrawLatex(0.47, 0.700, Form("p_{1} = %.3f +/- %.3f", myfit->GetParameter(1), myfit->GetParError(1)));
      latex.DrawLatex(0.47, 0.675, Form("p_{2} = %.3f +/- %.3f", myfit->GetParameter(2), myfit->GetParError(2)));
      latex.DrawLatex(0.47, 0.65, Form("#chi^{2} / ndf = %.3f / %d", myfit->GetChisquare(), int(fit.chisq_y.size()-3)));
    }
    c1.Print(ps_name);
    //c1.Print( "plots/LTF_chi2.pdf");

    c1.Print( (string(ps_name)+"]").c_str() );
}


// __________________________________________________________________________________ //
//!
//!
//!  Plot a LiTeFit object using ROOT
//!
//!  The binning needs to be provided to the plotting function,
//!  since this is not included in LTF::LiTeFit
//! 
void LTF_ROOTTools::plotLiTeFit_2D(const LTF::LiTeFit& fit, const vector<double> bins )
{

   gStyle->SetOptStat(0);
   gSystem->mkdir("plots");
   auto& M = fit.M;
   if ( M.cols() != 3 ) {cout<<"Error! only 2-dim plotting implemented."<<endl;exit(1);}
   Eigen::VectorXd reference_values1 = M.col(1);
   Eigen::VectorXd reference_values2 = M.col(2);

   map<double,TH1D*> templates;
   set<double> r1,r2;
   for ( int iref = 0 ; iref<reference_values1.size() ; iref++ ) {
      //double ref = reference_values1(iref);
      templates[iref] = MakeHistogram(fit.Y.col(iref),bins);
      r1.insert(reference_values1(iref));
      r2.insert(reference_values2(iref));
   }
   
   if ( r1.size()<=1 ) {
      cout<<"ERROR! at least two distinct reference points for dimension 1 must be given"<<endl;
      exit(1);
   }
   if ( r2.size()<=1 ) {
      cout<<"ERROR! at least two distinct reference points for dimension 2 must be given"<<endl;
      exit(1);
   }

   double xmin = 2.* (*r1.begin())  - (*(++r1.begin()))*0.9999;
   double xmax = 2.* (*r1.rbegin()) - (*(++r1.rbegin()))*1.00001;
   double ymin = 2.* (*r2.begin())  - (*(++r2.begin()))*0.9999;
   double ymax = 2.* (*r2.rbegin()) - (*(++r2.rbegin()))*1.00001;
   
   TH1D* data    = MakeHistogram(fit.Dt,bins,fit.Vs);
   TH1D* TheoFit = MakeHistogram(fit.TheoFit,bins);
   
   TCanvas c1("c1","LTF plots",800,800);
   c1.SetRightMargin(0.05);
   c1.SetLeftMargin(0.15);
   c1.SetTopMargin(0.08);
   // c1.SetRightMargin(0.02);

   const char* ps_name = "LTF2D_plots.ps";
   c1.Print( (string(ps_name)+"[").c_str() );

   // ---------------------------------------------- //
   // main plot
   // ---------------------------------------------- //
   vector<Color_t> colors = {kCyan-9, kGray+1, kYellow-6, kOrange+7, kOrange-7, kViolet+1, kGray, kRed+1, kOrange-3, kBlue+7};
   //int colors = {kP10Blue, kP10Red, kP10Yellow, kP10Gray, kP10Violet, kP10Brown, kP10Orange, kP10Green, kP10Ash, kP10Cyan};
   cout<<colors[0]<<endl;
   for ( int iref = 0 ; iref<reference_values1.size() ; iref++ ) {
      templates[iref]->SetLineWidth(2);
      if ( iref == 0 ) {
         templates[0]->SetTitle("Linear Template Fit;Observable [unit]; Value [unit]");
         templates[0]->SetLineColor(kRed+2);
         templates[0]->SetMinimum(0);
         templates[0]->SetMaximum(templates[0]->GetMaximum()*1.7);
         templates[0]->SetLineWidth(3);
         templates[0]->DrawClone("hist");
      }
      else if ( iref==reference_values1.size()-1) {
         templates[iref]->SetFillColorAlpha(kBlue,0.15);
         templates[iref]->SetLineColor(kBlue+2);
         templates[iref]->SetLineWidth(3);
         templates[iref]->Draw("histsame");
      }
      else {
         //int color = iref+1;
	 //if (color >= 10 ) color=(color-10)*2+28;
	 templates[iref]->SetLineColor(colors[iref-1]);
	 //templates[iref]->SetFillColorAlpha(color,0.08);
         templates[iref]->Draw("histsame");
      }
   }   
   templates[0]->SetFillColorAlpha(kRed,0.15);
   templates[0]->Draw("histsame");

   data->SetMarkerStyle(20);
   data->SetMarkerSize(1.4);
   data->SetLineColor(kBlack);
   data->Draw("e0same");

   TheoFit->SetLineColor(923);
   TheoFit->SetLineWidth(4);
   TheoFit->SetLineStyle(2);
   TheoFit->Draw("histsame");

   TLegend legend(0.18,0.70,0.94,0.92,"","NDC");
   legend.SetNColumns(3);
   legend.SetFillStyle(0);
   legend.SetBorderSize(0);
   legend.AddEntry(data,"Data","E0P");
   for ( int iref = 0 ; iref<reference_values1.size() ; iref++ ) {
      legend.AddEntry(templates[iref],Form("Tpl. #alpha_{0}=%5.1f, #alpha_{1}=%3.1f",reference_values1[iref],reference_values2[iref]),"FL");
   }
   legend.AddEntry(TheoFit,"Estimated best model","L");
   legend.Draw();

   c1.Print(ps_name);
   c1.Print("plots/LTF2D_plot.pdf");


   // ---------------------------------------------- //
   // print linear-functions in every bin
   // ---------------------------------------------- //
   c1.Clear();
   c1.SetRightMargin(0.02);
   c1.SetTopMargin(0.02);

   c1.SetLeftMargin(0.12);
   c1.SetBottomMargin(0.12);

   for ( int ibin = 0 ; ibin<fit.Dt.size() ; ibin++ ) {
      //TGraphErrors* data  = MakeTGraph(fit.ahat.row(0),fit.Dt.row(ibin));
      TGraph2DErrors* gdata = new TGraph2DErrors();
      gdata->SetPoint(0, fit.ahat(0), fit.ahat(1), fit.Dt(ibin));
      gdata->SetPointError(0,0,0,data->GetBinError(ibin+1));
      gdata->SetMarkerStyle(20);
      gdata->SetMarkerSize(2);

      TGraph2DErrors* graph = MakeTGraph2D(fit.M.col(1),fit.M.col(2),ibin,fit.Y,fit.SysY);
      TF2* f2 = new TF2(Form("Linear Template Fit (2-dim., bin %d)",ibin),"[0]+[1]*x+[2]*y",
                        xmin,xmax,ymin,ymax
         );
      graph->Fit(f2,"QW");

      TGraph2D* gtheo = new TGraph2D();
      gtheo->SetPoint(0, fit.ahat(0), fit.ahat(1), f2->Eval(fit.ahat(0), fit.ahat(1)));
      gtheo->SetMarkerStyle(20);
      gtheo->SetMarkerSize(0.4);
      gtheo->SetMarkerColor(kBlack);
      for ( int it = 0 ; it<fit.M.col(1).size(); it++ ) {
         gtheo->SetPoint(it+1, reference_values1(it), reference_values2(it),
                         f2->Eval(reference_values1(it), reference_values2(it)));
      }

      f2->SetLineColor(922);
      f2->SetLineStyle(3);
      f2->SetLineWidth(1);
      f2->SetMinimum(0);
      f2->Draw("");
      f2->GetHistogram()->GetXaxis()->SetNdivisions(r1.size()+1);
      f2->GetHistogram()->GetYaxis()->SetNdivisions(r2.size()+1);
      f2->GetHistogram()->GetXaxis()->SetTitleOffset(2.2);
      f2->GetHistogram()->GetYaxis()->SetTitleOffset(2.2);
      f2->GetHistogram()->GetZaxis()->SetTitleOffset(1.8);
      f2->GetHistogram()->SetTitle(";Reference value 0 (#alpha_{0}) [unit]  ;Reference value 1 (#alpha_{1}) [unit];Value [unit]");

      graph->SetMarkerStyle(25);
      //graph->SetMarkerSize(2);
      graph->SetMarkerColor(kRed+3);
      graph->SetLineColor(kRed+3);
      graph->SetTitle(";Reference value (#alpha) [unit];Value [unit]");

      //graph->SetMinimum(0);
      //graph->SetMaximum( max(gdata->GetY()[0],max(graph->GetY()[0],graph->GetY()[graph->GetN()-1]))*1.2);
      //graph->Draw("ap");
      //pol1w->Draw("Lsame");
      //graph->GetFunction("pol1")->Draw("Lsame");
      //graph->GetHistogram()->GetXaxis()->SetNdivisions(graph->GetN()+1);

      f2->Draw("surf1");
      gtheo->Draw("P0 same");
      graph->Draw("err p0 same");
      gdata->Draw("err P0 same");


      TF2* f2leg = (TF2*)f2->Clone("f2leg");
      f2leg->SetLineWidth(3);
      gtheo->SetMarkerStyle(24);
      gtheo->SetMarkerSize(0.4);
      gdata->SetMarkerStyle(24);
      graph->SetMarkerStyle(24);
      if ( ibin==1 ) {
         TLegend legend(0.42,0.18,0.80,0.43,"","NDC");
         //legend.SetNColumns(3);
         legend.SetFillStyle(0);
         legend.SetBorderSize(0);
         //legend.SetTextSize(0.03);
         legend.AddEntry(graph,"Templates","PE0");
         legend.AddEntry(f2leg,"Linearized model","L");
         legend.AddEntry(gtheo,"Projections onto model","P");
         legend.AddEntry(gdata,"Data","E0P");
         legend.DrawClone();
      }

      TLatex text;
      text.SetNDC();
      //text.SetTextFont(43);
      text.SetTextAlign(13);
      text.SetTextSize(12);
      //text.DrawLatex(0.75,0.25,Form("Bin %d",ibin));
      text.DrawLatex(0.02,0.97,Form("Bin %d",ibin));

      c1.Print(ps_name);
      c1.Print( Form("plots/LTF2D_bin_%02d.pdf",ibin));
      
   }
   
   c1.Print( (string(ps_name)+"]").c_str() );
   
}






// __________________________________________________________________________________ //
//!
//!
//!  Plot a LiTeFit object using ROOT
//!
//!  The binning needs to be provided to the plotting function,
//!  since this is not included in LTF::LiTeFit
//! 
void LTF_ROOTTools::plotLiTeFitPol2Test(const LTF::LiTeFit& fit, const vector<double>& bins, 
                 const string& yaxistitle,
                 const string& referencename,
                 const string& observablename) 
{
   gStyle->SetOptStat(0);
   gSystem->mkdir("plots");
   auto& M = fit.M;

   // sanity check
   if ( M.cols() != 2 ) {cout<<"Error! only 1-dim plotting implemented."<<endl;exit(1);}
   Eigen::VectorXd reference_values = M.col(1);
   
   map<double,TH1D*> templates;
   for ( int iref = 0 ; iref<reference_values.size() ; iref++ ) {
      templates[iref] = MakeHistogram(fit.Y.col(iref),bins);
   }

   TH1D* data    = MakeHistogram(fit.Dt,bins,fit.Vs);
   TH1D* TheoFit = MakeHistogram(fit.TheoFit,bins);
   
   TCanvas c1("c1","LTF plots",800,800);
   c1.SetRightMargin(0.05);
   c1.SetLeftMargin(0.15);
   c1.SetTopMargin(0.08);

   const char* ps_name = fit.GetLogNormal() ?
      "LTFlog_plots.ps" :
      "LTF_plots.ps";
   c1.Print( (string(ps_name)+"[").c_str() );

   // ---------------------------------------------- //
   // main plot
   // ---------------------------------------------- //
   for ( int iref = 0 ; iref<reference_values.size() ; iref++ ) {
      templates[iref]->SetLineWidth(2);
      if ( iref == 0 ) {
         templates[0]->SetTitle(("Quadratic Template Fit;"+observablename+";"+yaxistitle).c_str());
         templates[0]->SetLineColor(kRed+1);
         if ( templates[0]->GetMaximum()>0 )templates[0]->SetMinimum(0);
         if ( !fit.GetLogNormal() ) 
            templates[0]->SetMaximum(templates[0]->GetMaximum()*1.7);
         // else 
         //    templates[0]->SetMaximum(templates[0]->GetMaximum()*3.); 
         templates[0]->SetLineWidth(3);
         templates[0]->DrawClone("hist");
      }
      else if ( iref==reference_values.size()-1) {
         if ( templates[0]->GetMaximum()>0 ) templates[iref]->SetFillColorAlpha(kBlue,0.15);
         templates[iref]->SetLineColor(kBlue+2);
         templates[iref]->SetLineWidth(3);
         templates[iref]->Draw("histsame");
      }
      else {
         templates[iref]->SetLineColor(iref+2);
         templates[iref]->SetLineWidth(2);
         templates[iref]->Draw("histsame");
      }
   }   
   if ( templates[0]->GetMaximum()>0 )templates[0]->SetFillColorAlpha(kRed,0.15);
   templates[0]->Draw("histsame");

   data->SetMarkerStyle(20);
   data->SetMarkerSize(1.4);
   data->SetLineColor(kBlack);
   data->Draw("e0same");

   TheoFit->SetLineColor(923);
   TheoFit->SetLineWidth(4);
   TheoFit->SetLineStyle(2);
   TheoFit->Draw("histsame");

   TLegend legend(0.18,0.70,0.94,0.92,"","NDC");
   legend.SetNColumns(3);
   legend.SetFillStyle(0);
   legend.SetBorderSize(0);
   legend.AddEntry(data,"Data","E0P");
   for ( int iref = 0 ; iref<reference_values.size() ; iref++ ) {
      legend.AddEntry(templates[iref],Form("Template #alpha=%6.2f",reference_values[iref]),"FL");
   }
   legend.AddEntry(TheoFit,"Estimated best model","L");
   legend.Draw();

   c1.Print(ps_name);
   c1.Print("plots/LTF_plot.pdf");


   // ---------------------------------------------- //
   // print linear-functions in every bin
   // ---------------------------------------------- //
   c1.Clear();
   c1.SetRightMargin(0.02);
   c1.SetTopMargin(0.02);
   c1.SetLeftMargin(0.16);
   c1.SetBottomMargin(0.16);

   gPad->SetTicky(1);

   //gStyle->SetLabelSize(0.05,"XYZ");
   //gStyle->SetTitleSize(0.05,"XYZ");
   gStyle->SetTitleOffset(1.1,"X");
   gStyle->SetTitleOffset(1.6,"Y");
   gStyle->SetMarkerSize(2);


   map < string, vector<double> > input_table = read_input_table2("data/CMS_data.txt",32);

   for ( int ibin = 0 ; ibin<fit.Dt.size() ; ibin++ ) {
      //TGraphErrors* data  = MakeTGraph(fit.ahat.row(0),fit.Dt.row(ibin));
      TGraphErrors* gdata = new TGraphErrors();
      gdata->SetPoint(0,fit.ahat(0),fit.Dt(ibin));
      gdata->SetPointError(0,0,data->GetBinError(ibin+1));
      gdata->SetMarkerStyle(20);

      TGraphErrors* graph = MakeTGraph(fit.M.col(1),ibin,fit.Y,fit.SysY);
      graph->Fit("pol1","QW"); // in this function it is "QW" (unweighted
      TF1* pol1 = (TF1*)graph->GetFunction("pol1")->Clone("pol1");
      pol1->SetLineColor(kRed);
      //pol1w->SetLineStyle(3);
      pol1->SetLineWidth(2);

      TF1* f2= new TF1("f2","[0]+[1]*pow(x,2)", -FLT_MIN,FLT_MAX );
      graph->Fit(f2,"QW");
      TF1* f1= new TF1("f1","[0]+[1]*pow(x,[2])", -FLT_MIN,FLT_MAX );
      

      graph->Fit("pol1","QW"); // "W": Ignore all point errors when fitting a TGraphErrors 
      bool UseLTFOutput = true;
      if ( UseLTFOutput ) {
         Eigen::VectorXd ltfpol1param =(fit.Y*fit.Mc().transpose()).row(ibin);
         //cout<<"M+*Y:"<<endl<< ltfpol1param <<endl;
         graph->GetFunction("pol1")->SetParameter(0,ltfpol1param(0));
         graph->GetFunction("pol1")->SetParameter(1,ltfpol1param(1));
         f1->SetParameter(0,ltfpol1param(0));
         f1->SetParameter(1,ltfpol1param(1));
         f1->SetParameter(2,fit.Gamma[0]);
      }
      if ( fit.Gamma[0]!=1 )   cout<<"Warning! Plotting with gamma factor !=1 not correctly implmeneted!"<<endl;

      graph->GetFunction("pol1")->SetLineWidth(3);
      
      graph->Fit("pol2","QW");
      graph->GetFunction("pol2")->SetLineWidth(83);

      graph->SetMarkerStyle(47);
      graph->SetMarkerColor(kRed+3);
      graph->SetLineColor(kRed+3);
      if ( fit.GetLogNormal() ) 
         graph->SetTitle((";"+referencename+";log("+yaxistitle+")").c_str()); // log(value/unit)
      else
         graph->SetTitle((";"+referencename+";"+yaxistitle).c_str());


      TF1* f1log = NULL;
      if ( fit.GetLogNormal() ) {
         //f1log = new TF1("pol1log","[0]+[1]*exp(x)", -FLT_MIN,FLT_MAX );
         f1log = new TF1("pol1log","log([0]+[1]*pow(x,1))", -FLT_MIN,FLT_MAX );

         TGraph* gexp = new TGraph();
         for ( int i = 0 ; i<graph->GetN() ; i++ ) 
            gexp->SetPoint(i,graph->GetX()[i], exp(graph->GetY()[i]));
         gexp->Fit("pol1","QW");

         Eigen::VectorXd ltfpol1param =(fit.Y*fit.Mc().transpose()).row(ibin);
         graph->GetFunction("pol1")->SetParameter(0,ltfpol1param(0));

         f1log->SetParameter(0, gexp->GetFunction("pol1")->GetParameter(0));
         f1log->SetParameter(1, gexp->GetFunction("pol1")->GetParameter(1));
         if ( fit.Gamma[0]!=1 )   cout<<"Warning! Plotting with gamma factor !=1 not correctly implmeneted!"<<endl;

         // //f1 = new TF1("pol1","[0]+[1]*x", -FLT_MIN,FLT_MAX );
         // graph->Fit(f1log,"QW"); // "W": Ignore all point errors when fitting a TGraphErrors 
         f1log->SetLineColor(kBlue+1);
         f1log->SetLineStyle(7);
         f1log->SetLineWidth(2);
      }


      if ( gdata->GetY()[0] > 0 )
         graph->SetMinimum(0);

      graph->SetMaximum( max(gdata->GetY()[0],max(graph->GetY()[0],graph->GetY()[graph->GetN()-1]))*1.2);
      graph->Fit("pol2","QW");
      graph->GetFunction("pol2")->SetLineColor(kBlue);;

      TF1* tang = new TF1("tang","[0]+[1]*x", -FLT_MIN,FLT_MAX );
      double f0  = graph->GetFunction("pol2")->Eval(gdata->GetX()[0]);
      double fp0 = graph->GetFunction("pol2")->GetParameter(1) + 2.*graph->GetFunction("pol2")->GetParameter(2)*gdata->GetX()[0];
      tang->SetParameter(0, f0 - fp0*gdata->GetX()[0] );
      tang->SetParameter(1, fp0);
      tang->SetLineWidth(2);
      tang->SetLineColor(kTeal+4);
      tang->SetLineStyle(7);
      
      pol1->SetLineStyle(3);
      graph->Draw("ap");
      if ( fit.GetLogNormal() )
         f1log->Draw("Lsame");
      else
         pol1->Draw("Lsame");
      //graph->GetFunction("pol1")->Draw("Lsame");
      if ( reference_values.size()+1<= 8 ) 
         graph->GetHistogram()->GetXaxis()->SetNdivisions(graph->GetN()+1);
      else
         graph->GetHistogram()->GetXaxis()->SetNdivisions(int(graph->GetN()/2)+1+200);
      //f2->Draw("same");

      tang->Draw("Lsame");
      gdata->Draw("P");
      
      if ( ibin==0 ) {
         //double xmin = fit.GetLogNormal() ? 0.36 : 0.45;
         //TLegend legend(xmin,0.19,0.96,0.47,"","NDC");
         TLegend legend(0.19,0.19,0.56,0.47,"","NDC"); 
         //legend.SetNColumns(3);
         legend.SetFillStyle(0);
         legend.SetBorderSize(0);
         //legend.SetTextSize(0.045);
         legend.AddEntry(data,"Data","E0P");
         legend.AddEntry(graph,"Templates","PE0");
         if ( fit.GetLogNormal() ) {
            legend.AddEntry(graph->GetFunction("pol1"),"Linear log(model)","L");
            legend.AddEntry(f1log,"#scale[0.9]{Linearized model #scale[0.7]{(unused)}}","L");
         }
         else {
            legend.AddEntry(graph->GetFunction("pol2"),"Second-order model","L");
            legend.AddEntry(tang,"Linear model","L");
            legend.AddEntry(pol1,"Linear Template Fit","L");
         }
         legend.DrawClone();
      }

      TLatex text;
      text.SetNDC();
      text.SetTextAlign(11);
      text.SetTextAlign(31);
      text.DrawLatex(0.955,0.20,Form("Bin %d",ibin));

      // text.DrawLatex(0.20,0.30,"CMS inclusive jets");
      // text.SetTextSize(0.04);
      // text.DrawLatex(0.20,0.25,Form("%3.1f_{ }<_{ }|y|_{ }<_{ }%3.1f",input_table["ylow"][ibin],input_table["yhigh"][ibin]));
      // text.DrawLatex(0.20,0.20,Form("%3.0f_{ }<_{ }p_{T}_{ }<_{ }%3.0f_{ }GeV",input_table["ptlow"][ibin],input_table["pthigh"][ibin]));
      
      // text.SetTextSize(0.04);
      // text.DrawLatex(0.20,0.93,Form("%3.1f_{ }<_{ }|y|_{ }<_{ }%3.1f",input_table["ylow"][ibin],input_table["yhigh"][ibin]));
      // text.DrawLatex(0.20,0.88,Form("%3.0f_{ }<_{ }p_{T}_{ }<_{ }%3.0f_{ }GeV",input_table["ptlow"][ibin],input_table["pthigh"][ibin]));

      cout<<input_table["ylow"][ibin]<<"\t"
          <<input_table["yhigh"][ibin]<<"\t"
          <<input_table["ptlow"][ibin]<<"\t"
          <<input_table["pthigh"][ibin]<<endl;

      c1.Print(ps_name);
      if ( fit.GetLogNormal() )
         c1.Print( Form("plots/LTFlog_bin_%02d.pdf",ibin));
      else
         c1.Print( Form("plots/LTF_bin_%02d.pdf",ibin));
      
   }
   

   // ---------------------------------------------- //
   //   chisq plot
   // ---------------------------------------------- //
    TGraph* gChi2 = new TGraph();
    int ndf = (fit.Dt.rows()-(fit.M.cols()-1));
    for ( int itmpl = 0 ; itmpl<fit.chisq_y.size() ; itmpl++ )  {
       gChi2->SetPoint(itmpl, reference_values[itmpl], fit.chisq_y(itmpl)/ndf);
    }
    gChi2->Fit("pol2","QW");
    gChi2->SetMarkerStyle(20);
    
    TGraph* gChi2LTF = new TGraph();
    gChi2LTF->SetPoint(0, fit.ahat(0), fit.chisq/ndf);
    gChi2LTF->SetMarkerStyle(29);
    gChi2LTF->SetMarkerSize(3.1);
    gChi2LTF->SetMarkerColor(kViolet+2);

    TGraph* gChi2chk = new TGraph();
    gChi2chk->SetPoint(0, fit.achk(0), fit.achk_chisq/ndf);
    gChi2chk->SetMarkerSize(2.2);
    gChi2chk->SetMarkerStyle(24);
    gChi2chk->SetMarkerColor(kRed);
    
    //gChi2LTF->Print("all");
    
    //gChi2->SetTitle(";#alpha_{0} [unit];#chi^{2}/ndf");
    gChi2->SetTitle((";"+referencename+";#chi^{2}/ndf").c_str());
    gChi2->SetMinimum(0.);
    gChi2->SetMaximum(20.2);
    gChi2->Draw("apc");
    if ( reference_values.size()+1<= 8 ) 
       gChi2->GetHistogram()->SetNdivisions(reference_values.size()+1+500,"X");
    else
       gChi2->GetHistogram()->SetNdivisions(int(reference_values.size()/2)+1+500,"X");
    
    TLine line;
    line.SetLineColor(920);
    line.SetLineStyle(3);
    line.DrawLine( 
       gChi2->GetHistogram()->GetXaxis()->GetXmin(),1.,
       gChi2->GetHistogram()->GetXaxis()->GetXmax(), 1);

    line.DrawLine( 
       fit.achk_chisq/ndf+1,1.,
       fit.achk_chisq/ndf+1,1.);

    gChi2chk->Draw("Psame");
    gChi2LTF->Draw("Psame");

    {
       TLegend legend(0.18,0.75,0.56,0.97,"","NDC");
       //legend.SetNColumns(3);
       legend.SetFillStyle(0);
       legend.SetBorderSize(0);
       //legend.SetTextSize(0.045);
       legend.AddEntry(gChi2LTF,"#hat#chi^{2} of the Quadratic Template Fit","P");
       legend.AddEntry(gChi2,   "#chi^{2}_{#font[12]{j}} of the individual templates","PL");
       legend.AddEntry(gChi2->GetFunction("pol2"),"Parabola","L");
       //legend.AddEntry(gChi2chk,"Minimum of #chi^{2} parabola #scale[0.8]{(#check#chi^{2})}","P"); //  (#check#chi^{2})
       legend.DrawClone();
    }
       
    c1.Print(ps_name);
    c1.Print( "plots/LTF_chi2.pdf");

    c1.Print( (string(ps_name)+"]").c_str() );
   
}

void LTF_ROOTTools::plotFitComparison(){
   TFile* file = new TFile("fit_quality.root", "READ");
   if (file->IsZombie()) {
        return;
    }

   TH1D* h_prob_linear = file->Get<TH1D>("h_chisq_prob_linear");
   TH1D* h_prob_quadratic = file->Get<TH1D>("h_chisq_prob_quadratic");
   TH1D* h_prob_ratio = file->Get<TH1D>("h_chisq_prob_ratio");
   TH1D* h_prob_ratio_rel = file->Get<TH1D>("h_chisq_prob_ratio_rel");
   TH1D* h_chisq_ratio = file->Get<TH1D>("h_chisq_ratio");
   TH1D* h_cheb2_sign = file->Get<TH1D>("h_cheb2_sign");
   TH1D* h_cheb3_sign = file->Get<TH1D>("h_cheb3_sign");
   TH1D* h_cheb_sign_ratio = file->Get<TH1D>("h_cheb_sign_ratio");
   TH2D* h_cheb2_sign_cheb3_sign = file->Get<TH2D>("h_cheb2_sign_cheb3_sign");
   TH2D* h_chisq_ratio_cheb2_sign = file->Get<TH2D>("h_chisq_ratio_cheb2_sign");
   TCanvas c1("c1","c1",800,800);
   c1.SetRightMargin(0.05);
   c1.SetLeftMargin(0.15);
   c1.SetTopMargin(0.08);

   c1.Print("plots/fit_quality.ps[");
	     
   {
     h_prob_linear->GetYaxis()->CenterTitle();
     h_prob_linear->GetYaxis()->SetTitle("Counts [a.u.]");
     h_prob_linear->GetXaxis()->SetTitle("p-value ( #chi^{2} probability)");
     h_prob_linear->SetLineColor(kBlack);
     h_prob_linear->SetLineWidth(2);
     h_prob_linear->Draw("hist");

     TLatex text;
     text.SetNDC();
     text.SetTextAlign(31);
     text.DrawLatex(0.80,0.80,"Linear fit");
     c1.Print("plots/fit_quality.ps");
   }
   {
     h_prob_quadratic->GetYaxis()->CenterTitle();
     h_prob_quadratic->GetYaxis()->SetTitle("Counts [a.u.]");
     h_prob_quadratic->GetXaxis()->SetTitle("p-value ( #chi^{2} probability)");
     h_prob_quadratic->SetLineColor(kBlack);
     h_prob_quadratic->SetLineWidth(2);
     h_prob_quadratic->Draw("hist");

     TLatex text;
     text.SetNDC();
     text.SetTextAlign(31);
     text.DrawLatex(0.80,0.80,"Quadratic fit");
     c1.Print("plots/fit_quality.ps");
   }
   {
     h_prob_ratio->GetYaxis()->CenterTitle();
     h_prob_ratio->GetYaxis()->SetTitle("Counts [a.u.]");
     h_prob_ratio->GetXaxis()->SetTitle("p-value (quad.) / p-value(lin) ");
     h_prob_ratio->SetLineColor(kBlack);
     h_prob_ratio->SetLineWidth(2);
     h_prob_ratio->Draw("hist");

     TLatex text;
     text.SetNDC();
     text.SetTextAlign(31);
     text.DrawLatex(0.80,0.80,"Ratio of probabilities");
     c1.Print("plots/fit_quality.ps");
   }
   {
     h_prob_ratio_rel->GetYaxis()->CenterTitle();
     h_prob_ratio_rel->GetYaxis()->SetTitle("Counts [a.u.]");
     h_prob_ratio_rel->GetXaxis()->SetTitle("(p-value (quad.) - p-value(lin)) / p-value (quad.) ");
     h_prob_ratio_rel->SetLineColor(kBlack);
     h_prob_ratio_rel->SetLineWidth(2);
     h_prob_ratio_rel->Draw("hist");

     TLatex text;
     text.SetNDC();
     text.SetTextAlign(31);
     text.DrawLatex(0.80,0.80,"Relative ratio");
     text.DrawLatex(0.80,0.76,"of probabilities");
     c1.Print("plots/fit_quality.ps");
   }
   {
     h_chisq_ratio->GetYaxis()->CenterTitle();
     h_chisq_ratio->GetYaxis()->SetTitle("Counts [a.u.]");
     h_chisq_ratio->GetXaxis()->SetTitle("(#chi^{2}_{quad} / ndf_{quad}) / (#chi^{2}_{linear} / ndf_{linear})");
     h_chisq_ratio->SetLineColor(kBlack);
     h_chisq_ratio->SetLineWidth(2);
     h_chisq_ratio->Draw("hist");

     TLatex text;
     text.SetNDC();
     text.SetTextAlign(31);
     text.DrawLatex(0.80,0.80,"Ratio of #chi^{2} values");
     c1.Print("plots/fit_quality.ps");
   }
   {
     h_cheb2_sign->GetYaxis()->CenterTitle();
     h_cheb2_sign->GetYaxis()->SetTitle("Counts [a.u.]");
     h_cheb2_sign->GetXaxis()->SetTitle("Chebyshev: abs(p2/p2err)");
     h_cheb2_sign->SetLineColor(kBlack);
     h_cheb2_sign->SetLineWidth(2);
     h_cheb2_sign->Draw("hist");
     c1.Print("plots/fit_quality.ps");
   }
   {
     h_cheb3_sign->GetYaxis()->CenterTitle();
     h_cheb3_sign->GetYaxis()->SetTitle("Counts [a.u.]");
     h_cheb3_sign->GetXaxis()->SetTitle("Chebyshev: abs(p3/p3err)");
     h_cheb3_sign->SetLineColor(kBlack);
     h_cheb3_sign->SetLineWidth(2);
     h_cheb3_sign->Draw("hist");
     c1.Print("plots/fit_quality.ps");
   }
   {
     h_cheb_sign_ratio->GetYaxis()->CenterTitle();
     h_cheb_sign_ratio->GetYaxis()->SetTitle("Counts [a.u.]");
     h_cheb_sign_ratio->GetXaxis()->SetTitle("abs(p3/p3err) / abs(p2/p2err)");
     h_cheb_sign_ratio->SetLineColor(kBlack);
     h_cheb_sign_ratio->SetLineWidth(2);
     h_cheb_sign_ratio->Draw("hist");
     c1.Print("plots/fit_quality.ps");
   }
   {
     c1.SetRightMargin(0.1);
     h_cheb2_sign_cheb3_sign->GetXaxis()->SetTitle("Chebyshev significance p2/p2err");
     h_cheb2_sign_cheb3_sign->GetYaxis()->SetTitle("Chebyshev significance p3/p3err");
     h_cheb2_sign_cheb3_sign->Draw("colz");
     c1.Print("plots/fit_quality.ps");
   }
   {
     c1.SetRightMargin(0.1);
     h_chisq_ratio_cheb2_sign->GetYaxis()->SetTitle("Chebyshev significance p2/p2err");
     h_chisq_ratio_cheb2_sign->GetXaxis()->SetTitle("(#chi^{2}_{quad} / ndf_{quad}) / (#chi^{2}_{linear} / ndf_{linear})");
     h_chisq_ratio_cheb2_sign->Draw("colz");
     c1.Print("plots/fit_quality.ps");
   }
   c1.Print("plots/fit_quality.ps]");
   
   file->Close();
}
