#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <cstring>

#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include "TCanvas.h"
#include <TMarker.h>
#include <TString.h>
#include <TVirtualFitter.h>

#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"

#include "boundaries.h"


using namespace std;

const int digi=3;


// change centrality and pT bins as and when necessary

void  LoadStyle();
void drawText2(const char */*text*/, float /*xp*/, float /*yp*/, int /*size*/);
void makeMultiPanelCanvas(TCanvas*& /*canv*/,
                          const Int_t /*columns*/,
                          const Int_t /*rows*/,
                          const Float_t /*leftOffset*/,
                          const Float_t /*bottomOffset*/,
                          const Float_t /*leftMargin*/,
                          const Float_t /*bottomMargin*/,
                          const Float_t /*edge*/,const Float_t /*asyoffset*/); 
void MakeHist(TH1F *&/*hist*/,int /*istat*/, const char */*xname*/, const char */*yname*/);
void MakeHistRMS(TH1F *&/*hRMS*/);
void MakeHistMean(TH1F *&/*Mean*/);
void MakeZero(TH1F *&/*hist*/);
double roundoff(double /*val*/);
TLegend *getLegend(double /*x1*/, double /*y1*/, double /*x2*/, double /*y2*/);

void FitGauss(TH1F* inHist_p, double &mean, double &meanError, double &res, double &resError);


void fit_double_gaussian(TH1F *&hrsp);
int fit_dscb(TH1F *&hrsp,
             const double nsigma,
             const double jtptmin,
             const int niter,
             const string alg);

/// double sided crystal ball function definition  
double fnc_dscb(double*xx,double*pp);
template <class T>
bool from_string(T& t,
                 const std::string& s,
                 std::ios_base& (*f)(std::ios_base&))
{
  std::istringstream iss(s);
  return !(iss >> f >> t).fail();
}


void fit_gaussian(TH1F*& hrsp,
                  const double nsigma,
                  const double jtptmin,
                  const int niter);
void adjust_fitrange(TH1 *h,double& min,double& max);

void set_range_truncatedRMS(TH1F *&hist,float frac);
double fracRMS = 1.00;

void FitDist(TH1F *&/*hrsp*/, double &/*mean*/, double &/*emean*/, double &/*sig*/, double &/*esig*/);

int plot_JESClosure(int wJetID=1,
		    std::string jetType = "PF")
{
  int iSave=1;

  int id = wJetID;

  cout <<" npt : " << nbins_pt << endl;

  LoadStyle();

  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  std::string algo="ak";
  
  TVirtualFitter::SetDefaultFitter("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-04); 
  
  TFile *finPP = TFile::Open("ppJECCheck_merged.root");

  TH1F *hrsp_pp[knj][nbins_eta][nbins_pt];
  TH1F *hMean_pp[knj][nbins_eta], *hSigma_pp[knj][nbins_eta];
  TH1F *hrsp_hin[knj][nbins_eta][nbins_pt];
  TH1F *hMean_hin[knj][nbins_eta], *hSigma_hin[knj][nbins_eta];

  TF1 *fgaus=0;

  for(int nj=0; nj<knj; nj++){
    for(int ny = 0; ny < nbins_eta; ++ny){

      hMean_pp[nj][ny] = new TH1F(Form("hMean_pp_ptbin%d_etabin%d",nj, ny),"",nbins_pt,ptbins);
      hMean_pp[nj][ny]->SetMarkerColor(1);
      hMean_pp[nj][ny]->SetMarkerStyle(20);
      hMean_pp[nj][ny]->SetLineColor(1);
      hMean_pp[nj][ny]->SetMarkerSize(1.3);
      MakeHistMean(hMean_pp[nj][ny]); 
      //MakeHistMean(hMean[nj],0.052,0.934); 
    
      hSigma_pp[nj][ny] = new TH1F(Form("hSigma_pp_ptbin%d_etabin%d",nj,ny),"",nbins_pt,ptbins);
      hSigma_pp[nj][ny]->SetMarkerColor(1);
      hSigma_pp[nj][ny]->SetMarkerStyle(20);
      hSigma_pp[nj][ny]->SetLineColor(1);
      hSigma_pp[nj][ny]->SetMarkerSize(1.3);
      MakeHistRMS(hSigma_pp[nj][ny]); 

      hMean_hin[nj][ny] = new TH1F(Form("hMean_hin_ptbin%d_etabin%d",nj, ny),"",nbins_pt,ptbins);
      hMean_hin[nj][ny]->SetMarkerColor(2);
      hMean_hin[nj][ny]->SetMarkerStyle(25);
      hMean_hin[nj][ny]->SetLineColor(1);
      hMean_hin[nj][ny]->SetMarkerSize(1.3);
      MakeHistMean(hMean_hin[nj][ny]); 
      //MakeHistMean(hMean[nj],0.052,0.934); 
    
      hSigma_hin[nj][ny] = new TH1F(Form("hSigma_hin_ptbin%d_etabin%d",nj,ny),"",nbins_pt,ptbins);
      hSigma_hin[nj][ny]->SetMarkerColor(2);
      hSigma_hin[nj][ny]->SetMarkerStyle(25);
      hSigma_hin[nj][ny]->SetLineColor(1);
      hSigma_hin[nj][ny]->SetMarkerSize(1.3);
      MakeHistRMS(hSigma_hin[nj][ny]); 
      
      for(int ip=0; ip<nbins_pt; ip++){
	
	hrsp_pp[nj][ny][ip]      = (TH1F*)finPP->Get(Form("hJEC_ppJEC_ptbin%d_etabin%d",ip, ny));
	double norm=-9;
	double err=-9; 
	double mean=-9; 
	double emean=-9;
	double sig=-9;
	double esig=-9;
	
	int fitstatus = 0;

	//hrsp_pp[nj][ny][ip]->Print("base");
	if(hrsp_pp[nj][ny][ip]->GetEntries()>20){
	  hrsp_pp[nj][ny][ip]->Scale(1./hrsp_pp[nj][ny][ip]->Integral());
	  hrsp_pp[nj][ny][ip]->Rebin(2);
	  norm  = hrsp_pp[nj][ny][ip]->GetMaximumStored();
	  err   = hrsp_pp[nj][ny][ip]->GetStdDevError();
	  mean  = hrsp_pp[nj][ny][ip]->GetMean();
	  emean = hrsp_pp[nj][ny][ip]->GetMeanError();
	  sig   = hrsp_pp[nj][ny][ip]->GetStdDev()/mean;
	  esig  = (pow(1/mean,2)*pow(err,2))+(pow(-sig/pow(mean,2),2)*pow(emean,2));

	  fgaus = new TF1("fgaus","gaus", 0.5,1.5);
	  fgaus->SetParameters(1, 1.009);
	  fgaus->SetParameters(2, 0.5);
	
	  fitstatus = hrsp_pp[nj][ny][ip]->Fit(fgaus,"RQ");
	  cout << " fitstatus : "<< fitstatus << " Fit Error : " << fgaus->GetParError(1) << " Hist Error " << hrsp_pp[nj][ny][ip]->GetMeanError() << endl;

	  mean  = (fitstatus!=0) ? hrsp_pp[nj][ny][ip]->GetMean()     : fgaus->GetParameter(1);
	  emean = (fitstatus!=0) ? hrsp_pp[nj][ny][ip]->GetMeanError(): fgaus->GetParError(1);//hrsp[nj][ip]->GetMeanError();
	  sig   = (fitstatus!=0) ? hrsp_pp[nj][ny][ip]->GetStdDev()/mean : fgaus->GetParameter(2)/fgaus->GetParameter(1);
	  esig  = (fitstatus!=0) ? sqrt((pow(1/mean,2)*pow(hrsp_pp[nj][ny][ip]->GetStdDevError(),2))+(pow(-hrsp_pp[nj][ny][ip]->GetStdDev()/pow(mean,2),2)*pow(emean,2))) : sqrt((pow(1/mean,2)*pow(hrsp_pp[nj][ny][ip]->GetStdDevError(),2))+(pow(-hrsp_pp[nj][ny][ip]->GetStdDev()/pow(mean,2),2)*pow(emean,2))); 

	  hMean_pp[nj][ny]->SetBinContent (ip+1, mean);
	  hMean_pp[nj][ny]->SetBinError   (ip+1, emean);
	  hSigma_pp[nj][ny]->SetBinContent (ip+1, sig);
	  hSigma_pp[nj][ny]->SetBinError   (ip+1, esig);
	  cout<<"ppJEC Mean = "<<mean<<"  sig = "<<sig<<endl;
	  
	}
	//   hMean_pp[nj][ny]->SetBinContent (ip+1, -9);
	//   hMean_pp[nj][ny]->SetBinError   (ip+1, 0);
	//   hSigma_pp[nj][ny]->SetBinContent (ip+1, -9);
	//   hSigma_pp[nj][ny]->SetBinError   (ip+1, 0);
	// }
	//! HIN JEC 
	hrsp_hin[nj][ny][ip]      = (TH1F*)finPP->Get(Form("hJEC_hinJEC_ptbin%d_etabin%d",ip, ny));

	//hrsp_hin[nj][ny][ip]->Print("base");
	if(hrsp_hin[nj][ny][ip]->GetEntries()>20){
	  hrsp_hin[nj][ny][ip]->Scale(1./hrsp_hin[nj][ny][ip]->Integral());
	  hrsp_hin[nj][ny][ip]->Rebin(2);
	  norm  = hrsp_hin[nj][ny][ip]->GetMaximumStored();
	  err   = hrsp_hin[nj][ny][ip]->GetStdDevError();
	  mean  = hrsp_hin[nj][ny][ip]->GetMean();
	  emean = hrsp_hin[nj][ny][ip]->GetMeanError();
	  sig   = hrsp_hin[nj][ny][ip]->GetStdDev()/mean;
	  esig = (pow(1/mean,2)*pow(err,2))+(pow(-sig/pow(mean,2),2)*pow(emean,2));

	  fgaus = new TF1("fgaus","gaus", 0.5,1.5);
	  fgaus->SetParameters(1, 1.009);
	  fgaus->SetParameters(2, 0.5);
	
	  fitstatus = 0;
	  fitstatus = hrsp_hin[nj][ny][ip]->Fit(fgaus,"RQ");
	  cout << " fitstatus : "<< fitstatus << " Fit Error : " << fgaus->GetParError(1) << " Hist Error " << hrsp_hin[nj][ny][ip]->GetMeanError() << endl;

	  mean  = (fitstatus!=0) ? hrsp_hin[nj][ny][ip]->GetMean()     : fgaus->GetParameter(1);
	  emean = (fitstatus!=0) ? hrsp_hin[nj][ny][ip]->GetMeanError(): fgaus->GetParError(1);//hrsp[nj][ip]->GetMeanError();
	  sig   = (fitstatus!=0) ? hrsp_hin[nj][ny][ip]->GetStdDev()/mean : fgaus->GetParameter(2)/fgaus->GetParameter(1);
	  esig  = (fitstatus!=0) ? sqrt((pow(1/mean,2)*pow(hrsp_hin[nj][ny][ip]->GetStdDevError(),2))+(pow(-hrsp_hin[nj][ny][ip]->GetStdDev()/pow(mean,2),2)*pow(emean,2))) : sqrt((pow(1/mean,2)*pow(hrsp_hin[nj][ny][ip]->GetStdDevError(),2))+(pow(-hrsp_hin[nj][ny][ip]->GetStdDev()/pow(mean,2),2)*pow(emean,2))); 

	  hMean_hin[nj][ny]->SetBinContent (ip+1, mean);
	  hMean_hin[nj][ny]->SetBinError   (ip+1, emean);
	  hSigma_hin[nj][ny]->SetBinContent (ip+1, sig);
	  hSigma_hin[nj][ny]->SetBinError   (ip+1, esig);
	  cout<<"hinJEC Mean = "<<mean<<"  sig = "<<sig<<endl;
	}
	// else{
	//   hMean_hin[nj][ny]->SetBinContent (ip+1, -9);
	//   hMean_hin[nj][ny]->SetBinError   (ip+1, 0);
	//   hSigma_hin[nj][ny]->SetBinContent (ip+1, -9);
	//   hSigma_hin[nj][ny]->SetBinError   (ip+1, 0);
	// }
	
	cout<<"Finished running for ptbin "<<ip<<"  etabin "<<ny<<"   radiusbin "<< nj<<endl;
      }
    }
  }

  cout<<"Going to make plots "<<endl;
  //gPad->Close();

  int maxc=6;
  int maxr=2;
  int ipad=0;

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TLegend *leg0 = getLegend(0.1950635,0.4625641,0.8242194,0.7497436);
  leg0->SetTextSize(0.10);
  cout<<"test break1"<<endl;
  leg0->AddEntry(hMean_hin[0][0],Form("HIN JEC ak%s",(srad[0]+jetType).c_str()),"p");
  leg0->AddEntry(hMean_pp[0][0],Form("pp  JEC ak%s",(srad[0]+jetType).c_str()),"p");
  cout<<"test break2"<<endl;
  TLine *l0 = new TLine(xmin,1.00,xmax,1.00);
  l0->SetLineStyle(2);
  l0->SetLineWidth(2);
  TLine *l1 = new TLine(xmin,1.01,xmax,1.01);
  l1->SetLineStyle(2);
  TLine *l2 = new TLine(xmin,0.99,xmax,0.99);
  l2->SetLineStyle(2);

  cout<<"plotting the mean/res"<<endl;
  
  TCanvas *c11[knj];
  for(int nj=0; nj<knj; nj++){
    c11[nj] = new TCanvas(Form("c11_%d",nj),Form("%s",(srad[nj]+jetType).c_str()),1695,530);
    makeMultiPanelCanvas(c11[nj],maxc,maxr,0.0,0.0,0.22,0.22,0.02,0);
    cout<<"r = "<<nj+4<<endl;
  
    for(int ny = 0; ny < nbins_eta-2; ++ny){
      ipad++;
      cout<<"ipad = "<<ipad<<" etabin = "<<ny<<endl;
      c11[nj]->cd(ipad);
      gPad->SetLogx();
      gStyle->SetErrorX(0);

      hSigma_pp[nj][ny]->GetXaxis()->SetRangeUser(xmin,xmax);
      hSigma_pp[nj][ny]->GetYaxis()->SetRangeUser(0.001, 0.4);
      hSigma_pp[nj][ny]->Draw("p");

      hSigma_hin[nj][ny]->Draw("p same");
      
      if(ipad==2){
	drawText2("PYTHIA8 CUETP8M1",0.19,0.85,21);
	leg0->Draw();
      }
      if(ipad==3){
	drawText2("pp #sqrt{s} = 5.02 TeV",0.10,0.85,19);
      }
      drawText2(Form("%s", legetabins[ny].c_str()),0.28,0.70,21);

      c11[nj]->cd(ipad+maxc);
      gPad->SetLogx();
      cout<<"ipad+maxc = "<<ipad+maxc<<" etabin = "<<ny<<endl;
    
      // if(ipad!=1)
      hMean_pp[nj][ny]->GetXaxis()->SetRangeUser(xmin,xmax);
      hMean_pp[nj][ny]->GetYaxis()->SetRangeUser(0.9,1.1);
      hMean_pp[nj][ny]->Draw("p");
      
      hMean_hin[nj][ny]->Draw("p same");

      l1->Draw();
      l0->Draw();
      l2->Draw();
    }
    if(iSave){
      c11[nj]->SaveAs(Form("Plots/JESJER_%d_pp_%s.pdf",wJetID,(algo+srad[nj]+jetType).c_str()));
    }
  }
  //return 0;

  ipad=0;
  TCanvas *c99[knj][nbins_eta-1];
  for(int nj=0;nj<knj;nj++){
    for(int ny = 0; ny < nbins_eta-1; ++ny){
      c99[nj][ny] = new TCanvas(Form("c99%d_%d",nj,ny),Form("%s %s ",srad[nj].c_str(),legetabins[ny].c_str()),100,102,1399,942);
      c99[nj][ny]->Divide(5,5,0,0);
      ipad=0;
      for(int ip=0;ip<nbins_pt;ip++){      
	c99[nj][ny]->cd(++ipad);
	
	gPad->SetLeftMargin(0.15);
	gPad->SetBottomMargin(0.15);
	gPad->SetRightMargin(0.1);
	
	double ymax = TMath::Max(hrsp_pp[nj][ny][ip]->GetMaximum(),hrsp_pp[nj][nbins_eta-1][ip]->GetMaximum());
	//double ymax = hrsp[nj][ic][ip]->GetMaximum();
	hrsp_pp[nj][ny][ip]->SetMaximum(1.45*ymax);
	hrsp_pp[nj][ny][ip]->SetMinimum(0);
	hrsp_pp[nj][ny][ip]->SetTitle(0);
	hrsp_pp[nj][ny][ip]->SetAxisRange(0, 2, "X");
	hrsp_pp[nj][ny][ip]->GetXaxis()->SetTitle("<(reco jet p_{T} / gen jet p_{T}>");
	hrsp_pp[nj][ny][ip]->GetXaxis()->SetTitleFont(42);
	hrsp_pp[nj][ny][ip]->GetXaxis()->SetLabelFont(42);
	hrsp_pp[nj][ny][ip]->GetXaxis()->SetLabelSize(0.08);
	hrsp_pp[nj][ny][ip]->GetXaxis()->SetTitleSize(0.07);
	hrsp_pp[nj][ny][ip]->GetYaxis()->SetTitle("");
	hrsp_pp[nj][ny][ip]->GetYaxis()->SetTitleFont(42);
	hrsp_pp[nj][ny][ip]->GetYaxis()->SetLabelFont(42);
	hrsp_pp[nj][ny][ip]->GetYaxis()->SetLabelSize(0.08);
      
	hrsp_pp[nj][ny][ip]->SetMarkerStyle(20);
	hrsp_pp[nj][ny][ip]->SetMarkerColor(1);
	hrsp_pp[nj][ny][ip]->SetLineColor(1);
	hrsp_pp[nj][ny][ip]->SetMarkerSize(1.3);
	hrsp_pp[nj][ny][ip]->Draw("p");  

	hrsp_hin[nj][ny][ip]->SetMarkerStyle(24);
	hrsp_hin[nj][ny][ip]->SetMarkerColor(2);
	hrsp_hin[nj][ny][ip]->SetLineColor(1);
	hrsp_hin[nj][ny][ip]->SetMarkerSize(1.3);
	hrsp_hin[nj][ny][ip]->Draw("p same");  
	
	std::ostringstream strs; 
	strs << ptbins[ip] << "< p_{T} (GeV/c) <" << ptbins[ip+1];
	std::string spt = strs.str();

	if(ipad==1){
	  //drawText2(Form("%s",(srad[nj]+jetType).c_str()),0.28,0.90,21);      

	  drawText2(legetabins[ny].c_str(),0.75,0.87,21);	  
	  drawText2(spt.c_str(),0.22,0.80,21);		

	  if(wJetID)drawText2("w JetID",0.22,0.65,21);
	  
	} else drawText2(spt.c_str(),0.17,0.80,21);

      }
      if(iSave){
	c99[nj][ny]->SaveAs(Form("Plots/IndDist_pp_etabin%d_%s.pdf",ny,(algo+srad[nj]+jetType).c_str()));
      }
    }
  }
  cout<<"finished making all plots"<<endl;
  TFile *fout = new TFile(Form("JERJES_R%s_pp.root",srad[0].c_str()),"RECREATE");
  fout->cd();
  fout->mkdir(Form("ak%sJetAnalyzer",srad[0].c_str()));
  fout->cd(Form("ak%sJetAnalyzer",srad[0].c_str()));
  for(int nj=0; nj<knj; nj++){    
    for(int ny = 0; ny < nbins_eta-1; ++ny){
      hMean_pp[nj][ny]->SetName(Form("hMean_ppJEC_R%s_PF_etabin%d",srad[nj].c_str(),ny));
      hMean_pp[nj][ny]->SetTitle(Form("hMean_ppJEC_R%s_PF_etabin%d",srad[nj].c_str(),ny));
      hSigma_pp[nj][ny]->SetName(Form("hSigma_ppJEC_R%s_PF_etabin%d",srad[nj].c_str(),ny));
      hSigma_pp[nj][ny]->SetTitle(Form("hSigma_ppJEC_R%s_PF_etabin%d",srad[nj].c_str(),ny));
      hMean_pp[nj][ny]->Write();
      hSigma_pp[nj][ny]->Write();
      hMean_hin[nj][ny]->SetName(Form("hMean_hinJEC_R%s_PF_etabin%d",srad[nj].c_str(),ny));
      hMean_hin[nj][ny]->SetTitle(Form("hMean_hinJEC_R%s_PF_etabin%d",srad[nj].c_str(),ny));
      hSigma_hin[nj][ny]->SetName(Form("hSigma_hinJEC_R%s_PF_etabin%d",srad[nj].c_str(),ny));
      hSigma_hin[nj][ny]->SetTitle(Form("hSigma_hinJEC_R%s_PF_etabin%d",srad[nj].c_str(),ny));
      hMean_hin[nj][ny]->Write();
      hSigma_hin[nj][ny]->Write();
    }
  }
  fout->Close();
  return 0;
}
void LoadStyle()
{
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetErrorX(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetFillColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetStatColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetPadBorderSize(0);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetPalette(1);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
}
void drawText2(const char *text, float xp, float yp, int size){
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(43);
  tex->SetTextSize(size);
  tex->SetTextColor(kBlack);
  tex->SetLineWidth(1);
  tex->SetNDC();
  tex->Draw();
}
void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge,const Float_t asyoffset) {
  if (canv==0) {
    cout<<"makeMultiPanelCanvas :  Got null canvas."<<endl;
    return;
  }
  canv->Clear();
  
  TPad* pad[columns][rows];
  
  Float_t Xlow[columns];
  Float_t Xup [columns];
  Float_t Ylow[rows];
  Float_t Yup [rows];
  
  Float_t PadWidth =
    (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
		      (1.0/(1.0-edge))+(Float_t)columns-2.0);
  Float_t PadHeight =
    (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
			(1.0/(1.0-edge))+(Float_t)rows-2.0);
  
  Xlow[0] = leftOffset - asyoffset;
  Xup[0]  = leftOffset + PadWidth/(1.0-leftMargin);
  Xup[columns-1] = 1;
  Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);
  
  Yup[0] = 1;
  Ylow[0] = 1.0-PadHeight/(1.0-edge);
  Ylow[rows-1] = bottomOffset;
  Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);
  
  for(Int_t i=1;i<columns-1;i++) {
    Xlow[i] = Xup[0] + (i-1)*PadWidth;
    Xup[i] = Xup[0] + (i)*PadWidth;
  }
  Int_t ct = 0;
  for(Int_t i=rows-2;i>0;i--) {
    Ylow[i] = Yup[rows-1] + ct*PadHeight;
    Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
    ct++;
  }
  
  TString padName;
  for(Int_t i=0;i<columns;i++) {
    for(Int_t j=0;j<rows;j++) {
      canv->cd();
      padName = Form("p_%d_%d",i,j);
      pad[i][j] = new TPad(padName.Data(),padName.Data(),
			   Xlow[i],Ylow[j],Xup[i],Yup[j]);
      
      // this is hacked version to create aysmmetric pads around low 
      if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
      else pad[i][j]->SetLeftMargin(0);
      
      if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
      else pad[i][j]->SetRightMargin(0);
      
      if(j==0) pad[i][j]->SetTopMargin(edge);
      else pad[i][j]->SetTopMargin(0);
      
      if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
      else pad[i][j]->SetBottomMargin(0);
      
      pad[i][j]->Draw();
      pad[i][j]->cd();
      pad[i][j]->SetNumber(columns*j+i+1);
    }
  }
}
void MakeHist(TH1F *&histo,int istat,const char *xname, const char *yname)
{
  histo->SetStats(istat);
  histo->SetMarkerStyle(24);
  histo->SetMarkerColor(1);
  histo->SetLineColor(1);
  histo->SetLineStyle(1);
  histo->GetXaxis()->SetTitle(xname);
  histo->GetXaxis()->CenterTitle(true);
  histo->GetYaxis()->SetTitle(yname);
  histo->GetYaxis()->CenterTitle(true);
}
void MakeHistMean(TH1F *&h1)
{
  // h1->SetMaximum(ymax);
  // h1->SetMinimum(ymin);
  h1->SetTitle("");
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
  h1->GetXaxis()->CenterTitle(true);
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetNoExponent();
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetTitleSize(0.07);
  h1->GetXaxis()->SetTitleOffset(1.15);
  h1->GetXaxis()->SetLabelSize(0.07);
  h1->GetXaxis()->SetLabelOffset(0.005);
  h1->GetXaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetTitle("#mu");
  h1->GetYaxis()->CenterTitle(true);
  h1->GetYaxis()->SetTitleSize(0.07);
  h1->GetYaxis()->SetTitleOffset(1.50);
  h1->GetYaxis()->SetLabelSize(0.07);
  h1->GetYaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetDecimals(true);
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetLabelFont(42);
}
void MakeHistRMS(TH1F *&h1)
{

  h1->SetTitle("");
  // h1->SetMaximum(ymax);
  // h1->SetMinimum(ymin);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
  h1->GetXaxis()->CenterTitle(true);
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetNoExponent();
  h1->GetXaxis()->SetNdivisions(507);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetLabelOffset(0.01);
  h1->GetXaxis()->SetLabelSize(0.07);
  h1->GetXaxis()->SetTitleSize(0.07);
  h1->GetXaxis()->SetTitleOffset(1.15);
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetTitle("#sigma / #mu");
  h1->GetYaxis()->CenterTitle(true);
  h1->GetYaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetLabelFont(42);
  h1->GetYaxis()->SetLabelOffset(0.01);
  h1->GetYaxis()->SetLabelSize(0.09);
  h1->GetYaxis()->SetTitleSize(0.09);
  h1->GetYaxis()->SetTitleOffset(1.12);
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetDecimals(true);

}

void fit_gaussian(TH1F *&hrsp,
                  const double nsigma,
                  const double jtptmin,
                  const int niter)
{
  if (0==hrsp) {
    cout<<"ERROR: Empty pointer to fit_gaussian()"<<endl;return;
  }
  
  string histname = hrsp->GetName();
  double mean     = hrsp->GetMean();
  double rms      = hrsp->GetRMS();
  double ptRefMax(1.0),rspMax(0.0);

  //cout << " Mean : "  << mean << " \t  RMS  : " << rms << endl;

  double norm  = hrsp->GetMaximumStored();
  double peak  = mean;
  double sigma = rms;

  rspMax =jtptmin/ptRefMax;

  double xmin  = hrsp->GetXaxis()->GetXmin();
  double xmax  = hrsp->GetXaxis()->GetXmax();
  TF1* fitfnc(0); int fitstatus(-1);
  for (int iiter=0;iiter<niter;iiter++) {
    vector<double> vv;
    vv.push_back(rspMax);
    vv.push_back(xmin);
    vv.push_back(peak-nsigma*sigma);
    double fitrange_min = *std::max_element(vv.begin(),vv.end());
    double fitrange_max = std::min(xmax,peak+nsigma*sigma);

    adjust_fitrange(hrsp,fitrange_min,fitrange_max);
    fitfnc = new TF1("fgaus","gaus",fitrange_min,fitrange_max);
    fitfnc->SetParNames("N","#mu","#sigma");
    fitfnc->SetParameter(0,norm);
    fitfnc->SetParameter(1,peak);
    fitfnc->SetParameter(2,sigma);
    fitstatus = hrsp->Fit(fitfnc,"RQ0");
    delete fitfnc;
    fitfnc = hrsp->GetFunction("fgaus");
    //fitfnc->ResetBit(TF1::kNotDraw);
    norm  = fitfnc->GetParameter(0);
    peak  = fitfnc->GetParameter(1);
    sigma = fitfnc->GetParameter(2);
  }
  if(hrsp->GetFunction("fgaus")==0){
      cout << "No function recorded in histogram " << hrsp->GetName() << endl;
    }
  //   if (0==fitstatus){
  //     cout<<"fit_gaussian() to "<<hrsp->GetName()    <<"  sucessful : " <<endl;
  //   }          



  if (0!=fitstatus){
    cout<<"fit_gaussian() to "<<hrsp->GetName()
        <<" failed. Fitstatus: "<<fitstatus
        <<" - FNC deleted."<<endl;
    hrsp->GetListOfFunctions()->Delete();
  }
}
void adjust_fitrange(TH1 *h,double& min,double& max)
{
  int imin=1; while (h->GetBinLowEdge(imin)<min) imin++;
  int imax=1; while (h->GetBinLowEdge(imax)<max) imax++;
  while ((imax-imin)<8) {
    if (imin>1) {imin--; min = h->GetBinCenter(imin); }
    if (imax<h->GetNbinsX()-1) { imax++; max=h->GetBinCenter(imax); }
  }
}
void MakeZero(TH1F *&h1)
{
  h1->GetYaxis()->SetLabelSize(0);
  h1->GetYaxis()->SetTitleSize(0);
}
int fit_dscb(TH1F *&hrsp,
             const double nsigma,
             const double jtptmin,
             const int niter,
             const string alg)
{
  if (0==hrsp) {
    cout<<"ERROR: Empty pointer to fit_dscb()"<<endl;return -1;
  }


  // first use a gaussian to constrain crystal ball gaussian core
  fit_gaussian(hrsp, nsigma, jtptmin, niter);
  TF1* fgaus = hrsp->GetFunction("fgaus");
  if (0==fgaus) {
    hrsp->GetListOfFunctions()->Delete();
    return -1;
  }

  // implementation of the low pt bias threshold 
  string histname = hrsp->GetName();
  //double ptRefMax(1.0),rspMax(0.0);

  double fitrange_min(0.4);
  double fitrange_max(1.6);

  //cout <<" \t \t  xmin : "  << fitrange_min << "\t" << fitrange_max << endl;


  TF1* fdscb = new TF1("fdscb",fnc_dscb,fitrange_min,fitrange_max,7);
  fdscb->SetLineWidth(2);
  fdscb->SetLineStyle(2);

  double norm = fgaus->GetParameter(0);
  double mean = fgaus->GetParameter(1);
  double sigma= fgaus->GetParameter(2);

  double aone(2.0),atwo(2.0),pone(10.0),ptwo(10.0);
  TVirtualFitter::SetDefaultFitter("Minuit");
  
  int fitstatus(0);
  for (int i=0;i<niter;i++) {
    fdscb->SetParameter(0,norm); // N
    fdscb->SetParameter(1,mean); // mean
    fdscb->SetParameter(2,sigma);// sigma
    fdscb->SetParameter(3,aone); // a1
    fdscb->SetParameter(4,pone); // p1
    fdscb->SetParameter(5,atwo); // a2
    fdscb->SetParameter(6,ptwo); // p2                

    fdscb->FixParameter(1,mean);
    fdscb->FixParameter(2,sigma);
    
    if (i>0) fdscb->FixParameter(3,aone);
    else fdscb->SetParLimits(3,1.,50.);
    
    if (i>1) fdscb->FixParameter(5,atwo);
    else fdscb->SetParLimits(5,1.,50.);

    fdscb->SetParLimits(4,0.,100.);
    fdscb->SetParLimits(6,0.,100.);

    fitstatus = hrsp->Fit(fdscb,"RQB+");
    if (0==fitstatus) i=999;
    delete fdscb;
    fdscb = hrsp->GetFunction("fdscb");
    if (0==fdscb) return -1;

    norm  = fdscb->GetParameter(0);
    aone  = fdscb->GetParameter(3);
    pone  = fdscb->GetParameter(4);
    atwo  = fdscb->GetParameter(5);
    ptwo  = fdscb->GetParameter(6);

    // reset sigma and mean to gauss values...
    fdscb->SetParameter(1,fgaus->GetParameter(1));
    fdscb->SetParError (1,fgaus->GetParError(1));
    fdscb->SetParameter(2,fgaus->GetParameter(2));
    fdscb->SetParError (2,fgaus->GetParError(2));
  }
  if (0!=fitstatus){
    cout<<"fit_fdscb() to "<<hrsp->GetName()
        <<" failed. Fitstatus: "<<fitstatus<<endl;
    hrsp->GetFunction("fdscb")->Delete();
  }
  return fitstatus;
}
double fnc_dscb(double*xx,double*pp)
{
  double x   = xx[0];
  // gaussian core
  double N   = pp[0];//norm
  double mu  = pp[1];//mean
  double sig = pp[2];//variance

  // transition parameters
  double a1  = pp[3];
  double p1  = pp[4];
  double a2  = pp[5];
  double p2  = pp[6];

  double u   = (x-mu)/sig;
  double A1  = pow(p1/fabs(a1),p1)*exp(-a1*a1/2);
  double A2  = pow(p2/abs(a2),p2)*exp(-a2*a2/2);
  double B1  = p1/fabs(a1) - fabs(a1);
  double B2  = p2/fabs(a2) - fabs(a2);

  double result(N);
  if      (u<-a1) result *= A1*pow(B1-u,-p1);
  else if (u<a2)  result *= exp(-u*u/2);
  else            result *= A2*pow(B2+u,-p2);
  return result;
}
void fit_double_gaussian(TH1F *&hrsp)
{
  if (0==hrsp) {
    cout<<"ERROR: Empty pointer to fit_double_gaussian()"<<endl;return;
  }
  
  string histname = hrsp->GetName();
  hrsp->Scale(1./hrsp->Integral());
  double mean     = hrsp->GetMean();
  double rms      = hrsp->GetRMS();

  int maxbin    = hrsp->GetMaximumBin();
  double norm1  = hrsp->GetBinContent(maxbin);
  double peak1  = hrsp->GetBinCenter(maxbin);
  double sigma1 = 0.04;

  double norm2  = norm1/2.0;
  double peak2  = mean;
  double sigma2 = 2*rms;

  //cout << " Mean  : "  << mean  << " \t  RMS  : " << rms    << endl;
  //cout << " norm1 : "  << norm1 << " \t  norm2 : " << norm2 << endl;
  //cout << " peak1 : "  << peak1 << " \t  sig1 : " << sigma1 << endl;
  //cout << " peak2 : "  << peak2 << " \t  sig2 : " << sigma2 << endl;

  double fitrange_min = 0.4;
  double fitrange_max = 1.4;

  TF1* fitfnc(0); int fitstatus(-1);
  TF1 *fitg1(0), *fitg2(0);
  fitfnc = new TF1("fdgaus","gaus(0)+gaus(3)",fitrange_min,fitrange_max);
  fitfnc->SetLineColor(1);
  fitfnc->SetLineStyle(2);
  fitfnc->SetLineWidth(2);

  fitfnc->SetParNames("N_{1}", "#mu_{1}", "#sigma_{1}",
		      "N_{2}", "#mu_{2}", "#sigma_{2}");
  fitfnc->SetParameters(norm1, peak1, sigma1, 
   			norm2, peak2, sigma2); 
  fitstatus = hrsp->Fit(fitfnc,"RQ");

  fitfnc->SetParLimits(0,0.01,50*norm1);
  // fitfnc->SetParLimits(1,0.7,1.2);
  // fitfnc->SetParLimits(2,0.01,5.0);

  //fitfnc->SetParLimits(3,0.0,2*norm2);
  // fitfnc->SetParLimits(4,0.2,1.7);
  // fitfnc->SetParLimits(5,1.0,10.0);

  //fitfnc->SetParLimits(4,peak2-3.0*sigma2,peak2+3.0*sigma2);
  //fitfnc->SetParLimits(5,0.10,2.0*sigma2);


  // if (0!=fitstatus){
  //   fitfnc->SetParLimits(4,0.2,1.7);
  //   fitfnc->SetParLimits(5,2.5,20.0);
  //   //cout <<" Not able to Fit this pt bin " << hrsp->GetName() << endl;
  // }
  fitstatus = hrsp->Fit(fitfnc,"RQ");
  hrsp->SetMaximum(norm1+0.2*norm1);
  fitg1 = new TF1("fg1","gaus(0)",fitrange_min,fitrange_max);
  fitg1->SetParameters(fitfnc->GetParameter(0),
		       fitfnc->GetParameter(1),
		       fitfnc->GetParameter(2));
  fitg1->SetLineColor(2);
  fitg1->SetLineStyle(2);
  hrsp->GetListOfFunctions()->Add(fitg1);

  fitg2 = new TF1("fg2","gaus(0)",fitrange_min,fitrange_max);
  fitg2->SetParameters(fitfnc->GetParameter(3),
		       fitfnc->GetParameter(4),
		       fitfnc->GetParameter(5));
  fitg2->SetLineColor(4);
  fitg2->SetLineStyle(4);
  hrsp->GetListOfFunctions()->Add(fitg2);

  // if(hrsp->GetFunction("fdgaus")==0){
  //   cout << "No function recorded in histogram " << hrsp->GetName() << endl;
  // }
  // if (0!=fitstatus){
  //   cout<<"fit_double_gaussian() to "<<hrsp->GetName()
  //       <<" failed. Fitstatus: "<<fitstatus
  //       <<" - FNC deleted."<<endl;
  //   hrsp->GetListOfFunctions()->Delete();
  // }
}
void set_range_truncatedRMS(TH1F *&hist,float frac)
{
  if (0==hist) return;

  const float nevts = hist->Integral(); if (0==nevts) return;
  const int   nbins = hist->GetNbinsX();

  if (frac<=0.0 || frac==1.) return;

  for (int ibin=1;ibin<nbins;++ibin) {
    int binx1   = ibin;
    int binx2   = nbins+1-ibin;
    float ievts = hist->Integral(binx1,binx2);

    if ( (ievts/nevts)>frac ) continue;
    else { hist->GetXaxis()->SetRange(binx1,binx2); break; }
  }
  return;
}
double roundoff(double val)
{
  return ceil( ( val * pow( 10, digi ) ) - 0.49 ) / pow( 10, digi );
}
TLegend *getLegend(double x1, double y1, double x2, double y2)
{
  TLegend *leg = new TLegend(x1,y1,x2,y2,NULL,"BRNDC");
  leg->SetHeader("");
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.06);
  leg->SetLineColor(1);
  leg->SetLineStyle(1);
  leg->SetLineWidth(1);
  leg->SetFillColor(10);
  leg->SetFillStyle(1001);
  return leg;
}

void FitDist(TH1F *&hrsp, double &mean, double &emean, double &sig, double &esig)
{
  //hrsp->Rebin(2);
  set_range_truncatedRMS(hrsp, fracRMS);

  double norm  = hrsp->GetMaximumStored();
  mean  = hrsp->GetMean();
  emean = hrsp->GetMeanError();
  double rms   = hrsp->GetRMS();
  double err   = hrsp->GetRMSError();
  sig   = hrsp->GetRMS()/mean;
  esig = (pow(1/mean,2)*pow(err,2))+(pow(-sig/pow(mean,2),2)*pow(emean,2));
  esig = sqrt(esig);
      
  //TF1 *fgaus = new TF1("fgaus","gaus", mean - 1.50*rms, mean + 1.50*rms);
  TF1 *fgaus = new TF1("fgaus","gaus", 0.65, 1.30);
  fgaus->SetParameters(norm, 1.0, 0.5);
  fgaus->SetParLimits(1,0.38,2.00);
  fgaus->SetParLimits(2,0.2,1.00);
  int fitstatus  = hrsp->Fit(fgaus,"RMLQ");
  //fitstatus=-1;
  // mean  = (fitstatus!=4000) ? hrsp->GetMean()     :  0.5*(fgaus->GetParameter(1) + hrsp->GetMean());
  // emean = (fitstatus!=4000) ? hrsp->GetMeanError(): sqrt(pow(fgaus->GetParError(1),2) + pow(hrsp->GetMeanError(),2));
  mean  = (fitstatus!=4000) ? hrsp->GetMean()     : fgaus->GetParameter(1);
  emean = (fitstatus!=4000) ? hrsp->GetMeanError(): fgaus->GetParError(1);
  sig   = (fitstatus!=4000) ? hrsp->GetRMS()/mean :  fgaus->GetParameter(2)/fgaus->GetParameter(1);
  esig  = (fitstatus!=4000) ? sqrt((pow(1/mean,2)*pow(hrsp->GetRMSError(),2))+(pow(-hrsp->GetRMS()/pow(mean,2),2)*pow(emean,2))) : sqrt((pow(1./fgaus->GetParameter(1),2)*pow(fgaus->GetParError(2),2))+pow(-fgaus->GetParameter(2)/pow(fgaus->GetParameter(1),2),2)*pow(fgaus->GetParError(1),2));

  hrsp->GetFunction("fgaus")->SetBit(TF1::kNotDraw);
  delete fgaus;

  // mean  = (fitstatus!=0) ? hrsp->GetMean()     :  ((fgaus->GetParameter(1)+hrsp->GetMean())*0.5);
  // emean = (fitstatus!=0) ? hrsp->GetMeanError():   0.5*sqrt(pow(hrsp->GetMeanError(),2)+pow(fgaus->GetParError(1),2));//fgaus->GetParError(1);
  // sig   = (fitstatus!=0) ? hrsp->GetRMS()/mean :  fgaus->GetParameter(2)/fgaus->GetParameter(1);
  // esig  = (fitstatus!=0) ? sqrt((pow(1/mean,2)*pow(hrsp->GetRMSError(),2))+(pow(-hrsp->GetRMS()/pow(mean,2),2)*pow(emean,2))) : sqrt((pow(1./fgaus->GetParameter(1),2)*pow(fgaus->GetParError(2),2))+pow(-fgaus->GetParameter(2)/pow(fgaus->GetParameter(1),2),2)*pow(fgaus->GetParError(1),2));

  //cout <<" fit status : " << fitstatus << endl;

  // fit_gaussian(hrsp, 1.0, 1.0, 5);      
  // TF1*  fgaus = (TF1*) hrsp->GetListOfFunctions()->Last();
  // mean  = (fgaus==0) ? hrsp->GetMean()     :  fgaus->GetParameter(1);      
  // emean = (fgaus==0) ? hrsp->GetMeanError():  fgaus->GetParError(1);
  // sig   = (fgaus==0) ? hrsp->GetRMS()/mean :  fgaus->GetParameter(2)/mean;
  // esig  = (fgaus==0) ? hrsp->GetRMSError() :  fgaus->GetParError(2);
	
  //mean  = (fgaus==0) ? hrsp->GetMean()     :  ((fgaus->GetParameter(1)+hrsp->GetMean())*0.5);
  //emean = (fgaus==0) ? hrsp->GetMeanError():   0.5*sqrt(pow(hrsp->GetMeanError(),2)+pow(fgaus->GetParError(1),2));//fgaus->GetParError(1);
  //esig  = (fgaus==0) ? sqrt((pow(1/mean,2)*pow(hrsp->GetRMSError(),2))+(pow(-hrsp->GetRMS()/pow(mean,2),2)*pow(emean,2))) : sqrt((pow(1./fgaus->GetParameter(1),2)*pow(fgaus->GetParError(2),2))+pow(-fgaus->GetParameter(2)/pow(fgaus->GetParameter(1),2),2)*pow(fgaus->GetParError(1),2));
	
	
  //! Double sided crystal ball
  // fit_dscb(hrsp, 1.0, 1.0, 15, Form("ak%d%s",nj+1,algname[nj].c_str()));
  // TF1*  frelrsp = (TF1*) hrsp->GetListOfFunctions()->Last();
  // //mean  = (frelrsp==0) ? hrsp->GetMean()     :  frelrsp->GetParameter(1);
  // //emean = (frelrsp==0) ? hrsp->GetMeanError():  frelrsp->GetParError(1);
  // mean  = (frelrsp==0) ? hrsp->GetMean()     :   ((frelrsp->GetParameter(1)+hrsp->GetMean())*0.5);
  // emean = (frelrsp==0) ? hrsp->GetMeanError():   0.5*sqrt(pow(hrsp->GetMeanError(),2)+pow(frelrsp->GetParError(1),2));
  // sig   = (frelrsp==0) ? hrsp->GetRMS()/mean :  frelrsp->GetParameter(2)/frelrsp->GetParameter(1);
  // esig  = (frelrsp==0) ? sqrt((pow(1/mean,2)*pow(hrsp->GetRMSError(),2))+(pow(-hrsp->GetRMS()/pow(mean,2),2)*pow(emean,2))) 
  // 	: sqrt((pow(1./frelrsp->GetParameter(1),2)*pow(frelrsp->GetParError(2),2))+pow(-frelrsp->GetParameter(2)/pow(frelrsp->GetParameter(1),2),2)*pow(frelrsp->GetParError(1),2));
	
  //sig   = (frelrsp==0) ? hrsp->GetRMS()/mean :  frelrsp->GetParameter(2);
  //esig  = (frelrsp==0) ? hrsp->GetRMSError() :  frelrsp->GetParError(2);
	
  //peak    =(frelrsp==0) ? hrsp[i][j]->GetMean()     : ((frelrsp->GetParameter(1)+hrsp[i][j]->GetMean())*0.5);
  //epeak   =(frelrsp==0) ? hrsp[i][j]->GetMeanError(): 0.5*sqrt(pow(hrsp[i][j]->GetMeanError(),2)+pow(frelrsp->GetParError(1),2));
	
	
  //fit_double_gaussian(hrsp);
  // TF1*  frelrsp = (TF1*) hrsp->GetFunction("fgaus");
  // mean         = (frelrsp==0) ? hrsp->GetMean()     :  frelrsp->GetParameter(1);
  // double emean = (frelrsp==0) ? hrsp->GetMeanError():  frelrsp->GetParError(1);
  // double sig   = (frelrsp==0) ? hrsp->GetRMS()      :  frelrsp->GetParameter(2);
  // double esig  = (frelrsp==0) ? hrsp->GetRMSError() :  frelrsp->GetParError(2);


}
void FitGauss(TH1F* inHist_p, double &mean, double &meanError, double &res, double &resError)
{
  inHist_p->Fit("gaus", "Q L M", "");
  mean = inHist_p->GetFunction("gaus")->GetParameter(1);
  meanError = inHist_p->GetFunction("gaus")->GetParError(1);
  res = inHist_p->GetFunction("gaus")->GetParameter(2);
  resError = inHist_p->GetFunction("gaus")->GetParError(2);
  Float_t prob = inHist_p->GetFunction("gaus")->GetProb();
  // if(TMath::Abs(1.00 - mean) < .01) return;
  Int_t meanBin = inHist_p->FindBin(mean);
  Float_t meanRMS = -1;
  Float_t total = inHist_p->Integral();
  for(Int_t iter = 0; iter < inHist_p->GetNbinsX(); iter++){
    Int_t lowBound = 0;
    if(meanBin - iter > 0) lowBound = meanBin - iter;
    if(inHist_p->Integral(lowBound, meanBin + iter)/total > .95 || lowBound == 0 || inHist_p->GetBinContent(lowBound) < .01){
      meanRMS = inHist_p->GetBinCenter(meanBin + iter) - inHist_p->GetBinCenter(meanBin);
      break;
    }
  }
  double minPt = inHist_p->GetBinCenter(meanBin) - meanRMS;
  double maxPt = inHist_p->GetBinCenter(meanBin) + meanRMS;
  minPt = std::max(std::max(minPt, 0.0), inHist_p->GetXaxis()->GetXmin());
  maxPt = std::min(maxPt, inHist_p->GetXaxis()->GetXmax());
  inHist_p->Fit("gaus", "Q L M", "", minPt, maxPt);
  if(TMath::Abs(1.00 - mean) < TMath::Abs(1.00 - inHist_p->GetFunction("gaus")->GetParameter(1)) && prob > 0.0001){
    inHist_p->Fit("gaus", "Q L M", "");
    return;
  }
  mean = inHist_p->GetFunction("gaus")->GetParameter(1);
  meanError = inHist_p->GetFunction("gaus")->GetParError(1);
  res = inHist_p->GetFunction("gaus")->GetParameter(2);
  resError = inHist_p->GetFunction("gaus")->GetParError(2);
  return;
}
