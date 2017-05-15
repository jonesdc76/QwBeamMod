#include "TString.h"
#include "TLine.h"
#include "TProfile.h"
#include "TSystem.h"
#include <iostream>
#include <stdio.h>
#include <fstream>
#include <cstdio>
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TChain.h"
#include "TMath.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"

Int_t plotDegPerMPS(Int_t runStart = 13846, Int_t runEnd = 19000){
  gStyle->SetOptFit(1111);
  char degrees[255], per[255], mps[255], deg[255];
  Int_t nRuns = 0, n= 0;
  string str;
  Double_t deg_per_mps[5000];
  Double_t runs[5000];
  for(int i = runStart;i<runEnd+1;i++){
    Bool_t endOfFile = 0;
    ifstream file(Form("%s_8Coil/diagnostics/diagnostic_%i_ChiSqMin.set0.dat",
		       gSystem->Getenv("BMOD_OUT"), i));
    if(file.is_open()&&file.good()){
      //skip to end of file where fit results are stored
      nRuns ++;
      getline(file,str);
      while(str != "Ramp information"){
	getline(file,str);
	if(file.eof()){
	  endOfFile = 1;
	  break;
	}
      }
      for(int j=0;j<4;++j)
	getline(file,str);
      if(!endOfFile){
	file>>degrees>>per>>mps>>deg;
	deg_per_mps[n] = atof(deg);
	runs[n] = i;
	getline(file,str);
	//	  if(periodErr[j][n]>1000)
	cout<<runs[n]<<": "<<deg_per_mps[n]<<endl;
      }
      ++n;
    }
    file.close();
  }
  cout<<nRuns<<" total runs found.\n";
  cout<<n<<" runs with fit results found.\n";
  TCanvas *cDeg_Per_Mps = new TCanvas("cDeg_Per_Mps","cDeg_Per_Mps",0,0,800,500);
  TGraph *grDeg_Per_Mps;
  grDeg_Per_Mps = new TGraphErrors(n-1,runs,deg_per_mps);
  grDeg_Per_Mps->Draw("ap");
  grDeg_Per_Mps->SetMarkerColor(kBlue);
  grDeg_Per_Mps->SetLineColor(kBlue);
  grDeg_Per_Mps->SetMarkerStyle(7);
  grDeg_Per_Mps->SetTitle(Form("Deg_Per_Mps Fit Result"));
  grDeg_Per_Mps->GetXaxis()->SetTitle("Run");
  grDeg_Per_Mps->GetYaxis()->SetTitle(Form("Deg_Per_Mps Fit Result"));
  grDeg_Per_Mps->Draw("ap");
  TF1 *f = new TF1("f","pol0",0,360);
  grDeg_Per_Mps->Fit(f);
  gPad->Update();
  return 0;
}
