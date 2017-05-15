#include "TMath.h"
#include "TFile.h"
#include "TLine.h"
#include "TChain.h"
#include "TPaveText.h"
#include "TString.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TROOT.h"
#include "TH1.h"
#include "TEventList.h"
#include "TCanvas.h"
#include "TStopwatch.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TObject.h"
#include "TObjArray.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
using namespace std;
const int nMON = 6;
TString monitor[nMON] = {"diff_qwk_targetX","diff_qwk_targetY","diff_qwk_targetXSlope","diff_qwk_targetYSlope","diff_qwk_bpm3c12X", "asym_qwk_charge"};
const TString units[nMON] = {"nm","nm","nrad","nrad","nm","ppb"};

void GetResiduals(int run, char *scheme, double *slope, double *err){
  gStyle->SetOptFit(1111);
  gStyle->SetTitleW(0.9);
  TString monitor[nMON] = {"diff_qwk_targetX","diff_qwk_targetY","diff_qwk_targetXSlope","diff_qwk_targetYSlope","diff_qwk_bpm3c12X", "asym_qwk_charge"};
  TString units[nMON] = {"nm","nm","nrad","nrad","nm","ppb"};
  double convert[nMON] = {1.0,1.0,1000,1000,1.0,1000};
  double resolution[nMON][2] = {{980,1720},{980,1720},{130,210},{130,210},{1000,1000},{85000,42500}};
  int lSlug = (run==1 ? 40: 137), hSlug = (run==1 ? 136: 321);
  char *det = "mdallbars", *mon = "diff_qwk_targetY";
  vector<double>sl , sle;
  vector<vector<double> >x, xe, y, ye;
  x.resize(nMON);
  xe.resize(nMON);
  y.resize(nMON);
  ye.resize(nMON);
  x.resize(nMON);
  TChain *ch = new TChain("reduced_tree");
  ch->Add(Form("%s/run%i/hydrogen_cell_reduced_tree.root", gSystem->Getenv("DB_ROOTFILES"), run));
  ch->AddFriend("slopes", Form("%s/run%i/hydrogen_cell_bmod_slopes_tree%s.root", gSystem->Getenv("DB_ROOTFILES"), run, scheme));
  ch->AddFriend("corrected_tree", Form("%s/run%i/hydrogen_cell_corrected_tree.root", gSystem->Getenv("DB_ROOTFILES"), run));
  ch->AddFriend("tree", Form("%s/run%i/HYDROGEN-CELL_off_tree.root", gSystem->Getenv("DB_ROOTFILES"), run));

  int islug = lSlug, n=0, N=0;
  TH1D *h[nMON][200], *hc[nMON][200], *hn[nMON][200];

  for(int imon=0;imon<nMON;++imon){
    islug = lSlug;
    N=0;
    while(islug<=hSlug){
      if(islug%5==0) cout<<"Processing slug "<<islug<<endl;
      if(ch->GetEntries(Form("(good&&slopes_exist&&slug==%i)",islug))>10&&islug!=49&&islug!=98&&islug!=99&&islug!=107){
	ch->Draw(Form("sign_correction*(asym_qwk_mdallbars*1000-(mdallbars_targetX*diff_qwk_targetX+mdallbars_targetXSlope*diff_qwk_targetXSlope*1000+mdallbars_targetY*diff_qwk_targetY+mdallbars_targetYSlope*diff_qwk_targetYSlope*1000+mdallbars_bpm3c12X*diff_qwk_bpm3c12X))>>h%i",n),Form("(good&&slopes_exist&&slug==%i)*pow(corrected_asym_qwk_mdallbars.err,-2)",islug),"goff");
	h[imon][N] = (TH1D*)gDirectory->Get(Form("h%i",n));
	ch->Draw(Form("sign_correction*%s*%f>>hc%i", monitor[imon].Data(), convert[imon], n),Form("(good&&slopes_exist&&slug==%i)*pow(corrected_asym_qwk_mdallbars.err,-2)",islug),"goff");
	hc[imon][N] = (TH1D*)gDirectory->Get(Form("hc%i",n));
	ch->Draw(Form("asym_qwk_mdallbars.n>>hn%i", n),Form("(good&&slopes_exist&&slug==%i)",islug),"goff");
	hn[imon][N] = (TH1D*)gDirectory->Get(Form("hn%i",n));
	double error = resolution[imon][run]/sqrt(hn[imon][N]->GetMean()
						  * hn[imon][N]->GetEntries());
	sl.push_back(islug);
	sle.push_back(0);
	x[imon].push_back(hc[imon][N]->GetMean());
	xe[imon].push_back(error);
	y[imon].push_back(h[imon][N]->GetMean());
	ye[imon].push_back(h[imon][N]->GetMeanError());
	++N;
      }
      ++islug;
      ++n;
    }
  } 
  TGraphErrors *gr[nMON];
  TCanvas *c = new TCanvas("c","c",0,0,1500,990);
  TF1 *f = new TF1("f","pol1",-1,1);
  c->Divide(3,2);
  for(int imon=0;imon<nMON;++imon){
    c->cd(imon+1)->SetGrid();
    gr[imon] = new TGraphErrors(N,&x[imon][0],&y[imon][0],&xe[imon][0],&ye[imon][0]);
    gr[imon]->SetLineColor(kBlue);
    gr[imon]->SetMarkerColor(kBlue);
    gr[imon]->SetTitle(Form("Dither-Corrected(%s) asym_%s vs %s ", scheme, det, monitor[imon].Data()));
    gr[imon]->Draw("ap");
    gr[imon]->GetYaxis()->SetTitle(Form("asym_%s (ppb)",det));
    gr[imon]->GetXaxis()->SetTitle(Form("%s (%s)",monitor[imon].Data(),units[imon].Data()));
    f->SetParameters(0,0);
    gr[imon]->Fit(f);
    slope[imon] = f->GetParameter(1);
    err[imon] = f->GetParError(1);
  }
}

int MakeResidualCorrelationsTable(int run=1){
  const int nSCHEME = 2;
  TString scheme[nSCHEME] = {"_omit_coil0","_omit_coil5"};//"_omit_coils3and8","_omit_coils1and6","_omit_coils1and9","_omit_coil0","_omit_coil3","_omit_coil4","_omit_coil5"};
  TString title[nSCHEME] = {"Omit 0","Omit 5"};//,"Omit 3,8","Omit 1,6","Omit 1,9","Omit 0","Omit 3","Omit 4","Omit 5"};
  double slope[nSCHEME][nMON], err[nSCHEME][nMON]; 
  for(int i=0;i<nSCHEME;++i){
    GetResiduals(run , (char*)scheme[i].Data(), slope[i], err[i]);
  }
  cout<<"Dithering Scheme | ";
  for(int imon=0;imon<nMON;++imon){
    char *c = "";
    cout<<monitor[imon].Data()<<" | ";
    if(monitor[imon].Contains("_qwk"))
      monitor[imon].Replace(0, 9, c);
  }
  cout<<endl;
  for(int i=0;i<nSCHEME;++i){
    cout<<title[i].Data()<<" | ";
      for(int imon=0;imon<nMON;++imon)
	cout<<slope[i][imon]<<" +/- "<<err[i][imon]<<" | ";
    cout<<endl;
  }

  cout<<"\\begin{tabular}[h]{|l|c|c|c|c|c|c|}\\hline"<<endl;
  cout<<"Dithering";
  for(int imon=0;imon<nMON;++imon)
    cout<<"&"<<monitor[imon].Data();
  cout<<"\\\\"<<endl;
  cout<<"~Scheme";
  for(int imon=0;imon<nMON;++imon)
    cout<<"&(ppb/"<<units[imon].Data()<<")";
  cout<<"\\\\\\hline"<<endl;
  for(int i=0;i<nSCHEME;++i){
    cout<<title[i].Data();
    for(int imon=0;imon<nMON;++imon)
      printf("& %0.2f$\\pm$%0.2f",slope[i][imon], err[i][imon]);
    cout<<"\\\\\\hline"<<endl;
  }
    cout<<"\\end{tabular}\n";

  return 0;
}
