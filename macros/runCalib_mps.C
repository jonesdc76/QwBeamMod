#include "Calib_mps.C"
#include <iostream>
#include <string>
#include "TString.h"
#include "TProfile.h"
#include "TSystem.h"
#include <vector>
#include <fstream>
#include <cstdio>
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TChain.h"
#include "TMath.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TFile.h"
#include "TTree.h"

using namespace std;

int runCalib_mps(int run, const char *stem = ""){

  char *dir = (char*)gSystem->Getenv("MPS_ONLY_ROOTFILES");
  //  TFile *file = new TFile(Form("%s/Calib_mps_%i.root", dir, run),"recreate");
  TChain *ch = new TChain("mps_slug");
  int n=0;
  n += ch->Add(Form("%s/mps_only_%i_0*.root", dir, run));
  n += ch->Add(Form("%s/mps_only_%i_5*.root", dir, run));
  n += ch->Add(Form("%s/mps_only_%i_1*.root", dir, run));
  n += ch->Add(Form("%s/mps_only_%i_2*.root", dir, run));
  n += ch->Add(Form("%s/mps_only_%i_3*.root", dir, run));
  n += ch->Add(Form("%s/mps_only_%i_ful*.root", dir, run));
  if(n>0){
    cout<<n<<" files added to chain. Adding friend tree.\n";
    if(!ch->AddFriend("mps_slug",Form("%s/mps_only_friend_%i.root", dir, run))){
      cout<<"Unable to add friend tree. Exiting.\n";
      return -1;
    }
  }
  cout<<ch->GetEntries()<<" total entries in chain\n";
  Calib_mps cal = Calib_mps(ch);
  cal.Loop();
  ifstream f(Form("%s/MonitorList%s.dat",dir,stem));
  if(f.is_open()&&f.good()){
    cout<<"Monitor list file found.\n";
  }else{
    cout<<"Monitor list file not found. Exiting.\n";
    return -1;
  }
  vector<TString>vMonitor;
  while(!f.eof()){
    string str;
    TString t;
    f>>t;
    if(!t.Contains("#"))
    vMonitor.push_back(t);
    getline(f,str);
    f.peek();
  }
  for(int i=0;i<(int)vMonitor.size();++i){
    cal.AddMonitor(vMonitor[i].Data());
    cout<<vMonitor[i].Data()<<" added to monitor list.\n";
  }
  printf("getting monitor covariances\n");
  cal.GetCovariances();
  printf("checking controls...\n");
  cal.CheckControls();
  printf("getting nomitor covariances\n");
  cal.GetCovariances();
  printf("generating corrected tree...\n");
  cal.ApplyCorrections();

  return 0;
}
