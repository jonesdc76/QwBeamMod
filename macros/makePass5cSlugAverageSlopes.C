#include <TChain.h>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include "TSystem.h"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TFile.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TTree.h"
#include "TEventList.h"
#include "TGraphErrors.h"
#include "TH1D.h"

const int nDET = 60, nMON = 5, nMOD = 5, nSLUG = 400, MAX_PER_SLUG = 500;
typedef struct leaf_t{
  int slug_lo;
  int slug_hi;
  int run_lo;
  int run_hi;
};

Int_t GetMonitorAndDetectorLists(TString *monitorList, TString *detectorList,
				 Bool_t trunc, char* config){
  string line;
  TChain *ch = new TChain("mps_slug");
  ch->Add(Form("%s/mps_only_13993*.root",gSystem->Getenv("MPS_ONLY_ROOTFILES")));
  ch->AddFriend("mps_slug", Form("%s/mps_only_friend_13993.root",
				 gSystem->Getenv("MPS_ONLY_ROOTFILES")));

  ifstream file(config);

  for(int i=0;i<nMON;i++){
    char id[4] = "   ", monitr[20] = "                   ";
    file>>id>>monitr;
    getline(file, line);
    monitorList[i] = TString(monitr);
    if(!ch->GetBranch(monitorList[i].Data())){
      cout<<monitorList[i].Data()<<" missing. Exiting.\n";
      return -1;
    }else{
      if( trunc && monitorList[i].Contains("qwk_"))
	monitorList[i].Replace(0,4,"");
      cout<<monitorList[i].Data()<<"\n";
    }
  }
  getline(file, line);
  Int_t nDet = 0;
  for(int i=0;i<nDET;i++){
    char id[4] = "   ", detectr[20] = "                   ";
    file>>id>>detectr;
    getline(file, line);
    detectorList[nDet] = TString(detectr);
    if(!ch->GetBranch(detectorList[nDet].Data())){
      cout<<detectorList[nDet].Data()<<" missing.\n";
    }else{
      if( trunc && detectorList[nDet].Contains("qwk_"))
	detectorList[nDet].Replace(0,4,"");
      cout<<detectorList[nDet].Data()<<endl;
      nDet++;
    }
    file.peek();
    if(file.eof())break;
  }
  delete ch;
  return nDet;
}

void makePass5cSlugAveragedSlopes(char* stem ="")
{
  TString MonitorList[nMON], MonitorListFull[nMON];
  TString DetectorList[nDET], DetectorListFull[nDET];

  char* configFile = Form("%s/config/setup_mpsonly%s.config",
			  gSystem->Getenv("BMOD_SRC"), stem);
  GetMonitorAndDetectorLists(MonitorListFull, DetectorListFull, 0, configFile);
  Int_t nDet = GetMonitorAndDetectorLists(MonitorList, DetectorList,1,configFile);
  if(nDet == -1){
    cout<<"Detector list not found. Exiting.\n"<<endl;
    return;
  }


  TChain *ch = new TChain("slug");
  ch->Add("/net/data1/paschkedata1/Set11_dithering_slopes.root");
  double slug;
  ch->SetBranchAddress("slug", &slug);
  double det[nDET][nMON];
  for(int idet=0;idet<nDet;++idet)
    for(int imon=0;imon<nMON;++imon){
      TString name = Form("slope_%s_%s", DetectorListFull[idet].Data(), 
		    MonitorList[imon].Data());
      cout<<"Setting address for "<<name.Data()<<endl;
      ch->SetBranchAddress(name.Data(), &det[idet][imon]);
    }
  for(int i=0;i<ch->GetEntries();++i){
    ch->GetEntry(i);
    cout<<"Slug: "<<(int)slug<<" "<<det[0][0]<<endl;
    for(int idet=0;idet<nDet;++idet){
      TString fileName = Form("/net/data1/paschkedata1/bmod_out%s/slopelists/slopeList_asym_%s_slug%i.dat", stem, DetectorListFull[idet].Data(), (int)slug);
      ofstream file(fileName.Data());
      if(!(file.good()&&file.is_open())){
	cout<<"Not able to open slope list file. Exiting.\n";
	return;
      }
      for(int imon=0;imon<nMON;++imon)
	file<<Form("%s  %15.10e   0.00000", MonitorList[imon].Data(),
		   (MonitorList[imon].Contains("Slope") ? 1.0 : 1000.0) 
		   * det[idet][imon])<<endl;
      file.close();
    }
  }

  return;
}
