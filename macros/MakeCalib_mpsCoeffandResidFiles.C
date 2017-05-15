#include <sstream>
#include <map>
#include <list>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdio>
#include "TString.h"
#include "TProfile.h"
#include "TSystem.h"
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
const int nMOD = 5, nCOIL = 10;


TH1D* lookup_histo(std::string name, int pattern, char comp)
{
   // Private convenience method for looking up the 1D histogram
   // containing fit coefficients for:
   //  *) name = signal name
   //  *) pattern = 11,12,13,14,15 for specific modulation patterns,
   //               or 0 for all
   //  *) comp = 'r' for cosine term, 'i' for sine term, 'c' for constant term

   std::stringstream sname;
   sname << name;
   if (pattern)
      sname << "_" << pattern;
   sname << comp;
   TH1D *hist = (TH1D*)gROOT->FindObject(sname.str().c_str());
   if (hist == 0) {
      std::cerr << "Calib_mps::lookup_histo error - "
                << "histogram named " << sname.str()
                << " not found, cannot continue."
                << std::endl;
      exit(1);
   }
   return hist;
}



int MakeCalib_mpsCoeffandResidFiles(char* stem, int runlo = 10000, 
				    int runhi = 19000)
{
  char *dir = (char*)gSystem->Getenv("MPS_ONLY_ROOTFILES");


  //Get list of monitor and detectors
  //////////////////////////////////////////////////
  ifstream fMonList(Form("%s/MonitorList%s.dat",dir, stem));
  if(fMonList.is_open()&&fMonList.good()){
    cout<<"Monitor list file found.\n";
  }else{
    cout<<"Monitor list file not found. Exiting.\n";
    return -1;
  }
  vector<TString>vMonitor;
  while(!fMonList.eof()){
    string str;
    TString t;
    fMonList>>t;
    if(!t.Contains("#"))
    vMonitor.push_back(t);
    getline(fMonList,str);
    fMonList.peek();
  }
  fMonList.close();
  ifstream fDetList(Form("%s/DetectorList%s.dat",dir, stem));
  if(fDetList.is_open()&&fDetList.good()){
    cout<<"Detector list file found.\n";
  }else{
    cout<<"Detector list file not found. Exiting.\n";
    return -1;
  }
  vector<TString>vDetector;
  while(!fDetList.eof()){
    string str;
    TString t;
    fDetList>>t;
    if(!t.Contains("#"))
    vDetector.push_back(t);
    getline(fDetList,str);
    fDetList.peek();
  }
  fDetList.close();
  int nMon = (int)vMonitor.size();
  int nDet = (int)vDetector.size();
  ////////////////////////////////////////////////////

  TFile f = TFile(Form("%s/../MacrocycleSlopesTree_ChiSqMin.set0.root",dir));
  TTree *slTree = (TTree*)f.Get("slopes");
  cout<<slTree->GetEntries()<<" entries in slopes tree.\n";
  slTree->SetBranchStatus("*",0);
  int good;
  double cyclet[nMOD], hEntry, lEntry, Id, macrocycle, run;
  for(int i=0;i<nMOD;++i){
    slTree->SetBranchStatus(Form("cyclet%i",i),1);
    slTree->SetBranchAddress(Form("cyclet%i",i), &cyclet[i]);
  }
  slTree->SetBranchStatus("run", 1);
  slTree->SetBranchAddress("run",&run);
  //  slTree->SetBranchStatus("good", 1);
  //  slTree->SetBranchAddress("good",&good);
  slTree->SetBranchStatus("lEntry", 1);
  slTree->SetBranchAddress("lEntry",&lEntry);
  slTree->SetBranchStatus("hEntry", 1);
  slTree->SetBranchAddress("hEntry",&hEntry);
  slTree->SetBranchStatus("macrocycle", 1);
  slTree->SetBranchAddress("macrocycle",&macrocycle);
  slTree->SetBranchStatus("hId", 1);
  slTree->SetBranchAddress("hId",&Id);
  int range = runhi-runlo+1, prev_run=0;
  TTree *tr;
  double *mon = new double[nMon];
  double *det = new double[nDet];
  std::map<std::string, double*>mMonitor, mDetector;
  std::map<std::string, TH1D*>mMonHisto, mDetHisto;
  for(int idet=0;idet<nDet;++idet){
    mDetector[vDetector[idet].Data()] = &det[idet];
  }
  for(int imon=0;imon<nMon;++imon){
    mMonitor[vMonitor[imon].Data()] = &mon[imon];
  }
  for(int i=0;i<range;++i){
    slTree->GetEntry(i);
    if(i%1000==0)cout<<"Processing entry "<<i<<" of "<<slTree->GetEntries()<<endl;
    //whenever run number changes load new root tree
    if(run!=prev_run){
      //      cout<<run<<"\n";
      TFile fileCurrentRun = TFile(Form("%s/Calib_mps_%i%s.root",dir, (int)run, 
					stem));
      if(!fileCurrentRun.IsOpen()){
	cout<<"Calib_mps rootfile not found.\n";
	continue;
      }
      tr = (TTree*)fileCurrentRun.Get("corrected");
      if(!tr){
	printf("Cannot find corrected tree in rootfile for run %i\n", int(run));
	continue;
      }
      cout<<tr->GetEntries()<<" entries.\n";
      tr->ResetBranchAddresses();
      tr->SetBranchStatus("*",0);
      std::map<std::string, Double_t*>::iterator iter;
      tr->SetBranchStatus("*",0);
      for(iter = mDetector.begin();iter != mDetector.end(); ++iter){
	std::string name(iter->first);
	//	cout<<name<<endl;
	tr->SetBranchStatus(name.c_str(), 1);
	tr->SetBranchAddress(name.c_str(), iter->second);
      }
      for(iter = mMonitor.begin();iter != mMonitor.end(); ++iter){
	std::string name(iter->first);
       	cout<<name<<endl;
	tr->SetBranchStatus(name.c_str(), 1);
	tr->SetBranchAddress(name.c_str(), iter->second);
      }
//       for(int ent=(int)lEntry;ent<=(int)hEntry;++ent){
// 	tr->GetEntry(ent);
// 	cout<<vDetector[0].Data()<<" "<<*mDetector[vDetector[0].Data()]<<endl;
// 	break;

//       }
//       for(int imod=0;imod<nMOD;++imod){
// 	for(int idet=0;idet<nDet;++idet){
// 	  //constant(offset), real(cosine) and imaginary(sine) terms
// 	  string name = vDetector[idet].Data();
// 	  mDetHisto[Form("%s_%ic", name.c_str(), imod+11)] = 
// 	    lookup_histo(name, imod+11, 'c');
// 	  mDetHisto[Form("%s_%ii",name.c_str(), imod+11)] = 
// 	    lookup_histo(name, imod+11, 'i');
// 	  mDetHisto[Form("%s_%ir",name.c_str(), imod+11)] = 
// 	    lookup_histo(name, imod+11, 'r');
// 	}
// 	for(int imon=0;imon<nMon;++imon){
// 	  //constant(offset), real(cosine) and imaginary(sine) terms
// 	  string name = vMonitor[imon].Data();
// 	  mMonHisto[Form("%s_%ic", name.c_str(), imod+11)] = 
// 	    lookup_histo(name, imod+11, 'c');
// 	  mMonHisto[Form("%s_%ii", name.c_str(), imod+11)] = 
// 	    lookup_histo(name, imod+11, 'i');
// 	  mMonHisto[Form("%s_%ir",name.c_str(), imod+11)] = 
// 	    lookup_histo(name, imod+11, 'r');
// 	}
//       }
    }
    //    mMonHisto["dslumi_even_11_r"]->Draw();
    prev_run = (int)run;
  }

  return 0;
}
