#include<iostream>
#include<cstdio>
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TString.h"
#include "TSystem.h"
#include <map>
#include <vector>

const int nBPM = 8;

void makeMonopoleTree(int run){
  TFile *oldfile = new TFile(Form("/net/data1/paschkedata1/pass5b_bmod_mpsslug/mps_only_%i_full.root", run));
  TTree *oldtree = (TTree*)oldfile->Get("mps_slug");
  cout<<oldtree->GetEntries()<<" entries"<<endl;
  TString bpmname[nBPM] = {"qwk_bpm3h07aX","qwk_bpm3h07aY", "qwk_bpm3h07bX","qwk_bpm3h07bY", "qwk_bpm3h07cX","qwk_bpm3h07cY","qwk_bpm3h09X","qwk_bpm3h09Y"};

  double pbpm[nBPM], mbpm[nBPM];

  TFile *newfile = new TFile(Form("/net/data1/paschkedata1/pass5b_bmod_mpsslug/mps_only_friend2_%i.root", run),"recreate");

  TTree *tree = new TTree("mps_slug","mps_slug");
  int n7X = 0, n7Y = 0, n7 = 0, n9X = 0,  n9Y = 0, n9 = 0;
  double avg3h07X, avg3h07Y, avg3h07, avg3h09X, avg3h09Y, avg3h09, avgX, avgY, 
    avgAll;
  TBranch *br7X = tree->Branch("avg3h07X",&avg3h07X,"avg3h07X/D");
  TBranch *br7Y = tree->Branch("avg3h07Y",&avg3h07Y,"avg3h07Y/D");
  TBranch *br7 = tree->Branch("avg3h07",&avg3h07,"avg3h07/D");
  TBranch *br9X = tree->Branch("avg3h09X",&avg3h09X,"avg3h09X/D");
  TBranch *br9Y = tree->Branch("avg3h09Y",&avg3h09Y,"avg3h09Y/D");
  TBranch *br9 = tree->Branch("avg3h09",&avg3h09,"avg3h09/D");
  TBranch *brAll = tree->Branch("avgAll",&avgAll,"avgAll/D");
  oldtree->SetBranchStatus("*",0);
  for(int i=0;i<nBPM;++i){
    char *nm = (char*)bpmname[i].Data();
    if(bpmname[i].Contains("7") && bpmname[i].Contains("X")){
      ++n7X;
      ++n7;
    }
    if(bpmname[i].Contains("7") && bpmname[i].Contains("Y")){
      ++n7Y;
      ++n7;
    }
    if(bpmname[i].Contains("9X")){
      ++n9X;
      ++n9;
    }
    if(bpmname[i].Contains("9Y")){
      ++n9Y;
      ++n9;
    }
    oldtree->SetBranchStatus(Form("%sM",nm),1);
    oldtree->SetBranchAddress(Form("%sM",nm), &mbpm[i]);
    oldtree->SetBranchStatus(Form("%sP",nm),1);
    oldtree->SetBranchAddress(Form("%sP",nm), &pbpm[i]);
  }

  for(int i=0;i<oldtree->GetEntries();++i){
    oldtree->GetEntry(i);
    avg3h07X = avg3h07Y = avg3h07 = avg3h09X = avg3h09Y = avg3h09 = avgAll = 0;
    for(int ibpm=0;ibpm<nBPM;++ibpm){
      if(bpmname[ibpm].Contains("7") && bpmname[ibpm].Contains("X")){
	avg3h07X += (pbpm[ibpm] + mbpm[ibpm]) / (2.0 * n7X);
	avg3h07 += (pbpm[ibpm] + mbpm[ibpm]) / (2.0 * n7);
	avgAll += (pbpm[ibpm] + mbpm[ibpm]) / (2.0 * nBPM);
      }
      if(bpmname[ibpm].Contains("7") && bpmname[ibpm].Contains("Y")){
	avg3h07Y += (pbpm[ibpm] + mbpm[ibpm]) / (2.0 * n7Y);
	avg3h07 += (pbpm[ibpm] + mbpm[ibpm]) / (2.0 * n7);
	avgAll += (pbpm[ibpm] + mbpm[ibpm]) / (2.0 * nBPM);
      }
  
      if(bpmname[ibpm].Contains("9X")){
	avg3h09X += (pbpm[ibpm] + mbpm[ibpm]) / (2.0 * n9X);
	avg3h09 += (pbpm[ibpm] + mbpm[ibpm]) / (2.0 * n9);
	avgAll += (pbpm[ibpm] + mbpm[ibpm]) / (2.0 * nBPM);
      }
      if(bpmname[ibpm].Contains("9Y")){
	avg3h09Y += (pbpm[ibpm] + mbpm[ibpm]) / (2.0 * n9Y);
	avg3h09 += (pbpm[ibpm] + mbpm[ibpm]) / (2.0 * n9);
	avgAll += (pbpm[ibpm] + mbpm[ibpm]) / (2.0 *nBPM);
      }
    }
    tree->Fill();
  }

  tree->Write("",TObject::kOverwrite);
  newfile->Close();
  return;
}
