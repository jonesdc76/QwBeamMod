#include "../include/headers.h"
#include "../include/QwMpsOnly.hh"
#include <iostream>

int AnalyzeBeamOffData(int run){

  TString filename;
  TChain *mps_tree = new TChain("mps_slug");
  QwMpsOnly *mps_only = new QwMpsOnly(mps_tree);

  filename = Form("%s/../beam_off_runs/mps_only_beam_off_%i_full.root", 
		  gSystem->Getenv("MPS_ONLY_ROOTFILES"), run);
  std::cout<<"File: "<<filename.Data()<<std::endl;
  if(mps_only->LoadRootFile(filename, mps_tree, 0)){
    mps_only->SetFileName(filename);
  }

  std::cout << "Creating friend tree with new monitors and ramp_filled leaves."
	    << std::endl;
  if(mps_only->MakeFriendTree(1)){
    std::cout<<"Failed to make friend tree. Exiting.\n";
    return -1;
  }
  std::cout<<"Friend tree successfully created.\n";

  return 0;

}
