#include "TChain.h"
#include "TString.h"
#include "TROOT.h"
#include "TSystem.h"

#include <iostream>
#include <fstream>
#include <cstdio>


using namespace std;
const int nDET = 100, nMON = 5, nMOD = 5, nCOIL = 10;
Int_t GetMonitorAndDetectorLists(TChain *ch, TString *monitorList, 
				 TString *detectorList, Bool_t trunc,
				 char *config){
  string str;
  ifstream file(config);

  for(int i=0;i<nMON;i++){
    char id[4] = "   ", monitr[20] = "                   ";
    file>>id>>monitr;
    getline(file, str);
    monitorList[i] = TString(monitr);
    if(!ch->GetBranch(monitorList[i].Data())){
      cout<<monitorList[i].Data()<<" "<<ch->GetEntries()<<" missing. Exiting.\n";
      return -1;
    }else{
      if( trunc && monitorList[i].Contains("qwk_"))
	monitorList[i].Replace(0,4,"");
      //      cout<<monitorList[i].Data()<<"\n";
    }
  }
  getline(file, str);
  Int_t nDet = 0;
  for(int i=0;i<nDET&&!file.eof();i++){
    char id[4] = "   ", detectr[20] = "                   ";
    file>>id>>detectr;
    getline(file, str);
    detectorList[nDet] = TString(detectr);
    if(!ch->GetBranch(detectorList[nDet].Data())){
      cout<<detectorList[nDet].Data()<<" missing.\n";
    }else{
      if( trunc && detectorList[nDet].Contains("qwk_"))
	detectorList[nDet].Replace(0,4,"");
      //      cout<<detectorList[nDet].Data()<<endl;
      nDet++;
    }
    file.peek();
  }
  file.close();
  return nDet;
}

int copy_coeff_files(int run, char *stem=""){

  TChain *ch = new TChain("mps_slug");
  ch->Add(Form("%s/mps_only_%i*root",gSystem->Getenv("MPS_ONLY_ROOTFILES"), run));
  //Get lists of monitors and detectors
  /////////////////////////////////////////////////////////
  int nDet;
  TString MonitorList[nMON], DetectorList[nDET];
  char * configFile = Form("%s/config/setup_mpsonly%s.config",
			   gSystem->Getenv("BMOD_SRC"), stem);
  nDet = GetMonitorAndDetectorLists(ch, MonitorList, DetectorList,0, configFile);
  //  GetMonitorAndDetectorLists(ch, mList, dList, 1, configFile);
  if(nDet == -1){
    cout<<"Detector list not found. Exiting.\n"<<endl;
    return -1;
  }
  /////////////////////////////////////////////////////////

  //Open coefficients file in $BMOD_OUT/macrocycles directory
  /////////////////////////////////////////////////////////
  ifstream coeffFile(Form("%s/macrocycle_slopes/coil_coeffs_%i_ChiSqMin.set0.dat",
			  gSystem->Getenv("BMOD_OUT"), run));
  if(!(coeffFile.good() && coeffFile.is_open())){
    cout<<"Coefficients file not found. Exiting.\n";
    return -1;
  }
  /////////////////////////////////////////////////////////

  //Open new coefficients file for writing
  /////////////////////////////////////////////////////////
  ofstream newCoeffFile(Form("%s%s/macrocycle_slopes/coil_coeffs_%i_ChiSqMin."
			     "set0.dat", gSystem->Getenv("BMOD_OUT"), stem, run));
  if(!(newCoeffFile.good() && newCoeffFile.is_open())){
    cout<<"Coefficients file failed to open in new directory. Exiting.\n";
    return -1;
  }
  /////////////////////////////////////////////////////////


  //Transfer desired data to new file
  /////////////////////////////////////////////////////////
  string line;
  getline(coeffFile, line);
  if(line.size()==0){
    cout<<"Empty coefficients file for run "<<run<<". Exiting.\n";
    return -1;
  }
  while(!coeffFile.eof()&& line.find("Macrocycle",0)!=string::npos){
    newCoeffFile<<line<<endl;
    for(int imon=0;imon<nMON;++imon){
      getline(coeffFile, line);
      newCoeffFile<<line<<endl;
      for(int icoil = 0; icoil<nCOIL;++icoil){
	getline(coeffFile, line);
	newCoeffFile<<line<<endl;
      }
    }
    for(int idet=0;idet<nDet;++idet){
      getline(coeffFile, line);
      newCoeffFile<<line<<endl;
      for(int icoil = 0; icoil<nCOIL;++icoil){
	getline(coeffFile, line);
	newCoeffFile<<line<<endl;
      }
    }
    while(!coeffFile.eof()&& line.compare(0,5,"Macro")!=0)
      getline(coeffFile, line);
  }
  /////////////////////////////////////////////////////////




  coeffFile.close();
  newCoeffFile.close();

  return 0;
}
///////////////////////////////////////////////////////////
//To compile at the command line                         //
//g++ -Wall `root-config --ldflags --libs --cflags` -O0  //
//copy_coeff_files.C -o copy_coeff_files                 //
///////////////////////////////////////////////////////////

Int_t main(int argc, char *argv[]){
  if (argc < 1){
    std::cout<<"Usage:copy_coeff_files(int runNumber, char *stem)"
	<<std::endl;
    return 0;
  }else if(argc==3){
  copy_coeff_files(atoi(argv[1]),argv[2]);
  }else copy_coeff_files(atoi(argv[1]));
  return 1;
}
