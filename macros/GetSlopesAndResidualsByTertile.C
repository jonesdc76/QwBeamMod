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
#include "TProfile.h"
#include "TPaveStats.h"
#include "TVectorD.h"
#include "TMatrixD.h"
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

#include <time.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

using namespace std;

const int nMON = 5, nCOIL = 10, nMOD = 5, nDET = 100;
double MatrixMultiply(TMatrixD a, TMatrixD b, int row, int col){
  double val = 0;
  for(int icol=0;icol<nMON;++icol){
    val += a(row,icol)*b(icol,col);
  }
  
  return val;
}

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
      cout<<monitorList[i].Data()<<"\n";
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
      cout<<detectorList[nDet].Data()<<endl;
      nDet++;
    }
    file.peek();
  }
  file.close();
  return nDet;
}

int GetSlopesAndResidualsByTertile( int run, char *stem = "", int tertile =1,
				    int coil_to_tweak = 0, 
				    double tweak_factor = 1.0, bool fast = 1){
  gStyle->SetOptFit(1111);
  char *chisquare = "_ChiSqMin"; 
  int debug = 0;
  clock_t time = clock();

  if(tweak_factor!=1.0){
    TString nm = stem;
    if(!nm.Contains(Form("%0.1f",tweak_factor))){
      cout<<"File stem does not contain tweak factor. Continue anyway?(y or n)\n";
      cin>>nm;
      if(!nm.Contains("y"))return -1;
    }
  }
  double  lBound = 0, hBound = 360;
  bool fractional_detector_slopes = 1, error_weighted_chi_square = 0;
  TChain *ch = new TChain("mps_slug");
  ch->Add(Form("%s/mps_only_%i_0*.root",gSystem->Getenv("MPS_ONLY_ROOTFILES"), 
	       run));
  ch->Add(Form("%s/mps_only_%i_5*.root",gSystem->Getenv("MPS_ONLY_ROOTFILES"), 
	       run));
  ch->Add(Form("%s/mps_only_%i_1*.root",gSystem->Getenv("MPS_ONLY_ROOTFILES"), 
	       run));
  ch->Add(Form("%s/mps_only_%i_2*.root",gSystem->Getenv("MPS_ONLY_ROOTFILES"), 
	       run));
  ch->Add(Form("%s/mps_only_%i_3*.root",gSystem->Getenv("MPS_ONLY_ROOTFILES"), 
	       run));
  ch->Add(Form("%s/mps_only_%i_ful*.root",gSystem->Getenv("MPS_ONLY_ROOTFILES"), 
	       run));
  if(!ch->AddFriend("mps_slug", Form("%s/mps_only_friend_%i.root",
				     gSystem->Getenv("MPS_ONLY_ROOTFILES"), run))){
    cout<<"Failed to add friend tree in run "<<run<<". Exiting\n";
    return -1;
  }
  if(ch->GetEntries()<1000){
    cout<<"Insufficient entries in run "<<run<<". Exiting\n";
    return -1;
  }


  //Get lists of monitors and detectors
  /////////////////////////////////////////////////////////
  int nDet;
  TString MonitorList[nMON], DetectorList[nDET], mList[nMON], dList[nDET];
  char * configFile = Form("%s/config/setup_mpsonly%s.config",
			   gSystem->Getenv("BMOD_SRC"), stem);
  nDet = GetMonitorAndDetectorLists(ch, MonitorList, DetectorList,0, configFile);
  GetMonitorAndDetectorLists(ch, mList, dList, 1, configFile);
  if(nDet == -1){
    cout<<"Detector list not found. Exiting.\n"<<endl;
    return -1;
  }
  /////////////////////////////////////////////////////////


  //Get coefficients and macrocycle boundaries
  /////////////////////////////////////////////////////////
  char *eq = "[0]+[1]*sin(x*TMath::Pi()/180.0)+[2]*cos(x*TMath::Pi()/180.0)";
  TF1 *f = new TF1("f", eq, lBound, hBound);
  vector<vector<vector<double> > >vDetCoeff, vMonCoeff, vDetMean, vMonMean;
  vDetCoeff.resize(nDet);
  vDetMean.resize(nDet);
  for(int i=0;i<nDet;++i){
    vDetCoeff[i].resize(nCOIL);
    vDetMean[i].resize(nMOD);
  }
  vMonCoeff.resize(nMON);
  vMonMean.resize(nMON);
  for(int i=0;i<nMON;++i){
    vMonCoeff[i].resize(nCOIL);
    vMonMean[i].resize(nMOD);
  }
  vector<int>vLEntry, vHEntry, vId, vMacrocycle;
  vector<vector<int> >vCyclet;

  vCyclet.resize(nMOD);

  TChain *chSlopes = new TChain("slopes");
  chSlopes->Add(Form("%s/../MacrocycleSlopesTree_ChiSqMin.set0_tertile%i.root", 
		     gSystem->Getenv("MPS_ONLY_ROOTFILES"), tertile));

//   TChain *chSlopes2 = new TChain("slopes");
//   chSlopes2->Add(Form("%s/../MacrocycleSlopesTree_ChiSqMin.set0_monitor_test.root",		     gSystem->Getenv("MPS_ONLY_ROOTFILES")));
//   chSlopes->AddFriend(chSlopes2);
  int good;
  double run_number;
  chSlopes->SetBranchStatus("*",0);
  chSlopes->SetBranchStatus("run",1);
  chSlopes->SetBranchAddress("run",&run_number);
  chSlopes->SetBranchStatus("good", 1);
  chSlopes->SetBranchAddress("good", &good);


  int n = 0;
  for(int i=0;i<chSlopes->GetEntries();++i){
    chSlopes->GetEntry(i);
    if(run_number>=run)break;
    n=i;
  }
  chSlopes->ResetBranchAddresses();
  //  ch->AddFriend(chSlopes2);
  chSlopes->SetBranchStatus("*",1);
  //  chSlopes2->SetBranchStatus("*",1);
  chSlopes->SetBranchAddress("run",&run_number);
  chSlopes->SetBranchAddress("good", &good);
  for(int i=n;i<chSlopes->GetEntries();++i){
    chSlopes->GetEntry(i);
    if(run_number>run)break;
    if(run_number != run)continue;
    //    if(!good)continue;
    vMacrocycle.push_back((int)chSlopes->GetLeaf(Form("macrocycle"))->GetValue());
    if(debug == 1) cout<<vMacrocycle.back()<<" ";
    for(int imod=0;imod<nMOD;++imod){
      vCyclet[imod].push_back((int)chSlopes->GetLeaf(Form("cyclet%i",imod))->GetValue());
      if(debug == 1) cout<<vCyclet[imod].back()<<" ";
    }
    vLEntry.push_back((int)chSlopes->GetLeaf("lEntry")->GetValue());
    vHEntry.push_back((int)chSlopes->GetLeaf("hEntry")->GetValue());
    vId.push_back((int)chSlopes->GetLeaf("Id")->GetValue());
    if(debug ==1)
      cout<<vLEntry.back()<<" "<<vHEntry.back()<<" "<<vId.back()<<endl;
    for(int icoil=0;icoil<nCOIL;++icoil){
      for(int imon=0;imon<nMON;++imon){
	if(imon<5){
	  TString nm = Form("%s_Coil%i",mList[imon].Data(),icoil);
	  //	  cout<< Form("%s_Coil%i\n",mList[imon].Data(),icoil);
	  vMonCoeff[imon][icoil].push_back(chSlopes->GetLeaf(nm.Data())->GetValue());
	  if(icoil<nMOD){
	    nm = Form("Coil%i_%s_mean",icoil, mList[imon].Data());
	    double mean = chSlopes->GetLeaf(nm.Data())->GetValue();
	    if(nm.Contains("Slope"))
	      mean *= 1.0e-6;
	    else if(nm.Contains("bpm")||nm.Contains("target"))
	      mean *= 1.0e-3;
	    vMonMean[imon][icoil].push_back(mean);
	  }
	}else if(imon==nMON-1 && icoil<nMOD){
	  ch->Draw(Form("%s:ramp_filled>>m%i%i",MonitorList[imon].Data(),imon,
			icoil),Form("(ErrorFlag == 67207296||ErrorFlag == 0)"
				    " && ramp_good && bm_pattern_number==%i",
				    11+icoil), "prof");
	  TProfile *pr = (TProfile*)gDirectory->Get(Form("m%i%i",imon,icoil));
	  pr->Fit(f,"QRME");
	  vMonMean[imon][icoil].push_back(f->GetParameter(0));
	  vMonCoeff[imon][icoil].push_back(f->GetParameter(1)*1000);
	  vMonCoeff[imon][icoil+nMOD].push_back(f->GetParameter(2)*1000);
	}
	printf("%s_Coil%i %f   %f", mList[imon].Data(),icoil, 
	       vMonCoeff[imon][icoil].back(), vMonMean[imon][icoil%nMOD].back());
	cout<<endl;
      }

      for(int idet=0;idet<nDet;++idet){
	TString nm = Form("%s_Coil%i",dList[idet].Data(),icoil);
	vDetCoeff[idet][icoil].push_back(chSlopes->GetLeaf(nm.Data())->GetValue());
	printf("%s_Coil%i %f   ", dList[idet].Data(),icoil, 
	       vDetCoeff[idet][icoil].back());
	if(icoil<nMOD){
	  nm = Form("Coil%i_%s_mean",icoil, dList[idet].Data());
	  vDetMean[idet][icoil].push_back(chSlopes->GetLeaf(nm.Data())->GetValue()*1.0e-6);
	  printf("%f", vDetMean[idet][icoil].back());
	}
	cout<<endl;
      }
    }
  }
  delete chSlopes;
  /////////////////////////////////////////////////////////


  //Tweak det coefficients and invert to find new slopes
  /////////////////////////////////////////////////////////
  int nCycle = (int)vId.size();
  cout<<nCycle<<" macrocycles found.\n";
  for(int idet = 0;idet<nDet;++idet)
    for(int icycle = 0;icycle<nCycle;++icycle)
      vDetCoeff[idet][coil_to_tweak][icycle] *= tweak_factor;

  double ***detSlopes = new double **[nDet];
  for(int idet=0;idet<nDet;++idet){
    detSlopes[idet] = new double *[nMON];
    for(int imon=0;imon<nMON;++imon){
	detSlopes[idet][imon] = new double[nCycle];
    }
  }
  for(int icycle=0;icycle<nCycle;++icycle){
    //Find monitor system matrix
    TMatrixD R(nMON,nMON);
    for(int irow=0;irow<nMON;++irow){
      for(int icol=irow;icol<nMON;++icol){
	double sum = 0;
	//	cout<<irow<<" "<<icol<<endl;
	for(int icoil=0;icoil<nCOIL;++icoil){
	  sum += vMonCoeff[irow][icoil][icycle] * vMonCoeff[icol][icoil][icycle];
	}
	R(irow, icol) = sum;
	if(irow != icol) R(icol, irow) = sum;
      }
    }
    cout<<"Rmatrix macrocycle "<<icycle<<endl;
    R.Print();
    TMatrixD Rinv = R;
    double determinant;
    Rinv.Invert(&determinant);
    Rinv.Print();
    TMatrixD Id(nMON,nMON);
    Id.Mult(R,Rinv);
    cout<<"Identity:\n";
    Id.Print();
    //Find detector vectors
    for(int idet=0;idet<nDet;++idet){
      TVectorD A(nMON);
      for(int irow=0;irow<nMON;++irow){
	double val = 0;
	for(int icoil=0;icoil<nCOIL;++icoil){
	  val += vDetCoeff[idet][icoil][icycle] * vMonCoeff[irow][icoil][icycle];
	}
	A(irow) = val;
      }
      double slope[nMON] = {0,0,0,0,0};
      for(int irow=0;irow<nMON;++irow){
	double unit_conversion = (mList[irow].Contains("Slope")? 1.0 : 0.001);
	for(int icol=0;icol<nMON;++icol){
	  slope[irow] += Rinv(irow,icol) * A(icol);
	}
	detSlopes[idet][irow][icycle] = slope[irow] * unit_conversion;
      }
    }
  }
  /////////////////////////////////////////////////////////


  //Save new slopes
  /////////////////////////////////////////////////////////
//   std::fstream macrocycle_slopes;
//   macrocycle_slopes.open(Form("%s%s/macrocycle_slopes/"
// 			      "macrocycle_slopes_%i_ChiSqMin.set0.dat", 
// 			      gSystem->Getenv("BMOD_OUT"), stem, run), fstream::out);
  
//   if(macrocycle_slopes.is_open() && macrocycle_slopes.good()){
//     for(int icycle=0;icycle<nCycle;++icycle){
//       macrocycle_slopes<<"Macrocycle "<<vMacrocycle[icycle]<<" \t";
//       for(Int_t imod = 0;imod<nMOD;imod++){
// 	macrocycle_slopes<<"Mod"<<imod<<" "<<vCyclet[imod][icycle]<<" \t";
//       }
//       macrocycle_slopes<<"Id "<<run<<vMacrocycle[icycle]<<"  first_entry "<<
// 	vLEntry[icycle]<<"  last_entry "<<vHEntry[icycle]<<std::endl;
//       for(Int_t idet = 0; idet < nDet; idet++){
// 	macrocycle_slopes << "det " << DetectorList[idet].Data() << std::endl;
// 	for(Int_t imon = 0; imon < nMON; imon++){
// 	  Double_t err = 0;
// 	  macrocycle_slopes << Form("%+14.7e\t%12.5e\n",detSlopes[idet][imon][icycle],err);
// 	}
//       }
//     }
//   }else{
//     std::cout<<"Macrocycle slopes file not open.Exiting.\n";
//     return -1;
//   }
  /////////////////////////////////////////////////////////

  //Find new residuals.
  /////////////////////////////////////////////////////////

  vector<vector<vector<double> > >v_cos_res, v_sin_res, 
    v_cos_res_err, v_sin_res_err;
  v_cos_res.resize(nDet);
  v_sin_res.resize(nDet);
  v_cos_res_err.resize(nDet);
  v_sin_res_err.resize(nDet);
  for(int i=0;i<nDet;++i){
    v_cos_res[i].resize(nMOD);
    v_sin_res[i].resize(nMOD);
    v_cos_res_err[i].resize(nMOD);
    v_sin_res_err[i].resize(nMOD);
  }
  //Using matrices (fast but no error bars calculated).
  /////////////////////////////////////////////////////////
  if(fast){
    for(int icycle=0;icycle<nCycle;++icycle){
      for(int idet=0;idet<nDet;++idet){
	TVectorD v(nMON), DetCoeff(nCOIL);
	for(int imon=0;imon<nMON;++imon){
	  v(imon) = detSlopes[idet][imon][icycle];
	}
	for(int imod=0;imod<nMOD;++imod){
	  DetCoeff(imod) = vDetCoeff[idet][imod][icycle];
	  DetCoeff(imod+nMOD) = vDetCoeff[idet][imod+nMOD][icycle];
	}
	TMatrixD MonCoeffMatrix(nCOIL,nMON);
	for(int irow=0;irow<nCOIL;++irow){
	  for(int icol=0;icol<nMON;++icol){
	    MonCoeffMatrix(irow,icol) = vMonCoeff[icol][irow][icycle];
	  }
	}
	cout<<"Macrocycle "<<icycle<<" "<<dList[0].Data()<<" residuals:\n";
	//    MonCoeffMatrix.Print();
	TVectorD v2(nCOIL);
	for(int icoil=0;icoil<nCOIL;++icoil){
	  double val=0;
	  for(int imon=0;imon<nMON;++imon){
	    double unit_conversion = (mList[imon].Contains("Slope")? 1.0:1000);
	    val += MonCoeffMatrix(icoil,imon) * 
	      detSlopes[idet][imon][icycle]*unit_conversion ;
	  }
	  v2(icoil) = val;
	}
	DetCoeff.Print();
	v2.Print();
	DetCoeff -= v2;
	DetCoeff *= 1e-6;
	for(int imod=0;imod<nMOD;++imod){
	  v_sin_res[idet][imod].push_back(DetCoeff(imod));
	  v_sin_res_err[idet][imod].push_back(0);
	  v_cos_res[idet][imod].push_back(DetCoeff(imod+nMOD));
	  v_cos_res_err[idet][imod].push_back(0);
	}
	DetCoeff.Print();
      }
    }
  }
  //By fitting (slow but with error bars)
  /////////////////////////////////////////////////////////
  if(!fast){
    double ramp_filled, det[nDET], mon[nMON], pat, eFlag;
    int ramp_good;
    ch->SetBranchStatus("*",0);
    ch->SetBranchStatus("ErrorFlag",1);
    ch->SetBranchAddress("ErrorFlag",&eFlag);
    ch->SetBranchStatus("bm_pattern_number",1);
    ch->SetBranchAddress("bm_pattern_number",&pat);
    if(ch->GetLeaf("ramp_good")==0){
      cout<<"Failed to find ramp_good leaf. Exiting.\n";
      return -1;
    }
    ch->SetBranchStatus("ramp_good",1);
    ch->SetBranchAddress("ramp_good",&ramp_good);
    ch->SetBranchStatus("ramp_filled",1);
    ch->SetBranchAddress("ramp_filled",&ramp_filled);
    for(int i=0;i<nDet;++i){
      ch->SetBranchStatus(DetectorList[i].Data(),1);
      ch->SetBranchAddress(DetectorList[i].Data(), &det[i]);
    }
    for(int i=0;i<nMON;++i){
      //    cout<<"Setting status for monitor "<<MonitorList[i].Data()<<endl;
      ch->SetBranchStatus(MonitorList[i].Data(), 1);

      ch->SetBranchAddress(MonitorList[i].Data(), &mon[i]);
    }
    for(int icycle=0;icycle<nCycle;++icycle){
      cout<<"Processing cycle "<<icycle+1<<" of "<<nCycle<<endl;

      vector<vector<double> >v_ramp_filled;
      //det[detector][bm_pattern_number][entry]
      vector<vector<vector<double> > >v_detector;
      v_ramp_filled.resize(nMOD);
      v_detector.resize(nDet);
      for(int i=0;i<nDet;++i){
	v_detector[i].resize(nMOD);
      }

      //Fill vectors with one macrocycle of data
      for(int i=vLEntry[icycle];i<=vHEntry[icycle];++i){
	ch->GetEntry(i);
	if(eFlag != 67207296 && eFlag != 0 || !ramp_good) continue;
	int coil = (int)(pat < 10 ? pat : pat - 11);
	v_ramp_filled[coil].push_back(ramp_filled);
	for(int idet=0;idet<nDet;++idet){
	  double val = det[idet] / abs(vDetMean[idet][coil][icycle]);
	  //	  cout<<dList[idet].Data()<<" "<<val<<endl;
	  double correction = 1.0;
// 	  if(coil==0) {
// 	    cout<<"correction "<<dList[idet].Data()<<"_Coil"<<coil<<"  "
// 		<<vLEntry[icycle]<<" "<<vHEntry[icycle]<<endl;
// 	    cout<<"1 - (";
// 	  }
	  for(int imon=0;imon<nMON;++imon){
// 	    if(coil==0) 
// 	      printf("(%f-%f)*(%f)%s",mon[imon],vMonMean[imon][coil]
// 		     [icycle], detSlopes[idet][imon][icycle],
// 		     (imon==nMON-1 ? ")" : " + "));
	    correction -= (mon[imon] - vMonMean[imon][coil][icycle]) 
	      * detSlopes[idet][imon][icycle] / val;
	  }
// 	  if(coil==0) 
// 	    cout<<"/("<<det[idet]<<"/"<<
// 	      abs(vDetMean[idet][coil][icycle])<<")\n";
       	  v_detector[idet][coil].push_back(val*correction);
	}
      }

      //Graph and fit residuals
      for(int idet=0;idet<nDet;++idet){
	for(int imod=0;imod<nMOD;++imod){
	  TGraph gr = TGraph((int)v_ramp_filled[imod].size(), 
			     &v_ramp_filled[imod][0], 
			     &v_detector[idet][imod][0]);
	  gr.Draw("ap");
	  f->SetParameters(1,0,0);
	  gr.Fit(f,"QRME");
	  if(idet==0)
	    cout<<"Entries: "<<gr.GetN()<<endl;
	  gPad->Update();
	  v_sin_res[idet][imod].push_back(f->GetParameter(1)/abs(f->GetParameter(0)));
	  v_cos_res[idet][imod].push_back(f->GetParameter(2)/abs(f->GetParameter(0)));
	  v_sin_res_err[idet][imod].push_back(f->GetParError(1)/abs(f->GetParameter(0)));
	  v_cos_res_err[idet][imod].push_back(f->GetParError(2)/abs(f->GetParameter(0)));
	}
      }
    }
  }
  /////////////////////////////////////////////////////////


  //Save new residuals
  /////////////////////////////////////////////////////////
  fstream sineData, cosineData;
  sineData.open(Form("%s%s/macrocycle_slopes/run%iSineAmpl%s.set0.dat",
		     gSystem->Getenv("BMOD_OUT"), stem, run, chisquare),fstream::out);
  cosineData.open(Form("%s%s/macrocycle_slopes/run%iCosineAmpl%s.set0.dat",
		     gSystem->Getenv("BMOD_OUT"), stem, run, chisquare),fstream::out);
  if(sineData.is_open() && sineData.good() && cosineData.is_open() && cosineData.good()){

    //Sine residuals
    for(int icycle=0;icycle<nCycle;++icycle){
      printf("Run %i macrocycle %i corrected detector residual SINE component from"
	     " fit.\n", run, vMacrocycle[icycle]);

    printf("                |         Xmod          |         Ymod          |        Emod"
	   "           |        XPmod          |        YPmod          |\n"
	   "***************************************************************************"
	   "**************************************************************\n");
      sineData<<"Macrocycle "<<vMacrocycle[icycle]<<" \t";
      for(Int_t imod = 0;imod<nMOD;imod++){
	sineData<<"Mod"<<imod<<" "<<vCyclet[imod][icycle]<<" \t";
      }
      sineData<<"Id "<<run<<vMacrocycle[icycle]<<"  first_entry "<<
	vLEntry[icycle]<<"  last_entry "<<vHEntry[icycle]<<std::endl;
      for(Int_t idet = 0; idet < nDet; idet++){
	TString detec = DetectorList[idet];
	detec.Resize(16);
	printf("%s%s|", ANSI_COLOR_RESET, detec.Data());
	sineData << DetectorList[idet].Data();
	for(Int_t imod=0;imod<nMOD;imod++){
	  Int_t good = 0;
	  if(abs(v_sin_res[idet][imod][icycle])<abs(v_sin_res_err[idet][imod][icycle]))
	    good = 1;
	  
	  //print to screen
	  printf("%s %+8.2e",(!good ? ANSI_COLOR_RED: ANSI_COLOR_RESET),
		 v_sin_res[idet][imod][icycle]);
	  printf("%s +/- %6.1e %s|",(!good ? ANSI_COLOR_RED:ANSI_COLOR_RESET),
		 v_sin_res_err[idet][imod][icycle],ANSI_COLOR_RESET);
	
	  //save to file
	  sineData << Form("   %+8.3e  %7.3e",v_sin_res[idet][imod][icycle],
			   v_sin_res_err[idet][imod][icycle]);
	}
	printf("\n");
	sineData<<endl;
      }
    }
    //Cosine residuals
    for(int icycle=0;icycle<nCycle;++icycle){
      printf("Run %i macrocycle %i corrected detector residual COSINE component from"
	     " fit.\n", run, vMacrocycle[icycle]);

      printf("                |         Xmod          |         Ymod          |        Emod"
	     "           |        XPmod          |        YPmod          |\n"
	     "***************************************************************************"
	     "**************************************************************\n");
      cosineData<<"Macrocycle "<<vMacrocycle[icycle]<<" \t";
      for(Int_t imod = 0;imod<nMOD;imod++){
	cosineData<<"Mod"<<imod<<" "<<vCyclet[imod][icycle]<<" \t";
      }
      cosineData<<"Id "<<run<<vMacrocycle[icycle]<<"  first_entry "<<
	vLEntry[icycle]<<"  last_entry "<<vHEntry[icycle]<<std::endl;
      for(Int_t idet = 0; idet < nDet; idet++){
	TString detec = DetectorList[idet];
	detec.Resize(16);
	printf("%s%s|", ANSI_COLOR_RESET, detec.Data());
	cosineData << DetectorList[idet].Data();
	for(Int_t imod=0;imod<nMOD;imod++){
	  Int_t good = 0;
	  if(abs(v_cos_res[idet][imod][icycle])<abs(v_cos_res_err[idet][imod][icycle]))
	    good = 1;
	  
	  //print to screen
	  printf("%s %+8.2e",(!good ? ANSI_COLOR_RED: ANSI_COLOR_RESET),
		 v_cos_res[idet][imod][icycle]);
	  printf("%s +/- %6.1e %s|",(!good ? ANSI_COLOR_RED:ANSI_COLOR_RESET),
		 v_cos_res_err[idet][imod][icycle],ANSI_COLOR_RESET);
	
	  //save to file
	  cosineData << Form("   %+8.3e  %7.3e",v_cos_res[idet][imod][icycle],
			   v_cos_res_err[idet][imod][icycle]);
	}
	printf("\n");
	cosineData<<endl;
      }
    }

  }else{
    perror ("Error opening cosine/sine file(s). Exiting.\n");
    return -1;
  }
  /////////////////////////////////////////////////////////

  sineData.close();
  cosineData.close();


  cout<<"Time: "<<double(clock()-time)/CLOCKS_PER_SEC<<" seconds"<<endl;
  return 0;
}
