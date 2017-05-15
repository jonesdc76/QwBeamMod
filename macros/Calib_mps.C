//
// Calib_mps : ROOT tree analyzer class for scanning Qweak ROOT trees
//             for driven beam modulation periods, and extracting the
//             sensitivity of each detector signal to the modulation.
//
// Author: richard.t.jones at uconn.edu
// Version: march 24, 2015
//
// For examples of how to use this class,
// see python script mps.py and functions therein.
//

//#define INDIVIDUAL_CYCLE_GAINS_EVALUATION 1
//#define PRINT_MCONTROL_EIGENVECTORS 1
#define PRINT_MON2MD_EIGENVECTORS 1

#define MIN_CALIBRATION_CYCLE 1000
#define EIGENVALUE_PRECISION 1e-2
#define MAX_MCONTROLS_RANK 5

#include <iostream>
#include <sstream>
#include <list>
#include <math.h>

#define Calib_mps_cxx
#include "Calib_mps.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMatrixDSym.h>
#include <TDecompSVD.h>

void Calib_mps::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L Calib_mps.C
//      Root > Calib_mps t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//    b_branchname->GetEntry(ientry); //read only this branch

   if (fChain == 0)
      return;

   TFile fout(fWorkfile.c_str(), "recreate");
   CreateOutputs();

   // in the first pass, extract the sensitivity coefficients
 
   double X00, X10, X01, X11, X20, X02;
   X00 = X10 = X01 = X11 = X20 = X02 = 0;
   std::map<std::string, double> Y00, Y10, Y01;
   std::map<std::string, double> Z00;
   int pcount = 0;
   int pattern_number = 0;
   int pattern_counter[20] = {0};
   Long64_t nentries = fChain->GetEntries();

   std::cout << "Starting first pass, " << nentries << " rows..." << std::endl;

   for (Long64_t jentry=0; jentry <= nentries; jentry++) {
      if (jentry < nentries)
         fChain->GetEntry(jentry);
      else
         bm_pattern_number = 0;

      // at the end of a pattern, compute the sensitivities
      // and store them and their errors in histograms
 
      if (bm_pattern_number != pattern_number) {
         ++pattern_counter[pattern_number];
         if (pcount > MIN_CALIBRATION_CYCLE) {
            TMatrixDSym M(3);
            M(0,0) = X20;
            M(0,1) = M(1,0) = X11;
            M(0,2) = M(2,0) = X10;
            M(1,1) = X02;
            M(1,2) = M(2,1) = X01;
            M(2,2) = X00;
            TMatrixDSym MI(M);
            double determ;
            MI.Invert(&determ);
            std::map<std::string, Double_t*>::iterator iter;
            for (iter = fSignames.begin(); iter != fSignames.end(); ++iter) {
               double B[3];
               B[0] = Y10[iter->first];
               B[1] = Y01[iter->first];
               B[2] = Y00[iter->first];
               double areal = MI(0,0) * B[0] + MI(0,1) * B[1] + MI(0,2) * B[2];
               double aimag = MI(1,0) * B[0] + MI(1,1) * B[1] + MI(1,2) * B[2];
               double acons = MI(2,0) * B[0] + MI(2,1) * B[1] + MI(2,2) * B[2];
               double chisqr = Z00[iter->first]
                             - 2 * (areal * B[0] + aimag * B[1] + acons * B[2])
                             + areal * areal * M(0,0) + aimag * aimag * M(1,1)
                             + acons * acons * M(2,2)
                             + 2 * (areal * aimag * M(0,1) 
                                    + areal * acons * M(0,2)
                                    + aimag * acons * M(1,2));
               double varresid = chisqr / (X00 - 3);
               double areal_err = sqrt(MI(0,0) * varresid);
               double aimag_err = sqrt(MI(1,1) * varresid);
               double acons_err = sqrt(MI(2,2) * varresid);
               TH1D *hreal = lookup_histo(iter->first, pattern_number, 'r');
               hreal->SetBinContent(pattern_counter[pattern_number], areal);
               hreal->SetBinError(pattern_counter[pattern_number], areal_err);
               TH1D *himag = lookup_histo(iter->first, pattern_number, 'i');
               himag->SetBinContent(pattern_counter[pattern_number], aimag);
               himag->SetBinError(pattern_counter[pattern_number], aimag_err);
               TH1D *hcons = lookup_histo(iter->first, pattern_number, 'c');
               hcons->SetBinContent(pattern_counter[pattern_number], acons);
               hcons->SetBinError(pattern_counter[pattern_number], acons_err);
            }
         }
         pcount = 0;
         pattern_number = bm_pattern_number;
         X00 = X10 = X01 = X11 = X20 = X02 = 0;
         std::map<std::string, Double_t*>::iterator iter;
         for (iter = fSignames.begin(); iter != fSignames.end(); ++iter) {
            Y00[iter->first] = 0;
            Y01[iter->first] = 0;
            Y10[iter->first] = 0;
            Z00[iter->first] = 0;
         }
      }
      if (jentry == nentries)
         break;

      // during modulation cycles, accummulate statistics on each signal

      if (ramp_good) {
         double costheta = cos(ramp_filled * M_PI / 180);
         double sintheta = sin(ramp_filled * M_PI / 180);
         X00 += 1;
         X10 += costheta;
         X01 += sintheta;
         X11 += costheta * sintheta;
         X20 += costheta * costheta;
         X02 += sintheta * sintheta;
         std::map<std::string, Double_t*>::iterator iter;
         for (iter = fSignames.begin(); iter != fSignames.end(); ++iter) {
            double y = *iter->second;
            Y00[iter->first] += y;
            Y10[iter->first] += y * costheta;
            Y01[iter->first] += y * sintheta;
            Z00[iter->first] += y * y;
         }
         ++pcount;
      }
   }
 
   // do a second pass to compute the residuals

   pcount = 0;
   pattern_number = 0;
   for (int pat=11; pat < 16; ++pat)
      pattern_counter[pat] = 0;
   TTree *resid = (TTree*)gROOT->FindObject("resid");

   std::cout << "Starting second pass..." << std::endl;

   std::map<std::string, double> areal, aimag;
   for (Long64_t jentry=0; jentry < nentries; jentry++) {
      fChain->GetEntry(jentry);
      if (bm_pattern_number != pattern_number) {
         pattern_number = bm_pattern_number;
         ++pattern_counter[pattern_number];
         std::map<std::string, Double_t*>::iterator iter;
         for (iter = fSignames.begin(); iter != fSignames.end(); ++iter) {
            FindFitCoeff(iter->first, 'r', pattern_number,
                         pattern_counter[pattern_number], &areal[iter->first]);
            FindFitCoeff(iter->first, 'i', pattern_number,
                         pattern_counter[pattern_number], &aimag[iter->first]);
         }
      }
      std::map<std::string, Double_t*>::iterator iter;
      for (iter = fSignames.begin(); iter != fSignames.end(); ++iter) {
         double costheta = cos(ramp_filled * M_PI / 180);
         double sintheta = sin(ramp_filled * M_PI / 180);
         double correction = areal[iter->first] * costheta 
                           + aimag[iter->first] * sintheta;
         *iter->second -= correction;
      }
      resid->Fill();
   }
   fout.Write();
   fout.Close();

   // If anyone has this workfile open, close and reopen it
   // for them so they can see the updates. When I reopen it
   // I do so in readonly mode because it is an i/o collision
   // if it was open for writing when this method has opened
   // it independently for writing using a private descriptor.

   TFile *fin = (TFile*)gROOT->FindObject(fWorkfile.c_str());
   if (fin && fin->IsOpen()) {
      fin->Close();
      fin->Open(fWorkfile.c_str());
   }
}

Int_t Calib_mps::CreateOutputs()
{
   // Two kinds of outputs are generated when Loop() runs:
   //  1) modulation sensitivities for each detector signal (1D histograms)
   //  2) residual detector signals with driven modulation removed (tree)
   // CreateOutputs allocates empty objects to receive and store outputs
   // for every signal variable in the mps_slug input tree.

   if (fChain == 0) 
      return 0;

   int period_counter[20] = {20 * 0};
   int pattern_number = 0;
   Long64_t nentries = fChain->GetEntries();
   for (Long64_t jentry=0; jentry < nentries; jentry++) {
      b_bm_pattern_number->GetEntry(fChain->LoadTree(jentry));
      if (bm_pattern_number != pattern_number) {
         pattern_number = bm_pattern_number;
         ++period_counter[pattern_number];
      }
   }

   std::string title(fChain->GetTitle());
   title += " residuals after modulation removal";
   TTree *resid = new TTree("resid", title.c_str());

   std::map<std::string, Double_t*>::iterator iter;
   for (iter = fSignames.begin(); iter != fSignames.end(); ++iter) {
      for (int pat = 11; pat < 16; ++pat) {
         int reps = period_counter[pat];
         make_histo("average offset", iter->first, pat, 'r', reps);
         make_histo("average offset", iter->first, pat, 'i', reps);
         make_histo("average offset", iter->first, pat, 'c', reps);
      }
      std::string bname(iter->first);
      resid->Branch(bname.c_str(), iter->second);
   }
   for (iter = fAuxnames.begin(); iter != fAuxnames.end(); ++iter) {
      std::string bname(iter->first);
      resid->Branch(bname.c_str(), iter->second);
   }
   return 0;
}

TH1D *Calib_mps::lookup_histo(std::string name, int pattern, char comp)
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

TH1D *Calib_mps::make_histo(std::string type, std::string name,
                            int pattern, char comp, int range)
{
   // Private convenience method to consolidate creation of 1D histograms
   // to hold coefficients from fits to the oscillator waveform.

   std::stringstream sname, stitle;
   sname << name;
   if (pattern)
      sname << "_" << pattern;
   sname << comp;
   stitle << type << " of " << name;
   if (pattern)
      stitle << " in modulation mode " << pattern;
   if (comp == 'r')
      stitle << ", cosine term";
   else if (comp == 'i')
      stitle << ", sine term";
   else
      stitle << ", constant term";
   TH1D *hist = new TH1D(sname.str().c_str(), stitle.str().c_str(),
                         range, 0, range);
   hist->SetStats(0);
   hist->GetXaxis()->SetTitle("modulation cycle");
   hist->GetYaxis()->SetTitle("gain coefficient");
   hist->GetYaxis()->SetTitleOffset(1.5);
   return hist;
}

Int_t Calib_mps::CorrelateMonitors()
{
   // Pass through the signal data in the input tree and compute the
   // covariance matrix of the signals designated as monitors. Only the
   // residue (data minus fit) of the signals are used in the computation
   // of the covariance matrix elements. The result is saved in a new
   // TMatrixDSym pointed to by fMcovar. If the fVevalues/fMevectors
   // pointers are loaded then the renormalized monitors are evaluated
   // and used, otherwise the monitor basis is the list of signals in
   // the order they are specified in fMonnames.

   TFile *fin = (TFile*)gROOT->FindObject(fWorkfile.c_str());
   if (fin == 0) {
      fin = new TFile(fWorkfile.c_str());
   }
   fin->cd();

   int N = fMonnames.size();
   int NxN = N * (N+1) / 2;
   double Cov[NxN];
   for (int nn = 0; nn < NxN; ++nn)
      Cov[nn] = 0;

   double a[N], b[N], c[N];
   int pattern_number = 0;
   int pattern_counter[20] = {0};
   int nentries = fChain->GetEntries();
   for (Long64_t jentry=0; jentry < nentries; jentry++) {
      fChain->GetEntry(jentry);
      if (bm_pattern_number != pattern_number) {
         pattern_number = bm_pattern_number;
         ++pattern_counter[pattern_number];
         for (unsigned int n = 0; n < fMonnames.size(); ++n) {
            FindFitCoeff(fMonnames[n], 'r', pattern_number,
                         pattern_counter[pattern_number], &a[n]);
            FindFitCoeff(fMonnames[n], 'i', pattern_number,
                         pattern_counter[pattern_number], &b[n]);
            FindFitCoeff(fMonnames[n], 'c', pattern_number,
                         pattern_counter[pattern_number], &c[n]);
         }
      }

      // subtract off the fit value from the measured one
 
      double *mvalue = new double[N];
      for (unsigned int n = 0; n < fMonnames.size(); ++n) {
         mvalue[n] = *fMonvalues[n];
         double costheta = cos(ramp_filled * M_PI / 180);
         double sintheta = sin(ramp_filled * M_PI / 180);
         mvalue[n] -= a[n] * costheta + b[n] * sintheta + c[n];
      }
 
      // apply the monitor normalization, if available
 
      double *nvalue = mvalue;
      if (fMnormal.GetNrows() == N) {
         nvalue = new double[N];
         for (int n = 0; n < N; ++n) {
            nvalue[n] = 0;
            for (int m = 0; m < N; ++m) {
               nvalue[n] += fMnormal(n, m) * mvalue[m];
            }
         }
         delete [] mvalue;
      }

      // accummulate the covariance statistics, assuming zero mean
 
      int index = 0;
      for (unsigned int n = 0; n < fMonnames.size(); ++n) {
         for (unsigned int m = 0; m <= n; ++m, ++index) {
            Cov[index] += nvalue[m] * nvalue[n];
         }
      }
      delete [] nvalue;
   }
   
   // save the result in the fMcovar matrix

   fMcovar.ResizeTo(N, N);
   int index = 0;
   for (unsigned int n = 0; n < fMonnames.size(); ++n) {
      for (unsigned int m = 0; m <= n; ++m, ++index) {
         fMcovar(m,n) = fMcovar(n,m) = Cov[index] / nentries;
      }
   }
   return N;
}

Int_t Calib_mps::NormalizeMonitors()
{
   // Diagonalize the monitors covariance matrix in fMcovar and store
   // the list of eigenvalues[eigenvectors] in fVevalues[fMevectors].
   // If fMcovar is null then CorrelateMonitors() is called to fill it.

   int N = fMonnames.size();
   if (N == 0)
      return 0;
   if (fMcovar.GetNrows() != N)
      CorrelateMonitors();

   // before trying to diagonalize the matrix, first balance the matrix

   TMatrixDSym Sigma(N);
   Sigma *= 0;
   for (int n = 0; n < N; ++n) {
      Sigma(n, n) = 1 / sqrt(fMcovar(n, n));
   }
   TMatrixDSym Correl(fMcovar);
   Correl.Similarity(Sigma);
   fMevectors.ResizeTo(N, N);
   fMevectors = (TMatrixD)Correl.EigenVectors(fVevalues);

   // compute the normalization transform matrix

   fMnormal.ResizeTo(N, N);
   fMnormal *= 0;
   int rank;
   for (rank = 0; rank < fVevalues.GetNrows(); ++rank) {
      if (fVevalues[rank] > 1e-30)
         fMnormal(rank, rank) = 1 / sqrt(fVevalues[rank]);
      else
         fMnormal(rank, rank) = 0;
   }
   fMnormal.Similarity(fMevectors);

   // now restore the normalization factors taken out earlier

   for (int n = 0; n < N; ++n) {
      for (int m = 0; m < N; ++m) {
         fMnormal(n, m) *= Sigma(m, m);
      }
   }
   return rank;
}

Int_t Calib_mps::EvaluateControls(Int_t pattern_number, Int_t pattern_counter[])
{
   // Evaluates the matrix that converts a vector of monitor values to an
   // estimate for the control vector, and stores it in fMcontrol. The
   // inputs to this calculation are the monitor normalization matrix stored
   // in fMnormal and the fit coefficients saved in ROOT histograms found
   // in the current working directory. The pattern_counter vector argument
   // selects which set of fit coefficients to use. The components of
   // pattern counter should be set to the index of the nearest modulation
   // period for each of the 5 valid types, indexed as bm_pattern_number.
   // As a by-product of this operation, the individual offsets for each
   // monitor signal are stored in fVoffsets. These offsets must be added
   // to the raw monitor signals before the control matrix is applied.

   TFile *fin = (TFile*)gROOT->FindObject(fWorkfile.c_str());
   if (fin == 0) {
      fin = new TFile(fWorkfile.c_str());
   }
   fin->cd();

   int Ncon = 10;
   int Nmon = fMnormal.GetNrows();
   if (Nmon < Ncon) {
      std::cerr << "Calib_mps::EvaluateControls error - "
                << "insufficient monitors have been declared, "
                << "or monitor normalization has not been done."
                << std::endl;
      exit(1);
   }
   fVoffsets.ResizeTo(Nmon);

   // pull the gains and offsets from stored histograms

   TMatrixD gains(Nmon, Ncon);
   for (int mon = 0; mon < Nmon; ++mon) {
      for (int pat = 11; pat < 16; ++pat) {
         int con = pat - 11;
         FindFitCoeff(fMonnames[mon], 'r', pat,
                      pattern_counter[pat], &gains(mon, con));
         FindFitCoeff(fMonnames[mon], 'i', pat, 
                      pattern_counter[pat], &gains(mon, con + 5));
      }
      FindFitCoeff(fMonnames[mon], 'c', pattern_number,
                   pattern_counter[pattern_number], &fVoffsets[mon]);
   }

   // apply the normalization transform and compute the SVD

   TMatrixD normalized_gains(Nmon, Ncon);
   normalized_gains.Mult(fMnormal, gains);
   TDecompSVD svd(normalized_gains);

   // apply the pseudo-inversion

   TMatrixD minv(Nmon, Nmon);
   minv.Transpose(svd.GetU());
   minv *= fMnormal;
   TVectorD sig(svd.GetSig());
   TMatrixD sinv(Ncon, Nmon);
   sinv *= 0;
   int rank = 0;
   for (int con = 0; con < Ncon; ++con) {
      if ((sig[con] / sig[0]) < EIGENVALUE_PRECISION)
         break;
      sinv(con, con) = 1 / sig[con];
      ++rank;
   }
   fMcontrols.ResizeTo(Ncon, Nmon);
   fMcontrols.Mult(svd.GetV(), sinv *= minv);

   TMatrixD Ctranspose(Nmon, Ncon);
   Ctranspose.Transpose(fMcontrols);
   TDecompSVD Csvd(Ctranspose);

# if PRINT_MCONTROL_EIGENVECTORS
   std::cout << "EvaluateControls at pattern_counters";
   for (int pat=11; pat < 16; ++pat)
      std::cout << " " << pattern_counter[pat];
   std::cout << std::endl;
   std::cout << "eigenvectors of Mcontrols on the monitors side" << std::endl;
   std::cout << " eigenvalue    ";
   for (int n = 0; n < MAX_MCONTROLS_RANK; ++n) {
      char s[99];
      sprintf(s, " %10.4g", Csvd.GetSig()[n]);
      std::cout << s;
   }
   std::cout << std::endl;
   for (int mon = 0; mon < (int)fMonnames.size(); ++mon) {
      char s[99];
      sprintf(s, "%15s", fMonnames[mon].c_str());
      std::cout << s;
      for (int n = 0; n < MAX_MCONTROLS_RANK; ++n) {
         sprintf(s, " %10.4g", Csvd.GetU()(mon, n));
         std::cout << s;
      }
      std::cout << std::endl;
   }
   std::cout << "eigenvectors of Mcontrols on the modulators side" << std::endl;
   std::cout << " eigenvalue    ";
   for (int n = 0; n < MAX_MCONTROLS_RANK; ++n) {
      char s[99];
      sprintf(s, " %10.4g", Csvd.GetSig()[n]);
      std::cout << s;
   }
   std::cout << std::endl;
   for (int con = 0; con < 5; ++con) {
      char s[99];
      sprintf(s, " D[%d] cos(ramp)", con);
      std::cout << s;
      for (int n = 0; n < MAX_MCONTROLS_RANK; ++n) {
         sprintf(s, " %10.4g", Csvd.GetV()(con, n));
         std::cout << s;
      }
      std::cout << std::endl;
      sprintf(s, " D[%d] sin(ramp)", con);
      std::cout << s;
      for (int n = 0; n < MAX_MCONTROLS_RANK; ++n) {
         sprintf(s, " %10.4g", Csvd.GetV()(con + 5, n));
         std::cout << s;
      }
      std::cout << std::endl;
   }
   std::cout << std::endl;
# endif

#if MAX_MCONTROLS_RANK < 10
   TMatrixD Cdiag(Ncon, Nmon);
   Cdiag *= 0;
   for (int ii = 0; ii < Ncon; ++ii) {
      if (ii < MAX_MCONTROLS_RANK)
         Cdiag(ii, ii) = Csvd.GetSig()[ii];
      else
         break;
   }
   TMatrixD Utranspose(Nmon, Nmon);
   Utranspose.Transpose(Csvd.GetU());
   fMcontrols.Mult(Csvd.GetV(), Cdiag *= Utranspose);
#endif

# if PRINT_MON2MD_EIGENVECTORS

   std::vector<std::string> mdsignals;
   for (int md = 1; md <= 8; ++md) {
      std::stringstream mdstr;
      mdstr << "qwk_md" << md << "pos";
      mdsignals.push_back(mdstr.str());
      mdstr.str("");
      mdstr << "qwk_md" << md << "neg";
      mdsignals.push_back(mdstr.str());
   }
   int Nmdsig = mdsignals.size();
   TMatrixD mdgains(Nmdsig, Ncon);
   for (int md = 0; md < Nmdsig; ++md) {
      for (int pat = 11; pat < 16; ++pat) {
         int con = pat - 11;
         FindFitCoeff(mdsignals[md], 'r', pat,
                      pattern_counter[pat], &mdgains(md, con));
         FindFitCoeff(mdsignals[md], 'i', pat,
                      pattern_counter[pat], &mdgains(md, con + 5));
      }
   }
   TMatrixD mon2md(Nmdsig, Nmon);
   mon2md.Mult(mdgains, fMcontrols);
   TMatrixD mon2mdTranspose(Nmon, Nmdsig);
   mon2mdTranspose.Transpose(mon2md);
   TDecompSVD Gsvd(mon2mdTranspose, 1e-14);
   std::cout << "EvaluateControls at pattern_counters";
   for (int pat=11; pat < 16; ++pat)
      std::cout << " " << pattern_counter[pat];
   std::cout << std::endl;
   std::cout << "eigenvectors of mon2md on the monitors side" << std::endl;
   std::cout << " eigenvalue    ";
   for (int n = 0; n < Gsvd.GetSig().GetNrows(); ++n) {
      if (Gsvd.GetSig()[n] < 1e-14)
         break;
      char s[99];
      sprintf(s, " %10.4g", Gsvd.GetSig()[n]);
      std::cout << s;
   }
   std::cout << std::endl;
   for (int mon = 0; mon < (int)fMonnames.size(); ++mon) {
      char s[99];
      sprintf(s, "%15s", fMonnames[mon].c_str());
      std::cout << s;
      for (int n = 0; n < Gsvd.GetSig().GetNrows(); ++n) {
         if (Gsvd.GetSig()[n] < 1e-14)
            break;
         sprintf(s, " %10.4g", Gsvd.GetU()(mon, n));
         std::cout << s;
      }
      std::cout << std::endl;
   }
   std::cout << "eigenvectors of mon2md on the md side" << std::endl;
   std::cout << " eigenvalue    ";
   for (int n = 0; n < Gsvd.GetSig().GetNrows(); ++n) {
      if (Gsvd.GetSig()[n] < 1e-14)
         break;
      char s[99];
      sprintf(s, " %10.4g", Gsvd.GetSig()[n]);
      std::cout << s;
   }
   std::cout << std::endl;
   for (int md = 0; md < Nmdsig; ++md) {
      char s[99];
      sprintf(s, "%15s", mdsignals[md].c_str());
      std::cout << s;
      for (int n = 0; n < Gsvd.GetSig().GetNrows(); ++n) {
         if (Gsvd.GetSig()[n] < 1e-14)
            break;
         sprintf(s, " %10.4g", Gsvd.GetV()(md, n));
         std::cout << s;
      }
      std::cout << std::endl;
   }
   std::cout << std::endl;
#endif

   return 0;
}

Int_t Calib_mps::GetControls()
{
   // Apply the transformation to the monitors vector and use it to
   // estimate the controls vector, storing it in fVcontrols. This
   // method is called on every event so it is optimized for minimum
   // overhead. The EvaluateControls method is assumed to have been
   // invoked already, but this is not checked.

   int Ncon = fMcontrols.GetNrows();
   int Nmon = fMonvalues.size();
   fVcontrols.ResizeTo(Nmon);
   for (int mon = 0; mon < Nmon; ++mon) {
      fVcontrols[mon] = *fMonvalues[mon] - fVoffsets[mon];
   }
   fVcontrols *= fMcontrols;
   return 0;
}

void Calib_mps::CheckControls()
{
   // Loop over all events in the input tree and apply the controls
   // transformation to the monitor signals, generating monitor-driven
   // estimates for the control inputs. A new tree "controls" is saved
   // with the event-by-event values for the controls vector. At the
   // same time, the components of the controls vector are analyzed
   // in a similar manner to what is done to the original signals in 
   // Loop(), and the fit coefficients to the sin + cos + constant terms
   // are saved in output histograms indexed by modulation interval.
   // Results are written to the same output file used by Loop().
 
   if (fChain == 0)
      return;

   if (fMnormal.GetNrows() != (int)fMonnames.size())
      NormalizeMonitors();

   TFile fout(fWorkfile.c_str(), "update");

   int Ncon = 10;
   int period = 0;
   double controls[Ncon];
   TTree contree("controls", "controls vector");
   contree.Branch("controls", controls, "controls[10]/D");
   contree.Branch("period", &period, "period/I");

   int Nperiods = 0;
   for (int pat = 11; pat < 16; ++pat) {
      TH1D *h = lookup_histo("qwk_charge", pat, 'r');
      Nperiods += h->GetNbinsX();
   }

   std::vector<TH1D*> control_r;
   std::vector<TH1D*> control_i;
   std::vector<TH1D*> control_c;
   for (int con = 0; con < Ncon; ++con) {
      std::stringstream sname;
      sname << "control_" << con;
      const char *name = sname.str().c_str();
      TH1D *hr = make_histo("control component", name, 0, 'r', Nperiods);
      TH1D *hi = make_histo("control component", name, 0, 'i', Nperiods);
      TH1D *hc = make_histo("control component", name, 0, 'c', Nperiods);
      control_r.push_back(hr);
      control_i.push_back(hi);
      control_c.push_back(hc);
   }

   // pass through the input tree
 
   double X00, X10, X01, X11, X20, X02;
   X00 = X10 = X01 = X11 = X20 = X02 = 0;
   double Y00[Ncon], Y10[Ncon], Y01[Ncon];
   double Z00[Ncon];
   int pcount = 0;
   int pattern_number = 0;
   int pattern_counter[20] = {0};
   Long64_t nentries = fChain->GetEntries();
   for (Long64_t jentry=0; jentry <= nentries; jentry++) {
      if (jentry < nentries)
         fChain->GetEntry(jentry);
      else
         bm_pattern_number = 0;

      // at the end of a pattern, compute the fit parameters
      // and store them and their errors in histograms
 
      if (bm_pattern_number != pattern_number) {
         if (pcount > MIN_CALIBRATION_CYCLE) {
            TMatrixDSym M(3);
            M(0,0) = X20;
            M(0,1) = M(1,0) = X11;
            M(0,2) = M(2,0) = X10;
            M(1,1) = X02;
            M(1,2) = M(2,1) = X01;
            M(2,2) = X00;
            TMatrixDSym MI(M);
            double determ;
            MI.Invert(&determ);
            for (int con = 0; con < Ncon; ++con) {
               double B[3];
               B[0] = Y10[con];
               B[1] = Y01[con];
               B[2] = Y00[con];
               double areal = MI(0,0) * B[0] + MI(0,1) * B[1] + MI(0,2) * B[2];
               double aimag = MI(1,0) * B[0] + MI(1,1) * B[1] + MI(1,2) * B[2];
               double acons = MI(2,0) * B[0] + MI(2,1) * B[1] + MI(2,2) * B[2];
               double chisqr = Z00[con]
                             - 2 * (areal * B[0] + aimag * B[1] + acons * B[2])
                             + areal * areal * M(0,0) + aimag * aimag * M(1,1)
                             + acons * acons * M(2,2)
                             + 2 * (areal * aimag * M(0,1) 
                                    + areal * acons * M(0,2)
                                    + aimag * acons * M(1,2));
               double varresid = chisqr / (X00 - 3);
               double areal_err = sqrt(MI(0,0) * varresid);
               double aimag_err = sqrt(MI(1,1) * varresid);
               double acons_err = sqrt(MI(2,2) * varresid);
               control_r[con]->SetBinContent(period, areal);
               control_r[con]->SetBinError(period, areal_err);
               control_i[con]->SetBinContent(period, aimag);
               control_i[con]->SetBinError(period, aimag_err);
               control_c[con]->SetBinContent(period, acons);
               control_c[con]->SetBinError(period, acons_err);
            }
         }
         ++period;
         pcount = 0;
         pattern_number = bm_pattern_number;
         ++pattern_counter[pattern_number];
         if (pattern_number != 0)
            EvaluateControls(pattern_number, pattern_counter);
         X00 = X10 = X01 = X11 = X20 = X02 = 0;
         for (int con = 0; con < Ncon; ++con) {
            Y00[con] = 0;
            Y01[con] = 0;
            Y10[con] = 0;
            Z00[con] = 0;
         }
      }
      if (jentry == nentries)
         break;

      // during modulation cycles, accummulate statistics on each control

      if (ramp_good) {
         double costheta = cos(ramp_filled * M_PI / 180);
         double sintheta = sin(ramp_filled * M_PI / 180);
         X00 += 1;
         X10 += costheta;
         X01 += sintheta;
         X11 += costheta * sintheta;
         X20 += costheta * costheta;
         X02 += sintheta * sintheta;
         GetControls();
         for (int con = 0; con < Ncon; ++con) {
            double y = fVcontrols[con];
            Y00[con] += y;
            Y10[con] += y * costheta;
            Y01[con] += y * sintheta;
            Z00[con] += y * y;
            controls[con] = y;
         }
         contree.Fill();
         ++pcount;
      }
   }
 
   fout.cd();
   contree.Write();
   for (int con = 0; con < Ncon; ++con) {
      control_r[con]->Write();
      control_i[con]->Write();
      control_c[con]->Write();
   }
   fout.Close();

   // If anyone has this workfile open, close and reopen it
   // for them so they can see the updates. When I reopen it
   // I do so in readonly mode because it is an i/o collision
   // if it was open for writing when this method has opened
   // it independently for writing using a private descriptor.

   TFile *fin = (TFile*)gROOT->FindObject(fWorkfile.c_str());
   if (fin && fin->IsOpen()) {
      fin->Close();
      fin->Open(fWorkfile.c_str());
   }
}

TH2D *Calib_mps::GetCovariances()
{
   // Extract the covariance matrix of the monitor residuals, or the
   // nomitors if the normalization has already been applied, and
   // return the results in a 2D histogram. A copy of the histogram
   // is automatically written to the workfile.
 
   int Nmon = fMonnames.size();
   if (Nmon == 0)
      return 0;
   std::string histo_name, histo_title;
   if (fMnormal.GetNrows() == (int)fMonnames.size()) {
      histo_name = "nom_covars";
      histo_title = "covariance matrix of the nomitor signals";
   }
   else {
      histo_name = "mon_covars";
      histo_title = "covariance matrix of the monitor signals";
   }

   CorrelateMonitors();
   TH2D *covar = new TH2D(histo_name.c_str(), histo_title.c_str(),
                          Nmon, 0, Nmon, Nmon, 0, Nmon); 
   for (int i = 0; i < Nmon; ++i) {
      for (int j = 0; j < Nmon; ++j) {
         covar->SetBinContent(i + 1, j + 1, fMcovar(i, j));
      }
   }
   covar->SetStats(0);

   TFile fout(fWorkfile.c_str(), "update");
   covar->Write();
   fout.Close();

   TFile *fin = (TFile*)gROOT->FindObject(fWorkfile.c_str());
   if (fin && fin->IsOpen()) {
      fin->Close();
      fin->Open(fWorkfile.c_str());
   }
   return (TH2D*)gROOT->FindObject(histo_name.c_str());
}

int Calib_mps::FindFitCoeff(std::string signame, char comp,
                            int pattern, int pattern_counter,
                            double *value, double *error)
{
   // Look up the coefficient of the cosine term, the sine term, and the
   // constant term for the named signal at the location in the tree
   // specified in pattern and pattern_counter.
   //  * signame = one of the names of a signal, eg. an input tree var
   //  * comp = 'r' for cos term, 'i' for sin term, 'c' for const term
   //  * pattern = 11, 12, 13, 14, or 15
   //  * pattern_counter = cycle number for this pattern type, starting
   //      from 1 and counting up as one moves down the rows of the tree.
   // If the stated modulation period exists and has a valid fit then
   // a return value 1 is returned and the fit coefficient and its error
   // are returned through the value and error argument references. If
   // the stated modulation period does not exist, or if the fit did not
   // succeed (eg. too few data points, large dc shifts in the middle,
   // etc.) then the value 0 is returned and the value,error arguments
   // are unchanged.

   TH1D *hist = lookup_histo(signame, pattern, comp);

#if INDIVIDUAL_CYCLE_GAINS_EVALUATION

   int pcnt = pattern_counter;
   int max_pcnt = hist->GetNbinsX();
   if (pcnt < 0 || pcnt > max_pcnt)
      return 0;
   else if (pcnt == 0 && max_pcnt > 0)
      pcnt = 1;
   double val = hist->GetBinContent(pcnt);
   double err = hist->GetBinError(pcnt);
   if (err == 0)
      return 0;
   else if (error != 0)
      *error = err;
   *value = hist->GetBinContent(pcnt);

#else 

   // compute just one global gains matrix for all cycles

   double S0, S1, S2;
   S0 = S1 = S2 = 0;
   int nbins = hist->GetNbinsX();
   for (int bin = 1; bin <= nbins; ++bin) {
      double err = hist->GetBinError(bin);
      if (err == 0)
         continue;
      double wgt = 1 / (err * err);
      S0 += wgt;
      double val = hist->GetBinContent(bin);
      S1 += wgt * val;
      S2 += wgt * val * val;
   }
   if (S0 == 0)
      return 0;
   double sigma2 = 1 / S0;
   double mean = S1 * sigma2;
   double chisqr = S2 - S1 * mean;
   int ndof = nbins - 1;
   if (error)
      *error = sqrt(sigma2);
   *value = mean;

#endif

   return 1;
}

void Calib_mps::ApplyCorrections()
{
   // Loop over all events in the input tree and apply the controls
   // transformation to the monitor signals, then use the computed
   // controls vector to compute a correction for each signal in the
   // input tree based upon the precomputed sensitivities. A new tree
   // "corrected" is saved containing a copy of the input tree, with
   // signal values overwritten with their corrected values. Results
   // are written to the same output file used by Loop().
 
   if (fChain == 0)
      return;

   if (fMnormal.GetNrows() != (int)fMonnames.size())
      NormalizeMonitors();

   TFile fout(fWorkfile.c_str(), "update");

   std::string title(fChain->GetTitle());
   title += " corrected by Calib_mps";
   TTree corrected("corrected", title.c_str());
   std::map<std::string, Double_t*>::iterator iter;
   for (iter = fSignames.begin(); iter != fSignames.end(); ++iter) {
      std::string bname(iter->first);
      corrected.Branch(bname.c_str(), iter->second);
   }
   for (iter = fAuxnames.begin(); iter != fAuxnames.end(); ++iter) {
      std::string bname(iter->first);
      corrected.Branch(bname.c_str(), iter->second);
   }

   // create vectors to hold the sensitivities for each signal,
   // updated during each modulation period in the event stream
 
   int Ncon = 10;
   std::map<std::string, TVectorD> gains;
   for (iter = fSignames.begin(); iter != fSignames.end(); ++iter) {
      gains[iter->first].ResizeTo(Ncon);
   }

   // pass through the input tree
 
   int pattern_number = 0;
   int pattern_counter[20] = {0};
   Long64_t nentries = fChain->GetEntries();
   for (Long64_t jentry=0; jentry < nentries; jentry++) {
      fChain->GetEntry(jentry);
      if (bm_pattern_number != pattern_number) {
         pattern_number = bm_pattern_number;
         ++pattern_counter[pattern_number];
         EvaluateControls(pattern_number, pattern_counter);
         for (iter = fSignames.begin(); iter != fSignames.end(); ++iter) {
            for (int pat = 11; pat < 16; ++pat) {
               int con = pat - 11;
               FindFitCoeff(iter->first, 'r', pat,
                            pattern_counter[pat], &gains[iter->first][con]);
               FindFitCoeff(iter->first, 'i', pat,
                            pattern_counter[pat], &gains[iter->first][con + 5]);
            }
         }
      }
      if (ramp_good) {
         GetControls();
         std::map<std::string, Double_t*>::iterator iter;
         for (iter = fSignames.begin(); iter != fSignames.end(); ++iter) {
            double correction = 0;
            for (int con = 0; con < Ncon; ++con) {
               correction += gains[iter->first][con] * fVcontrols[con];
            }
            *iter->second -= correction;
         }
         corrected.Fill();
      }
   }
   fout.cd();
   corrected.Write();
   fout.Close();

   // If anyone has this workfile open, close and reopen it
   // for them so they can see the updates. When I reopen it
   // I do so in readonly mode because it is an i/o collision
   // if it was open for writing when this method has opened
   // it independently for writing using a private descriptor.

   TFile *fin = (TFile*)gROOT->FindObject(fWorkfile.c_str());
   if (fin && fin->IsOpen()) {
      fin->Close();
      fin->Open(fWorkfile.c_str());
   }
}
