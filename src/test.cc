#include "printme.hh"
#include "math.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TPad.h"

int main(int argc, char** argv)
{ 
  printme pm;
  pm.setname("Donald");
  for(int i=0; i<5; i++)
    pm.printname();
  TCanvas *c = new TCanvas("c","c",0,0,600,600);
  TH1F *h = new TH1F("h","h",100,-1.2,1.2);
  for(int i=0;i<1e6;i++)
    h->Fill(sin(i));
  h->Draw();
  gPad->Update();
  c->SaveAs("~/test.png");
  return 0;
}
