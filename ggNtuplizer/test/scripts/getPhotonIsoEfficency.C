#define getPhotonIsoEfficency_cxx
#include "getPhotonIsoEfficency.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void getPhotonIsoEfficency::Loop()
{
   if (fChain == 0) return;

   TFile* f1 = new TFile("hist_PhotonIsoEfficiency.root", "recreate");
   TH1F* h_Num = new TH1F("h_Num", "photon pT", 100, 0, 200);
   TH1F* h_Den = new TH1F("h_Den", "photon pT", 100, 0, 200);
   TH1F* h_Eff = new TH1F("h_Eff", "photon pT", 100, 0, 200);

   
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
   
   
   }
}
