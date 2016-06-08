#define makeDiPhotonMass_cxx
#include "makeDiPhotonMass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void makeDiPhotonMass::Loop()
{
   if (fChain == 0) return;
   TH1F* hMass = new TH1F("hMass", "", 100, 0, 200);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      vector <int> iphotons;
      for(int ipho = 0; ipho < nPho; ++ipho){
	 if((*phoEt)[ipho] < 20) continue;

	 if(fabs((*phoSCEta)[ipho])>2.5) continue;
	 if(fabs((*phoSCEta)[ipho])>1.4442 && fabs((*phoSCEta)[ipho]) < 1.566) continue;

	 int mc_truth = MCTruthMatch(ipho);

	 if(phoCut(ipho)) iphotons.push_back(ipho);
      }
      if(iphotons.size() < 2) continue;

      TLorentzVector p1, p2;
      int i = iphotons[0];
      p1.SetPtEtaPhiM((*phoEt)[i],(*phoSCEta)[i],(*phoSCPhi)[i],0);
      int j = iphotons[1];
      p2.SetPtEtaPhiM((*phoEt)[j],(*phoSCEta)[j],(*phoSCPhi)[j],0);

      double mass = (p1 + p2).M();

      hMass->Fill(mass);
   }
   TFile* f1 = new TFile("hist_DiPhotonMass.root", "recreate");
   hMass->Write();
   f1->Write();
   f1->Close();

}

double makeDiPhotonMass::dR(double eta1, double phi1, double eta2, double phi2){
   double dphi = phi1 - phi2;
   double deta = eta1 - eta2;
   return sqrt(dphi*dphi + deta*deta);
}

int makeDiPhotonMass::MCTruthMatch(int jpho){
   int phoInd = -1;

   for(int imc = 0; imc < nMC; ++imc){
      if( mcPID->at(imc) != 22) continue;
      if( mcPt->at(imc) < 20) continue;
      if( !((*mcStatusFlag)[imc]>>1&1) ) continue;
      bool match_gen = dR((*mcEta)[imc], (*mcPhi)[imc], (*phoSCEta)[jpho], (*phoSCPhi)[jpho]) < 0.05;
      if(match_gen && phoInd < 0) phoInd = imc;
   }

   if(phoInd >= 0)
      return 1;
   else
      return 2;
}

bool makeDiPhotonMass::phoCut(int i){

   bool pass = false;
   if (fabs((*phoSCEta)[i]) < 1.4442) {
      if((*phoSigmaIEtaIEtaFull5x5)[i] < 0.0102  &&
	    (*phoHoverE)[i] < 0.05	    
	) pass = true;
   } else {
      if((*phoSigmaIEtaIEtaFull5x5)[i] < 0.0268 &&
	    (*phoHoverE)[i] < 0.05	    
	) pass = true;
   }
   return pass;
}


