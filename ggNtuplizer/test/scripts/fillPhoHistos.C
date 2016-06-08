#define fillPhoHistos_cxx
#include "fillPhoHistos.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void fillPhoHistos::Loop()
{

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   TFile* f1 = new TFile("hist_Photon.root", "recreate");
   TH1F* hPt = new TH1F("hPt", "photon pT", 100, 0, 200);
   TH1F* hEta = new TH1F("hEta", "photon eta", 100,-3, 3);
   TH1F* hSigmaIetaIeta_barrel = new TH1F("hSigmaIetaIeta_barrel", "", 100, 0, 0.05);
   TH1F* hSigmaIetaIeta_endcap = new TH1F("hSigmaIetaIeta_endcap", "", 100, 0, 0.05);


   Long64_t nbytes = 0, nb = 0;
   //for (Long64_t jentry=0; jentry<nentries;jentry++) {
   for (Long64_t jentry=0; jentry<10;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
  
      for(int ipho = 0; ipho < nPho; ++ipho){
	 if((*phoEt)[ipho] < 15) continue;
	 int mc_truth = MCTruthMatch(ipho);
	 ////for signal prompt photon when reconstructed photon match to gen photon
	 if(mc_truth != 1 ) continue;
	 //    //for bkgs, 
	 //if( !(mc_truth==2 || mc_truth==3 ) continue;
	 // // Fill histograms
	 hPt->Fill( (*phoEt)[ipho] );
	 hEta->Fill( (*phoSCEta)[ipho] );
	 bool isBarrel = fabs( (*phoSCEta)[ipho] ) < 1.479 ? true : false;
	 if( isBarrel ) {
	    hSigmaIetaIeta_barrel->Fill((*phoSigmaIEtaIEtaFull5x5)[ipho]);
	 }else{
	    hSigmaIetaIeta_endcap->Fill((*phoSigmaIEtaIEtaFull5x5)[ipho]);
	 }
      } // end loop over for photons
   }
   f1->Write();
   f1->Close();
}

double fillPhoHistos::dR(double eta1, double phi1, double eta2, double phi2){
   double dphi = phi1 - phi2;
   double deta = eta1 - eta2;
   return sqrt(dphi*dphi + deta*deta);
}

int fillPhoHistos::MCTruthMatch(int jpho){
   int phoInd = -1;

   for(int imc = 0; imc < nMC; ++imc){
      if( mcPID->at(imc) != 22) continue;
      if( mcPt->at(imc) < 20) continue;
      bool match_gen = dR((*mcEta)[imc], (*mcPhi)[imc], (*phoSCEta)[jpho], (*phoSCPhi)[jpho]) < 0.05;
      if(match_gen && phoInd < 0) phoInd = imc;
   }

   if(phoInd >= 0){
      if(((*mcParentage)[phoInd]&4)==0) return 1;
      else
	 return 2;
   } else {
      return 3;
   }

}


