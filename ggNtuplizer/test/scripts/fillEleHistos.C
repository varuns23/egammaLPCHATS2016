#define fillEleHistos_cxx
#include "fillEleHistos.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

void fillEleHistos::Loop()
{
   if (fChain == 0) return;
   TFile* f1 = new TFile("hist_Electrons_TT.root", "recreate");
   TH1F* hPt = new TH1F("hPt", "electron pT", 100, 0, 200);
   TH1F* hEta = new TH1F("hEta", "electron eta", 100,-3, 3);
   TH1F* hSigmaIetaIeta_barrel = new TH1F("hSigmaIetaIeta_barrel", "", 100, 0, 0.05);
   TH1F* hSigmaIetaIeta_endcap = new TH1F("hSigmaIetaIeta_endcap", "", 100, 0, 0.05);


   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   
      for(int iele = 0; iele < nEle; ++iele){
	 if((*elePt)[iele] < 15) continue;
	 int mc_truth = MCTruthMatch(iele);
	 ////for signal prompt electrons when reconstructed electrons match to gen electrons
	 if(mc_truth != 1) continue;
	 //    //for bkgs, 
	 //if( !(mc_truth==2 || mc_truth==3) ) continue;
	 // // Fill histograms
	 hPt->Fill( (*elePt)[iele] );
	 hEta->Fill( (*eleSCEta)[iele] );
	 bool isBarrel = fabs( (*eleSCEta)[iele] ) < 1.479 ? true : false;
	 if( isBarrel ) {
	    hSigmaIetaIeta_barrel->Fill((*eleSigmaIEtaIEtaFull5x5)[iele]);
	 }else{
	    hSigmaIetaIeta_endcap->Fill((*eleSigmaIEtaIEtaFull5x5)[iele]);
	 }
      } // end loop over for electrons
   }
   f1->Write();
   f1->Close();
}

double fillEleHistos::dR(double eta1, double phi1, double eta2, double phi2){
   double dphi = phi1 - phi2;
   double deta = eta1 - eta2;
   return sqrt(dphi*dphi + deta*deta);
}

int fillEleHistos::MCTruthMatch(int jele){
   int eleInd = -1;

   for(int imc = 0; imc < nMC; ++imc){
      if(fabs((*mcPID)[imc]) != 11) continue;
      if((*mcPt)[imc] < 10) continue;
      bool match_gen = dR((*mcEta)[imc], (*mcPhi)[imc], (*eleSCEta)[jele], (*eleSCPhi)[jele]) < 0.05;
      if(match_gen && eleInd < 0) eleInd = imc;
   }

   if(eleInd >= 0){
      if(((*mcParentage)[eleInd]& 4)==0) return 1;
      else
	 return 2;
   } else {
      return 3;
   }
}


