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

      for(int ipho = 0; ipho < nPho; ++ipho){
	 if((*phoEt)[ipho] < 15) continue;
	 int mc_truth = MCTruthMatch(ipho);
	 ////for signal prompt photon when reconstructed photon match to gen photon
	 if(mc_truth != 1) continue;


	 // // Fill histograms

	 bool isBarrel = fabs( (*phoSCEta)[ipho] ) < 1.479 ? true : false;
	 if( isBarrel && (*phoEt)[ipho] > 15.0 && getID(ipho)) {
	    h_Den->Fill( (*phoEt)[ipho] );

	    if( getIso(ipho)){
	       h_Num->Fill( (*phoEt)[ipho] );
	    }
	 }
      } // end loop over for photons
   }

   h_Eff = (TH1F*) h_Num->Clone("h_Eff");
   h_Eff->Divide(h_Den);

   f1->Write();
   f1->Close();

}

double getPhotonIsoEfficency::dR(double eta1, double phi1, double eta2, double phi2){
   double dphi = phi1 - phi2;
   double deta = eta1 - eta2;
   return sqrt(dphi*dphi + deta*deta);
}

int getPhotonIsoEfficency::MCTruthMatch(int jpho){
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


bool getPhotonIsoEfficency::getID(int jpho){

   bool pass = false;
   if (fabs((*phoSCEta)[jpho]) < 1.4442) {
      if((*phoSigmaIEtaIEtaFull5x5)[jpho] < 0.0102  &&
	    (*phoHoverE)[jpho] < 0.05
	) pass = true;
   } else {
      if((*phoSigmaIEtaIEtaFull5x5)[jpho] < 0.0268 &&
	    (*phoHoverE)[jpho] < 0.05
	) pass = true;
   }
   return pass;
}

bool getPhotonIsoEfficency::getIso(int jpho){
   bool pass = false;

   if( fabs((*phoSCEta)[jpho]) < 1.4442){ // EB
      pass = ( 
	    ((*phoPFChIso)[jpho]               <  1.37   )  &&
	    (TMath::Max(((*phoPFNeuIso)[jpho] - rho*EAneutral((*phoSCEta)[jpho])), 0.0) < (1.06 + (0.014 * (*phoEt)[jpho]) + (0.000019 * pow((*phoEt)[jpho], 2.0))) )  &&
	    (TMath::Max(((*phoPFPhoIso)[jpho] - rho*EAphoton( (*phoSCEta)[jpho])), 0.0) < (0.28 + (0.0053 * (*phoEt)[jpho])) ) ); 
   }
   if( fabs((*phoSCEta)[jpho]) > 1.4442){ // EE
      pass = ( 
	    ((*phoPFChIso)[jpho]               <  1.10   ) &&
	    (TMath::Max(((*phoPFNeuIso)[jpho] - rho*EAneutral((*phoSCEta)[jpho])), 0.0) < (2.69 + (0.0139 * (*phoEt)[jpho]) + (0.000025 * pow((*phoEt)[jpho], 2.0))) )  &&
	    (TMath::Max(((*phoPFPhoIso)[jpho] - rho*EAphoton( (*phoSCEta)[jpho])), 0.0) < (0.39 + (0.0034 * (*phoEt)[jpho])) ) ); 
   }   

   return pass;

}



