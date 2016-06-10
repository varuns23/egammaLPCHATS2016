#define makeDiElectronMass_cxx                                                                                                                                                     
#include "makeDiElectronMass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
 
using namespace std;
 
void makeDiElectronMass::Loop()
{
   if (fChain == 0) return;
 
   TFile* f1 = new TFile("hist_DiElectronMass.root", "recreate");
 
   TH1F* hMass_reco = new TH1F("hMass_reco", "", 100, 0, 200);
   TH1F* hMass_prompt = new TH1F("hMass_prompt", "", 100, 0, 200);
   TH1F* hMass_cleaned = new TH1F("hMass_cleaned", "", 100, 0, 200);
 
   Long64_t nentries = fChain->GetEntriesFast();
 
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
 
      vector <int> reco_electrons;
      vector <int> prompt_electrons;
      vector <int> cleaned_electrons;
      for(int iele = 0; iele < nEle; ++iele){
         if((*elePt)[iele] < 15) continue;
 
         if(fabs((*eleSCEta)[iele])>2.5) continue;
         if(fabs((*eleSCEta)[iele])>1.4442 && fabs((*eleSCEta)[iele]) < 1.566) continue;
 
         reco_electrons.push_back(iele);
 
         int mc_truth = MCTruthMatch(iele);
         if(mc_truth != 1 ) continue;
         prompt_electrons.push_back(iele);
 
         if(eleCut(iele))
            cleaned_electrons.push_back(iele);
      }
      if(reco_electrons.size() >= 2) 
         hMass_reco->Fill( getMass(reco_electrons[0], reco_electrons[1]) );
       
      if(prompt_electrons.size() >= 2) 
         hMass_prompt->Fill( getMass(prompt_electrons[0], prompt_electrons[1]) );
       
      if(cleaned_electrons.size() >= 2) 
         hMass_cleaned->Fill( getMass(cleaned_electrons[0], cleaned_electrons[1]) );
       
   }   
       
   f1->Write();
   f1->Close();
       
       
}      
       
       
double makeDiElectronMass::dR(double eta1, double phi1, double eta2, double phi2){
   double dphi = phi1 - phi2;
   double deta = eta1 - eta2;
   return sqrt(dphi*dphi + deta*deta);
}      
       
int makeDiElectronMass::MCTruthMatch(int jele){
   int eleInd = -1;
       
   for(int imc = 0; imc < nMC; ++imc){
      if(fabs((*mcPID)[imc]) != 11) continue;
      if((*mcPt)[imc] < 10) continue;
      if(!((*mcStatusFlag)[imc]>>1&1)) continue;
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
bool makeDiElectronMass::eleCut(int i){
 
   bool pass = false;
   if (fabs((*eleSCEta)[i]) < 1.4442) {
      if((*eleSigmaIEtaIEtaFull5x5)[i] < 0.0114 &&  (*eleHoverE)[i] < 0.181 ) pass = true;
      } else {
      if((*eleSigmaIEtaIEtaFull5x5)[i] < 0.0352 &&  (*eleHoverE)[i] < 0.116 ) pass = true;
   }
   return pass;
}
 
 
double makeDiElectronMass::getMass(int i, int j){
      double mass = 0.0;
 
      TLorentzVector e1, e2;
      e1.SetPtEtaPhiM((*elePt)[i],(*eleEta)[i],(*elePhi)[i],0);
      e2.SetPtEtaPhiM((*elePt)[j],(*eleEta)[j],(*elePhi)[j],0);
 
      mass = (e1 + e2).M();
 
      return mass;
}
