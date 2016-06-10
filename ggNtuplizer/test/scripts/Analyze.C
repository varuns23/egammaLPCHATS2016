#define Analyze_cxx
#include "Analyze.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#include<iostream> // for std // cout etc..
#include <TMath.h>  // for abs, max

using namespace std;

float zero = 0.0;

void Analyze::Loop()
{
//

   TFile *MyFile = new TFile("ELE_FakeRate.root","RECREATE");
   MyFile->cd();
   
   TH1F *h_EEMass50 = new TH1F("E-E Invariant Mass50","E-E Invariant Mass50",100,60,120);
   TH1F *h_EEMass5060 = new TH1F("E-E Invariant Mass5060","E-E Invariant Mass5060",100,60,120);
   TH1F *h_EEMass6080 = new TH1F("E-E Invariant Mass6080","E-E Invariant Mass6080",100,60,120);
   TH1F *h_EEMass80100 = new TH1F("E-E Invariant Mass80100","E-E Invariant Mass80100",100,60,120);
   TH1F *h_EEMass100 = new TH1F("E-E Invariant Mass100","E-E Invariant Mass100",100,60,120);
   TH1F *h_EPMass50 = new TH1F("E-P Invariant Mass50","E-P Invariant Mass50",100,60,120);
   TH1F *h_EPMass5060 = new TH1F("E-P Invariant Mass5060","E-P Invariant Mass5060",100,60,120);
   TH1F *h_EPMass6080 = new TH1F("E-P Invariant Mass6080","E-P Invariant Mass6080",100,60,120);
   TH1F *h_EPMass80100 = new TH1F("E-P Invariant Mass80100","E-P Invariant Mass80100",100,60,120);
   TH1F *h_EPMass100 = new TH1F("E-P Invariant Mass100","E-P Invariant Mass100",100,60,120);
   
   //   TH1F *h_passPPt = new TH1F("Pass probe Pt","Pass probe Pt",80,100,420);
   //   TH1F *h_failPPt = new TH1F("Fail probe Pt","Fail probe Pt",80,100,420);
   //   TH1F *h_EoP = new TH1F("ElectronEoP","ElectronEoP",100,0,1.9);
   //   TH1F *h_counter = new TH1F("counter","counter",10,0,10);
   
   int tagcount = 0;
   int eprobecount = 0;
   int pprobecount = 0;
   int nevents = 0;

   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntries();
   Long64_t nbytes = 0, nb = 0;
   std::cout<< " Will Analyze = " << nentries << " events" <<std::endl;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     Long64_t ientry = LoadTree(jentry);
     if (ientry < 0) break;
     nb = fChain->GetEntry(jentry);   nbytes += nb;
     nevents++;

     if (!( HLTEleMuX>>6&1)) continue;

     for ( int j=0; j<nEle; j++ ){
       TLorentzVector etag;
       if ((TightElectron(j)  && (elePt->at(j) > 30))){
       	   tagcount++;
	   etag.SetPtEtaPhiE(elePt->at(j), eleEta->at(j), elePhi->at(j), eleEn->at(j) );
	   
       }
       for( int i=0; i<nPho; i++ ){
	 float dR = dRCalc(phoSCEta->at(i),phoSCPhi->at(i), eleSCEta->at(j),eleSCPhi->at(j));	
	 if (dR < 0.5) continue;
	 if ( MediumPhoton(i) && Collision(i) && (phohasPixelSeed->at(i)) ) {
	   if ((phoEt->at(i) > 10) && (phoEt->at(i) < 50)){
	     TLorentzVector elecpr50;	
	     elecpr50.SetPtEtaPhiE(phoEt->at(i), phoEta->at(i), phoPhi->at(i), phoE->at(i) );
	     double InvMass50  = ( etag + elecpr50).M();
	     if ((InvMass50 > 60) && (InvMass50 < 120)){
	       h_EEMass50->Fill(InvMass50);
	     }
	   }
	   if ((phoEt->at(i) > 50) && (phoEt->at(i) < 60)){
	     TLorentzVector elecpr5060;	
	     elecpr5060.SetPtEtaPhiE(phoEt->at(i), phoEta->at(i), phoPhi->at(i), phoE->at(i) );
	     double InvMass5060  = ( etag + elecpr5060).M();
	     if ((InvMass5060 > 60) && (InvMass5060 < 120)){
	       h_EEMass5060->Fill(InvMass5060);
	     }
	   }
	   if ((phoEt->at(i) > 60) && (phoEt->at(i) < 80)){
	     TLorentzVector elecpr6080;	
	     elecpr6080.SetPtEtaPhiE(phoEt->at(i), phoEta->at(i), phoPhi->at(i), phoE->at(i) );
	     double InvMass6080  = ( etag + elecpr6080).M();
	     if ((InvMass6080 > 60) && (InvMass6080 < 120)){
	       h_EEMass6080->Fill(InvMass6080);
	     }
	   }
	   if ((phoEt->at(i) > 80) && (phoEt->at(i) < 100)){
	     TLorentzVector elecpr80100;	
	     elecpr80100.SetPtEtaPhiE(phoEt->at(i), phoEta->at(i), phoPhi->at(i), phoE->at(i) );
	     double InvMass80100  = ( etag + elecpr80100).M();
	     if ((InvMass80100 > 60) && (InvMass80100 < 120)){
	       h_EEMass80100->Fill(InvMass80100);
	     }
	   }
	   if ((phoEt->at(i) > 100)){
	     TLorentzVector elecpr100;	
	     elecpr100.SetPtEtaPhiE(phoEt->at(i), phoEta->at(i), phoPhi->at(i), phoE->at(i) );
	     double InvMass100  = ( etag + elecpr100).M();
	     if ((InvMass100 > 60) && (InvMass100 < 120)){
	       h_EEMass100->Fill(InvMass100);
	     }
	   }
	 }//Electron probe
	 
	 if (MediumPhoton(i) && Collision(i) && (phohasPixelSeed->at(i) == 0)) {
	   if ((phoEt->at(i) > 10) && (phoEt->at(i) < 50) ){
	     TLorentzVector phopr50;
	     phopr50.SetPtEtaPhiE(phoEt->at(i), phoEta->at(i), phoPhi->at(i), phoE->at(i) );
	     double InvMass50  = ( etag + phopr50 ).M();
	     if ((InvMass50 > 60) && (InvMass50 < 120)){
	       h_EPMass50->Fill(InvMass50);
	     }
	   }
	   if ((phoEt->at(i) > 50) && (phoEt->at(i) < 60) ){
	     TLorentzVector phopr5060;
	     phopr5060.SetPtEtaPhiE(phoEt->at(i), phoEta->at(i), phoPhi->at(i), phoE->at(i) );
	     double InvMass5060  = ( etag + phopr5060 ).M();
	     if ((InvMass5060 > 60) && (InvMass5060 < 120)){
	       h_EPMass5060->Fill(InvMass5060);
	     }
	   }
	   if ((phoEt->at(i) > 60) && (phoEt->at(i) < 80) ){
	     TLorentzVector phopr6080;
	     phopr6080.SetPtEtaPhiE(phoEt->at(i), phoEta->at(i), phoPhi->at(i), phoE->at(i) );
	     double InvMass6080  = ( etag + phopr6080 ).M();
	     if ((InvMass6080 > 60) && (InvMass6080 < 120)){
	       h_EPMass6080->Fill(InvMass6080);
	     }
	   }
	   if ((phoEt->at(i) > 80) && (phoEt->at(i) < 100) ){
	     TLorentzVector phopr80100;
	     phopr80100.SetPtEtaPhiE(phoEt->at(i), phoEta->at(i), phoPhi->at(i), phoE->at(i) );
	     double InvMass80100  = ( etag + phopr80100 ).M();
	     if ((InvMass80100 > 60) && (InvMass80100 < 120)){
	       h_EPMass80100->Fill(InvMass80100);
	     }
	   }
	   if ((phoEt->at(i) > 100) ){
	     TLorentzVector phopr100;
	     phopr100.SetPtEtaPhiE(phoEt->at(i), phoEta->at(i), phoPhi->at(i), phoE->at(i) );
	     double InvMass100  = ( etag + phopr100 ).M();
	     if ((InvMass100 > 60) && (InvMass100 < 120)){
	       h_EPMass100->Fill(InvMass100);
	     }
	   }
	 }//Photon probe
       }//for( int i=0; i<nPho; i++ ){
     }	 //for ( int i=0; i<nEle; i++ ){
   }//for (Long64_t jentry=0; jentry<nentries;jentry++)
   
     cout<<"Total events : "<<nevents<<endl;  
   //   cout<<"total tag electrons  :"<<tagcount<<" total  probe electrons   :"<<eprobecount<<" total probe photons     :"<<pprobecount<<endl;
  
   MyFile->Write();
   MyFile->Close();
}

bool Analyze::TightElectron(int& i){
  
  //  float pfChIso      = TMath::Max( phoPFChIso->at(i)  - rho * EaChargedHadrons(phoSCEta->at(i)), zero );
  //  float pfNeuIso     = TMath::Max( phoPFNeuIso->at(i) - rho * EaNeutralHadrons(phoSCEta->at(i)), zero );
  //  float pfPhoIso     = TMath::Max( phoPFPhoIso->at(i) - rho * EaPhotons       (phoSCEta->at(i)), zero );
  //  double eleooEmooP = (1/) ;
  double neu = elePFPhoIso->at(i) + elePFNeuIso->at(i) - effArea(eleSCEta->at(i)) * rho;
  double combinedRelIso = (elePFChIso->at(i) + max(neu, 0.)) / (elePt->at(i));

  
  bool isTightElectron = false;
  if (fabs(eleSCEta->at(i)) < 1.479){
    isTightElectron  = 
      eleHoverE->at(i)                <   0.0597 && 
      eleSigmaIEtaIEtaFull5x5->at(i)  <   0.0101 && 
      eledEtaAtVtx->at(i)             <   0.00926 &&		                         
      eledPhiAtVtx->at(i)             <   0.0336  &&  
      fabs(eleD0->at(i))               <   0.0111  &&   
      fabs(eleDz->at(i))               <   0.0466  &&
      eleMissHits->at(i)              <=  2 &&
      eleConvVeto->at(i)               == 1 &&
      eleEoverPInv->at(i)                <   0.012 &&
      combinedRelIso                  <   0.0354;
  }
  
  
  return isTightElectron;
}
float Analyze::effArea(float eta)  {

  float ea = 0.0;                                
                                                 
  if           (fabs(eta) < 1.0)                        ea = 0.1752;        
  else if      (fabs(eta) > 1.0   && fabs(eta) < 1.479) ea = 0.1862;        
  else if      (fabs(eta) > 1.479 && fabs(eta) < 2.0)   ea = 0.1411;        
  else if      (fabs(eta) > 2.0   && fabs(eta) < 2.2)   ea = 0.1534;        
  else if      (fabs(eta) > 2.2   && fabs(eta) < 2.3)   ea = 0.1903;        
  else if      (fabs(eta) > 2.4   && fabs(eta) < 2.4)   ea = 0.2243;         
  else if      (                     fabs(eta) > 2.4)   ea = 0.2687;         
                                                 
  return ea;                                     


}


   
bool Analyze::MediumPhoton(int& i){
     
  float pfChIso      = phoPFChIso->at(i);
  float pfNeuIso     = TMath::Max( phoPFNeuIso->at(i) - rho * EaNeutralHadrons(phoSCEta->at(i)), zero );
  float pfPhoIso     = TMath::Max( phoPFPhoIso->at(i) - rho * EaPhotons       (phoSCEta->at(i)), zero );
  
  bool isMediumPhoton = false;
  if (fabs(phoSCEta->at(i)) < 1.4442){
    isMediumPhoton  =  
      phoHoverE->at(i)                <   0.05 && 
      phoSigmaIEtaIEtaFull5x5->at(i)  <   0.0102 && 
      //		       phohasPixelSeed->at (i)  ==  0 &&
      pfChIso                         <   1.37 &&
      pfNeuIso                        <   1.06 + 0.014 * phoEt->at(i) + 0.000019 * (phoEt->at(i) * phoEt->at(i)) &&
      pfPhoIso                        <   0.28 + 0.0053 * phoEt->at(i); 
  }
  
  return isMediumPhoton;
}

float Analyze::EaNeutralHadrons(float eta)  
{                                                
  float ea = 0.0;                                
                                                 
  if           (fabs(eta) < 1.0)                        ea = 0.0599;        
  else if      (fabs(eta) > 1.0   && fabs(eta) < 1.479) ea = 0.0819;        
  else if      (fabs(eta) > 1.479 && fabs(eta) < 2.0)   ea = 0.0696;        
  else if      (fabs(eta) > 2.0   && fabs(eta) < 2.2)   ea = 0.0360;        
  else if      (fabs(eta) > 2.2   && fabs(eta) < 2.3)   ea = 0.0360;        
  else if      (fabs(eta) > 2.4   && fabs(eta) < 2.4)   ea = 0.0462;         
  else if      (                     fabs(eta) > 2.4)   ea = 0.0656;             
                                          
  return ea;                                     
}                                                                                                 
float Analyze::EaPhotons(float eta)         
{                                                
  float ea = 0.0;                                
                                                 
  if           (fabs(eta) < 1.0)                        ea = 0.1271;        
  else if      (fabs(eta) > 1.0   && fabs(eta) < 1.479) ea = 0.1101;        
  else if      (fabs(eta) > 1.479 && fabs(eta) < 2.0)   ea = 0.0756;        
  else if      (fabs(eta) > 2.0   && fabs(eta) < 2.2)   ea = 0.1175;        
  else if      (fabs(eta) > 2.2   && fabs(eta) < 2.3)   ea = 0.1498;        
  else if      (fabs(eta) > 2.4   && fabs(eta) < 2.4)   ea = 0.1857;         
  else if      (                     fabs(eta) > 2.4)   ea = 0.2183;         
                                                 
  return ea;                                     
}                                                

float Analyze::dRCalc(float etaEle, float phiEle, float etaPho, float phiPho){
    
  float dphi = fabs(phiEle - phiPho);
  if (dphi > TMath::Pi()) dphi = TMath::Pi()*2. - dphi;
  float deta = fabs(etaEle - etaPho);
  float dR = sqrt(deta*deta + dphi*dphi);
  return dR;
}

bool Analyze::Collision(int i){
   bool isnonColl = false;  
   isnonColl = (  phoSigmaIEtaIEtaFull5x5->at(i) > 0.001 && 
   phoSigmaIPhiIPhiFull5x5->at(i) > 0.001);
   
   return isnonColl;
}
