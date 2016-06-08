#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "egammaLPCHATS2016/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;
using namespace reco;

// (local) variables associated with tree branches
Int_t          nEle_;
vector<int>    eleCharge_;
vector<float>  eleEn_;
vector<float>  eleSCEn_;
vector<float>  eleD0_;
vector<float>  eleDz_;
vector<float>  elePt_;
vector<float>  eleEta_;
vector<float>  elePhi_;
vector<float>  eleSCEta_;
vector<float>  eleSCPhi_;
vector<float>  eleSCRawEn_;
vector<float>  eleHoverE_;
vector<float>  eleEoverP_;
vector<float>  eleEoverPout_;
vector<float>  eleEoverPInv_;
vector<float>  eledEtaAtVtx_;
vector<float>  eledPhiAtVtx_;
vector<float>  eleSigmaIEtaIEtaFull5x5_;
vector<float>  eleSigmaIPhiIPhiFull5x5_;
vector<int>    eleConvVeto_;
vector<int>    eleMissHits_;
vector<float>  elePFChIso_;
vector<float>  elePFPhoIso_;
vector<float>  elePFNeuIso_;
vector<float>  elePFPUIso_;
vector<float>  eleIDMVANonTrg_;
vector<float>  eleIDMVATrg_;
vector<float>  eleR9Full5x5_;
vector<int>    eleEcalDrivenSeed_;

vector<UShort_t> eleIDbit_;


void ggNtuplizer::branchesElectrons(TTree* tree) {

  tree->Branch("nEle",                    &nEle_);
  tree->Branch("eleCharge",               &eleCharge_);
  tree->Branch("eleEn",                   &eleEn_);
  tree->Branch("eleSCEn",                 &eleSCEn_);
  tree->Branch("eleD0",                   &eleD0_);
  tree->Branch("eleDz",                   &eleDz_);
  tree->Branch("elePt",                   &elePt_);
  tree->Branch("eleEta",                  &eleEta_);
  tree->Branch("elePhi",                  &elePhi_);
  tree->Branch("eleSCEta",                &eleSCEta_);
  tree->Branch("eleSCPhi",                &eleSCPhi_);
  tree->Branch("eleSCRawEn",              &eleSCRawEn_);
  tree->Branch("eleHoverE",               &eleHoverE_);
  tree->Branch("eleEoverP",               &eleEoverP_);
  tree->Branch("eleEoverPout",            &eleEoverPout_);
  tree->Branch("eleEoverPInv",            &eleEoverPInv_);
  tree->Branch("eledEtaAtVtx",            &eledEtaAtVtx_);
  tree->Branch("eledPhiAtVtx",            &eledPhiAtVtx_);
  tree->Branch("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5_);
  tree->Branch("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5_);
  tree->Branch("eleConvVeto",             &eleConvVeto_);
  tree->Branch("eleMissHits",             &eleMissHits_);
  tree->Branch("elePFChIso",              &elePFChIso_);
  tree->Branch("elePFPhoIso",             &elePFPhoIso_);
  tree->Branch("elePFNeuIso",             &elePFNeuIso_);
  tree->Branch("elePFPUIso",              &elePFPUIso_);
  tree->Branch("eleIDMVANonTrg",          &eleIDMVANonTrg_);
  tree->Branch("eleIDMVATrg",             &eleIDMVATrg_);
  tree->Branch("eleR9Full5x5",                &eleR9Full5x5_);
  tree->Branch("eleEcalDrivenSeed",           &eleEcalDrivenSeed_);

  if (runeleIDVID_) tree->Branch("eleIDbit",  &eleIDbit_);
  
}

void ggNtuplizer::fillElectrons(const edm::Event &e, const edm::EventSetup &es, math::XYZPoint &pv) {
    
  // cleanup from previous execution
  eleCharge_                  .clear();
  eleEn_                      .clear();
  eleSCEn_                    .clear();
  eleD0_                      .clear();
  eleDz_                      .clear();
  elePt_                      .clear();
  eleEta_                     .clear();
  elePhi_                     .clear();
  eleSCEta_                   .clear();
  eleSCPhi_                   .clear();
  eleSCRawEn_                 .clear();
  eleHoverE_                  .clear();
  eleEoverP_                  .clear();
  eleEoverPout_               .clear();
  eleEoverPInv_               .clear();
  eledEtaAtVtx_               .clear();
  eledPhiAtVtx_               .clear();
  eleSigmaIEtaIEtaFull5x5_    .clear();
  eleSigmaIPhiIPhiFull5x5_    .clear();
  eleConvVeto_                .clear();
  eleMissHits_                .clear();
  elePFChIso_                 .clear();
  elePFPhoIso_                .clear();
  elePFNeuIso_                .clear();
  elePFPUIso_                 .clear();
  eleIDMVANonTrg_             .clear();
  eleIDMVATrg_                .clear();
  eleEcalDrivenSeed_          .clear();
  eleR9Full5x5_               .clear();
  eleIDbit_                   .clear();

  nEle_ = 0;

  edm::Handle<edm::View<pat::Electron> > electronHandle;
  e.getByToken(electronCollection_, electronHandle);

  if (!electronHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Electrons in event";
    return;
  }

  edm::Handle<edm::ValueMap<bool> >  veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  tight_id_decisions; 
  edm::Handle<edm::ValueMap<bool> >  heep_id_decisions;
  edm::Handle<edm::ValueMap<float> > eleNonTrgMVAValues;
  edm::Handle<edm::ValueMap<float> > eleTrgMVAValues;

  if (runeleIDVID_) {
    e.getByToken(eleVetoIdMapToken_ ,         veto_id_decisions);
    e.getByToken(eleLooseIdMapToken_ ,        loose_id_decisions);
    e.getByToken(eleMediumIdMapToken_,        medium_id_decisions);
    e.getByToken(eleTightIdMapToken_,         tight_id_decisions);
    e.getByToken(eleHEEPIdMapToken_ ,         heep_id_decisions);
    e.getByToken(eleNonTrgMVAValuesMapToken_, eleNonTrgMVAValues);
    e.getByToken(eleTrgMVAValuesMapToken_,    eleTrgMVAValues);
  }

  edm::Handle<reco::VertexCollection> recVtxs;
  e.getByToken(vtxLabel_, recVtxs);
  
  EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

  for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {

    eleCharge_          .push_back(iEle->charge());
    eleEn_              .push_back(iEle->energy());
    eleD0_              .push_back(iEle->gsfTrack()->dxy(pv));
    eleDz_              .push_back(iEle->gsfTrack()->dz(pv));
    elePt_              .push_back(iEle->pt());
    eleEta_             .push_back(iEle->eta());
    elePhi_             .push_back(iEle->phi());
    eleSCEn_            .push_back(iEle->superCluster()->energy());
    eleSCEta_           .push_back(iEle->superCluster()->eta());
    eleSCPhi_           .push_back(iEle->superCluster()->phi());
    eleSCRawEn_         .push_back(iEle->superCluster()->rawEnergy());
    eleHoverE_          .push_back(iEle->hcalOverEcal());

    ///https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d8/dac/GsfElectron_8h_source.html
    eleEoverP_          .push_back(iEle->eSuperClusterOverP());
    eleEoverPout_       .push_back(iEle->eEleClusterOverPout());
    eledEtaAtVtx_       .push_back(iEle->deltaEtaSuperClusterTrackAtVtx());
    eledPhiAtVtx_       .push_back(iEle->deltaPhiSuperClusterTrackAtVtx());
    eleConvVeto_        .push_back((Int_t)iEle->passConversionVeto()); // ConvVtxFit || missHit == 0
    eleMissHits_        .push_back(iEle->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));

    // VID calculation of (1/E - 1/p)
    if (iEle->ecalEnergy() == 0)   eleEoverPInv_.push_back(1e30);
    else if (!std::isfinite(iEle->ecalEnergy()))  eleEoverPInv_.push_back(1e30);
    else  eleEoverPInv_.push_back((1.0 - iEle->eSuperClusterOverP())/iEle->ecalEnergy());


    reco::GsfElectron::PflowIsolationVariables pfIso = iEle->pfIsolationVariables();
    elePFChIso_         .push_back(pfIso.sumChargedHadronPt);
    elePFPhoIso_        .push_back(pfIso.sumPhotonEt);
    elePFNeuIso_        .push_back(pfIso.sumNeutralHadronEt);
    elePFPUIso_         .push_back(pfIso.sumPUPt);

    /////quantities which were used for Run1 - these do not
    ///calculated through PF (meaning no energy is subtracted
    ///using PF)
    ///https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d9/d44/ElectronIDValueMapProducer_8cc_source.html
    ///line 120

    eleSigmaIEtaIEtaFull5x5_    .push_back(iEle->full5x5_sigmaIetaIeta());
    eleSigmaIPhiIPhiFull5x5_    .push_back(iEle->full5x5_sigmaIphiIphi());
    eleR9Full5x5_               .push_back(iEle->full5x5_r9());

    ///For HEEP ID
    eleEcalDrivenSeed_          .push_back(iEle->ecalDrivenSeed());



    //
    // Look up and save the ID decisions
    // 
    if (runeleIDVID_) {
      //edm::Ptr<reco::GsfElectron> recoEl(iEle);
      
      //const auto el = electrons->ptrAt(nEle_);
      const auto el = electronHandle->ptrAt(nEle_);
      
      UShort_t tmpeleIDbit = 0;

        ///el->electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto") also works
        bool isPassVeto  = (*veto_id_decisions)[el];
        if(isPassVeto) setbit(tmpeleIDbit, 0);
        //cout<<"isVeto: "<<isPassVeto<<endl;

        bool isPassLoose  = (*loose_id_decisions)[el];
        if(isPassLoose) setbit(tmpeleIDbit, 1);
        //cout<<"isLoose: "<<isPassLoose<<endl;

        bool isPassMedium = (*medium_id_decisions)[el];
        if(isPassMedium) setbit(tmpeleIDbit, 2);
        //cout<<"isMedium: "<<isPassMedium<<endl;

        bool isPassTight  = (*tight_id_decisions)[el];
        if(isPassTight) setbit(tmpeleIDbit, 3);
        //cout<<"isTight: "<<isPassTight<<endl;
      
        bool isPassHEEP = (*heep_id_decisions)[el];
        if(isPassHEEP) setbit(tmpeleIDbit, 4);
        //cout<<"isHeep: "<<isPassHEEP<<endl;

        eleIDMVANonTrg_.push_back((*eleNonTrgMVAValues)[el]);
        eleIDMVATrg_   .push_back((*eleTrgMVAValues)[el]);

      eleIDbit_.push_back(tmpeleIDbit);
      
      //cout<<"tmpele : eleIDbit: "<<tmpeleIDbit<<":"<<eleIDbit_[nEle_]<<endl;
    }//if(runeleIDVID_)

    nEle_++;
  }

}
