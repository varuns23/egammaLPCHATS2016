#include <TString.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "egammaLPCHATS2016/ggNtuplizer/interface/GEDPhoIDTools.h"
#include "egammaLPCHATS2016/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

Int_t          nPho_;
vector<float>  phoE_;
vector<float>  phoEt_;
vector<float>  phoEta_;
vector<float>  phoPhi_;
vector<float>  phoSCE_;
vector<float>  phoSCRawE_;
vector<float>  phoSCEta_;
vector<float>  phoSCPhi_;
vector<float>  phoSCEtaWidth_;
vector<float>  phoSCPhiWidth_;
vector<int>    phohasPixelSeed_;
vector<int>    phoEleVeto_;
vector<float>  phoR9_;
vector<float>  phoHoverE_;
vector<float>  phoSigmaIEtaIEtaFull5x5_;
vector<float>  phoSigmaIEtaIPhiFull5x5_;
vector<float>  phoSigmaIPhiIPhiFull5x5_;
vector<float>  phoR9Full5x5_;
vector<float>  phoPFChIso_;
vector<float>  phoPFPhoIso_;
vector<float>  phoPFNeuIso_;
vector<float>  phoIDMVA_;

//Necessary for the Photon Footprint removal
template <class T, class U>
bool isInFootprint(const T& thefootprint, const U& theCandidate) {
  for ( auto itr = thefootprint.begin(); itr != thefootprint.end(); ++itr ) {

    if( itr.key() == theCandidate.key() ) return true;
    
  }
  return false;
}

void ggNtuplizer::branchesPhotons(TTree* tree) {
  
  tree->Branch("nPho",                    &nPho_);
  tree->Branch("phoE",                    &phoE_);
  tree->Branch("phoEt",                   &phoEt_);
  tree->Branch("phoEta",                  &phoEta_);
  tree->Branch("phoPhi",                  &phoPhi_);
  tree->Branch("phoSCE",                  &phoSCE_);
  tree->Branch("phoSCRawE",               &phoSCRawE_);
  tree->Branch("phoSCEta",                &phoSCEta_);
  tree->Branch("phoSCPhi",                &phoSCPhi_);
  tree->Branch("phoSCEtaWidth",           &phoSCEtaWidth_);
  tree->Branch("phoSCPhiWidth",           &phoSCPhiWidth_);
  tree->Branch("phohasPixelSeed",         &phohasPixelSeed_);
  tree->Branch("phoEleVeto",              &phoEleVeto_);
  tree->Branch("phoR9",                   &phoR9_);
  tree->Branch("phoHoverE",               &phoHoverE_);
  tree->Branch("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5_);
  tree->Branch("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5_);
  tree->Branch("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5_);
  tree->Branch("phoR9Full5x5",            &phoR9Full5x5_);
  tree->Branch("phoPFChIso",              &phoPFChIso_);
  tree->Branch("phoPFPhoIso",             &phoPFPhoIso_);
  tree->Branch("phoPFNeuIso",             &phoPFNeuIso_);
  tree->Branch("phoIDMVA",                        &phoIDMVA_);

}

void ggNtuplizer::fillPhotons(const edm::Event& e, const edm::EventSetup& es) {
  
  // cleanup from previous execution
  phoE_                 .clear();
  phoEt_                .clear();
  phoEta_               .clear();
  phoPhi_               .clear();
  phoSCE_               .clear();
  phoSCRawE_            .clear();
  phoSCEta_             .clear();
  phoSCPhi_             .clear();
  phoSCEtaWidth_        .clear();
  phoSCPhiWidth_        .clear();
  phohasPixelSeed_      .clear();
  phoEleVeto_           .clear();
  phoR9_                .clear();
  phoHoverE_            .clear();
  phoSigmaIEtaIEtaFull5x5_.clear();
  phoSigmaIEtaIPhiFull5x5_.clear();
  phoSigmaIPhiIPhiFull5x5_.clear();
  phoR9Full5x5_         .clear();
  phoPFChIso_           .clear();
  phoPFPhoIso_          .clear();
  phoPFNeuIso_          .clear();
  phoIDMVA_             .clear();

  nPho_ = 0;

  edm::Handle<edm::View<pat::Photon> > photonHandle;
  e.getByToken(photonCollection_, photonHandle);


  if (!photonHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Photons in event";
    return;
  }


  edm::Handle<reco::PhotonCollection> recoPhotonHandle;
  e.getByToken(recophotonCollection_, recoPhotonHandle);

  ///Photon ID in VID framwork - 11th may, 2015
  // Get the photon ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  // If you need more info, check with the EGM group.
  edm::Handle<edm::ValueMap<bool> >  loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  tight_id_decisions;
  edm::Handle<edm::ValueMap<float> > mvaValues;
  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
  
  if (runphoIDVID_) {
    e.getByToken(phoLooseIdMapToken_ ,  loose_id_decisions);
    e.getByToken(phoMediumIdMapToken_,  medium_id_decisions);
    e.getByToken(phoTightIdMapToken_ ,  tight_id_decisions);
    e.getByToken(phoMVAValuesMapToken_, mvaValues);

    e.getByToken(phoChargedIsolationToken_,       phoChargedIsolationMap);
    e.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
    e.getByToken(phoPhotonIsolationToken_,        phoPhotonIsolationMap);
  }

  EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

  edm::Handle<reco::VertexCollection> recVtxs;
  e.getByToken(vtxLabel_, recVtxs);

  edm::Handle<reco::VertexCollection> recVtxsBS;
  e.getByToken(vtxBSLabel_, recVtxsBS);

  edm::Handle<reco::TrackCollection> tracksHandle;
  e.getByToken(tracklabel_, tracksHandle);

  edm::Handle<reco::GsfElectronCollection> gsfElectronHandle;
  e.getByToken(gsfElectronlabel_, gsfElectronHandle);


  for (edm::View<pat::Photon>::const_iterator iPho = photonHandle->begin(); iPho != photonHandle->end(); ++iPho) {

    phoE_             .push_back(iPho->energy());
    phoEt_            .push_back(iPho->et());
    phoEta_           .push_back(iPho->eta());
    phoPhi_           .push_back(iPho->phi());
    phoSCE_           .push_back((*iPho).superCluster()->energy());
    phoSCRawE_        .push_back((*iPho).superCluster()->rawEnergy());
    phoSCEta_         .push_back((*iPho).superCluster()->eta());
    phoSCPhi_         .push_back((*iPho).superCluster()->phi());
    phoSCEtaWidth_    .push_back((*iPho).superCluster()->etaWidth());
    phoSCPhiWidth_    .push_back((*iPho).superCluster()->phiWidth());
    phohasPixelSeed_  .push_back((Int_t)iPho->hasPixelSeed());
    phoEleVeto_       .push_back((Int_t)iPho->passElectronVeto());
    phoR9_            .push_back(iPho->r9());
    phoHoverE_        .push_back(iPho->hadTowOverEm());

    std::vector<float> vCov = lazyToolnoZS.localCovariances( *((*iPho).superCluster()->seed()) );
    //const float see = (isnan(vCov[0]) ? 0. : sqrt(vCov[0]));
    const float spp = (isnan(vCov[2]) ? 0. : sqrt(vCov[2]));
    const float sep = vCov[1];

    phoSigmaIEtaIEtaFull5x5_ .push_back(iPho->full5x5_sigmaIetaIeta());
    phoSigmaIEtaIPhiFull5x5_ .push_back(sep);
    phoSigmaIPhiIPhiFull5x5_ .push_back(spp);
    phoR9Full5x5_            .push_back(iPho->full5x5_r9());

       

    if (runphoIDVID_) {
      //Photon ID in VID framwork - 11th may, 2015
      // Look up and save the ID decisions
      // 
      const auto pho = photonHandle->ptrAt(nPho_);
      
      UShort_t tmpphoIDbit = 0;

        phoPFChIso_              .push_back((*phoChargedIsolationMap)[pho]);
        phoPFPhoIso_             .push_back((*phoPhotonIsolationMap)[pho]);
        phoPFNeuIso_             .push_back((*phoNeutralHadronIsolationMap)[pho]);

        //cout<<"Photons "<<endl;
        bool isPassLoose  = (*loose_id_decisions)[pho];
        if(isPassLoose) setbit(tmpphoIDbit, 0);
        //cout<<"isPassLoose "<<isPassLoose<<endl;

        bool isPassMedium = (*medium_id_decisions)[pho];
        if(isPassMedium) setbit(tmpphoIDbit, 1);
        //cout<<"isPassMedium "<<isPassMedium<<endl;

        bool isPassTight  = (*tight_id_decisions)[pho];
        if(isPassTight) setbit(tmpphoIDbit, 2);
        //cout<<"isPassTight "<<isPassTight<<endl;

        phoIDMVA_.push_back((*mvaValues)[pho]);
      }

    nPho_++;
  }

}

void ggNtuplizer::cleanupPhotons() {

}
