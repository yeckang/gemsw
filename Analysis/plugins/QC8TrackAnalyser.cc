// cd /cms/ldap_home/iawatson/scratch/GEM/CMSSW_10_1_5/src/ && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// cd ../../.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10
// system include files
#include <memory>
#include <cmath>
#include <iostream>
#include <map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
// GEM
#include "DataFormats/GEMDigi/interface/GEMDigiCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/GEMRecHit/interface/GEMSegmentCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/GEMStripTopology.h"
// Muon
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackingRecHit/interface/KfComponentsHolder.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TTree.h"

using namespace std;

typedef std::tuple<int, int> Key2;
typedef std::tuple<int, int, int> Key3;

class QC8TrackAnalyzer : public edm::EDAnalyzer {
public:
  explicit QC8TrackAnalyzer(const edm::ParameterSet&);
  ~QC8TrackAnalyzer();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------
  edm::Service<TFileService> fs;
  edm::EDGetTokenT<reco::TrackCollection> tracks_;
  edm::EDGetTokenT<vector<Trajectory>> trajs_;
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;

  MuonServiceProxy* theService_;

  TH1D* trackChi2_;

  std::map<Key2, TH2D*> track_occ_;
  std::map<Key2, TH2D*> rechit_occ_;
};

QC8TrackAnalyzer::QC8TrackAnalyzer(const edm::ParameterSet& iConfig)
{ 
  tracks_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
  trajs_ = consumes<vector<Trajectory>>(iConfig.getParameter<edm::InputTag>("trajs"));
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHitLabel"));
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters, consumesCollector(), MuonServiceProxy::UseEventSetupIn::RunAndEvent);
}

QC8TrackAnalyzer::~QC8TrackAnalyzer(){}

void QC8TrackAnalyzer::beginRun(const edm::Run& run, const edm::EventSetup& iSetup) {
  edm::ESHandle<GEMGeometry> hGEMGeom;
  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

  trackChi2_ = fs->make<TH1D>("track_chi2", "Normalized Track Chi2", 100, 0, 10);

  for (auto etaPart : GEMGeometry_->etaPartitions()){    
    auto etaPartId = etaPart->id();
    auto ch = etaPartId.chamber();
    auto ieta = etaPartId.ieta();

    Key2 key2(ch, ieta);
    track_occ_[key2] = fs->make<TH2D>(Form("track_occ_ch%d_ieta%d", ch, ieta),
                                      Form("track_occ_ch%d_ieta%d", ch, ieta),
                                      200, -100, 100,
                                      200, -100, 100);
    rechit_occ_[key2] = fs->make<TH2D>(Form("rechit_occ_ch%d_ieta%d", ch, ieta),
                                       Form("rechit_occ_ch%d_ieta%d", ch, ieta),
                                       200, -100, 100,
                                       200, -100, 100);
  }
}

void
QC8TrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  /* GEM Geometry */
  edm::ESHandle<GEMGeometry> hGEMGeom;
  iSetup.get<MuonGeometryRecord>().get(hGEMGeom);
  const GEMGeometry* GEMGeometry_ = &*hGEMGeom;

  edm::Handle<vector<reco::Track>> tracks;
  iEvent.getByToken(tracks_, tracks);

  edm::Handle<vector<Trajectory>> trajs;
  iEvent.getByToken(trajs_, trajs);
  
  edm::Handle<GEMRecHitCollection> gemRecHits;
  iEvent.getByToken(gemRecHits_,gemRecHits);

  theService_->update(iSetup);

  std::map<GEMDetId, TrajectoryStateOnSurface> tsosMap;
  for (std::vector<reco::Track>::const_iterator track = tracks->begin(); track != tracks->end(); ++track)
  {
    trackChi2_->Fill(track->normalizedChi2());
    auto traj = trajs->begin();
    for (auto trajMeas : traj->measurements()) {
      auto tsos = trajMeas.predictedState();
      auto rechit = trajMeas.recHit();
      auto gemId = GEMDetId(rechit->geographicalId());
      tsosMap[gemId] = tsos;
    }

    for (auto hit : track->recHits()) {
      auto etaPart = GEMGeometry_->etaPartition(hit->geographicalId());
      auto etaPartId = etaPart->id();

      //if (!hit->isValid()) continue;
      if (tsosMap.find(etaPartId) == tsosMap.end()) continue;
      auto tsos = tsosMap[etaPartId];

      auto lp_track = tsos.localPosition();

      auto range = gemRecHits->get(etaPartId);

      bool hasHit = false;
      float maxR = 500000;

      int chamber = etaPartId.chamber();
      int ieta = etaPartId.ieta();

      Key2 key2(chamber, ieta);
    
      track_occ_[key2]->Fill(lp_track.x(), lp_track.y());

      for (auto rechit = range.first; rechit != range.second; ++rechit) {
        auto lp_rechit = rechit->localPosition();

        auto deltaR = (lp_rechit - lp_track).mag();
        if (deltaR < maxR) {
          hasHit = true;
          maxR = deltaR; 
        }
      }
      if (hasHit) {
        rechit_occ_[key2]->Fill(lp_track.x(), lp_track.y());
      }
    }
  }
}

void QC8TrackAnalyzer::beginJob(){}
void QC8TrackAnalyzer::endJob(){}

void QC8TrackAnalyzer::endRun(edm::Run const&, edm::EventSetup const&){}
                   
//define this as a plug-in
DEFINE_FWK_MODULE(QC8TrackAnalyzer);
