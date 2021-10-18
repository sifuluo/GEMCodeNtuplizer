//////////////////////////////////////////////////////////////////////
//               Making Ntuples for GEM CSC analysis                //
//               Author: Sifu Luo                                   //
//               sifuluo@tamu.edu                                   //
//////////////////////////////////////////////////////////////////////

// CMSSW
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

// GEMCode
#include "GEMCode/GEMValidation/interface/MatcherManager.h"

// NtupleMaker tools
#include "GEMCode/GEMValidation/test/TreeDigi.cc"

// Other tools
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Muons
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "L1Trigger/L1TMuon/interface/MuonRawDigiTranslator.h"
#include "L1Trigger/L1TMuon/interface/RegionalMuonRawDigiTranslator.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCandFwd.h"
#include "DataFormats/L1TMuon/interface/RegionalMuonCand.h"
#include "L1Trigger/L1TMuon/interface/MicroGMTConfiguration.h"

// KBMTF
#include "DataFormats/L1TMuon/interface/L1MuKBMTrack.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"

// ROOT
#include <TROOT.h>
#include <TTree.h>

// STD
#include <iomanip>
#include <sstream>
#include <iostream>
#include <memory>
#include <string>
#include <math.h>
#include <bitset>
#include <fstream>
#include <cmath>
#include <vector>
#include <set>
#include <sstream>
#include <stdexcept>

// Unclear headers
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTCluster.h"
#include "DataFormats/L1TrackTrigger/interface/TTStub.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTClusterAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTStubAssociationMap.h"
#include "SimTracker/TrackTriggerAssociation/interface/TTTrackAssociationMap.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"

#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"
#include "Geometry/CommonDetUnit/interface/GeomDetType.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetType.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelTopologyBuilder.h"
#include "Geometry/Records/interface/StackedTrackerGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"

#include "DataFormats/CSCDigi/interface/CSCALCTDigi.h"
#include "DataFormats/CSCDigi/interface/CSCCLCTDigi.h"

using namespace std;
using namespace edm;

class NtupleMaker : public edm::one::EDAnalyzer<> {
public:
  explicit NtupleMaker(const edm::ParameterSet& iConfig);
  virtual ~NtupleMaker();

  virtual void beginJob();
  virtual void endJob();
  virtual void analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup);

  virtual void PrintHits(int iset, int idigi = -1);
  virtual int SaveHitMatrix(std::vector< std::vector<unsigned short> > hits, std::vector<int>* b_hit, std::vector<int>* b_pos, bool doprint = false, bool isclct = false);
  virtual std::vector<int> IntsToBinary(int n);
  virtual GlobalPoint getGlobalPointDigi(unsigned int rawId, const GEMDigi& d);

protected:

private:
  // ParameterSet passed from python configuration
  edm::ParameterSet config;

  int MyProcess;
  bool DebugMode;
  double TP_minPt;
  double TP_maxEta;
  double TP_maxZ0;
  bool IsRun4;
  bool useGEMs;
  bool Print_matchCscStubs;
  bool Print_allCscStubs;
  bool Print_all;
  bool Print_ALCT;
  bool Print_CLCT;

  int nEventMultiHitLayer;

  edm::InputTag TrackingParticleInputTag;

  std::string getFloatPointDataWord(const l1t::RegionalMuonCand& l1mu) const;
  std::string getGlobalPhi(const l1t::RegionalMuonCand& l1mu) const;

  edm::EDGetToken m_emtfToken;
  // edm::EDGetToken m_bmtfToken;
  // edm::EDGetToken m_omtfToken;

  edm::EDGetTokenT< std::vector< TrackingParticle > > TrackingParticleToken_;
  edm::EDGetTokenT<edm::SimVertexContainer> simVertexInput_;

  edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> lctToken_;
  edm::EDGetTokenT<CSCALCTDigiCollection> alctToken_;
  edm::EDGetTokenT<CSCCLCTDigiCollection> clctToken_;
  edm::EDGetTokenT<GEMDigiCollection> gemDigiToken_;
  edm::EDGetTokenT<GEMPadDigiCollection> gemPadDigiToken_;
  edm::EDGetTokenT<GEMPadDigiClusterCollection> gemPadDigiClusterToken_;

  // Ntuple
  TTree* eventTree;

  // Regional Muon candidates
  // std::vector<int>*   m_EMTF_muon_n;
  // std::vector<float>* m_EMTF_muon_pt;
  // std::vector<float>* m_EMTF_muon_eta;
  // std::vector<float>* m_EMTF_muon_phi;
  // std::vector<int>*   m_EMTF_muon_c;

  // std::vector<int>*   m_OMTF_muon_n;
  // std::vector<float>* m_OMTF_muon_pt;
  // std::vector<float>* m_OMTF_muon_eta;
  // std::vector<float>* m_OMTF_muon_phi;
  // std::vector<int>*   m_OMTF_muon_c;
  //
  // std::vector<int>*   m_BMTF_muon_n;
  // std::vector<float>* m_BMTF_muon_pt;
  // std::vector<float>* m_BMTF_muon_eta;
  // std::vector<float>* m_BMTF_muon_phi;
  // std::vector<int>*   m_BMTF_muon_c;

  std::vector<float>* m_matchmuon_pt;
  std::vector<float>* m_matchmuon_eta;
  std::vector<float>* m_matchmuon_phi;
  std::vector<int>*   m_matchmuon_charge;
  std::vector<int>*   m_matchmuon_type;
  std::vector<int>*   m_matchmuon_quality;

  TreeDigi *tp, *cscSimHit, *gemSimHit;
  TreeDigi *allCscStubsLCT, *allCscStubsALCT, *allCscStubsCLCT;
  TreeDigi *allALCT, *allCLCT, *allGemDigi;
  TreeDigi *matchCscStubsLCT, *matchCscStubsALCT, *matchCscStubsCLCT, *matchGemDigi;
  TreeDigi *matchCscGEM1, *matchCscGEM2, *allCscGEM1, *allCscGEM2, *gemPadDigi;
  TreeDigi *matchGemPadDigiCluster, *allGemPadDigiCluster;

  std::unique_ptr<MatcherManager> match;

  const TrackingGeometry* geometry_;

  edm::ESGetToken<GEMGeometry, MuonGeometryRecord> geomToken_;
  const GEMGeometry* gemGeometry_;
};


NtupleMaker::NtupleMaker(edm::ParameterSet const& iConfig) :
config(iConfig)
{
  MyProcess = iConfig.getParameter< int >("MyProcess");
  DebugMode        = iConfig.getParameter< bool >("DebugMode");
  TP_minPt         = iConfig.getParameter< double >("TP_minPt");
  TP_maxEta        = iConfig.getParameter< double >("TP_maxEta");
  TP_maxZ0         = iConfig.getParameter< double >("TP_maxZ0");
  IsRun4           = iConfig.getParameter< bool >("IsRun4");
  Print_matchCscStubs = iConfig.getParameter< bool >("Print_matchCscStubs");
  Print_allCscStubs   = iConfig.getParameter< bool >("Print_allCscStubs");
  Print_all         = iConfig.getParameter< bool >("Print_all");
  Print_ALCT        = iConfig.getParameter< bool >("Print_ALCT");
  Print_CLCT        = iConfig.getParameter< bool >("Print_CLCT");
  useGEMs           = iConfig.getParameter< bool >("useGEMs");

  TrackingParticleInputTag = iConfig.getParameter<edm::InputTag>("TrackingParticleInputTag");
  TrackingParticleToken_ = consumes< std::vector< TrackingParticle > >(TrackingParticleInputTag);

  // m_emtfToken = consumes<l1t::RegionalMuonCandBxCollection>(edm::InputTag("simEmtfDigis","EMTF"));
  // m_bmtfToken = consumes<l1t::RegionalMuonCandBxCollection>(edm::InputTag(IsRun4 ? "simBmtfDigis": "gmtStage2Digis","BMTF"));
  // m_omtfToken = consumes<l1t::RegionalMuonCandBxCollection>(edm::InputTag(IsRun4 ? "simOmtfDigis": "gmtStage2Digis","OMTF"));

  const auto& simVertex = iConfig.getParameter<edm::ParameterSet>("simVertex");
  simVertexInput_ = consumes<edm::SimVertexContainer>(simVertex.getParameter<edm::InputTag>("inputTag"));

  const auto& P_cscLCT = iConfig.getParameter<edm::ParameterSet>("cscLCT");
  lctToken_ = consumes<CSCCorrelatedLCTDigiCollection>(P_cscLCT.getParameter<edm::InputTag>("inputTag"));
  cout <<endl<< P_cscLCT <<endl;

  const auto& P_cscALCT = iConfig.getParameter<edm::ParameterSet>("cscALCT");
  alctToken_ = consumes<CSCALCTDigiCollection>(P_cscALCT.getParameter<edm::InputTag>("inputTag"));

  const auto& P_cscCLCT = iConfig.getParameter<edm::ParameterSet>("cscCLCT");
  clctToken_ = consumes<CSCCLCTDigiCollection>(P_cscCLCT.getParameter<edm::InputTag>("inputTag"));

  const auto& P_gemDigi = iConfig.getParameter<edm::ParameterSet>("gemStripDigi");
  gemDigiToken_ = consumes<GEMDigiCollection>(P_gemDigi.getParameter<edm::InputTag>("inputTag"));

  const auto& P_gemPadDigi = iConfig.getParameter<edm::ParameterSet>("gemPadDigi");
  gemPadDigiToken_ = consumes<GEMPadDigiCollection>(P_gemPadDigi.getParameter<edm::InputTag>("inputTag"));

  const auto& P_gemPadDigiCluster = iConfig.getParameter<edm::ParameterSet>("gemPadCluster");
  gemPadDigiClusterToken_ = consumes<GEMPadDigiClusterCollection>(P_gemPadDigiCluster.getParameter<edm::InputTag>("inputTag"));

  geomToken_ = esConsumes<GEMGeometry, MuonGeometryRecord>();

  match.reset(new MatcherManager(iConfig, consumesCollector()));
}

NtupleMaker::~NtupleMaker()
{
}

void NtupleMaker::endJob()
{
  cerr << "Number of event with layers having multiple hits: " << nEventMultiHitLayer <<endl;
  cerr << "NtupleMaker::endJob" << endl;
}

void NtupleMaker::beginJob()
{
  cerr << "NtupleMaker::beginJob" << endl;

  // if (true) {
  //   cout << " Types of LCTs:" <<endl;
  //   CSCCorrelatedLCTDigi tmplct;
  //   tmplct.setType(CSCCorrelatedLCTDigi::CLCTALCT);    0 //CLCTALCT
  //   cout << tmplct.getType()<<" , ";
  //   tmplct.setType(CSCCorrelatedLCTDigi::ALCTCLCT);    1 //ALCTCLCT
  //   cout << tmplct.getType()<<" , ";
  //   tmplct.setType(CSCCorrelatedLCTDigi::ALCTCLCTGEM); 2 //ALCTCLCTGEM
  //   cout << tmplct.getType()<<" , ";
  //   tmplct.setType(CSCCorrelatedLCTDigi::ALCTCLCT2GEM);3 //ALCTCLCT2GEM
  //   cout << tmplct.getType()<<" , ";
  //   tmplct.setType(CSCCorrelatedLCTDigi::ALCT2GEM);    4 //ALCT2GEM
  //   cout << tmplct.getType()<<" , ";
  //   tmplct.setType(CSCCorrelatedLCTDigi::CLCT2GEM);    5 //CLCT2GEM
  //   cout << tmplct.getType()<<" , ";
  //   tmplct.setType(CSCCorrelatedLCTDigi::CLCTONLY);    6 //CLCTONLY
  //   cout << tmplct.getType()<<" , ";
  //   tmplct.setType(CSCCorrelatedLCTDigi::ALCTONLY);    7 //ALCTONLY
  //   cout << tmplct.getType()<<endl;
  // }

  nEventMultiHitLayer = 0;

  edm::Service<TFileService> fs;

  // m_EMTF_muon_n   = new std::vector<int>;
  // m_EMTF_muon_pt  = new std::vector<float>;
  // m_EMTF_muon_eta = new std::vector<float>;
  // m_EMTF_muon_phi = new std::vector<float>;
  // m_EMTF_muon_c   = new std::vector<int>;

  // m_OMTF_muon_n = new std::vector<int>;
  // m_OMTF_muon_pt = new std::vector<float>;
  // m_OMTF_muon_eta = new std::vector<float>;
  // m_OMTF_muon_phi = new std::vector<float>;
  // m_OMTF_muon_c = new std::vector<int>;
  //
  // m_BMTF_muon_n = new std::vector<int>;
  // m_BMTF_muon_pt = new std::vector<float>;
  // m_BMTF_muon_eta = new std::vector<float>;
  // m_BMTF_muon_phi = new std::vector<float>;
  // m_BMTF_muon_c = new std::vector<int>;

  m_matchmuon_pt      = new std::vector<float>;
  m_matchmuon_eta     = new std::vector<float>;
  m_matchmuon_phi     = new std::vector<float>;
  m_matchmuon_charge  = new std::vector<int>;
  m_matchmuon_type    = new std::vector<int>;
  m_matchmuon_quality = new std::vector<int>;

  eventTree = fs->make<TTree>("eventTree", "Event tree");

  // eventTree->Branch("EMTF_muon_n",	 &m_EMTF_muon_n);
  // eventTree->Branch("EMTF_muon_pt",  	 &m_EMTF_muon_pt);
  // eventTree->Branch("EMTF_muon_eta", 	 &m_EMTF_muon_eta);
  // eventTree->Branch("EMTF_muon_phi", 	 &m_EMTF_muon_phi);
  // eventTree->Branch("EMTF_muon_c", 	 &m_EMTF_muon_c);

  // eventTree->Branch("OMTF_muon_n",	 &m_OMTF_muon_n);
  // eventTree->Branch("OMTF_muon_pt", 	 &m_OMTF_muon_pt);
  // eventTree->Branch("OMTF_muon_eta", 	 &m_OMTF_muon_eta);
  // eventTree->Branch("OMTF_muon_phi", 	 &m_OMTF_muon_phi);
  // eventTree->Branch("OMTF_muon_c", 	 &m_OMTF_muon_c);
  //
  // eventTree->Branch("BMTF_muon_n",	 &m_BMTF_muon_n);
  // eventTree->Branch("BMTF_muon_pt", 	 &m_BMTF_muon_pt);
  // eventTree->Branch("BMTF_muon_eta", 	 &m_BMTF_muon_eta);
  // eventTree->Branch("BMTF_muon_phi", 	 &m_BMTF_muon_phi);
  // eventTree->Branch("BMTF_muon_c", 	 &m_BMTF_muon_c);

  eventTree->Branch("matchmuon_pt", &m_matchmuon_pt);
  eventTree->Branch("matchmuon_eta", &m_matchmuon_eta);
  eventTree->Branch("matchmuon_phi", &m_matchmuon_phi);
  eventTree->Branch("matchmuon_charge",&m_matchmuon_charge);
  eventTree->Branch("matchmuon_type",&m_matchmuon_type);
  eventTree->Branch("matchmuon_quality",&m_matchmuon_quality);

  tp                     = new TreeDigi();
  cscSimHit              = new TreeDigi();
  gemSimHit              = new TreeDigi();
  allCscStubsLCT         = new TreeDigi();
  allCscStubsALCT        = new TreeDigi();
  allCscStubsCLCT        = new TreeDigi();
  allALCT                = new TreeDigi();
  allCLCT                = new TreeDigi();
  allGemDigi             = new TreeDigi();
  matchCscStubsLCT       = new TreeDigi();
  matchCscStubsALCT      = new TreeDigi();
  matchCscStubsCLCT      = new TreeDigi();
  matchGemDigi           = new TreeDigi();
  matchCscGEM1           = new TreeDigi();
  matchCscGEM2           = new TreeDigi();
  allCscGEM1             = new TreeDigi();
  allCscGEM2             = new TreeDigi();
  gemPadDigi             = new TreeDigi();
  matchGemPadDigiCluster = new TreeDigi();
  allGemPadDigiCluster   = new TreeDigi();

  tp                    ->Init(eventTree,"tp"                    ,"TP"               ,false);
  cscSimHit             ->Init(eventTree,"cscSimHit"             ,"SimHit"           ,true);
  gemSimHit             ->Init(eventTree,"gemSimHit"             ,"SimHit"           ,true);
  allCscStubsLCT        ->Init(eventTree,"allCscStubsLCT"        ,"LCT"              ,false);
  allCscStubsALCT       ->Init(eventTree,"allCscStubsALCT"       ,"ALCT"             ,false);
  allCscStubsCLCT       ->Init(eventTree,"allCscStubsCLCT"       ,"CLCT"             ,false);
  allALCT               ->Init(eventTree,"allALCT"               ,"ALCT"             ,false);
  allCLCT               ->Init(eventTree,"allCLCT"               ,"CLCT"             ,false);
  allGemDigi            ->Init(eventTree,"allGemDigi"            ,"GEM"              ,false);
  matchCscStubsLCT      ->Init(eventTree,"matchCscStubsLCT"      ,"LCT"              ,true);
  matchCscStubsALCT     ->Init(eventTree,"matchCscStubsALCT"     ,"ALCT"             ,false);
  matchCscStubsCLCT     ->Init(eventTree,"matchCscStubsCLCT"     ,"CLCT"             ,false);
  matchGemDigi          ->Init(eventTree,"matchGemDigi"          ,"GEM"              ,true);
  matchCscGEM1          ->Init(eventTree,"matchCscGEM1"          ,"GEMPad"           ,true);
  matchCscGEM2          ->Init(eventTree,"matchCscGEM2"          ,"GEMPad"           ,true);
  allCscGEM1            ->Init(eventTree,"allCscGEM1"            ,"GEMPad"           ,true);
  allCscGEM2            ->Init(eventTree,"allCscGEM2"            ,"GEMPad"           ,true);
  gemPadDigi            ->Init(eventTree,"gemPadDigi"            ,"GEMPad"           ,false);
  matchGemPadDigiCluster->Init(eventTree,"matchGemPadDigiCluster","GEMPadDigiCluster",true);
  allGemPadDigiCluster  ->Init(eventTree,"allGemPadDigiCluster"  ,"GEMPadDigiCluster",false);
}

void NtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if (DebugMode) cout <<endl<<endl<< "Starting a new event-------------------------------------------------------------------------------------------------------" << endl;
  if (!(MyProcess==13 || MyProcess==11 || MyProcess==211 || MyProcess==6 || MyProcess==15 || MyProcess==1)) {
    cout << "The specified MyProcess is invalid! Exiting..." << endl;
    return;
  }

  // m_EMTF_muon_n->clear();
  // m_EMTF_muon_pt->clear();
  // m_EMTF_muon_eta->clear();
  // m_EMTF_muon_phi->clear();
  // m_EMTF_muon_c->clear();

  // m_OMTF_muon_n->clear();
  // m_OMTF_muon_pt->clear();
  // m_OMTF_muon_eta->clear();
  // m_OMTF_muon_phi->clear();
  // m_OMTF_muon_c->clear();
  //
  // m_BMTF_muon_n->clear();
  // m_BMTF_muon_pt->clear();
  // m_BMTF_muon_eta->clear();
  // m_BMTF_muon_phi->clear();
  // m_BMTF_muon_c->clear();

  m_matchmuon_pt->clear();
  m_matchmuon_eta->clear();
  m_matchmuon_phi->clear();
  m_matchmuon_charge->clear();
  m_matchmuon_type->clear();
  m_matchmuon_quality->clear();

  tp->Reset();
  cscSimHit->Reset();
  gemSimHit->Reset();
  allCscStubsLCT->Reset();
  allCscStubsALCT->Reset();
  allCscStubsCLCT->Reset();
  allALCT->Reset();
  allCLCT->Reset();
  allGemDigi->Reset();
  matchCscStubsLCT->Reset();
  matchCscStubsALCT->Reset();
  matchCscStubsCLCT->Reset();
  matchGemDigi->Reset();
  matchCscGEM1->Reset();
  matchCscGEM2->Reset();
  allCscGEM1->Reset();
  allCscGEM2->Reset();
  gemPadDigi->Reset();
  matchGemPadDigiCluster->Reset();
  allGemPadDigiCluster->Reset();

  if (DebugMode) cout << "Finished branch initialization" << endl;

  edm::Handle< std::vector< TrackingParticle > > TrackingParticleHandle;
  iEvent.getByToken(TrackingParticleToken_, TrackingParticleHandle);

  edm::ESHandle<TrackerGeometry> geometryHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(geometryHandle);

  edm::ESHandle<TrackerTopology> tTopoHandle;
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);

  edm::ESHandle<TrackerGeometry> tGeomHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(tGeomHandle);
  match->init(iEvent,iSetup);

  edm::Handle<edm::SimVertexContainer> sim_vertices;
  iEvent.getByToken(simVertexInput_, sim_vertices);
  const edm::SimVertexContainer & sim_vert = *sim_vertices.product();

  gemGeometry_ = &iSetup.getData(geomToken_);

  int tp_index = 0;
  std::vector< TrackingParticle >::const_iterator iterTP;
  if (DebugMode) cout << "Started TP iteration" << endl;
  for (iterTP = TrackingParticleHandle->begin(); iterTP != TrackingParticleHandle->end(); ++iterTP) {
    if (DebugMode) cout << "Started "<< tp_index << " TP information Collecting" << endl;
    edm::Ptr< TrackingParticle > tp_ptr(TrackingParticleHandle, tp_index);

    int tmp_eventid = iterTP->eventId().event();
    if (MyProcess != 1 && tmp_eventid > 0) continue; //only care about tracking particles from the primary interaction (except for MyProcess==1, i.e. looking at all TPs)

    float tmp_tp_pt  = iterTP->pt();
    float tmp_tp_eta = iterTP->eta();
    float tmp_tp_phi = iterTP->phi();
    float tmp_tp_vz  = iterTP->vz();
    float tmp_tp_vx  = iterTP->vx();
    float tmp_tp_vy  = iterTP->vy();
    int tmp_tp_pdgid = iterTP->pdgId();
    float tmp_tp_z0_prod = tmp_tp_vz;
    float tmp_tp_d0_prod = -tmp_tp_vx*sin(tmp_tp_phi) + tmp_tp_vy*cos(tmp_tp_phi);

    if (MyProcess==13 && abs(tmp_tp_pdgid) != 13) continue;
    if (MyProcess==11 && abs(tmp_tp_pdgid) != 11) continue;
    if ((MyProcess==6 || MyProcess==15 || MyProcess==211) && abs(tmp_tp_pdgid) != 211) continue;

    if (tmp_tp_pt < TP_minPt) continue;
    if (fabs(tmp_tp_eta) > TP_maxEta) continue;

    // Calculation of d0 and z0

    float tmp_tp_t = tan(2.0*atan(1.0)-2.0*atan(exp(-tmp_tp_eta)));

    float delx = -tmp_tp_vx;
    float dely = -tmp_tp_vy;

    float A = 0.01*0.5696;
    float Kmagnitude = A / tmp_tp_pt;

    float tmp_tp_charge = tp_ptr->charge();
    float K = Kmagnitude * tmp_tp_charge;
    float d = 0;

    float tmp_tp_x0p = delx - (d + 1./(2. * K)*sin(tmp_tp_phi));
    float tmp_tp_y0p = dely + (d + 1./(2. * K)*cos(tmp_tp_phi));
    float tmp_tp_rp = sqrt(tmp_tp_x0p*tmp_tp_x0p + tmp_tp_y0p*tmp_tp_y0p);
    float tmp_tp_d0 = tmp_tp_charge*tmp_tp_rp - (1. / (2. * K));

    tmp_tp_d0 = tmp_tp_d0*(-1); //fix d0 sign

    static double pi = 4.0*atan(1.0);
    float delphi = tmp_tp_phi-atan2(-K*tmp_tp_x0p,K*tmp_tp_y0p);
    if (delphi<-pi) delphi+=2.0*pi;
    if (delphi>pi) delphi-=2.0*pi;
    float tmp_tp_z0 = tmp_tp_vz+tmp_tp_t*delphi/(2.0*K);
    // ----------------------------------------------------------------------------------------------

    if (fabs(tmp_tp_z0) > TP_maxZ0) continue;

    float dxy = sqrt(tmp_tp_vx*tmp_tp_vx + tmp_tp_vy*tmp_tp_vy);
    float tmp_tp_dxy = dxy;
    if (MyProcess==6 && (dxy > 1.0)) continue;

    tp->pt->push_back(tmp_tp_pt);
    tp->eta->push_back(tmp_tp_eta);
    tp->phi->push_back(tmp_tp_phi);
    tp->dxy->push_back(tmp_tp_dxy);
    tp->z0->push_back(tmp_tp_z0);
    tp->d0->push_back(tmp_tp_d0);
    tp->z0_prod->push_back(tmp_tp_z0_prod);
    tp->d0_prod->push_back(tmp_tp_d0_prod);
    tp->pdgid->push_back(tmp_tp_pdgid);
    tp->eventid->push_back(tmp_eventid);
    tp->charge->push_back(tmp_tp_charge);

    const SimTrack& t(tp_ptr->g4Tracks()[0]);
    if(abs(t.type())==13){
      if (DebugMode) cout << "Finished " << tp_index << " TP information and started matcher as this TP is a muon" <<endl;
      const SimVertex v = (t.vertIndex() < int(sim_vert.size())) ? sim_vert[t.vertIndex()] : SimVertex();


      match->match(t, v);

      // cout << "Mathcer initialized" <<endl;
      std::shared_ptr<CSCSimHitMatcher> cscsimhits = match->cscSimHits();
      std::shared_ptr<GEMSimHitMatcher> gemsimhits;
      gemsimhits = match->gemSimHits();
      // cout << "GEMSimHits initialized" <<endl;
      for (int istation = 1; istation < 5; ++istation) {
        const auto& cscIds = cscsimhits->chamberIdsStation(istation);
        for (const auto& p1 : cscIds) {
          const auto& hits = cscsimhits->hitsInChamber(p1);
          for (auto& hit : hits) {
            PSimHitContainer hitc;
            hitc.push_back(hit);
            GlobalPoint gp = cscsimhits->simHitsMeanPosition(hitc);
            cscSimHit->FillGP(gp);
            cscSimHit->FillSimHit(istation, tp_index);
          }
        }
      }
      if (DebugMode) cout << "cscSimHits Finished, starting gemSimHits" << endl;
      const auto& gemIds = gemsimhits->detIds();
      for (const auto&p1 :gemIds) {
        GEMDetId id1(p1);
        int istation = id1.station();
        const auto& hits = gemsimhits->hitsInDetId(p1);
        for (auto& hit : hits) {
          PSimHitContainer hitc;
          hitc.push_back(hit);
          GlobalPoint gp = gemsimhits->simHitsMeanPosition(hits);
          gemSimHit->FillGP(gp);
          gemSimHit->FillSimHit(istation, tp_index);
        }
      }
      if (DebugMode) cout << "gemSimhits Finished, starting muonCandidate" <<endl;
      auto muonCandidate = match->l1Muons()->emtfCand();
      if(muonCandidate){
        //std::cout<<"found match, extract information"<<std::endl;
        m_matchmuon_pt->push_back(muonCandidate->pt());
        m_matchmuon_eta->push_back(muonCandidate->eta());
        m_matchmuon_phi->push_back(muonCandidate->phi());
        m_matchmuon_charge->push_back(muonCandidate->charge());
        m_matchmuon_quality->push_back(muonCandidate->quality());
        int type=-1;
        m_matchmuon_type->push_back(type);
      }
      else{
        m_matchmuon_pt->push_back(-999.);
        m_matchmuon_eta->push_back(-999.);
        m_matchmuon_phi->push_back(-999.);
        m_matchmuon_charge->push_back(-999);
        m_matchmuon_type->push_back(-999);
        m_matchmuon_quality->push_back(-999);
      }

      //CSCStubMatcher
      auto cscStubs = match->cscStubs();
      if (DebugMode) cout << "Started Loop over matchCscStubsA/CLCT (match->cscStubs()) " << endl;
      int digicount = 0;
      for (int detid_int : cscStubs->chamberIdsLCT(0)) {
        CSCDetId detid_(detid_int);
        int digi_index = 0;
        bool doprinta = Print_matchCscStubs && Print_ALCT;
        bool doprintc = Print_matchCscStubs && Print_CLCT;
        for (auto digi_ : cscStubs->lctsInChamber(detid_) ){
          if (DebugMode) cout << "" << endl;
          auto gp = cscStubs->getGlobalPosition(detid_int,digi_);
          // auto gp2 = cscStubs->getGlobalPosition2(detid_int,digi_);
          // if (gp2.eta() != gp.eta() || gp2.phi() != gp.phi()) {
          //   cout << endl << "matched CscStubsLCT GlobalPoint before change: eta = " << gp2.eta()<< " , phi = " << gp2.phi() << endl;
          //   cout         << "matched CscStubsLCT GlobalPoint after  change: eta = " << gp.eta() << " , phi = " << gp.phi()  << endl;
          // }
          matchCscStubsLCT->FillGP(gp);
          matchCscStubsLCT->FillLCT(digi_,detid_.rawId(),tp_index);
          if (Print_matchCscStubs) cout << "detid_int = " <<detid_int<<", detid = " << int(detid_) << ", rawId = " << detid_.rawId() << endl;

          // matchCscStubsALCT
          const auto& alctDigi = digi_.getALCT();
          std::vector< std::vector<unsigned short> > alcthits = alctDigi.getHits();
          if (doprinta) cout << "For matchCscStubsALCTs; DetId: "<< detid_.rawId() <<", Digi Index: " << digi_index << ", keywire: "<< alctDigi.getKeyWG()<<" , lctDigiWire = "<<digi_.getKeyWG()<<endl;
          if (alctDigi.getKeyWG() != digi_.getKeyWG()) cout << "Inconsistent keywire in matchCscStubsLCT: alctDigiWire = "<< alctDigi.getKeyWG() <<" , lctDigiWire = "<<digi_.getKeyWG()<<endl;
          // int alctmultihit = SaveHitMatrix(alcthits, m_matchCscStubsALCT_hit, m_matchCscStubsALCT_position,doprinta,false);
          int alctmultihit = SaveHitMatrix(alcthits, matchCscStubsALCT->hit, matchCscStubsALCT->position,doprinta,false);
          if (doprinta && alctmultihit) cout << "ALCT Multihit found for matchCscStubsALCT" << endl;
          matchCscStubsALCT->FillALCT(alctDigi,detid_.rawId());

          // matchCscStubsCLCT
          const auto& clctDigi = digi_.getCLCT();
          std::vector< std::vector<unsigned short> > clcthits = clctDigi.getHits();
          if (doprintc) cout << "For matchCscStubsCLCTs; DetId: "<< detid_.rawId() <<", Digi Index: " << digi_index << ", strip: "<< clctDigi.getStrip() <<" , lctDigiStrip = "<< digi_.getStrip()<<endl;
          // if (clctDigi.getStrip() != digi_.getStrip()) cout << "Inconsistent strip: clctDigiStrip = "<< clctDigi.getStrip() <<" , lctDigiStrip = "<< digi_.getStrip() <<endl;
          if (clctDigi.getSlope() != digi_.getSlope()) cout << "Inconsistent Slope in matchCscStubsLCT: clctDigiSlope = "<< clctDigi.getSlope() <<" , lctDigiSlope = "<< digi_.getSlope() <<endl;
          // int clctmultihit = SaveHitMatrix(clcthits, m_matchCscStubsCLCT_hit, m_matchCscStubsCLCT_position,doprintc,true);
          int clctmultihit = SaveHitMatrix(clcthits, matchCscStubsCLCT->hit, matchCscStubsCLCT->position,doprintc,true);
          if (doprintc && clctmultihit) cout << "CLCT Multihit found for matchCscStubsCLCT" << endl;
          matchCscStubsCLCT->FillCLCT(clctDigi,detid_.rawId());

          // GEMPad for matchCscStubs
          if (DebugMode) cout << " Finishing a matchCscStub, starting matchedCscStub GEMPads" <<endl;
          if (detid_.ring() == 1 and (detid_.station() == 1 or detid_.station() == 2)) {
            bool matchl1(false), matchl2(false), printgem(false);
            if (digi_.getGEM1().pad() != 255 || digi_.getGEM2().pad() != 255 || digi_.getGEM1().nPartitions() != 8 || digi_.getGEM2().nPartitions() != 8) {
              printgem = true;
              cout << "In station " << detid_.station() << " , ring " << detid_.ring() << " , eta = " << gp.eta() << " , matchCscStubsLCT rawId: " <<detid_.rawId() << endl;
            }
            const GEMDetId gemDetIdL1(detid_.zendcap(), 1, detid_.station(), 1, detid_.chamber(), 0);
            for (const auto& p : match->gemDigis()->padsInChamber(gemDetIdL1.rawId())) {
              if (p == digi_.getGEM1() && p.isValid()) {
                cout << "GEM1 is Valid" <<endl;
                // auto gp1 = match->gemDigis()->getGlobalPointPad(gemDetIdL1.rawId(),p);
                // matchCscGEM1->FillGP(gp1);
                matchCscGEM1->FillGEMPad(p,gemDetIdL1.rawId(),digicount);
                matchl1 = true;
                break;
              }
              else if (p == digi_.getGEM1() && !(p.isValid())){
                cout <<"GEM1 is Invalid" <<endl;
              }
            }
            if (!matchl1) {
              matchCscGEM1->FillGP0();
              matchCscGEM1->FillGEMPad0(gemDetIdL1.rawId(),digicount);
            }

            // Check if matched to an GEM pad L2
            const GEMDetId gemDetIdL2(detid_.zendcap(), 1, detid_.station(), 2, detid_.chamber(), 0);
            for (const auto& p : match->gemDigis()->padsInChamber(gemDetIdL2.rawId())) {
              if (p == digi_.getGEM2() && p.isValid()) {
                cout << "GEM2 is Valid" <<endl;
                // auto gp2 = match->gemDigis()->getGlobalPointPad(gemDetIdL2.rawId(),p);
                // matchCscGEM2->FillGP(gp2);
                matchCscGEM2->FillGEMPad(p,gemDetIdL2.rawId(),digicount);
                matchl2 = true;
                break;
              }
              else if (p == digi_.getGEM2() && !(p.isValid())){
                cout <<"GEM2 is Invalid" <<endl;
              }
            }
            if (!matchl2) {
              matchCscGEM2->FillGP0();
              matchCscGEM2->FillGEMPad0(gemDetIdL2.rawId(),digicount);
            }

            if (printgem) {

              cout << "GEM1 pad = " << digi_.getGEM1().pad() << ", part = " << digi_.getGEM1().nPartitions() << (digi_.getGEM1().isValid()? ", Valid " : ", Invalid ") << (matchl1 ? ", Filled " : ", Not Filled") << endl;
              cout << "GEM2 pad = " << digi_.getGEM2().pad() << ", part = " << digi_.getGEM2().nPartitions() << (digi_.getGEM2().isValid()? ", Valid " : ", Invalid ") << (matchl2 ? ", Filled " : ", Not Filled") << endl;
            }
            // const auto& gem1 = digi_.getGEM1();
            // const auto& gem2 = digi_.getGEM2();
            // cout << " gem1 " << (gem1.isValid() ? "Valid" : "Invalid") <<endl;
            // cout << " gem2 " << (gem2.isValid() ? "Valid" : "Invalid") <<endl;

            // cout << "Read" <<endl;
            // const GEMDetId gemDetIdL1(detid_.zendcap(), 1, detid_.station(), 1, detid_.chamber(), gem1.nPartitions());
            // const GEMDetId gemDetIdL2(detid_.zendcap(), 1, detid_.station(), 2, detid_.chamber(), gem2.nPartitions());
            // cout << "DetId = " << detid_.rawId() << ", gemDetIdL1 = " << gemDetIdL1.rawId() << ", gemDetIdL2 = " << gemDetIdL2.rawId() <<endl;
            // auto gp1 = match->gemDigis()->getGlobalPointPad(gemDetIdL1.rawId(), gem1);
            // auto gp2 = match->gemDigis()->getGlobalPointPad(gemDetIdL2.rawId(), gem2);
            // cout << "GP" << endl;
            // matchCscGEM1->FillGP(gp1);
            // matchCscGEM2->FillGP(gp2);
            // matchCscGEM1->FillGEMPad(gem1, gemDetIdL1.rawId(), digicount);
            // matchCscGEM2->FillGEMPad(gem2, gemDetIdL2.rawId(), digicount);
          }
          digicount++;
        }
      }
      if (DebugMode) cout << "Finished Loop over matchCscStubsA/CLCT (match->cscStubs()) " << endl;
      auto gemDigis_ = match->gemDigis();
      const auto& detidsDigi = gemDigis_->detIdsDigi();
      // if (detidsDigi.size() == 0) cout << "!!!! detidsDigi is empty" <<endl;
      // else cout << "~~~~ detidsDigi is not empty" <<endl;
      for (const auto& id : detidsDigi) {
        // if (gemDigis_->digisInDetId(id).size() == 0) cout  << "!!!! digisInDetId is empty while there are " << detidsDigi.size() << " detIds"<< endl;
        for (auto gemdigi : gemDigis_->digisInDetId(id) ){
          if (DebugMode) cout << "Starting processing a matchgemdigi" <<endl;
          auto gp = gemDigis_->getGlobalPointDigi(id, gemdigi);
          matchGemDigi->FillGP(gp);
          if (DebugMode) cout << "gp filled" <<endl;
          matchGemDigi->FillGEM(gemdigi,id,tp_index);
          if (DebugMode) cout << "gem saved" <<endl;
        }
      }
      if (DebugMode) cout << "Finished Loop over matchGemDigi (match->gemDigis()) " << endl;

      // matchGemPadDigiCluster
      auto gemPadDigiClusters_ = match->gemDigis();
      const auto& detIdsCluster = gemDigis_->detIdsCluster();
      for (const auto& id : detIdsCluster) {
        for (auto cl : gemDigis_->clustersInDetId(id)) {
          GEMPadDigi mid = GEMPadDigi(cl.pads()[cl.pads().size() / 2], cl.bx(), cl.station(), cl.nPartitions());
          auto gp = gemDigis_->getGlobalPointPad(id, mid);
          matchGemPadDigiCluster->FillGP(gp);
          matchGemPadDigiCluster->FillGEMPadDigiCluster(cl, id, tp_index);
        }
      }

      if (DebugMode) cout << "Finished Loop over matchGemPadDigiCluster" <<endl;
    } // End of Muon Loop

    else{
      m_matchmuon_pt->push_back(-999.);
      m_matchmuon_eta->push_back(-999.);
      m_matchmuon_phi->push_back(-999.);
      m_matchmuon_charge->push_back(-999);
      m_matchmuon_type->push_back(-999);
      m_matchmuon_quality->push_back(-999);
    }

    tp_index++;
  } // End of Tracking Particle Loop

  if (DebugMode) cout << "Finished TP Iteration" << endl<< endl;

  // Handle< BXVector<l1t::RegionalMuonCand> > emtfs;
  // Handle< BXVector<l1t::RegionalMuonCand> > omtfs;
  // Handle< BXVector<l1t::RegionalMuonCand> > bmtfs;

  // iEvent.getByToken(m_emtfToken,emtfs);
  // iEvent.getByToken(m_omtfToken,omtfs);
  // iEvent.getByToken(m_bmtfToken,bmtfs);

  // int nEMTF=0;
  // for (auto it = emtfs->begin(0); it != emtfs->end(0); it++){
  //   m_EMTF_muon_eta->push_back(it->hwEta()*0.010875);
  //   int globPhi=l1t::MicroGMTConfiguration::calcGlobalPhi(it->hwPhi(), it->trackFinderType(), it->processor());
  //   m_EMTF_muon_phi->push_back(globPhi*2*M_PI/576.);
  //   m_EMTF_muon_pt->push_back(it->hwPt()*0.5);
  //   if(!it->hwSignValid()) m_EMTF_muon_c->push_back(0);
  //   else{
  //     if(it->hwSign())   m_EMTF_muon_c->push_back(-1);
  //     else 		   m_EMTF_muon_c->push_back(1);
  //   }
  //   ++nEMTF;
  // }
  // m_EMTF_muon_n->push_back(nEMTF);

  // int nOMTF=0;
  // for (auto it = omtfs->begin(0); it != omtfs->end(0); it++){
  //   m_OMTF_muon_eta->push_back(it->hwEta()*0.010875);
  //   m_OMTF_muon_phi->push_back(l1t::MicroGMTConfiguration::calcGlobalPhi(it->hwPhi(), it->trackFinderType(), it->processor())*2*M_PI/576.);
  //   m_OMTF_muon_pt->push_back(it->hwPt()*0.5);
  //   if(!it->hwSignValid()) m_OMTF_muon_c->push_back(0);
  //   else{
  //     if(it->hwSign())   m_OMTF_muon_c->push_back(-1);
  //     else 		   m_OMTF_muon_c->push_back(1);
  //   }
  //   ++nOMTF;
  // }
  // m_OMTF_muon_n->push_back(nOMTF);
  //
  // int nBMTF=0;
  // for (auto it = bmtfs->begin(0); it != bmtfs->end(0); it++){
  //   m_BMTF_muon_eta->push_back(it->hwEta()*0.010875);
  //   m_BMTF_muon_phi->push_back(l1t::MicroGMTConfiguration::calcGlobalPhi(it->hwPhi(), it->trackFinderType(), it->processor())*2*M_PI/576.);
  //   m_BMTF_muon_pt->push_back(it->hwPt()*0.5);
  //   if(!it->hwSignValid()) m_BMTF_muon_c->push_back(0);
  //   else{
  //     if(it->hwSign())   m_BMTF_muon_c->push_back(-1);
  //     else 		   m_BMTF_muon_c->push_back(1);
  //   }
  //   ++nBMTF;
  // }
  // m_BMTF_muon_n->push_back(nBMTF);

  if (DebugMode) cout << "Finished regional Muons started allCscStubs (CSCCorrelatedLCTDigiCollection)" << endl;

  edm::Handle<CSCCorrelatedLCTDigiCollection> lctsH_;
  iEvent.getByToken(lctToken_, lctsH_);
  const CSCCorrelatedLCTDigiCollection& lcts = *lctsH_.product();
  int allCscStubs_index = 0;
  int digicount = 0;
  for (auto it = lcts.begin(); it != lcts.end(); ++it) {
    const auto& digivec = (*it).second;
    const CSCDetId& detid = (*it).first;
    int digi_index = 0;
    bool doprinta = Print_allCscStubs && Print_ALCT;
    bool doprintc = Print_allCscStubs && Print_CLCT;
    for (auto itdigi = digivec.first; itdigi != digivec.second; ++itdigi) {
      if (DebugMode) cout << " Started "<< digi_index << "th Digi for "<< allCscStubs_index << "th allCscStubsLCTs" <<endl;
      // cout << "allCSCStubs in endc stat ring:" << detid.endcap()<< " " << detid.station() << " " << detid.ring() << ", ";
      // cout << "Slope: " << (*itdigi).getSlope();
      if ((*itdigi).getQuality() == 6) cout << "Found a 6!" <<endl;
      auto gp = match->cscStubs()->getGlobalPosition(detid.rawId(), *itdigi);
      // auto gp2 = match->cscStubs()->getGlobalPosition2(detid.rawId(), *itdigi);
      // if (gp2.eta() != gp.eta() || gp2.phi() != gp.phi()) {
      //   cout << endl << "all CscStubsLCT GlobalPoint before change: eta = " << gp2.eta()<< " , phi = " << gp2.phi() << endl;
      //   cout         << "all CscStubsLCT GlobalPoint after  change: eta = " << gp.eta() << " , phi = " << gp.phi()  << endl;
      // }
      allCscStubsLCT->FillGP(gp);
      allCscStubsLCT->FillLCT(*itdigi,detid.rawId());
      // cout << ", Filled is: " << allCscStubsLCT->slope->back() <<endl;
      // allCscStubsALCTs
      if (DebugMode) cout << " Started allCSCStubsALCT"<< endl;
      const auto& alctDigi = (*itdigi).getALCT();
      std::vector< std::vector<unsigned short> > alcthits = alctDigi.getHits();
      if (doprinta) cout << "For allCscStubsALCTs; DetId: "<< detid.rawId() <<", Digi Index: " << digi_index << ", keywire: "<< alctDigi.getKeyWG()<<" , lctDigiWire = "<<(*itdigi).getKeyWG() << endl;
      if (alctDigi.getKeyWG() != (*itdigi).getKeyWG()) cout << "Inconsistence keywire in allCscStubsLCT: alctDigiWire = "<< alctDigi.getKeyWG() <<" , lctDigiWire = "<<(*itdigi).getKeyWG()<<endl;
      // int alctmultihit = SaveHitMatrix(alcthits, m_allCscStubsALCT_hit, m_allCscStubsALCT_position,doprinta,false);
      int alctmultihit = SaveHitMatrix(alcthits, allCscStubsALCT->hit, allCscStubsALCT->position,doprinta,false);
      if (doprinta && alctmultihit) cout << "ALCT Multihit found for allCscStubsALCT" << endl;
      allCscStubsALCT->FillALCT(alctDigi,detid.rawId());

      // allCscStubsCLCTs
      if (DebugMode) cout << " Started allCSCStubsCLCT"<< endl;
      const auto& clctDigi = (*itdigi).getCLCT();
      std::vector< std::vector<unsigned short> > clcthits = clctDigi.getHits();
      if (doprintc) cout << "For allCscStubsCLCTs; DetId: "<< detid.rawId() <<", Digi Index:" << digi_index << ", strip: "<< clctDigi.getStrip() <<" , lctDigiStrip = "<< (*itdigi).getStrip()<<endl;
      // if (clctDigi.getStrip() != (*itdigi).getStrip()) cout << "Inconsistence strip: clctDigiStrip = "<< clctDigi.getStrip() <<" , lctDigiStrip = "<< (*itdigi).getStrip() <<endl;
      if (clctDigi.getSlope() != (*itdigi).getSlope()) cout << "Inconsistence slope in allCscStubsLCT: clctDigiSlope = "<< clctDigi.getSlope() <<" , lctDigiSlope = "<< (*itdigi).getSlope() <<endl;
      // int clctmultihit = SaveHitMatrix(clcthits, m_allCscStubsCLCT_hit, m_allCscStubsCLCT_position,doprintc,true);
      int clctmultihit = SaveHitMatrix(clcthits, allCscStubsCLCT->hit, allCscStubsCLCT->position,doprintc,true);
      if (doprintc && clctmultihit) cout << "CLCT Multihit found for allCscStubsCLCT" << endl;
      allCscStubsCLCT->FillCLCT(clctDigi,detid.rawId());
      // GEMPad for allCscStubs
      if (DebugMode) cout << " Finishing a allCscStub, starting allCscStub GEMPads" <<endl;
      if (detid.ring() == 1 and (detid.station() == 1 or detid.station() == 2)) {
        bool matchl1(false), matchl2(false), printgem(false);
        if ((*itdigi).getGEM1().pad() != 255 || (*itdigi).getGEM2().pad() != 255 || (*itdigi).getGEM1().nPartitions() != 8 || (*itdigi).getGEM2().nPartitions() != 8) {
          printgem = true;
          cout << "In station " << detid.station() << " , ring " << detid.ring() << " , eta = " << gp.eta() << " ,  allCscStubsLCT rawId: " << detid.rawId() << endl;
        }
        const GEMDetId gemDetIdL1(detid.zendcap(), 1, detid.station(), 1, detid.chamber(), 0);
        for (const auto& p : match->gemDigis()->padsInChamber(gemDetIdL1.rawId())) {
          if (p == (*itdigi).getGEM1() && p.isValid()) {
            cout << "GEM1 is Valid" <<endl;
            // auto gp1 = match->gemDigis()->getGlobalPointPad(gemDetIdL1.rawId(),p);
            // allCscGEM1->FillGP(gp1);
            allCscGEM1->FillGEMPad(p,gemDetIdL1.rawId(),digicount);
            matchl1 = true;
            break;
          }
        }
        if (!matchl1) {
          allCscGEM1->FillGP0();
          allCscGEM1->FillGEMPad0(gemDetIdL1.rawId(),digicount);
        }
        // Check if matched to an GEM pad L2
        const GEMDetId gemDetIdL2(detid.zendcap(), 1, detid.station(), 2, detid.chamber(), 0);
        for (const auto& p : match->gemDigis()->padsInChamber(gemDetIdL2.rawId())) {
          if (p == (*itdigi).getGEM2() && p.isValid()) {
            cout << "GEM2 is Valid" <<endl;
            // auto gp2 = match->gemDigis()->getGlobalPointPad(gemDetIdL2.rawId(),p);
            // allCscGEM2->FillGP(gp2);
            allCscGEM2->FillGEMPad(p,gemDetIdL2.rawId(),digicount);
            matchl2 = true;
            break;
          }
        }
        if (!matchl2) {
          allCscGEM2->FillGP0();
          allCscGEM2->FillGEMPad0(gemDetIdL2.rawId(),digicount);
        }

        if (printgem) {
          cout << "GEM1 pad = " << (*itdigi).getGEM1().pad() << " , part = " << (*itdigi).getGEM1().nPartitions() << ((*itdigi).getGEM1().isValid()? ", Valid " : ", Invalid ") << (matchl1 ? ", Filled " : ", Not Filled") <<endl;
          cout << "GEM2 pad = " << (*itdigi).getGEM2().pad() << " , part = " << (*itdigi).getGEM2().nPartitions() << ((*itdigi).getGEM2().isValid()? ", Valid " : ", Invalid ") << (matchl2 ? ", Filled " : ", Not Filled") <<endl;
        }
        // const auto& gem1 = digi_.getGEM1();
        // const auto& gem2 = digi_.getGEM2();
        // cout << " gem1 " << (gem1.isValid() ? "Valid" : "Invalid") <<endl;
        // cout << " gem2 " << (gem2.isValid() ? "Valid" : "Invalid") <<endl;

        // cout << "Read" <<endl;
        // const GEMDetId gemDetIdL1(detid_.zendcap(), 1, detid_.station(), 1, detid_.chamber(), gem1.nPartitions());
        // const GEMDetId gemDetIdL2(detid_.zendcap(), 1, detid_.station(), 2, detid_.chamber(), gem2.nPartitions());
        // cout << "DetId = " << detid_.rawId() << ", gemDetIdL1 = " << gemDetIdL1.rawId() << ", gemDetIdL2 = " << gemDetIdL2.rawId() <<endl;
        // auto gp1 = match->gemDigis()->getGlobalPointPad(gemDetIdL1.rawId(), gem1);
        // auto gp2 = match->gemDigis()->getGlobalPointPad(gemDetIdL2.rawId(), gem2);
        // cout << "GP" << endl;
        // matchCscGEM1->FillGP(gp1);
        // matchCscGEM2->FillGP(gp2);
        // matchCscGEM1->FillGEMPad(gem1, gemDetIdL1.rawId(), digicount);
        // matchCscGEM2->FillGEMPad(gem2, gemDetIdL2.rawId(), digicount);
      }
      digi_index++;
      digicount++;
    }
    allCscStubs_index++;
  }

  if (DebugMode) cout <<endl<< "Finished allCscStubs, started allALCT (CSCALCTDigiCollection)" << endl;

  bool multihit = false;
  int digi_index = 0;
  // allALCTs
  edm::Handle<CSCALCTDigiCollection> alctsH_;
  iEvent.getByToken(alctToken_, alctsH_);
  const CSCALCTDigiCollection& alcts = *alctsH_.product();
  for (auto it = alcts.begin(); it != alcts.end(); ++it) {
    const auto& digivec = (*it).second;
    const auto& detid = (*it).first;
    for (auto itdigi = digivec.first; itdigi != digivec.second; ++itdigi) {
      std::vector< std::vector<unsigned short> > alcthits = ((*itdigi).getHits());
      if (Print_all && Print_ALCT) cout << "For allALCTs; DetId: "<< detid.rawId() <<", Digi Index:" << digi_index << ", keywire: "<< (*itdigi).getKeyWG()<< endl;
      // int alctmultihit = SaveHitMatrix(alcthits, m_allALCT_hit, m_allALCT_position,Print_all&&Print_ALCT,false);
      int alctmultihit = SaveHitMatrix(alcthits, allALCT->hit, allALCT->position,Print_all&&Print_ALCT,false);
      if (alctmultihit) multihit = true;
      if (Print_all && Print_ALCT && alctmultihit) cout << "ALCT Multihit found for allALCT" << endl;
      allALCT->FillALCT(*itdigi,detid.rawId());
      ++digi_index;
    }
  }

  if (DebugMode) cout << "Finished allALCT, started allCLCT" << endl;
  // allCLCTs
  digi_index = 0;
  edm::Handle<CSCCLCTDigiCollection> clctsH_;
  iEvent.getByToken(clctToken_, clctsH_);
  const CSCCLCTDigiCollection& clcts = *clctsH_.product();
  for (auto it = clcts.begin(); it != clcts.end(); ++it) {
    const auto& digivec = (*it).second;
    const auto& detid = (*it).first;
    for (auto itdigi = digivec.first; itdigi != digivec.second; ++itdigi) {
      std::vector< std::vector<unsigned short> > clcthits = ((*itdigi).getHits());
      if (Print_all && Print_CLCT) cout << "For allCLCTs; DetId: "<< detid.rawId() <<", Digi Index:" << digi_index<< ", strip: "<< (*itdigi).getStrip() << endl;
      // int clctmultihit = SaveHitMatrix(clcthits, m_allCLCT_hit, m_allCLCT_position,Print_all&&Print_CLCT,true);
      int clctmultihit = SaveHitMatrix(clcthits, allCLCT->hit, allCLCT->position,Print_all&&Print_CLCT,true);
      if (clctmultihit) multihit = true;
      if (Print_all && Print_CLCT && clctmultihit) cout << "CLCT Multihit found for allCLCT" << endl;
      allCLCT->FillCLCT(*itdigi,detid.rawId());
      ++digi_index;
    }
  }
  if (multihit) nEventMultiHitLayer++;
  // if (floor(double(m_allCLCT_hit->size()) / 6.0 ) != m_allCLCT_detId->size()) {
  //   throw std::runtime_error("Hit size mismatch Digi size");
  // }

  // 0 for allA/CLCT, 1 for allCscStubsA/CLCT, 2 for matchCscStubsA/CLCT
  if (Print_all) PrintHits(0);
  if (Print_allCscStubs) PrintHits(1);
  if (Print_matchCscStubs) PrintHits(2);

  if (DebugMode) cout << "Finished allCLCT, started GEMDigis" << endl;
  // All GEMDigis
  edm::Handle<GEMDigiCollection> gemDigisH_;
  iEvent.getByToken(gemDigiToken_,gemDigisH_);
  const GEMDigiCollection& gems = *gemDigisH_.product();
  for (auto it = gems.begin(); it != gems.end(); ++it) {
    const auto& digivec = (*it).second;
    const GEMDetId& detid = (*it).first;
    for (auto itdigi = digivec.first; itdigi != digivec.second; ++itdigi) {
      if (DebugMode) cout << "Starting processing a allgemdigi" <<endl;
      auto gp = match->gemDigis()->getGlobalPointDigi(detid, *itdigi);
      if (DebugMode) cout << "gp obtained" <<endl;
      allGemDigi->FillGP(gp);
      if (DebugMode) cout << "gp filled" <<endl;
      allGemDigi->FillGEM(*itdigi,detid.rawId());
      if (DebugMode) cout << "gem saved" <<endl;
    }
  }
  if (DebugMode) cout << "Finished GEMDigis, started GEMPadDigi" << endl;

  // GEMPadDigi
  edm::Handle<GEMPadDigiCollection> gemPadDigisH_;
  iEvent.getByToken(gemPadDigiToken_, gemPadDigisH_);
  const GEMPadDigiCollection& gempads = *gemPadDigisH_.product();
  for (auto it = gempads.begin(); it != gempads.end(); ++it) {
    const auto& digivec = (*it).second;
    const GEMDetId& detid = (*it).first;
    for (auto itdigi = digivec.first; itdigi != digivec.second; ++itdigi) {
      auto gp = match->gemDigis()->getGlobalPointPad(detid, *itdigi);
      gemPadDigi->FillGP(gp);
      gemPadDigi->FillGEMPad(*itdigi, detid.rawId());
    }
  }
  if (DebugMode) cout << "Finished GEMPadDigis, started GEMPadDigiClusters" <<endl;

  // GEMPadDigiCluster
  edm::Handle<GEMPadDigiClusterCollection> gemPadDigiClustersH_;
  iEvent.getByToken(gemPadDigiClusterToken_,gemPadDigiClustersH_);
  const GEMPadDigiClusterCollection& cls = *gemPadDigiClustersH_.product();
  for (auto it = cls.begin(); it != cls.end(); ++it) {
    const auto& clvec = (*it).second;
    const GEMDetId& detid = (*it).first;
    for (auto itcl = clvec.first; itcl != clvec.second; ++itcl) {
      GEMPadDigi mid = GEMPadDigi(itcl->pads()[itcl->pads().size() / 2], itcl->bx(), itcl->station(), itcl->nPartitions());
      auto gp = match->gemDigis()->getGlobalPointPad(detid, mid);
      allGemPadDigiCluster->FillGP(gp);
      allGemPadDigiCluster->FillGEMPadDigiCluster(*itcl, detid.rawId());
    }
  }
  if (DebugMode) cout << "Finished GEMPadDigiClusters, started filling the tree" <<endl;

  eventTree->Fill();
} // End of analyze


void NtupleMaker::PrintHits(int iset, int idigi) {
  vector<vector<int> > SavedallCscStubs{*(allCscStubsALCT->hit), *(allCscStubsALCT->position), *(allCscStubsALCT->detId), *(allCscStubsALCT->keywire), *(allCscStubsCLCT->hit), *(allCscStubsCLCT->position), *(allCscStubsCLCT->detId), *(allCscStubsCLCT->strip), *(allCscStubsLCT->detId), *(allCscStubsLCT->keywire), *(allCscStubsLCT->strip)};
  vector<vector<int> > SavedmatchCscStubs{*(matchCscStubsALCT->hit), *(matchCscStubsALCT->position), *(matchCscStubsALCT->detId), *(matchCscStubsALCT->keywire), *(matchCscStubsCLCT->hit), *(matchCscStubsCLCT->position), *(matchCscStubsCLCT->detId), *(matchCscStubsCLCT->strip), *(matchCscStubsLCT->detId), *(matchCscStubsLCT->keywire), *(matchCscStubsLCT->strip)};
  vector<vector<int> > Savedall{*(allALCT->hit), *(allALCT->position), *(allALCT->detId), *(allALCT->keywire), *(allCLCT->hit), *(allCLCT->position), *(allCLCT->detId), *(allCLCT->strip)};
  vector<string> printtitle{"hits","positions","detids","keywires","hits","positions","detids","strips","LCT Digi detids","LCT Digi keywires","LCT Digi strips"};
  vector<vector<int> > printinfo;
  cout <<endl;
  if (iset == 0) {
    cout << "Saved info in allALCT and allCLCT" <<endl;
    printinfo = Savedall;
  }
  else if (iset == 1){
    cout << "Saved info in allCscStubsALCT and allCscStubsCLCT" <<endl;
    printinfo = SavedallCscStubs;
  }
  else if (iset == 2){
    cout << "Saved info in matchCscStubsALCT and matchCscStubsCLCT" <<endl;
    printinfo = SavedmatchCscStubs;
  }

  if (idigi == -1) {
    for (unsigned iprint = 0; iprint < printinfo.size(); ++iprint) {
      if ((!Print_ALCT) && iprint < 4) continue;
      if ((!Print_CLCT) && iprint >= 4 && iprint < 8) continue;
      if (iprint == 0) cout << "Saved ALCT hits format:" <<endl;
      if (iprint == 4) cout << "Saved CLCT hits format:" <<endl;
      if (iprint == 8) cout << "Saved LCT format:" <<endl;
      cout << printtitle.at(iprint) << ": {";
      vector<int>& vals = printinfo.at(iprint);
      bool printbracket = false;
      if (printtitle.at(iprint) == "hits" || printtitle.at(iprint) == "positions") printbracket = true;
      for (unsigned ival = 0; ival < vals.size(); ++ival) {
        if (ival % 6 == 0 && printbracket) cout << "(";
        cout << vals.at(ival);
        if (ival % 6 == 5 && printbracket) cout << ")";
        if (ival != vals.size() - 1) cout <<", ";
      }
      cout <<"}"<<endl;
    }
  }

}

int NtupleMaker::SaveHitMatrix(std::vector< std::vector<unsigned short> > hits, std::vector<int>* b_hit, std::vector<int>* b_pos, bool doprint, bool isclct) {
  std::vector<int> tmp_hit;
  std::vector<int> tmp_pos;
  int multihit = 0; // number of layers that have more than 1 hits
  if (hits.size() != 6) {
    cout << "nLayer != 6" <<endl;
    int ExceptionCode = -3;
    for (unsigned itlayer = 0; itlayer < 6; ++itlayer) {
      b_hit->push_back(ExceptionCode);
      b_pos->push_back(ExceptionCode);
      if (doprint) {
        tmp_hit.push_back(ExceptionCode);
        tmp_pos.push_back(ExceptionCode);
      }
    }
    return ExceptionCode;
  }
  if (doprint) cout << "--- Printing Hit Matrix( "<< hits.size() << " * " << hits.at(0).size() << " ) ---" << endl;
  for (unsigned itlayer = 0; itlayer < hits.size() ; ++itlayer) {
    if (doprint) cout << "||layer: " << itlayer << " : ";
    int ncount = 0;
    bool haszero = false;
    bool printzero = false;
    for (unsigned itpos = 0; itpos < hits.at(itlayer).size(); ++itpos) {
      int hitval = hits.at(itlayer).at(itpos);
      if (doprint) cout << " " << hitval << ",";
      if (hitval == 0 && isclct) haszero = true;
      if (hitval != 65535 && hitval != 0) {
        if (ncount == 0) {
          b_hit->push_back(hitval);
          b_pos->push_back(itpos);
          if (doprint) {
            tmp_hit.push_back(hitval);
            tmp_pos.push_back(itpos);
          }
        }
        // Just keep the first hit if there are multiple
        // else {
        //
        //   int ExceptionCode = -2;
        //   b_hit->back() = ExceptionCode;
        //   b_pos->back() = ExceptionCode;
        //   if (doprint) {
        //     tmp_hit.back() = ExceptionCode;
        //     tmp_pos.back() = ExceptionCode;
        //   }
        // }
        ncount++;
      }
    }
    if (haszero && printzero) cout << "Emergency!!!!!! CLCT hits contains 0!!!!!!";
    if (ncount > 1) {
      ++multihit;
      // cout << " Layer with more than 1 hits: " << ncount + 1 << " hits in " <<endl;
      if (DebugMode) cout << " -- Layer: " << itlayer << " contains " << ncount << " hits" <<endl;
    }
    if (ncount == 0) {
      b_hit->push_back(-1);
      b_pos->push_back(-1);
      if (doprint) {
        tmp_hit.push_back(-1);
        tmp_pos.push_back(-1);
      }
    }
    if (doprint) cout << "||" << endl;
  }

  if (doprint) {
    cout << "--- End of Hit Matrix ---" <<endl;
    cout << "Saved Hit for this Digi:";
    for (unsigned iv = 0; iv < tmp_hit.size(); ++ iv) {
      cout <<" "<<tmp_hit.at(iv)<<",";
    }
    if (tmp_hit.size() != 6) cout << "This Digi does not has 6 layers!!!";
    cout <<endl;
    cout << "Saved Position for this Digi:";
    for (unsigned iv = 0; iv < tmp_pos.size(); ++ iv) {
      cout <<" "<<tmp_pos.at(iv)<<",";
    }
    if (tmp_hit.size() != 6) cout << "This Digi does not has 6 layers!!!";
    cout <<endl;
  }
  return multihit;
}

std::vector<int> NtupleMaker::IntsToBinary(int n) {
  std::vector<int> binary_;
  if ( n / 2 != 0) binary_ = IntsToBinary(n/2);
  binary_.push_back(n%2);
  return binary_;
}

GlobalPoint NtupleMaker::getGlobalPointDigi(unsigned int rawId, const GEMDigi& d) {
  GEMDetId gem_id(rawId);
  const LocalPoint& gem_lp = gemGeometry_->etaPartition(gem_id)->centreOfStrip(d.strip());
  const GlobalPoint& gem_gp = gemGeometry_->idToDet(gem_id)->surface().toGlobal(gem_lp);
  return gem_gp;
}

DEFINE_FWK_MODULE(NtupleMaker);
