//////////////////////////////////////////////////////////////////////
//               Digi handler for GEM CSC analysis                //
//               Author: Sifu Luo                                   //
//               sifuluo@tamu.edu                                   //
//////////////////////////////////////////////////////////////////////

#ifndef TREEDIGI_CC
#define TREEDIGI_CC

// STD
#include <vector>
#include <cmath>

// ROOT
#include <TROOT.h>
#include <TTree.h>
#include <TString.h>

// GEMCode
#include "GEMCode/GEMValidation/interface/MatcherManager.h"
#include "GEMCode/GEMValidation/test/GEMTool.cc"

using namespace std;

class TreeDigi{
public:
  TreeDigi() {

  }

  void Init(TTree* evttree, TString name_, int data_type_, GEMTool* t, bool match_ = false, bool reader_ = false) {
    name = name_;
    data_type = data_type_; // 0 for LCT, 1 for ALCT, 2 for CLCT, 3 for GemDigi
    IsMatched = match_;
    reader = reader_;
    GEMConverter = t;
    // LCT
    if (data_type == 0) {
      InitGP(evttree);
      InitDet(evttree);
      bend = new std::vector<int>;
      pattern = new std::vector<int>;
      slope = new std::vector<int>;
      quality = new std::vector<int>;
      keywire = new std::vector<int>;
      strip = new std::vector<int>;
      strip8 = new std::vector<int>;
      valid = new std::vector<bool>;
      type = new std::vector<int>;
      GEM1pad = new std::vector<int>;
      GEM1strip = new std::vector<int>;
      GEM1strip8 = new std::vector<int>;
      GEM1strip_me1a = new std::vector<int>;
      GEM1strip8_me1a = new std::vector<int>;
      GEM1keywire_min = new std::vector<int>;
      GEM1keywire_max = new std::vector<int>;
      GEM1roll = new std::vector<int>;
      GEM1part = new std::vector<int>;
      GEM2pad = new std::vector<int>;
      GEM2strip = new std::vector<int>;
      GEM2strip8 = new std::vector<int>;
      GEM2strip_me1a = new std::vector<int>;
      GEM2strip8_me1a = new std::vector<int>;
      GEM2keywire_min = new std::vector<int>;
      GEM2keywire_max = new std::vector<int>;
      GEM2roll = new std::vector<int>;
      GEM2part = new std::vector<int>;
      evttree->Branch(name+"_bend", &bend);
      evttree->Branch(name+"_pattern", &pattern);
      evttree->Branch(name+"_slope", &slope);
      evttree->Branch(name+"_quality", &quality);
      evttree->Branch(name+"_keywire", &keywire);
      evttree->Branch(name+"_strip", &strip);
      evttree->Branch(name+"_strip8", &strip8);
      evttree->Branch(name+"_valid", &valid);
      evttree->Branch(name+"_type", &type);
      evttree->Branch(name+"_GEM1pad", &GEM1pad);
      evttree->Branch(name+"_GEM1strip", &GEM1strip);
      evttree->Branch(name+"_GEM1strip8", &GEM1strip8);
      evttree->Branch(name+"_GEM1strip_me1a", &GEM1strip_me1a);
      evttree->Branch(name+"_GEM1strip8_me1a", &GEM1strip8_me1a);
      evttree->Branch(name+"_GEM1keywire_min", &GEM1keywire_min);
      evttree->Branch(name+"_GEM1keywire_max", &GEM1keywire_max);
      evttree->Branch(name+"_GEM1roll", &GEM1roll);
      evttree->Branch(name+"_GEM1part", &GEM1part);
      evttree->Branch(name+"_GEM2pad", &GEM2pad);
      evttree->Branch(name+"_GEM2strip", &GEM2strip);
      evttree->Branch(name+"_GEM2strip8", &GEM2strip8);
      evttree->Branch(name+"_GEM2strip_me1a", &GEM2strip_me1a);
      evttree->Branch(name+"_GEM2strip8_me1a", &GEM2strip8_me1a);
      evttree->Branch(name+"_GEM2keywire_min", &GEM2keywire_min);
      evttree->Branch(name+"_GEM2keywire_max", &GEM2keywire_max);
      evttree->Branch(name+"_GEM2roll", &GEM2roll);
      evttree->Branch(name+"_GEM2part", &GEM2part);
      if (IsMatched) {
        matchIndex = new std::vector<int>;
        evttree->Branch(name+"_matchIndex", &matchIndex);
      }
    }

    // ALCT
    else if (data_type == 1) {
      InitDet(evttree);
      keywire = new std::vector<int>;
      hit = new std::vector<int>;
      position = new std::vector<int>;
      valid = new std::vector<bool>;
      evttree->Branch(name+"_keywire", &keywire);
      evttree->Branch(name+"_hit", &hit);
      evttree->Branch(name+"_position", &position);
      evttree->Branch(name+"_valid", &valid);
    }

    //CLCT
    else if (data_type == 2) {
      InitDet(evttree);
      strip = new std::vector<int>;
      strip8 = new std::vector<int>;
      hit = new std::vector<int>;
      position = new std::vector<int>;
      valid = new std::vector<bool>;
      bend = new std::vector<int>;
      pattern = new std::vector<int>;
      slope = new std::vector<int>;
      evttree->Branch(name+"_strip", &strip);
      evttree->Branch(name+"_strip8", &strip8);
      evttree->Branch(name+"_hit", &hit);
      evttree->Branch(name+"_position", &position);
      evttree->Branch(name+"_valid", &valid);
      evttree->Branch(name+"_bend", &bend);
      evttree->Branch(name+"_pattern", &pattern);
      evttree->Branch(name+"_slope", &slope);
    }

    // GEMDigi
    else if (data_type == 3) {
      InitGP(evttree);
      InitDet(evttree);
      strip = new std::vector<int>;
      evttree->Branch(name+"_strip", &strip);
      if (IsMatched) {
        matchIndex = new std::vector<int>;
        evttree->Branch(name+"_matchIndex", &matchIndex);
      }
    }

    // GEMPad
    else if (data_type == 4) {
      InitGP(evttree);
      InitDet(evttree);
      pad = new std::vector<int>;
      strip = new std::vector<int>;
      strip8 = new std::vector<int>;
      strip_me1a = new std::vector<int>;
      strip8_me1a = new std::vector<int>;
      keywire_min = new std::vector<int>;
      keywire_max = new std::vector<int>;
      part = new std::vector<int>;
      evttree->Branch(name+"_pad", &pad);
      evttree->Branch(name+"_strip", &strip);
      evttree->Branch(name+"_strip8", &strip8);
      evttree->Branch(name+"_strip_me1a", &strip_me1a);
      evttree->Branch(name+"_strip8_me1a", &strip8_me1a);
      evttree->Branch(name+"_keywire_min", &keywire_min);
      evttree->Branch(name+"_keywire_max", &keywire_max);
      evttree->Branch(name+"_part", &part);
      if (IsMatched) {
        matchIndex = new std::vector<int>;
        evttree->Branch(name+"_matchIndex", &matchIndex);
      }
    }

    //SimHit
    else if (data_type == 5) {
      if (!IsMatched) cout << name << " not matched? SimHits are always matched!" << endl;
      InitGP(evttree);
      InitDet(evttree);
      matchIndex = new std::vector<int>;
      evttree->Branch(name+"_matchIndex", &matchIndex);
    }

    // TP
    else if (data_type == 6) {
      if (IsMatched) cout << name << " Matched? TP can't match to anything" << endl;
      pt           = new std::vector<float>;
      eta          = new std::vector<float>;
      phi          = new std::vector<float>;
      dxy          = new std::vector<float>;
      d0           = new std::vector<float>;
      z0           = new std::vector<float>;
      d0_prod      = new std::vector<float>;
      z0_prod      = new std::vector<float>;
      pdgid        = new std::vector<int>;
      nmatch       = new std::vector<int>;
      nloosematch  = new std::vector<int>;
      nstub        = new std::vector<int>;
      eventid      = new std::vector<int>;
      charge       = new std::vector<int>;
      evttree->Branch(name+"_pt",          &pt);
      evttree->Branch(name+"_eta",         &eta);
      evttree->Branch(name+"_phi",         &phi);
      evttree->Branch(name+"_dxy",         &dxy);
      evttree->Branch(name+"_d0",          &d0);
      evttree->Branch(name+"_z0",          &z0);
      evttree->Branch(name+"_d0_prod",     &d0_prod);
      evttree->Branch(name+"_z0_prod",     &z0_prod);
      evttree->Branch(name+"_pdgid",       &pdgid);
      evttree->Branch(name+"_nmatch",      &nmatch);
      evttree->Branch(name+"_nloosematch", &nloosematch);
      evttree->Branch(name+"_nstub",       &nstub);
      evttree->Branch(name+"_eventid",     &eventid);
      evttree->Branch(name+"_charge",      &charge);
    }

    // GEMPadDigiCluster
    else if (data_type == 7) {
      InitGP(evttree);
      InitDet(evttree);
      pads = new std::vector<int>;
      pad = new std::vector<int>;
      strip = new std::vector<int>;
      strip8 = new std::vector<int>;
      strip_me1a = new std::vector<int>;
      strip8_me1a = new std::vector<int>;
      keywire_min = new std::vector<int>;
      keywire_max = new std::vector<int>;
      part = new std::vector<int>;
      len = new std::vector<int>;
      evttree->Branch(name+"_pads", &pads);
      evttree->Branch(name+"_pad", &pad);
      evttree->Branch(name+"_strip", &strip);
      evttree->Branch(name+"_strip8", &strip8);
      evttree->Branch(name+"_strip_me1a", &strip_me1a);
      evttree->Branch(name+"_strip8_me1a", &strip8_me1a);
      evttree->Branch(name+"_keywire_min", &keywire_min);
      evttree->Branch(name+"_keywire_max", &keywire_max);
      evttree->Branch(name+"_part", &part);
      evttree->Branch(name+"_len",&len);
      if (IsMatched) {
        matchIndex = new std::vector<int>;
        evttree->Branch(name+"_matchIndex", &matchIndex);
      }
    }
  }

  void Init(TTree* evttree, TString name_, TString data_type_st, GEMTool* t, bool match_ = false) {
    int data_type_ = -1;
    if (data_type_st == "LCT") data_type_ = 0;
    else if (data_type_st == "ALCT") data_type_ = 1;
    else if (data_type_st == "CLCT") data_type_ = 2;
    else if (data_type_st == "GEM") data_type_ = 3;
    else if (data_type_st == "GEMPad") data_type_ = 4;
    else if (data_type_st == "SimHit") data_type_ = 5;
    else if (data_type_st == "TP") data_type_ = 6;
    else if (data_type_st == "GEMPadDigiCluster") data_type_ = 7;


    if (data_type_ == -1) cout << "Wrong Data Type input for " << name_ << " as " << data_type_st << endl;
    else Init(evttree, name_, data_type_, t, match_);
  }

  void Reset() {
    if (data_type == 0) { //LCT
      ResetGP();
      ResetDet();
      bend->clear();
      pattern->clear();
      slope->clear();
      quality->clear();
      keywire->clear();
      strip->clear();
      strip8->clear();
      valid->clear();
      type->clear();
      GEM1pad->clear();
      GEM1strip->clear();
      GEM1strip8->clear();
      GEM1strip_me1a->clear();
      GEM1strip8_me1a->clear();
      GEM1keywire_min->clear();
      GEM1keywire_max->clear();
      GEM1roll->clear();
      GEM1part->clear();
      GEM2pad->clear();
      GEM2strip->clear();
      GEM2strip8->clear();
      GEM2strip_me1a->clear();
      GEM2strip8_me1a->clear();
      GEM2keywire_min->clear();
      GEM2keywire_max->clear();
      GEM2roll->clear();
      GEM2part->clear();
      if (IsMatched) {
        matchIndex->clear();
      }
    }
    else if (data_type == 1) { //ALCT
      ResetDet();
      keywire->clear();
      hit->clear();
      position->clear();
      valid->clear();
    }
    else if (data_type == 2) { //CLCT
      ResetDet();
      strip->clear();
      strip8->clear();
      hit->clear();
      position->clear();
      valid->clear();
      bend->clear();
      pattern->clear();
      slope->clear();
    }
    else if (data_type == 3) { //GEMDigi
      ResetGP();
      ResetDet();
      strip->clear();
      if (IsMatched) {
        matchIndex->clear();
      }
    }
    else if (data_type == 4) { //GEMPad
      ResetGP();
      ResetDet();
      pad->clear();
      strip->clear();
      strip8->clear();
      strip_me1a->clear();
      strip8_me1a->clear();
      keywire_min->clear();
      keywire_max->clear();
      part->clear();
      if (IsMatched) {
        matchIndex->clear();
      }
    }
    else if (data_type == 5) { //SimHit
      ResetGP();
      ResetDet();
      matchIndex->clear();
    }
    else if (data_type == 6) { // TP
      pt->clear();
      eta->clear();
      phi->clear();
      dxy->clear();
      d0->clear();
      z0->clear();
      d0_prod->clear();
      z0_prod->clear();
      pdgid->clear();
      nmatch->clear();
      nloosematch->clear();
      nstub->clear();
      eventid->clear();
      charge->clear();
    }
    else if (data_type == 7) { // GEMPadDigiCluster
      ResetGP();
      ResetDet();
      pads->clear();
      pad->clear();
      strip->clear();
      strip8->clear();
      strip_me1a->clear();
      strip8_me1a->clear();
      keywire_min->clear();
      keywire_max->clear();
      part->clear();
      len->clear();
      if (IsMatched) {
        matchIndex->clear();
      }
    }
  }

  void SetGEMConverter(GEMTool* t) {
    GEMConverter = t;
  }

  void ResetGP() {
    eta->clear();
    phi->clear();
    z->clear();
    r->clear();
  }

  void InitGP(TTree* evttree) {
    phi = new std::vector<float>;
    eta = new std::vector<float>;
    z = new std::vector<float>;
    r = new std::vector<float>;
    evttree->Branch(name+"_phi", &phi);
    evttree->Branch(name+"_eta", &eta);
    evttree->Branch(name+"_z", &z);
    evttree->Branch(name+"_r", &r);
  }

  void FillGP(GlobalPoint gp) {
    phi->push_back(gp.phi());
    eta->push_back(gp.eta());
    z->push_back(gp.z());
    r->push_back(gp.perp());
  }

  void FillGP0() {
    phi->push_back(0.);
    eta->push_back(0.);
    z->push_back(0.);
    r->push_back(0.);
  }

  void ResetDet() {
    detId->clear();
    zendcap->clear();
    ring->clear();
    station->clear();
    layer->clear();
    chamber->clear();
    roll->clear();
  }

  void InitDet(TTree* evttree) {
    detId = new std::vector<int>;
    zendcap = new std::vector<int>;
    ring = new std::vector<int>;
    station = new std::vector<int>;
    layer = new std::vector<int>;
    chamber = new std::vector<int>;
    roll = new std::vector<int>;
    evttree->Branch(name+"_detId", &detId);
    evttree->Branch(name+"_zendcap", &zendcap);
    evttree->Branch(name+"_ring", &ring);
    evttree->Branch(name+"_station", &station);
    evttree->Branch(name+"_layer", &layer);
    evttree->Branch(name+"_chamber", &chamber);
    evttree->Branch(name+"_roll", &roll);
  }

  void FillDet(CSCDetId id) {
    detId->push_back(id.rawId());
    zendcap->push_back(id.zendcap());
    ring->push_back(id.ring());
    station->push_back(id.station());
    layer->push_back(id.layer());
    chamber->push_back(id.chamber());
  }

  void FillDet(GEMDetId id) {
    detId->push_back(id.rawId());
    zendcap->push_back(id.region());
    ring->push_back(id.ring());
    station->push_back(id.station());
    layer->push_back(id.layer());
    chamber->push_back(id.chamber());
    roll->push_back(id.roll());
  }

  void FillLCT(CSCCorrelatedLCTDigi lct, CSCDetId id, int tp_index = -1) {
    FillDet(id);
    bend->push_back(lct.getBend());
    pattern->push_back(lct.getRun3Pattern());
    slope->push_back(lct.getSlope());
    // cout << "lct reads slope as " << lct.getSlope()<<endl;
    quality->push_back(lct.getQuality());
    keywire->push_back(lct.getKeyWG());
    strip->push_back(lct.getStrip());
    // strip8->push_back(lct.getFractionalStrip(8));
    strip8->push_back(lct.getStrip(8));
    valid->push_back(lct.isValid());
    type->push_back(lct.getType());
    if (tp_index != -1) matchIndex->push_back(tp_index);
  }

  void FillLCTGEMPads(CSCCorrelatedLCTDigi lct, GEMDetId id, TString info = "") {
    const GEMDetId id1(id.region(), id.ring(), id.station(), 1, id.chamber(), id.roll());
    const GEMDetId id2(id.region(), id.ring(), id.station(), 2, id.chamber(), id.roll());
    vector<int> gem1{-1,-1,-1,-1,-1,-1}, gem2{-1,-1,-1,-1,-1,-1};
    if (lct.getGEM1().isValid() || lct.getGEM2().isValid()) {
      gem1 = GEMConverter->ConvertPad(lct.getGEM1(), id1, info);
      gem2 = GEMConverter->ConvertPad(lct.getGEM2(), id2, info);
    }
    GEM1pad->push_back(lct.getGEM1().pad());
    GEM1strip->push_back(gem1[0]);
    GEM1strip8->push_back(gem1[1]);
    GEM1strip_me1a->push_back(gem1[2]);
    GEM1strip8_me1a->push_back(gem1[3]);
    GEM1keywire_min->push_back(gem1[4]);
    GEM1keywire_max->push_back(gem1[5]);
    GEM1roll->push_back(id.roll());
    GEM1part->push_back(lct.getGEM1().nPartitions());
    GEM2pad->push_back(lct.getGEM2().pad());
    GEM2strip->push_back(gem2[0]);
    GEM2strip8->push_back(gem2[1]);
    GEM2strip_me1a->push_back(gem2[2]);
    GEM2strip8_me1a->push_back(gem2[3]);
    GEM2keywire_min->push_back(gem2[4]);
    GEM2keywire_max->push_back(gem2[5]);
    GEM2roll->push_back(id.roll());
    GEM2part->push_back(lct.getGEM2().nPartitions());
  }

  void FillLCTGEMPads0(CSCCorrelatedLCTDigi lct){
    GEM1pad->push_back(255);
    GEM1strip->push_back(-1);
    GEM1strip8->push_back(-1);
    GEM1strip_me1a->push_back(-1);
    GEM1strip8_me1a->push_back(-1);
    GEM1keywire_min->push_back(-1);
    GEM1keywire_max->push_back(-1);
    GEM1roll->push_back(-1);
    GEM1part->push_back(0);
    GEM2pad->push_back(255);
    GEM2strip->push_back(-1);
    GEM2strip8->push_back(-1);
    GEM2strip_me1a->push_back(-1);
    GEM2strip8_me1a->push_back(-1);
    GEM2keywire_min->push_back(-1);
    GEM2keywire_max->push_back(-1);
    GEM2roll->push_back(-1);
    GEM2part->push_back(0);
  }

  void FillALCT(CSCALCTDigi alct, CSCDetId id, int tp_index = -1) {
    FillDet(id);
    keywire->push_back(alct.getKeyWG());
    valid->push_back(alct.isValid());
    if (tp_index != -1) matchIndex->push_back(tp_index);
  }

  void FillCLCT(CSCCLCTDigi clct, CSCDetId id, int tp_index = -1) {
    FillDet(id);
    strip->push_back(clct.getStrip());
    // strip8->push_back(clct.getFractionalStrip(8));
    strip8->push_back(clct.getKeyStrip(8));
    valid->push_back(clct.isValid());
    bend->push_back(clct.getBend());
    pattern->push_back(clct.getRun3Pattern());
    slope->push_back(clct.getSlope());
    if (tp_index != -1) matchIndex->push_back(tp_index);
  }

  void FillGEM(GEMDigi gemdigi, GEMDetId id, int tp_index = -1) {
    FillDet(id);
    strip->push_back(gemdigi.strip());
    if (tp_index != -1) matchIndex->push_back(tp_index);
  }

  void FillGEMPad(GEMPadDigi gempad, GEMDetId id, int csc_index = -1, TString info = "") {
    FillDet(id);
    vector<int> gem = GEMConverter->ConvertPad(gempad, id, info);
    pad->push_back(gempad.pad());
    strip->push_back(gem[0]);
    strip8->push_back(gem[1]);
    strip_me1a->push_back(gem[2]);
    strip8_me1a->push_back(gem[3]);
    keywire_min->push_back(gem[4]);
    keywire_max->push_back(gem[5]);
    part->push_back(gempad.nPartitions());
    if (csc_index != -1) matchIndex->push_back(csc_index);
  }

  void FillGEMPad0(GEMDetId id, int csc_index = -1) {
    FillDet(id);
    pad->push_back(255);
    strip->push_back(-1);
    strip8->push_back(-1);
    strip_me1a->push_back(-1);
    strip8_me1a->push_back(-1);
    keywire_min->push_back(-1);
    keywire_max->push_back(-1);
    part->push_back(0);
    if (csc_index != -1) matchIndex->push_back(csc_index);
  }

  void FillSimHit(GEMDetId id, int matchIndex_) {
    FillDet(id);
    matchIndex->push_back(matchIndex_);
  }

  void FillSimHit(CSCDetId id, int matchIndex_) {
    FillDet(id);
    matchIndex->push_back(matchIndex_);
  }

  void FillGEMPadDigiCluster(GEMPadDigiCluster cluster, GEMDetId id, int tp_index = -1, TString info = "") {
    FillDet(id);
    for (uint16_t pad_ : cluster.pads()) pads->push_back(pad_);
    int midpad = cluster.pads()[cluster.pads().size() / 2];
    vector<int> gem = GEMConverter->ConvertPad(midpad,id, info);
    len->push_back(cluster.pads().size());
    pad->push_back(midpad);
    strip->push_back(gem[0]);
    strip8->push_back(gem[1]);
    strip_me1a->push_back(gem[2]);
    strip8_me1a->push_back(gem[3]);
    keywire_min->push_back(gem[4]);
    keywire_max->push_back(gem[5]);
    part->push_back(cluster.nPartitions());
    if (tp_index != -1) matchIndex->push_back(tp_index);
  }

  void FillMatch(int tp_index) {
    matchIndex->push_back(tp_index);
  }

  std::vector<int>*   matchIndex;

  //Matrix extraction
  std::vector<int>*   hit;
  std::vector<int>*   position;

  // GP
  std::vector<float>* phi;
  std::vector<float>* eta;
  std::vector<float>* z;
  std::vector<float>* r;

  //Digi
  std::vector<int>*   bend;
  std::vector<int>*   pattern;
  std::vector<int>*   slope;
  std::vector<int>*   quality;
  std::vector<int>*   keywire;
  std::vector<int>*   strip;
  std::vector<int>*   strip8;
  std::vector<bool>*  valid;
  std::vector<int>*   type;
  std::vector<int>*   GEM1pad;
  std::vector<int>*   GEM1strip;
  std::vector<int>*   GEM1strip8;
  std::vector<int>*   GEM1strip_me1a;
  std::vector<int>*   GEM1strip8_me1a;
  std::vector<int>*   GEM1keywire_min;
  std::vector<int>*   GEM1keywire_max;
  std::vector<int>*   GEM1roll;
  std::vector<int>*   GEM1part;
  std::vector<int>*   GEM2pad;
  std::vector<int>*   GEM2strip;
  std::vector<int>*   GEM2strip8;
  std::vector<int>*   GEM2strip_me1a;
  std::vector<int>*   GEM2strip8_me1a;
  std::vector<int>*   GEM2keywire_min;
  std::vector<int>*   GEM2keywire_max;
  std::vector<int>*   GEM2roll;
  std::vector<int>*   GEM2part;
  std::vector<int>*   part;
  std::vector<int>*   pad;
  std::vector<int>*   strip_me1a;
  std::vector<int>*   strip8_me1a;
  std::vector<int>*   keywire_min;
  std::vector<int>*   keywire_max;
  std::vector<int>*   pads;
  std::vector<int>*   len;

  //TP
  std::vector<float>* pt;
  std::vector<float>* dxy;
  std::vector<float>* d0;
  std::vector<float>* z0;
  std::vector<float>* d0_prod;
  std::vector<float>* z0_prod;
  std::vector<int>* pdgid;
  std::vector<int>* nmatch;
  std::vector<int>* nloosematch;
  std::vector<int>* nstub;
  std::vector<int>* eventid;
  std::vector<int>* charge;

  //DetInfo
  // GEMDetId int region(+/-1), int ring, int station, int layer, int chamber, int ieta(roll)
  // CSCDetId int iendcap, int istation, int iring, int ichamber, int ilayer
  std::vector<int>* detId;
  std::vector<int>* zendcap;
  std::vector<int>* ring;
  std::vector<int>* station;
  std::vector<int>* layer;
  std::vector<int>* chamber;
  std::vector<int>* roll;


private:
  TString name;
  int data_type;
  bool IsMatched;
  bool reader;
  GEMTool* GEMConverter;
};

#endif
