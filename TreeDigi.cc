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

using namespace std;

class TreeDigi{
public:
  TreeDigi() {

  }

  void Init(TTree* evttree, TString name_, int data_type_, bool match_ = false, bool reader_ = false) {
    name = name_;
    data_type = data_type_; // 0 for LCT, 1 for ALCT, 2 for CLCT, 3 for GemDigi
    IsMatched = match_;
    reader = reader_;
    // LCT
    if (data_type == 0) {
      InitGP(evttree);
      bend = new std::vector<int>;
      pattern = new std::vector<int>;
      slope = new std::vector<int>;
      quality = new std::vector<int>;
      detId = new std::vector<int>;
      keywire = new std::vector<int>;
      strip = new std::vector<int>;
      strip8 = new std::vector<int>;
      valid = new std::vector<bool>;
      type = new std::vector<int>;
      evttree->Branch(name+"_bend", &bend);
      evttree->Branch(name+"_pattern", &pattern);
      evttree->Branch(name+"_slope", &slope);
      evttree->Branch(name+"_quality", &quality);
      evttree->Branch(name+"_detId", &detId);
      evttree->Branch(name+"_keywire", &keywire);
      evttree->Branch(name+"_strip", &strip);
      evttree->Branch(name+"_strip8", &strip8);
      evttree->Branch(name+"_valid", &valid);
      evttree->Branch(name+"_type", &type);
      if (IsMatched) {
        matchTp = new std::vector<int>;
        evttree->Branch(name+"_matchTp", &matchTp);
      }
    }

    // ALCT
    else if (data_type == 1) {
      detId = new std::vector<int>;
      keywire = new std::vector<int>;
      hit = new std::vector<int>;
      position = new std::vector<int>;
      valid = new std::vector<bool>;
      evttree->Branch(name+"_detId", &detId);
      evttree->Branch(name+"_keywire", &keywire);
      evttree->Branch(name+"_hit", &hit);
      evttree->Branch(name+"_position", &position);
      evttree->Branch(name+"_valid", &valid);
    }

    //CLCT
    else if (data_type == 2) {
      detId = new std::vector<int>;
      strip = new std::vector<int>;
      strip8 = new std::vector<int>;
      hit = new std::vector<int>;
      position = new std::vector<int>;
      valid = new std::vector<bool>;
      bend = new std::vector<int>;
      pattern = new std::vector<int>;
      slope = new std::vector<int>;
      evttree->Branch(name+"_detId", &detId);
      evttree->Branch(name+"_strip", &strip);
      evttree->Branch(name+"_strip8", &strip8);
      evttree->Branch(name+"_hit", &hit);
      evttree->Branch(name+"_position", &position);
      evttree->Branch(name+"_valid", &valid);
      evttree->Branch(name+"_bend", &bend);
      evttree->Branch(name+"_pattern", &pattern);
      evttree->Branch(name+"_slope", &slope);
    }

    // GEM
    else if (data_type == 3) {
      InitGP(evttree);
      detId = new std::vector<int>;
      strip = new std::vector<int>;
      evttree->Branch(name+"_detId", &detId);
      evttree->Branch(name+"_strip", &strip);
      if (IsMatched) {
        matchTp = new std::vector<int>;
        evttree->Branch(name+"_matchTp", &matchTp);
      }
    }

    // GEMPad
    else if (data_type == 4) {
      InitGP(evttree);
      detId = new std::vector<int>;
      pad = new std::vector<int>;
      part = new std::vector<int>;
      evttree->Branch(name+"_detId", &detId);
      evttree->Branch(name+"_pad", &pad);
      evttree->Branch(name+"_part", &part);
      if (IsMatched) {
        matchCSC = new std::vector<int>;
        evttree->Branch(name+"_matchCSC", &matchCSC);
      }
    }

    //SimHit
    else if (data_type == 5) {
      if (!IsMatched) cout << name << " not matched? SimHits are always matched!" << endl;
      InitGP(evttree);
      station = new std::vector<int>;
      matchTp = new std::vector<int>;
      evttree->Branch(name+"_station", &station);
      evttree->Branch(name+"_matchTp", &matchTp);
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
      detId = new std::vector<int>;
      pads = new std::vector<int>;
      part = new std::vector<int>;
      len = new std::vector<int>;
      evttree->Branch(name+"_detId", &detId);
      evttree->Branch(name+"_pads", &pads);
      evttree->Branch(name+"_part", &part);
      evttree->Branch(name+"_len",&len);
      if (IsMatched) {
        matchTp = new std::vector<int>;
        evttree->Branch(name+"_matchTp", &matchTp);
      }
    }
  }



  void Init(TTree* evttree, TString name_, TString data_type_st, bool match_ = false) {
    int data_type_ = -1;
    if (data_type_st == "LCT") data_type_ = 0;
    else if (data_type_st == "ALCT") data_type_ = 1;
    else if (data_type_st == "CLCT") data_type_ = 2;
    else if (data_type_st == "GEM") data_type_ = 3;
    else if (data_type_st == "GEMPad") data_type_ = 4;
    else if (data_type_st == "SimHit") data_type_ = 5;
    else if (data_type_st == "TP") data_type_ = 6;
    else if (data_type_st == "GEMPadDigiCluster") data_type_ = 7;
    // else if (data_type_st == "MatchMuon") data_type_ = 8;


    if (data_type_ == -1) cout << "Wrong Data Type input for " << name_ << " as " << data_type_st << endl;
    else Init(evttree, name_, data_type_, match_);
  }

  void Reset() {
    if (data_type == 0) {
      ResetGP();
      bend->clear();
      pattern->clear();
      slope->clear();
      quality->clear();
      detId->clear();
      keywire->clear();
      strip->clear();
      strip8->clear();
      valid->clear();
      type->clear();
      if (IsMatched) {
        matchTp->clear();
      }
    }
    else if (data_type == 1) {
      detId->clear();
      keywire->clear();
      hit->clear();
      position->clear();
      valid->clear();
    }
    else if (data_type == 2) {
      detId->clear();
      strip->clear();
      strip8->clear();
      hit->clear();
      position->clear();
      valid->clear();
      bend->clear();
      pattern->clear();
      slope->clear();
    }
    else if (data_type == 3) {
      ResetGP();
      detId->clear();
      strip->clear();
      if (IsMatched) {
        matchTp->clear();
      }
    }
    else if (data_type == 4) {
      ResetGP();
      detId->clear();
      pad->clear();
      part->clear();
      if (IsMatched) {
        matchCSC->clear();
      }
    }
    else if (data_type == 5) {
      ResetGP();
      station->clear();
      matchTp->clear();
    }
    else if (data_type == 6) {
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
    else if (data_type == 7) {
      ResetGP();
      detId->clear();
      pads->clear();
      part->clear();
      len->clear();
      if (IsMatched) {
        matchTp->clear();
      }
    }
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

  void FillLCT(CSCCorrelatedLCTDigi lct, int rawid, int tp_index = -1) {
    bend->push_back(lct.getBend());
    pattern->push_back(lct.getRun3Pattern());
    slope->push_back(lct.getSlope());
    quality->push_back(lct.getQuality());
    keywire->push_back(lct.getKeyWG());
    strip->push_back(lct.getStrip());
    strip8->push_back(lct.getFractionalStrip(8));
    valid->push_back(lct.isValid());
    type->push_back(lct.getType());
    detId->push_back(rawid);
    if (tp_index != -1) matchTp->push_back(tp_index);
  }

  void FillALCT(CSCALCTDigi alct, int rawid, int tp_index = -1) {
    keywire->push_back(alct.getKeyWG());
    valid->push_back(alct.isValid());
    detId->push_back(rawid);
    if (tp_index != -1) matchTp->push_back(tp_index);
  }

  void FillCLCT(CSCCLCTDigi clct, int rawid, int tp_index = -1) {
    strip->push_back(clct.getStrip());
    strip8->push_back(clct.getFractionalStrip(8));
    valid->push_back(clct.isValid());
    bend->push_back(clct.getBend());
    pattern->push_back(clct.getRun3Pattern());
    slope->push_back(clct.getSlope());
    detId->push_back(rawid);
    if (tp_index != -1) matchTp->push_back(tp_index);
  }

  void FillGEM(GEMDigi gemdigi, int rawid, int tp_index = -1) {
    strip->push_back(gemdigi.strip());
    detId->push_back(rawid);
    if (tp_index != -1) matchTp->push_back(tp_index);
  }

  void FillGEMPad(GEMPadDigi gempad, int rawid, int csc_index = -1) {
    pad->push_back(gempad.pad());
    part->push_back(gempad.nPartitions());
    detId->push_back(rawid);
    if (csc_index != -1) matchCSC->push_back(csc_index);
  }

  void FillGEMPad0(int rawid, int csc_index = -1) {
    pad->push_back(255);
    part->push_back(0);
    detId->push_back(rawid);
    if (csc_index != -1) matchCSC->push_back(csc_index);
  }

  void FillSimHit(int station_, int matchTp_) {
    station->push_back(station_);
    matchTp->push_back(matchTp_);
  }

  void FillGEMPadDigiCluster(GEMPadDigiCluster cluster, int rawid, int tp_index = -1) {
    // pad->insert(pad->end(), cluster.pads().begin(), cluster.pads().end());
    for (uint16_t pad_ : cluster.pads()) pads->push_back(pad_);
    pads->push_back(-1);
    len->push_back(cluster.pads().size());
    part->push_back(cluster.nPartitions());
    detId->push_back(rawid);
    if (tp_index != -1) matchTp->push_back(tp_index);
  }

  void FillDetId(int rawid) {
    detId->push_back(rawid);
  }
  void FillMatch(int tp_index) {
    matchTp->push_back(tp_index);
  }

  std::vector<int>*   detId;
  std::vector<int>*   matchTp;
  std::vector<int>*   matchCSC;

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
  std::vector<int>*   part;
  std::vector<int>*   pad;
  std::vector<int>*   pads;
  std::vector<int>*   len;

  //SimHit
  std::vector<int>* station;

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


private:
  TString name;
  int data_type;
  bool IsMatched;
  bool reader;
};

#endif
