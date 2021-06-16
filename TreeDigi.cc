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

  void Init(TTree* evttree, TString name_, int lct_type_, bool match_ = false) {
    name = name_;
    lct_type = lct_type_; // 0 for LCT, 1 for ALCT, 2 for CLCT, 3 for GemDigi
    IsMatched = match_;

    // LCT
    if (lct_type == 0) {
      phi = new std::vector<float>;
      eta = new std::vector<float>;
      z = new std::vector<float>;
      r = new std::vector<float>;
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
      evttree->Branch(name+"_phi", &phi);
      evttree->Branch(name+"_eta", &eta);
      evttree->Branch(name+"_z", &z);
      evttree->Branch(name+"_r", &r);
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
    else if (lct_type == 1) {
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
    else if (lct_type == 2) {
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
    else if (lct_type == 3) {
      phi = new std::vector<float>;
      eta = new std::vector<float>;
      z = new std::vector<float>;
      r = new std::vector<float>;
      detId = new std::vector<int>;
      strip = new std::vector<int>;
      evttree->Branch(name+"_phi", &phi);
      evttree->Branch(name+"_eta", &eta);
      evttree->Branch(name+"_z", &z);
      evttree->Branch(name+"_r", &r);
      evttree->Branch(name+"_detId", &detId);
      evttree->Branch(name+"_strip", &strip);

      if (IsMatched) {
        matchTp = new std::vector<int>;
        evttree->Branch(name+"_matchTp", &matchTp);
      }
    }

    // GEMPad
    else if (lct_type == 4) {
      phi = new std::vector<float>;
      eta = new std::vector<float>;
      z = new std::vector<float>;
      r = new std::vector<float>;
      detId = new std::vector<int>;
      pad = new std::vector<int>;
      part = new std::vector<int>;
      evttree->Branch(name+"_phi", &phi);
      evttree->Branch(name+"_eta", &eta);
      evttree->Branch(name+"_z", &z);
      evttree->Branch(name+"_r", &r);
      evttree->Branch(name+"_detId", &detId);
      evttree->Branch(name+"_pad", &pad);
      evttree->Branch(name+"_part", &part);

      if (IsMatched) {
        matchCSC = new std::vector<int>;
        evttree->Branch(name+"_matchCSC", &matchCSC);
      }
    }
  }

  void Reset() {
    if (lct_type == 0) {
      phi->clear();
      eta->clear();
      z->clear();
      r->clear();
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
    else if (lct_type == 1) {
      detId->clear();
      keywire->clear();
      hit->clear();
      position->clear();
      valid->clear();
    }
    else if (lct_type == 2) {
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
    else if (lct_type == 3) {
      phi->clear();
      eta->clear();
      z->clear();
      r->clear();
      detId->clear();
      strip->clear();
      if (IsMatched) {
        matchTp->clear();
      }
    }
    else if (lct_type == 4) {
      phi->clear();
      eta->clear();
      z->clear();
      r->clear();
      detId->clear();
      pad->clear();
      part->clear();
      if (IsMatched) {
        matchCSC->clear();
      }
    }
  }

  void FillGP(GlobalPoint gp) {
    phi->push_back(gp.phi());
    eta->push_back(gp.eta());
    z->push_back(gp.z());
    r->push_back(gp.perp());
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

private:
  TString name;
  int lct_type;
  bool IsMatched;
};

#endif
