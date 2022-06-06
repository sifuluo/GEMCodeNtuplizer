//////////////////////////////////////////////////////////////////////
//               GEMPad Strip Converter                             //
//               Author: Sifu Luo                                   //
//               sifuluo@tamu.edu                                   //
//////////////////////////////////////////////////////////////////////

#ifndef GEMPADTOSTRIPCONVERTER_CC
#define GEMPADTOSTRIPCONVERTER_CC

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/one/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include "CondFormats/DataRecord/interface/CSCL1TPLookupTableCCLUTRcd.h"
#include "CondFormats/DataRecord/interface/CSCL1TPLookupTableME11ILTRcd.h"
#include "CondFormats/DataRecord/interface/CSCL1TPLookupTableME21ILTRcd.h"
#include "CondFormats/CSCObjects/interface/CSCL1TPLookupTableCCLUT.h"
#include "CondFormats/CSCObjects/interface/CSCL1TPLookupTableME21ILT.h"
#include "CondFormats/CSCObjects/interface/CSCL1TPLookupTableME11ILT.h"

#include <vector>
#include <string>
#include <iostream>
#include "TString.h"

#include "GEMCode/GEMValidation/interface/MatcherManager.h"

using namespace std;

class GEMTool {
public:
  GEMTool() {
  };

  void Init(const CSCL1TPLookupTableME11ILT* conf1, const CSCL1TPLookupTableME21ILT* conf2){
    lookupTableME11ILT_ = conf1;
    lookupTableME21ILT_ = conf2;
  };

  vector<int> ConvertPad(GEMPadDigi digi, GEMDetId id, TString info = "") {
    if (!(digi.isValid())) return vector<int>{-1,-1,-1,-1,-1,-1};
    if (digi.isValid() && digi.pad() == 255) {
      cout << " Valid Digi with pad = 255" <<endl;
      return vector<int>{-1,-1,-1,-1,-1,-1};
    }
    return ConvertPad(digi.pad(),id,info);
  }

  vector<int> ConvertPad(int pad, GEMDetId id, TString info = "") {
    bool print = false;
    if (info != "") print = true;
    int station = id.station();
    int layer = id.layer();
    int chamber = id.chamber();
    bool isEven_ = chamber % 2;
    int roll = id.roll() - 1; // need to subtract 1 to use the LUTs
    int hs(-1), es(-1), hs_me1a(-1), es_me1a(-1), wg_min(-1), wg_max(-1);
    if (print) cout << info << ", Pad = " << pad << " is calling ConvertPad for GEMDetId: " << id.rawId() <<  ", " << (isEven_ ? "Even, " : "Odd, ") << id <<endl;
    if (station == 1) {
      if (isEven_) {
        es = lookupTableME11ILT_->GEM_pad_CSC_es_ME11b_even(pad);
        es_me1a = lookupTableME11ILT_->GEM_pad_CSC_es_ME11a_even(pad);
        if (roll >= 0) {
          if (layer == 1) {
            wg_min = lookupTableME11ILT_->GEM_roll_CSC_min_wg_ME11_even(roll);
            wg_max = lookupTableME11ILT_->GEM_roll_CSC_max_wg_ME11_even(roll);
          }
          else {
            wg_min = lookupTableME11ILT_->GEM_roll_CSC_min_wg_ME11_even(roll);
            wg_max = lookupTableME11ILT_->GEM_roll_CSC_max_wg_ME11_even(roll);
          }
        }
      }
      else {
        es = lookupTableME11ILT_->GEM_pad_CSC_es_ME11b_odd(pad);
        es_me1a = lookupTableME11ILT_->GEM_pad_CSC_es_ME11a_odd(pad);
        if (roll >= 0) {
          if (layer == 1) {
            wg_min = lookupTableME11ILT_->GEM_roll_CSC_min_wg_ME11_odd(roll);
            wg_max = lookupTableME11ILT_->GEM_roll_CSC_max_wg_ME11_odd(roll);
          }
          else {
            wg_min = lookupTableME11ILT_->GEM_roll_CSC_min_wg_ME11_odd(roll);
            wg_max = lookupTableME11ILT_->GEM_roll_CSC_max_wg_ME11_odd(roll);
          }
        }
      }
    }
    if (station == 2) {
      if (isEven_) {
        es = lookupTableME21ILT_->GEM_pad_CSC_es_ME21_even(pad);
        if (roll >=0) {
          if (layer == 1) {
            wg_min = lookupTableME21ILT_->GEM_roll_L1_CSC_min_wg_ME21_even(roll);
            wg_max = lookupTableME21ILT_->GEM_roll_L1_CSC_max_wg_ME21_even(roll);
          }
          else {
            wg_min = lookupTableME21ILT_->GEM_roll_L2_CSC_min_wg_ME21_even(roll);
            wg_max = lookupTableME21ILT_->GEM_roll_L2_CSC_max_wg_ME21_even(roll);
          }
        }
      }
      else {
        es = lookupTableME21ILT_->GEM_pad_CSC_es_ME21_odd(pad);
        if ( roll >=0 ) {
          if (layer == 1) {
            wg_min = lookupTableME21ILT_->GEM_roll_L1_CSC_min_wg_ME21_odd(roll);
            wg_max = lookupTableME21ILT_->GEM_roll_L1_CSC_max_wg_ME21_odd(roll);
          }
          else {
            wg_min = lookupTableME21ILT_->GEM_roll_L2_CSC_min_wg_ME21_odd(roll);
            wg_max = lookupTableME21ILT_->GEM_roll_L2_CSC_max_wg_ME21_odd(roll);
          }
        }
      }
    }
    vector<int> out{hs, es, hs_me1a, es_me1a, wg_min, wg_max};
    if (print) cout << "===>> Returning hs = " << hs << ", es = " << es << ", hs_me1a = " << hs_me1a<< ", es_me1a = " << es_me1a<< ", wg_min = " << wg_min<< ", wg_max = " << wg_max  <<endl;
    return out;
  }

private:
  const CSCL1TPLookupTableME11ILT* lookupTableME11ILT_;
  const CSCL1TPLookupTableME21ILT* lookupTableME21ILT_;
};

#endif
