//////////////////////////////////////////////////////////////////////
//               GEMPad Strip Converter                             //
//               Author: Sifu Luo                                   //
//               sifuluo@tamu.edu                                   //
//////////////////////////////////////////////////////////////////////

#ifndef GEMPADTOSTRIPCONVERTER_CC
#define GEMPADTOSTRIPCONVERTER_CC

#include "FWCore/ParameterSet/interface/ParameterSet.h"
// #include "DataFormats/GEMDigi/interface/GEMPadDigiCollection.h"
// #include "DataFormats/GEMDigi/interface/GEMPadDigiClusterCollection.h"
// #include "DataFormats/GEMDigi/interface/GEMCoPadDigi.h"
// #include "L1Trigger/CSCTriggerPrimitives/interface/GEMInternalCluster.h"
#include "L1Trigger/CSCTriggerPrimitives/interface/CSCLUTReader.h"

#include <vector>
#include <string>
#include <iostream>
#include "TString.h"

#include "GEMCode/GEMValidation/interface/MatcherManager.h"

using namespace std;

class GEMTool {
public:
  GEMTool(edm::ParameterSet const& conf) {
    Init(conf);
  };

  void Init(edm::ParameterSet const& conf){
    padToHsME1aFiles_ = conf.getParameter<vector<string>>("padToHsME1aFiles");
    padToHsME1bFiles_ = conf.getParameter<vector<string>>("padToHsME1bFiles");

    padToEsME1aFiles_ = conf.getParameter<vector<string>>("padToEsME1aFiles");
    padToEsME1bFiles_ = conf.getParameter<vector<string>>("padToEsME1bFiles");

    rollToMaxWgME11Files_ = conf.getParameter<vector<string>>("rollToMaxWgME11Files");
    rollToMinWgME11Files_ = conf.getParameter<vector<string>>("rollToMinWgME11Files");

    GEMCSCLUT_pad_hs_ME1a_even = make_unique<CSCLUTReader>(padToHsME1aFiles_[0]);
    GEMCSCLUT_pad_hs_ME1a_odd = make_unique<CSCLUTReader>(padToHsME1aFiles_[1]);
    GEMCSCLUT_pad_hs_ME1b_even = make_unique<CSCLUTReader>(padToHsME1bFiles_[0]);
    GEMCSCLUT_pad_hs_ME1b_odd = make_unique<CSCLUTReader>(padToHsME1bFiles_[1]);

    GEMCSCLUT_pad_es_ME1a_even = make_unique<CSCLUTReader>(padToEsME1aFiles_[0]);
    GEMCSCLUT_pad_es_ME1a_odd = make_unique<CSCLUTReader>(padToEsME1aFiles_[1]);
    GEMCSCLUT_pad_es_ME1b_even = make_unique<CSCLUTReader>(padToEsME1bFiles_[0]);
    GEMCSCLUT_pad_es_ME1b_odd = make_unique<CSCLUTReader>(padToEsME1bFiles_[1]);

    GEMCSCLUT_roll_l1_min_wg_ME11_even = make_unique<CSCLUTReader>(rollToMinWgME11Files_[0]);
    GEMCSCLUT_roll_l1_min_wg_ME11_odd = make_unique<CSCLUTReader>(rollToMinWgME11Files_[1]);
    GEMCSCLUT_roll_l2_min_wg_ME11_even = make_unique<CSCLUTReader>(rollToMinWgME11Files_[2]);
    GEMCSCLUT_roll_l2_min_wg_ME11_odd = make_unique<CSCLUTReader>(rollToMinWgME11Files_[3]);

    GEMCSCLUT_roll_l1_max_wg_ME11_even = make_unique<CSCLUTReader>(rollToMaxWgME11Files_[0]);
    GEMCSCLUT_roll_l1_max_wg_ME11_odd = make_unique<CSCLUTReader>(rollToMaxWgME11Files_[1]);
    GEMCSCLUT_roll_l2_max_wg_ME11_even = make_unique<CSCLUTReader>(rollToMaxWgME11Files_[2]);
    GEMCSCLUT_roll_l2_max_wg_ME11_odd = make_unique<CSCLUTReader>(rollToMaxWgME11Files_[3]);

    padToHsME21Files_ = conf.getParameter<vector<string>>("padToHsME21Files");
    padToEsME21Files_ = conf.getParameter<vector<string>>("padToEsME21Files");

    rollToMaxWgME21Files_ = conf.getParameter<vector<string>>("rollToMaxWgME21Files");
    rollToMinWgME21Files_ = conf.getParameter<vector<string>>("rollToMinWgME21Files");

    GEMCSCLUT_pad_hs_ME21_even = make_unique<CSCLUTReader>(padToHsME21Files_[0]);
    GEMCSCLUT_pad_hs_ME21_odd = make_unique<CSCLUTReader>(padToHsME21Files_[1]);
    GEMCSCLUT_pad_es_ME21_even = make_unique<CSCLUTReader>(padToEsME21Files_[0]);
    GEMCSCLUT_pad_es_ME21_odd = make_unique<CSCLUTReader>(padToEsME21Files_[1]);

    GEMCSCLUT_roll_l1_min_wg_ME21_even = make_unique<CSCLUTReader>(rollToMinWgME21Files_[0]);
    GEMCSCLUT_roll_l1_min_wg_ME21_odd = make_unique<CSCLUTReader>(rollToMinWgME21Files_[1]);
    GEMCSCLUT_roll_l2_min_wg_ME21_even = make_unique<CSCLUTReader>(rollToMinWgME21Files_[2]);
    GEMCSCLUT_roll_l2_min_wg_ME21_odd = make_unique<CSCLUTReader>(rollToMinWgME21Files_[3]);

    GEMCSCLUT_roll_l1_max_wg_ME21_even = make_unique<CSCLUTReader>(rollToMaxWgME21Files_[0]);
    GEMCSCLUT_roll_l1_max_wg_ME21_odd = make_unique<CSCLUTReader>(rollToMaxWgME21Files_[1]);
    GEMCSCLUT_roll_l2_max_wg_ME21_even = make_unique<CSCLUTReader>(rollToMaxWgME21Files_[2]);
    GEMCSCLUT_roll_l2_max_wg_ME21_odd = make_unique<CSCLUTReader>(rollToMaxWgME21Files_[3]);
    cout << rollToMaxWgME11Files_[0] <<endl;
  };

  vector<int> ConvertPad(GEMPadDigi digi, GEMDetId id, TString info = "") {
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
    // cout << "Start to look up pad = " << pad << " in station " << station <<", layer " <<layer << ", chamber " <<chamber<<"(" <<(isEven_ ? "Even" : "Odd") <<"), roll " << roll <<", " << id <<endl;
    if (station == 1) {
      if (isEven_) {
        hs = GEMCSCLUT_pad_hs_ME1b_even->lookup(pad);
        hs_me1a = GEMCSCLUT_pad_hs_ME1a_even->lookup(pad);
        es = GEMCSCLUT_pad_es_ME1b_even->lookup(pad);
        es_me1a = GEMCSCLUT_pad_es_ME1a_even->lookup(pad);
        if (roll >= 0) {
          if (layer == 1) {
            wg_min = GEMCSCLUT_roll_l1_min_wg_ME11_even->lookup(roll);
            wg_max = GEMCSCLUT_roll_l1_max_wg_ME11_even->lookup(roll);
          }
          else {
            wg_min = GEMCSCLUT_roll_l2_min_wg_ME11_even->lookup(roll);
            wg_max = GEMCSCLUT_roll_l2_max_wg_ME11_even->lookup(roll);
          }
        }
      }
      else {
        hs = GEMCSCLUT_pad_hs_ME1b_odd->lookup(pad);
        hs_me1a = GEMCSCLUT_pad_hs_ME1a_odd->lookup(pad);
        es = GEMCSCLUT_pad_es_ME1b_odd->lookup(pad);
        es_me1a = GEMCSCLUT_pad_es_ME1a_odd->lookup(pad);
        if (roll >= 0) {
          if (layer == 1) {
            wg_min = GEMCSCLUT_roll_l1_min_wg_ME11_odd->lookup(roll);
            wg_max = GEMCSCLUT_roll_l1_max_wg_ME11_odd->lookup(roll);
          }
          else {
            wg_min = GEMCSCLUT_roll_l2_min_wg_ME11_odd->lookup(roll);
            wg_max = GEMCSCLUT_roll_l2_max_wg_ME11_odd->lookup(roll);
          }
        }
      }
    }
    if (station == 2) {
      if (isEven_) {
        hs = GEMCSCLUT_pad_hs_ME21_even->lookup(pad);
        es = GEMCSCLUT_pad_es_ME21_even->lookup(pad);
        if (roll >=0) {
          if (layer == 1) {
            wg_min = GEMCSCLUT_roll_l1_min_wg_ME21_even->lookup(roll);
            wg_max = GEMCSCLUT_roll_l1_max_wg_ME21_even->lookup(roll);
          }
          else {
            wg_min = GEMCSCLUT_roll_l2_min_wg_ME21_even->lookup(roll);
            wg_max = GEMCSCLUT_roll_l2_max_wg_ME21_even->lookup(roll);
          }
        }
      }
      else {
        hs = GEMCSCLUT_pad_hs_ME21_odd->lookup(pad);
        es = GEMCSCLUT_pad_es_ME21_odd->lookup(pad);
        if ( roll >=0 ) {
          if (layer == 1) {
            wg_min = GEMCSCLUT_roll_l1_min_wg_ME21_odd->lookup(roll);
            wg_max = GEMCSCLUT_roll_l1_max_wg_ME21_odd->lookup(roll);
          }
          else {
            wg_min = GEMCSCLUT_roll_l2_min_wg_ME21_odd->lookup(roll);
            wg_max = GEMCSCLUT_roll_l2_max_wg_ME21_odd->lookup(roll);
          }
        }
      }
    }
    vector<int> out{hs, es, hs_me1a, es_me1a, wg_min, wg_max};
    if (print) cout << "===>> Returning hs = " << hs << ", es = " << es << ", hs_me1a = " << hs_me1a<< ", es_me1a = " << es_me1a<< ", wg_min = " << wg_min<< ", wg_max = " << wg_max  <<endl;
    return out;
  }

private:
  // strings to paths of LUTs
  vector<string> padToHsME1aFiles_;
  vector<string> padToHsME1bFiles_;
  vector<string> padToHsME21Files_;

  vector<string> padToEsME1aFiles_;
  vector<string> padToEsME1bFiles_;
  vector<string> padToEsME21Files_;

  vector<string> rollToMaxWgME11Files_;
  vector<string> rollToMinWgME11Files_;
  vector<string> rollToMaxWgME21Files_;
  vector<string> rollToMinWgME21Files_;

  // unique pointers to the luts
  unique_ptr<CSCLUTReader> GEMCSCLUT_pad_hs_ME1a_even;
  unique_ptr<CSCLUTReader> GEMCSCLUT_pad_hs_ME1a_odd;
  unique_ptr<CSCLUTReader> GEMCSCLUT_pad_hs_ME1b_even;
  unique_ptr<CSCLUTReader> GEMCSCLUT_pad_hs_ME1b_odd;
  unique_ptr<CSCLUTReader> GEMCSCLUT_pad_hs_ME21_even;
  unique_ptr<CSCLUTReader> GEMCSCLUT_pad_hs_ME21_odd;

  unique_ptr<CSCLUTReader> GEMCSCLUT_pad_es_ME1a_even;
  unique_ptr<CSCLUTReader> GEMCSCLUT_pad_es_ME1a_odd;
  unique_ptr<CSCLUTReader> GEMCSCLUT_pad_es_ME1b_even;
  unique_ptr<CSCLUTReader> GEMCSCLUT_pad_es_ME1b_odd;
  unique_ptr<CSCLUTReader> GEMCSCLUT_pad_es_ME21_even;
  unique_ptr<CSCLUTReader> GEMCSCLUT_pad_es_ME21_odd;

  unique_ptr<CSCLUTReader> GEMCSCLUT_roll_l1_max_wg_ME11_even;
  unique_ptr<CSCLUTReader> GEMCSCLUT_roll_l1_max_wg_ME11_odd;
  unique_ptr<CSCLUTReader> GEMCSCLUT_roll_l1_min_wg_ME11_even;
  unique_ptr<CSCLUTReader> GEMCSCLUT_roll_l1_min_wg_ME11_odd;
  unique_ptr<CSCLUTReader> GEMCSCLUT_roll_l1_max_wg_ME21_even;
  unique_ptr<CSCLUTReader> GEMCSCLUT_roll_l1_max_wg_ME21_odd;
  unique_ptr<CSCLUTReader> GEMCSCLUT_roll_l1_min_wg_ME21_even;
  unique_ptr<CSCLUTReader> GEMCSCLUT_roll_l1_min_wg_ME21_odd;

  unique_ptr<CSCLUTReader> GEMCSCLUT_roll_l2_max_wg_ME11_even;
  unique_ptr<CSCLUTReader> GEMCSCLUT_roll_l2_max_wg_ME11_odd;
  unique_ptr<CSCLUTReader> GEMCSCLUT_roll_l2_min_wg_ME11_even;
  unique_ptr<CSCLUTReader> GEMCSCLUT_roll_l2_min_wg_ME11_odd;
  unique_ptr<CSCLUTReader> GEMCSCLUT_roll_l2_max_wg_ME21_even;
  unique_ptr<CSCLUTReader> GEMCSCLUT_roll_l2_max_wg_ME21_odd;
  unique_ptr<CSCLUTReader> GEMCSCLUT_roll_l2_min_wg_ME21_even;
  unique_ptr<CSCLUTReader> GEMCSCLUT_roll_l2_min_wg_ME21_odd;
};

#endif
