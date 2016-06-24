#ifndef JECUNCERTAINTY_H
#define JECUNCERTAINTY_H
//-----------------------------------------------------------------------------
// File:        JECUncertainty.h
// Description: Return the jet energy correction uncertainty given the
//              corrected jet pT and eta. This code differs from 
//              JetCorrectionUncertainty in that this code reads from a
//              Table stored as a 2-D histogram. Therefore, only Root is
//              needed to use it.
// Usage:
//   string rooFileName("mytable.root");
//   JECUncertainty jecunc(filename);
//           :  :
//   double unc = jecunc(pt, eta);
//           :  :
//   double x = .... // +/-1 or randomly sampled from a TRandom3::Gaus()
//   double ptnew = pt * (1 + x*unc);           
//
// Created:     17-May-2016 Harrison B. Prosper   (Bari, Italy)   
//-----------------------------------------------------------------------------
#include <string>
#include "TFile.h"
#include "TH2F.h"
//-----------------------------------------------------------------------------
class JECUncertainty
{
 public:
  JECUncertainty() {}
  JECUncertainty(std::string rootFileName);
  ~JECUncertainty();
  double operator()(double pt, double eta);

  int    ptBins()  { return ptbins; }    
  double ptMin()   { return ptmin; }
  double ptMax()   { return ptmax; }

  int    etaBins() { return etabins; }    
  double etaMin()  { return etamin; }
  double etaMax()  { return etamax; }

  TH2F*  table()   { return table_; }

 private:
  TFile* rfile;
  TH2F*  table_;

  int    ptbins;
  double ptmin;
  double ptmax;
  double ptstep;

  int    etabins;
  double etamin;
  double etamax;
  double etastep;
};

#endif

