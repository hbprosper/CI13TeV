//-----------------------------------------------------------------------------
// File:        JECUncertainty.cc
// Description: see header
// Created:     17-May-2016 Harrison B. Prosper      
//-----------------------------------------------------------------------------
#include <iostream>
#include <cassert>
#include "JECUncertainty.h"
//-----------------------------------------------------------------------------
using namespace std;
//-----------------------------------------------------------------------------

JECUncertainty::JECUncertainty(string rootFileName)
  : rfile(new TFile(rootFileName.c_str())),
    table_(0)
{
  if ( ! rfile->IsOpen() )
    {
      cout << "JECUncertainty ** file " << rootFileName 
	   << " not opened " << endl
	   << "check file name"
	   << endl;
      exit(0);
    }
  table_ = dynamic_cast<TH2F*>(rfile->Get("htable"));
  assert(table_);

  TAxis* xaxis = table_->GetXaxis();
  ptbins = xaxis->GetNbins();
  ptmin  = xaxis->GetBinLowEdge(1);
  ptmax  = xaxis->GetBinUpEdge(ptbins);
  ptstep = (ptmax-ptmin)/ptbins;

  TAxis* yaxis = table_->GetYaxis();
  etabins= yaxis->GetNbins();
  etamin = yaxis->GetBinLowEdge(1);
  etamax = yaxis->GetBinUpEdge(etabins);
  etastep= (etamax-etamin)/etabins;
}
    
JECUncertainty::~JECUncertainty()

{ 
  if ( rfile ) rfile->Close();
}

double JECUncertainty::operator()(double pt, double eta)
{
  if ( pt  < ptmin ) return -1;
  if ( pt  > ptmax ) return -2;
  if ( eta < etamin) return -3;
  if ( eta > etamax) return -4;

  int ix = (int)((pt - ptmin)/ptstep) + 1;
  int iy = (int)((eta - etamin)/etastep) + 1;
  return table_->GetBinContent(ix, iy);
}


