#ifndef JETSPECTRUMSMEARED_H
#define JETSPECTRUMSMEARED_H
// ---------------------------------------------------------------------------
// file: JetSpectrumSmeared.h
// HBP 2012 - 2014
// Updated: 23-Jun-2016 HBP use JECUncertainty 
// ---------------------------------------------------------------------------
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include "TMath.h"
#include "TH1D.h"
#include "Math/Interpolator.h"
#include "Math/WrappedFunction.h"
#include "Math/Integrator.h"
// --------------------=-----------------------------------------------------
// Inclusive jet spectrum in |y| < 0.5
// --------------------------------------------------------------------------
class JetSpectrum;
class JECUncertainty;

class JetSpectrumSmeared
{
 public:
  JetSpectrumSmeared() {}

  JetSpectrumSmeared(JetSpectrum* spectrum_,
		     JECUncertainty* JESunc_=0,
		     double JERunc_=0.1,
		     double x_=0, double y_=0,
		     double pTmin_=500.0,
		     double pTmax_=2800.0,
		     int npT_=46);

  JetSpectrumSmeared(const JetSpectrumSmeared&);

  ~JetSpectrumSmeared()
    {
      if (f1) delete f1; if (Intf1) delete Intf1;
      if (f2) delete f2; if (Intf2) delete Intf2;
    }
  
  double operator()(double pT);
  double operator()(double pTlow, double pThigh);
  double response(double pTreco, double pT);
  double sigmapT(double pT);
  double integrand(double pT);
  void setxy(double x_, double y_) { x = x_; y = y_; }

  double pTreco() { return pTreco_; }

private:
  JetSpectrum* spectrum;
  JECUncertainty* JESunc;
  double JERunc;
  double x, y;
  double pTmin, pTmax;
  int npT;
  std::vector<double> pT;
  std::vector<double> xsection;
  ROOT::Math::Interpolator* interp;

  ROOT::Math::WrappedMemFunction<JetSpectrumSmeared, 
    double (JetSpectrumSmeared::*)(double)>* f1;
 
  ROOT::Math::Integrator* Intf1;

  ROOT::Math::WrappedMemFunction<JetSpectrumSmeared, 
    double (JetSpectrumSmeared::*)(double)>* f2;
 
  ROOT::Math::Integrator* Intf2;  
  
  double pTreco_;
  double applySmearing_(double pT);

};

#endif
