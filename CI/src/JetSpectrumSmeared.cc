// ---------------------------------------------------------------------------
// file: JetSpectrumSmeared.cc
// apply jet response function.
// HBP 2012 - 2014
// Updated: 11-October-2014 HBP - add non-perturbative correction
//          18-Jun-2016 HBP update NP corrections and smearing (SMP-15-007)
//          23-Jun-2016 HBP use JECUncertainty
//          19-Sep-2016 HBP fix integration bug introduced since June 2016!
// ---------------------------------------------------------------------------
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>

#include "TMath.h"
#include "TFile.h"
#include "TH1D.h"

#include "JetSpectrumSmeared.h"
#include "JetSpectrum.h"
#include "JECUncertainty.h"
// ---------------------------------------------------------------------------
using namespace std;
using namespace ROOT::Math;

namespace {
  const double RELTOL=1.e-4;
  const double PTMIN = 100;
  const double PTMAX =3100;
};
// ---------------------------------------------------------------------------
// Inclusive jet spectrum in |y| < 0.5
//----------------------------------------------------------------------------
JetSpectrumSmeared::JetSpectrumSmeared(JetSpectrum* spectrum_,
				       JECUncertainty* JESunc_,
				       double JERunc_,
				       double x_, double y_,
				       double pTmin_, double pTmax_,
				       int npT_)
  : spectrum(spectrum_), 
    JESunc(JESunc_),
    JERunc(JERunc_),
    x(x_), y(y_),
    pTmin(pTmin_), pTmax(pTmax_), npT(npT_),
    pT(std::vector<double>(npT+1,0)),
    xsection(std::vector<double>(npT+1,0)),
    
    interp(0),
    
    f1(new ROOT::Math::WrappedMemFunction<JetSpectrumSmeared, 
       double (JetSpectrumSmeared::*)(double)>
       (*this, &JetSpectrumSmeared::integrand)),
    
    Intf1(new ROOT::Math::Integrator(*f1)),

    f2(new ROOT::Math::WrappedMemFunction<JetSpectrumSmeared, 
       double (JetSpectrumSmeared::*)(double)>
       (*this, &JetSpectrumSmeared::operator())),
    
    Intf2(new ROOT::Math::Integrator(*f2))      
{
  Intf1->SetRelTolerance(RELTOL);
  Intf2->SetRelTolerance(RELTOL);    
  
  if ( spectrum->null() )
    interp = 0;
  else
    {
      double pTstep = (pTmax-pTmin)/npT;
      for(int c=0; c <= npT; c++)
	{
	  pT[c] = pTmin + c*pTstep;
	  if ( JESunc )
	    xsection[c] = applySmearing_(pT[c]);
	  else
	    xsection[c] = (*spectrum)(pT[c]);
	  
	  if ( spectrum->positive() )
	    xsection[c] = log(xsection[c]);
	}
      interp = new Interpolator(pT, xsection, Interpolation::kLINEAR);
    }
}

JetSpectrumSmeared::JetSpectrumSmeared(const JetSpectrumSmeared& o)
  : spectrum(o.spectrum), 
    JESunc(o.JESunc),
    JERunc(o.JERunc),
    x(o.x),
    y(o.y),
    pTmin(o.pTmin),
    pTmax(o.pTmax),
    npT(o.npT),
    pT(o.pT),
    xsection(o.xsection),
    interp(o.interp)
{
}

double JetSpectrumSmeared::operator()(double pt)
{
  if ( pt < pT.front() )  return 0;
  if ( pt > pT.back() )   return 0;
  if ( spectrum->null() ) return 0;
  double y = interp->Eval(pt);
  if ( y != y )
    {
      cout << "JetSpectrumSmeared::operator()(double) - NAN at pT = "
	   << pt << endl;
      exit(0);
    }  
  if ( spectrum->positive() )
    return exp(y);
  else
    return y;
}

double JetSpectrumSmeared::operator()(double pTlow, double pThigh)
{
  if ( pTlow  < pT.front() ) return 0;
  if ( pThigh > pT.back() )  return 0;
  if ( spectrum->null() )    return 0;
  return Intf2->Integral(pTlow, pThigh);
}

// Convolution of response function with NLO spectrum;
double JetSpectrumSmeared::applySmearing_(double pTreco)
{
  pTreco_ = pTreco; // NB: cache reco-level pT
  if ( spectrum->null() ) return 0;

  double offset = 5 * sigmapT(pTreco);
  double ptmin = TMath::Max( 100.0, pTreco - offset);
  double ptmax = TMath::Min(3100.0, pTreco + offset);
  double y = Intf1->Integral(ptmin, ptmax);
  return y;
}
 
double JetSpectrumSmeared::integrand(double pT)
{
  double y = response(pTreco_, pT) * (*spectrum)(pT);
  return y;
}

// Jet response function
double JetSpectrumSmeared::response(double pTreco, double pT)
{
  JECUncertainty& JECunc = *JESunc;
  double X = TMath::Max(1.e-3, 1.0 + x * JECunc(pTreco, 0));
  double Y = TMath::Max(1.e-3, 1.0 + y * JERunc);
  return TMath::Gaus(pTreco/X, pT, Y*sigmapT(pT), kTRUE); 
}

double JetSpectrumSmeared::sigmapT(double pT)
{    
  // SMP-15-007-PreApproval.pdf
  double a = 0.0257;
  double b = 1.091;
  double c = 0.5748;
  double d =-0.002826;
  return pT*(a + b / (pow(pT,c) + d*pT));
}
