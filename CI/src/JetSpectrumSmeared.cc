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

  // From Paolo Feb. 2018
  double gres_ofy0_1(double x){return  sqrt(1.992*abs(1.992)/(x*x)+0.5539*0.5539*pow(x,-0.7915)+0.02733*0.02733);}
  double gres_ofy0_2(double x){return  sqrt(3.301*abs(3.301)/(x*x)+0.6863*0.6863*pow(x,-0.8636)+0.02909*0.02909);}
  double gres_ofy0_3(double x){return  sqrt( 4.42*abs( 4.42)/(x*x)+0.7891*0.7891*pow(x,-0.9105)+0.03046*0.03046);}
  double gres_ofy0_4(double x){return  sqrt(5.566*abs(5.566)/(x*x)+0.8104*0.8104*pow(x,-0.9153)+0.03034*0.03034);}
  double gres_ofy0_5(double x){return  sqrt(6.782*abs(6.782)/(x*x)+0.8801*0.8801*pow(x,-0.9417)+0.03078*0.03078);}

  double wRho1=0.07858, wRho2=0.42870, wRho3=0.36185, wRho4=0.11092, wRho5=0.01995;
  
  double gres_ofy0_all(double x){return
      sqrt( wRho1*pow(gres_ofy0_1(x),2) +
	    wRho2*pow(gres_ofy0_2(x),2) +
	    wRho3*pow(gres_ofy0_3(x),2) +
	    wRho4* pow(gres_ofy0_4(x),2) +
	    wRho5*pow(gres_ofy0_5(x),2) );}
  
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
  
  // define range over which to integrate over true pT
  double offset = 5 * sigmapT(pTreco);
  double ptmin = TMath::Max( 100.0, pTreco - offset);
  double ptmax = TMath::Min(4200.0, pTreco + offset);
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

// double JetSpectrumSmeared::sigmapT(double pT)
// {    
//   // SMP-15-007-PreApproval.pdf
//   double a = 0.0257;
//   double b = 1.091;
//   double c = 0.5748;
//   double d =-0.002826;
//   return pT*(a + b / (pow(pT,c) + d*pT));
// }

double JetSpectrumSmeared::sigmapT(double pT)
{    
  return pT * gres_ofy0_all(pT);
}
