//--------------------------------------------------------------
// File: Bayes.cc
// Description: Implement standalone Bayes calculator
// 
// Created: 11 Jan 2011 HBP
// Updated: 13 Mar 2011 HBP - fix normalization of post.
//          06 Mar 2014 HBP
//          30 May 2015 HBP - implement direct RooFit interface.
//--------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <stdlib.h>

#include "Bayes.h"
#include "PDFWrapper.h"

#include "TMinuit.h"
#include "TMath.h"
#include "Math/WrappedFunction.h"
#include "Math/Integrator.h"
#include "Math/RootFinder.h"

using namespace std;
// ---------------------------------------------------------------------------

// function to be minimized
namespace {
  const int RULE=5;
  const double RELTOL=1.e-4;
  const int MAXITER=10000;
  const double TOLERANCE=1.e-5;
  Bayes* OBJ=0;
  void nlpFunc(int&    /*npar*/, 
	       double* /*grad*/, 
	       double& fval,
	       double* xval,
	       int     /*iflag*/)
  {
    double poi = xval[0];
    fval = -log(OBJ->posterior(poi));
  }
};

Bayes::Bayes(PDFunction& model,
	     double poimin,
	     double poimax,
	     double cl,
	     PriorFunction* prior_)
  : _pdf(&model),
    _poimin(poimin),
    _poimax(poimax),
    _cl(cl),
    _prior(prior_),
    _rfprior(0),
    _rfpoi(0),
    _normalize(true),
    _interp(0),
    _nsteps(200),
    _x(vector<double>()),
    _y(vector<double>()),
    _verbosity(-1)
{
  if ( getenv("limits_verbosity") != (char*)0 )
    _verbosity = atoi(getenv("limits_verbosity"));
  OBJ = this;
}

Bayes::Bayes(RooAbsPdf& pdf, RooArgSet& obs, RooRealVar& poi,
	     double cl,
	     RooAbsPdf* prior_)
  : _pdf(new PDFWrapper(pdf, obs, poi)),
    _poimin(poi.getMin()),
    _poimax(poi.getMax()),
    _cl(cl),
    _prior(0),
    _rfprior(prior_),
    _rfpoi(&poi),
    _normalize(true),
    _interp(0),
    _nsteps(200),
    _x(vector<double>()),
    _y(vector<double>()),
    _verbosity(-1)
{
  if ( getenv("limits_verbosity") != (char*)0 )
    _verbosity = atoi(getenv("limits_verbosity"));
  OBJ = this;
}


Bayes::~Bayes() 
{
  // if _rfpoi is non-zero, this means we are using the RooFit
  // interface
  if (_rfpoi) delete _pdf;
  if (_interp) delete _interp;
}

double 
Bayes::prior(double poi)
{
  if (_prior)
    {
      return (*_prior)(poi);
    }
  else if (_rfprior)
    {
      _rfpoi->setVal(poi);
      return _rfprior->getVal();
    }
  return 1;
}

void Bayes::setData(std::vector<double>& d) 
{ 
  _data = d;
  _pdf->setData(_data);
}

double 
Bayes::likelihood(double poi)
{
  return (*_pdf)(poi);
}

double 
Bayes::normalize()
{
  ROOT::Math::WrappedMemFunction<Bayes, double (Bayes::*)(double)> 
    fn(*this, &Bayes::_likeprior);
  ROOT::Math::Integrator ifn(fn);
  ifn.SetRelTolerance(RELTOL);
  _normalization = ifn.Integral(_poimin, _poimax);
  _normalize = false;
  
  // Compute cdf at several points
  _x.clear();
  _y.clear();  
  _x.push_back(_poimin);
  _y.push_back(0);
  double step = (_poimax - _poimin)/_nsteps;
  double ysum = 0;
  for(int i=1; i < _nsteps+1; i++)
    {
      double xx = _poimin + i*step;
      _x.push_back(xx);
      double y = _likeprior(xx);
      ysum += y;
      _y.push_back(y + _y.back());
    }
  for(int i=0; i < _nsteps+1; i++) _y[i] /= ysum;

  try
    {
      delete _interp;
    }
  catch (...)
    {}
  _interp = new ROOT::Math::Interpolator(_x.size(),
					 ROOT::Math::
					 Interpolation::kLINEAR);
  _interp->SetData(_x, _y);

  return _normalization;
}

double 
Bayes::posterior(double poi)
{
  if ( _normalize ) normalize();
  return _likeprior(poi) / _normalization;
}

double 
Bayes::cdf(double poi)
{
  if ( poi < _poimin ) 
    return 0;
  else if ( poi > _poimax )
    return 1;

  if ( _normalize ) normalize();
  // Compute CDF
  ROOT::Math::WrappedMemFunction<Bayes, double (Bayes::*)(double)> 
    fn(*this, &Bayes::posterior);
  ROOT::Math::Integrator ifn(fn);
  return ifn.Integral(_poimin, poi);
}

double 
Bayes::percentile(double p)
{
  if ( p > 0 ) _cl = p; // Credibility level
  if ( _normalize ) normalize();

  // function whose root is to be found
  ROOT::Math::WrappedMemFunction<Bayes, double (Bayes::*)(double)> 
    fn(*this, &Bayes::_q);
  ROOT::Math::RootFinder rootfinder;
  
  rootfinder.SetFunction(fn, _poimin, _poimax);
  int status = rootfinder.Solve();
  if ( status != 1 )
    {
      cout << "*** Bayes *** RootFinder failed to find quantile"
           << endl;
      return -1;
    }
  return rootfinder.Root();
}

pair<double, double>
Bayes::estimate(double guess)
{
  if ( _normalize ) normalize();
    
  pair<double, double> results(0, 0);
  
  TMinuit minuit(1);
  minuit.SetPrintLevel(_verbosity);
  minuit.SetFCN(nlpFunc);
  //double chi2 = TMath::NormQuantile((1.0+CL)/2);
  //chi2 = chi2*chi2/2;
  minuit.SetErrorDef(0.5);

  int status=0;
  if ( guess <= 0 ) guess = 0.8*_poimin + 0.2*_poimax;
  double stepsize = (_poimax-_poimin)/2000;
  minuit.mnparm(0, "poi", guess, stepsize, 
		_poimin, _poimax, status);
  double args[2] = {MAXITER, TOLERANCE};
  minuit.mnexcm("MIGRAD", args, 2, status);
  if ( status != 0 )
    {
      cout << "Bayes::MAP failed to find MAP" << endl;
      return results;
    }
  
  // get fit result
  minuit.GetParameter(0, results.first, results.second);

  // // get MINOS errors
  // double upper, lower, perror, gcc;
  // minuit.mnerrs(0, upper, lower, perror, gcc);
  // if ( (int)results.size() > 0 ) results[0] = poihat;
  // if ( (int)results.size() > 1 ) results[1] = poierr;
  // if ( (int)results.size() > 2 ) results[2] = lower;
  // if ( (int)results.size() > 3 ) results[3] = upper;
  // if ( (int)results.size() > 4 ) results[4] = gcc;
  return results;
}

double 
Bayes::_likeprior(double poi)
{
  if ( poi != poi )
    {
      cout << "*** Bayes - this is Baaaad! poi = " 
           << poi << endl;
      exit(0);
    }
  return likelihood(poi) * prior(poi);
}

double 
Bayes::_q(double poi)
{
  if ( poi < _poimin ) 
    return 0;
  else if ( poi > _poimax )
    return 1;
  else
    return _interp->Eval(poi) - _cl;
}

