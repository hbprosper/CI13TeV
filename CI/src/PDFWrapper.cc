//--------------------------------------------------------------
//
// Author: Harrison B. Prosper (harry@hep.fsu.edu)
// File: PDFWrapper.cc
// Description: Implements wrapper for RooAbsPdfs.
// 
// Created: June 11, 2010
//--------------------------------------------------------------
#include <iostream>
#include <cmath>
#include "RooDataSet.h"
#include "RooLinkedListIter.h"
#include "TMath.h"
#include "PDFWrapper.h"

ClassImp(PDFWrapper)

using namespace std;
//--------------------------------------------------------------
PDFWrapper::PDFWrapper()
  : PDFunction(),
    _data(vector<double>())
{}

PDFWrapper::PDFWrapper(RooAbsPdf& pdf, RooArgSet& obs, RooRealVar& poi) 
  : PDFunction(),
    _pdf(&pdf),
    _obs(&obs),
    _poi(&poi),
    _data(vector<double>(_obs->getSize()))
{}

PDFWrapper::PDFWrapper(const PDFWrapper& other)
  : PDFunction(),
    _pdf(other._pdf),
    _obs(other._obs),
    _poi(other._poi),
    _data(other._data)
{}

PDFWrapper::~PDFWrapper() 
{}

vector<double>&
PDFWrapper::generate(double poi)
{
  _poi->setVal(poi);
  const RooArgSet* row = _pdf->generate(*_obs, 1)->get();
  RooArgList d(*row);
  for(unsigned int i=0; i < _data.size(); i++)
    _data[i] = dynamic_cast<RooRealVar*>(&d[i])->getVal();
  delete row;
  return _data;
}

void 
PDFWrapper::setData(std::vector<double>& data)
{
  RooArgList list(*_obs);
  for(size_t c=0; c < (size_t)list.getSize(); c++)
  {
    _data[c] = data[c];
    RooAbsArg* arg = list.at(c);
    (dynamic_cast<RooRealVar*>(arg))->setVal(_data[c]);
  }
}

double 
PDFWrapper::operator() (double poi)
{
  _poi->setVal(poi);
  double p =  _pdf->getVal();
  return p;
}



