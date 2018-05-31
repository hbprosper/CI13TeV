#ifndef PDFWrapper_H
#define PDFWrapper_H
//--------------------------------------------------------------
//
// Author: Harrison B. Prosper (harry@hep.fsu.edu)
// File: PDFWrapper
// Description: Wrapper for RooAbsPdf
// 
// Created: 14-Mar-2013
// Updated: 06-Mar-2014
//
//--------------------------------------------------------------
#include <vector>
#include "TROOT.h"
#include "PDFunction.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooRealVar.h"

/**  A wrapper for RooAbsPdf objects.
*/
class PDFWrapper : public PDFunction
{
 public:
  
  /**
   */
  PDFWrapper();
  
  /** Constructor:  
      @param pdf - An object of type RooAbsPdf
      @param obs - Set of data objects
      @param poi - parameter of interest
   */
  PDFWrapper(RooAbsPdf& pdf, RooArgSet& obs, RooRealVar& poi);
  
  PDFWrapper(const PDFWrapper& other);

  /**
   */
  virtual ~PDFWrapper();
  
  /** Generate and cache data for one experiment.
      @param poi - parameter of interest
  */
  std::vector<double>& generate(double poi);
  
  /** Computes PDF.
      @param poi  - parameter of interest 
  */
  double operator() (double poi);
  virtual RooRealVar* getPoi() {return _poi;}
  virtual std::string GetTitle() {return _pdf->GetTitle(); }

  void setData(std::vector<double>& data);

 private:

  RooAbsPdf*  _pdf;
  RooArgSet*  _obs;
  RooRealVar* _poi;
  RooArgList  _list;

  std::vector<double> _data;
  
  ClassDef(PDFWrapper,1)
};

#endif

