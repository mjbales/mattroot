#ifndef MATTIO_H_INCLUDED
#define MATTIO_H_INCLUDED

#include <iostream>

#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TGraphErrors.h"

#include "MRHistDim.h"

void simpleHistOut(TH1* histOut, TString outputFileName, TString binUnitString, TString valueUnitString);
void convertTH1ToTXT(TH1* inpHist, TString outputFileName, bool printWidths=false);
void convertXYDataWErrorToTXT(int numData, double* x, double* y, double* eX, double* eY, TString titleString, TString outputFileName);
void convertTGraphErrorsToTXT(TGraphErrors* inpGraph, TString outputFileName);
void convertTH2ToTXT(TH2* inpHist, TString outputFilePath);
TH1D* getTabSeperatedHist(TString fileName, int numBins, double start, double end, char linedelim = '\n');
TH1D* getTabSeperatedHist(TString fileName, HistDim inpHistDim, char linedelim = '\n');
TH1D* getTabSeperatedHist(TString fileName, TH1D* histPrototype, char linedelim = '\n');
void convertTObjectPrint2StringStream(TObject* theObject, std::stringstream* theStream);
void convertTTree2Text(TString rootFileName, TString treeName, TString outputFileName, int numEntriesToRun);
TH1D* getCSVHist(TString fileName, int numBins, double start, double end, char endlineChar = '\n');
TH1D* getCSVHist2(TString fileName, int numBins, double start, double end);
TGraphErrors* getTabSeperatedTGraphErrors(TString filePath, char linedelim = '\n');
#endif // MATTIO_H_INCLUDED
