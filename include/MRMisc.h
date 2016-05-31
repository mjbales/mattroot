#ifndef MATTMISC_H_INCLUDED
#define MATTMISC_H_INCLUDED

#include <vector>

#include "TH1.h"
#include "TH3.h"
#include "TString.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TBranch.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "MRHistDim.h"

enum BASIC_MATH_OPERATION
{
	ADDITION, SUBTRACTION, MULTIPLICATION, DIVISION
};

TH3I* convertTwoTH1IToTH3I(TH1I* h1, TH1I* h2);
void makeHistNiceForLog(TH1* inpHist, int valToSet);
TH1D* convertTH1IToTH1D(TH1I* inpHist);
double getHistExtentHigh(TH1* inpHist);
double getHistExtentLow(TH1* inpHist);
void scaleHistogramsTogether(TH1* hist1, TH1* hist2, int lowBin = 1, int highBin = -1);  //Scale hist2 to have the same number within the range as hist1
TH1D* makeIntegratedHistogram(TH1D* histToIntegrate, bool doTopToBottom = false);
double getMaximumofHists(int numHists, TH1** theHists);
void scaleHistTo(TH1* theHist, double integralNumber);
TCanvas* plotHist2D(TH2D* inpHist, TString titleString);
double calcTotalResidualAndError(TH1* hist1, TH1* hist2, double& totalError, int binLow = 1, int binHigh = 0);  //(hist1-hist2)/hist2
void plotResidualsToImage(int numHistsA, TH1** histsA, int numHistsB, TH1** histsB, TString imageTitleString, TString imageFileName, double lowRange = -1., double highRange = 1., bool normalizedResidualsFlag = true);
void setHistNameTitle(TH1* inpHist, TString nameTitleString);
TH1* makeEmptyTH1Clone(TH1* inpHist, TString newName, TString newTitle = "");
int checkArrayForInt(int toCheckFor, const int numInArray, const int* theArray);
void plotResidualsAndComparisonToImage(int numHistsA, TH1** histsA, int numHistsB, TH1** histsB, TString imageTitleString, TString imageFileName, double lowRange = -1., double highRange = 1., bool normalizedResidualsFlag = true);
TString getBranchVarType(const TBranch* inpBranch);
void changeGainOfHist(TH1D* inpHist, double gainFactor);
bool checkTolerance(double target, double variable, double tolerance);  //Returns true if variable is within tolerance (in decimal percent) of target
int fitHistWithGaus(TH1D* inpHist, TF1* gaussFuncToFit); //Fits gaussian to histogram, then refits after centering range on mean. Returns fit status
double getBasicErrorProgation(double x, double y, double xErr, double yErr, BASIC_MATH_OPERATION inpOperation); //Outputs error
TH1* calcAvgHistWithErrorFromHists(const int numHists, TH1** inpHists, TString histName);
TString getCurrentDateString();
void createBlankThesisPlot(TString imageName, TString imageTitleString, TString legendList);
double GetIdealMaxForPlotting(TH1* inpHist, double sigmaLimit = 3);
double GetIdealMinForPlotting(TH1* inpHist, double sigmaLimit = 3);
double getCountMeanMinusOutliers(TH1* inpHist, double sigmaLimit, double& meanError);
TCanvas* getThesisPlotCanvas(TString name = "thesisCanvas");
TCanvas* getHalfPresentationPlot(TString name);
TCanvas* getHalfPresentationPlot2(TString name);
TH1* trimHist(TH1* inpHist, TString newName = "", int start = 0, int end = 0);
TH1D* makeTH1DFromDim(TString name, TString title, HistDim inpHistDim);
TH2D* makeTH2DFromDim(TString name, TString title, HistDim inpXHistDim, HistDim inpYHistDim);
TH3D* makeTH3DFromDim(TString name, TString title, HistDim inpXHistDim, HistDim inpYHistDim, HistDim inpZHistDim);
HistDim getDimFromHist(TH1* inpHist);
int getVarArrayFromHist(TH1* inpHist, double** outArray);
TH1D* adjustXScaleHist(TH1D* inpHist, double newStart, double newEnd, TString adjustedHistName = "");  //Creates a new histogram with the same data and bin count but with a new start and end
TH1D* recenterPeriodicHist(TH1D* inpHist, int numBinsShift, TString adjustedHistName = "");
double reducePeriodicNumber(double inp, double start, double end);
TVector3 getRandomPointInCylinder(double radius, double height);
TVector2 getRandomPointInCircle(double radius);
TVector3 getRandomDirection();
TVector2 getRandomDirection2D();
double meanArray(int numData, double* theData);
double stDevArray(int numData, double* theData, bool useBesselCorrection);
double absMinArray(int numData, double* theData);
double maxFromTGraphs(std::vector<TGraph*> theGraphs, bool includeError, bool forX = false);
double minFromTGraphs(std::vector<TGraph*> theGraphs, bool includeError, bool forX = false);
TH1I* histFromHist(TH1I* inpHist);
void scaleTGraphErrors(TGraphErrors* inpGraph, double scaleFactor);
double* makeEvenSampleOverLog(const int numPoints, const double begin, const double end);
double meanVector(const std::vector<double>& theData);
double carefulMeanVector(const std::vector<double>& theData);
double stDevVector(const std::vector<double>& theData, const bool useBesselCorrection);
double carefullStDevVector(const std::vector<double>& theData, const bool useBesselCorrection);
double minVector(const std::vector<double>& theData);
double maxVector(const std::vector<double>& theData);
double reducePeriodicToMeanInVector(std::vector<double>& theData);
TH1* makeVarHistScaledByWidth(TH1* inpHist, std::vector<double> binSpacing);
#endif // MATTMISC_H_INCLUDED
