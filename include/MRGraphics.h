#ifndef INCLUDE_MRGRAPHICS_H_
#define INCLUDE_MRGRAPHICS_H_

//////////////////////////////////////////////////
//MRGraphics.h - by Matthew Bales
//misc graphics related functions
//
///////////////////////////////////////////////////

#include <vector>

#include "TGraphErrors.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "MRHistDim.h"
#include "TString.h"

static const int MR_GRAPH_NUMCOLORS = 11;
static const int MR_GRAPH_COLOR_LIST[MR_GRAPH_NUMCOLORS] = { kBlue, kRed, kGreen + 2, kCyan, kMagenta, kPink, kOrange, kViolet, kAzure, kYellow };
static const int MR_MARKER_NUMSTYLES = 13;
static const int MR_MARKER_STYLE_LIST[MR_MARKER_NUMSTYLES] = { 20, 22, 23, 25, 24, 32, 33, 3, 5, 28, 26, 30, 29 };

TCanvas* plotTGraphs(std::vector<TGraph*> theGraphs, bool autoFormat, TString titleString, bool logX = false, bool logY = false);
void plotTGraphs(std::vector<TGraph*> theGraphs, bool autoFormat, TString titleString, TString imagePath, bool logX = false, bool logY = false);
void plotTGraphs(std::vector<TString> filePaths, TString titleString, TString imagePath, bool logX = false, bool logY = false);
void plotTGraph(TString filePath, TString titleString, TString imagePath, bool logX = false, bool logY = false);

TCanvas* plotExpVersusMC(int numExpHists, TH1** expHists, int numMCHists, TH1** mcHists, TString titleString);
void plotExpVersusMCToImage(int numExpHists, TH1** expHists, int numMCHists, TH1** mcHists, TString titleString, TString imageName, bool useLog = false);
void plotExpVersusMCToImage(TString fileNameExp, TString fileNameMC, HistDim inpDim, TString titleString, TString imageName, bool useLog = false);

#endif /* INCLUDE_MRGRAPHICS_H_ */
