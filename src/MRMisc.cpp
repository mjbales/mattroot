#include "MRMisc.h"

//Standard Libraries
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <limits>
#include <ctime>
#include <math.h>
#include <cfloat>
#include <algorithm>

//ROOT Libraries
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFrame.h"
#include "TCut.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TF1.h"
#include "TGraphErrors.h"

//Matt ROOT Libraries
#include "MRIO.h"
#include "MRPhys.h"
#include "MRHistDim.h"
#include "MRGraphics.h"
#include "MRText.h"

using namespace std;

void makeHistNiceForLog(TH1* inpHist, int valToSet)
{
	for (int i = 0; i <= inpHist->GetEntries() + 1; i++)
	{
		if(inpHist->GetBinContent(i) < 1) inpHist->SetBinContent(i, valToSet);
	}
}

TH1D* convertTH1IToTH1D(TH1I* inpHist)
{
	TString newName = inpHist->GetName();
	newName += "_db";
	int nBins = inpHist->GetNbinsX();
	TH1D* resultHist = new TH1D(newName, inpHist->GetTitle(), nBins, inpHist->GetBinLowEdge(1), inpHist->GetBinLowEdge(nBins) + inpHist->GetBinWidth(nBins));

	for (int i = 1; i <= nBins; i++)
	{
		resultHist->SetBinContent(i, inpHist->GetBinContent(i));
		resultHist->SetBinError(i, inpHist->GetBinError(i));
	}

	return resultHist;
}

double getHistExtentHigh(TH1* inpHist)
{
	return inpHist->GetBinLowEdge(inpHist->GetNbinsX()) + inpHist->GetBinWidth(inpHist->GetNbinsX());
}

double getHistExtentLow(TH1* inpHist)
{
	return inpHist->GetBinLowEdge(1);
}

void scaleHistogramsTogether(TH1* hist1, TH1* hist2, int lowBin, int highBin)  //Scale hist2 to have the same number within the range as hist1
{
	int newHighBin = highBin;
	if(highBin < 0) newHighBin = hist1->GetNbinsX();

	double scaleFactor = hist1->Integral(lowBin, newHighBin) / hist2->Integral(lowBin, newHighBin);
	hist2->Scale(scaleFactor);
}

TH1D* makeIntegratedHistogram(TH1D* histToIntegrate, bool doTopToBottom)
{
	int numBins = histToIntegrate->GetNbinsX();
	double lowEdge = getHistExtentLow(histToIntegrate);
	double highEdge = getHistExtentHigh(histToIntegrate);
	TString histName = histToIntegrate->GetName();
	histName += "_integrated";
	TString histTitle = histToIntegrate->GetTitle();
	histTitle += " (Integrated)";
	TH1D* integratedHist = new TH1D(histName, histTitle, numBins, lowEdge, highEdge);

	for (int i = 1; i <= numBins; i++)
	{
		if(doTopToBottom)
		{
			integratedHist->SetBinContent(i, histToIntegrate->Integral(1, i));
		}
		else
		{
			integratedHist->SetBinContent(i, histToIntegrate->Integral(i, numBins));
		}

	}

	return integratedHist;

}

double getMaximumofHists(int numHists, TH1** theHists)
{
	double maximum = -DBL_MAX;
	double value = 0;
	for (int i = 0; i < numHists; i++)
	{
		value = theHists[i]->GetBinContent(theHists[i]->GetMaximumBin());
		if(value > maximum) maximum = value;
	}
	return value;
}

void scaleHistTo(TH1* theHist, double integralNumber)
{
	double currentIntegral = theHist->Integral();
	theHist->Scale(integralNumber / currentIntegral);
}

TCanvas* plotHist2D(TH2D* inpHist, TString titleString)
{

	TCanvas* theCanvas = new TCanvas("2DComparisonCanvas", "2DComparisonCanvas", 10, 10, 1024, 800);
	gStyle->SetTitleColor(kWhite);
	// gStyle->SetPadRightMargin(2.5);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetFillColor(kWhite);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
//    gPad->SetLeftMargin(1.4);
//    theCanvas->SetRightMargin(2.4);
//    gPad->SetRightMargin(2.4);
//    gPad->SetTopMargin(1.6);
//    inpHist->GetYaxis()->SetTitleOffset(1.3);

	inpHist->Draw("COLZ");

	TPaveText *pt = new TPaveText(.05, .925, .95, .995, "NDC");
	pt->AddText(titleString);
	pt->SetFillColor(kWhite);
	pt->SetLineColor(kWhite);

	pt->SetShadowColor(kWhite);
	pt->SetFillStyle(4000);  //transparent
	pt->SetLineWidth(0);
	pt->Draw();

	theCanvas->Update();
	return theCanvas;
}

double calcTotalResidualAndError(TH1* hist1, TH1* hist2, double& totalError, int binLow, int binHigh)
{
	double error1, error2;
	double integral1, integral2, totalIntegral;
	int binHighUsed = binHigh;
	if(binHigh == 0) binHighUsed = hist1->GetNbinsX();

	integral1 = hist1->IntegralAndError(binLow, binHighUsed, error1);

	integral2 = hist2->IntegralAndError(binLow, binHighUsed, error2);

//    integralSub=integral1-integral2;
//
//    totalIntegral=integralSub/integral2;
//
//    errorSub=sqrt(error1*error1+error2*error2);
//
//    totalError=abs(totalIntegral)*sqrt(pow(errorSub/integralSub,2)+pow(error2/integral2,2));

	totalIntegral = integral1 / integral2 - 1.;

	totalError = abs(integral1 / integral2) * sqrt(pow(error1 / integral1, 2) + pow(error2 / integral2, 2));

	return totalIntegral;
}

//One of numHistsA and numHistsB must be 1
void plotResidualsToImage(int numHistsA, TH1** histsA, int numHistsB, TH1** histsB, TString imageTitleString, TString imagePath, double lowRange, double highRange, bool normalizedResidualsFlag)
{
	TH1** subHists;
	int numSubHists;

	if(numHistsB == 1)
	{
		numSubHists = numHistsA;
	}
	else if(numHistsA == 1 || numHistsA == numHistsB)
	{
		numSubHists = numHistsB;
	}
	else
	{
		return;
	}
	double* residualTotal = new double[numSubHists];
	double* residualError = new double[numSubHists];

	subHists = new TH1*[numSubHists];

	for (int i = 0; i < numSubHists; i++)
	{
		TString residHistName;
		residHistName.Form("residHist%d", i);
		subHists[i] = (TH1*) histsA[0]->Clone(residHistName);
		subHists[i]->Reset();
		subHists[i]->SetTitle("Residual");
		subHists[i]->Sumw2();
		subHists[i]->GetXaxis()->SetRangeUser(0, 1200);

		TString titleString;
		if(numHistsB == 1)
		{
			subHists[i]->Add(histsA[i], histsB[0], 1, -1);
			titleString = histsA[i]->GetTitle();
			if(normalizedResidualsFlag)
			{
				subHists[i]->Divide(histsB[0]);
				residualTotal[i] = calcTotalResidualAndError(histsA[i], histsB[0], residualError[i], histsB[0]->GetXaxis()->GetFirst(), histsB[0]->GetXaxis()->GetLast());
//                residualTotal[i]=calcTotalResidualAndError(histsA[i],histsB[0],residualError[i],1,histsB[0]->GetNbinsX());
			}
			else
			{
				residualTotal[i] = subHists[i]->IntegralAndError(histsB[0]->GetXaxis()->GetFirst(), histsB[0]->GetXaxis()->GetLast(), residualError[i]);
//                residualTotal[i]=subHists[i]->IntegralAndError(1,subHists[0]->GetNbinsX(),residualError[i]);

			}
			titleString += TString::Format(": %.3f #pm %.3f", residualTotal[i], residualError[i]);

			subHists[i]->SetTitle(titleString);
		}
		else if(numHistsA == 1)
		{
			subHists[i]->Add(histsA[0], histsB[i], 1, -1);
			titleString = histsB[i]->GetTitle();
			if(normalizedResidualsFlag)
			{
				subHists[i]->Divide(histsB[i]);
				residualTotal[i] = calcTotalResidualAndError(histsA[0], histsB[i], residualError[i], histsA[0]->GetXaxis()->GetFirst(), histsA[0]->GetXaxis()->GetLast());
//                residualTotal[i]=calcTotalResidualAndError(histsA[0],histsB[i],residualError[i],1,histsA[0]->GetNbinsX());
			}
			else
			{
				residualTotal[i] = subHists[i]->IntegralAndError(histsA[0]->GetXaxis()->GetFirst(), histsA[0]->GetXaxis()->GetLast(), residualError[i]);
//                residualTotal[i]=subHists[i]->IntegralAndError(1,subHists[0]->GetNbinsX(),residualError[i]);
			}
			titleString += TString::Format(": %.3f #pm %.3f", residualTotal[i], residualError[i]);

			subHists[i]->SetTitle(titleString);
		}
		else //They are equal
		{
			subHists[i]->Add(histsA[i], histsB[i], 1, -1);
			titleString = histsB[i]->GetTitle();
			if(normalizedResidualsFlag)
			{
				subHists[i]->Divide(histsB[i]);
				residualTotal[i] = calcTotalResidualAndError(histsA[i], histsB[i], residualError[i], histsA[i]->GetXaxis()->GetFirst(), histsA[i]->GetXaxis()->GetLast());
				titleString += TString::Format(": %.3f #pm %.3f", residualTotal[i], residualError[i]);
			}
			subHists[i]->SetTitle(titleString);
		}

		subHists[i]->GetYaxis()->SetRangeUser(lowRange, highRange);
		subHists[i]->GetXaxis()->SetRange(histsA[0]->GetXaxis()->GetFirst(), histsA[0]->GetXaxis()->GetLast());

	}

	TCanvas* theCanvas = plotExpVersusMC(numSubHists, subHists, 0, nullptr, imageTitleString);

//    theCanvas->SetLogx();

	double x[2] = { -1, 2001 };
	double y[2] = { 0, 0 };
	TGraph zeroAxis(2, x, y);
	zeroAxis.SetLineWidth(3);
	zeroAxis.Draw("same");

	TString testTitle = histsA[0]->GetXaxis()->GetTitle();
	subHists[0]->GetXaxis()->SetTitle(testTitle);

	if(normalizedResidualsFlag)
		subHists[0]->GetYaxis()->SetTitle("Normalized Residual");
	else
		subHists[0]->GetYaxis()->SetTitle(histsA[0]->GetYaxis()->GetTitle());

	theCanvas->Update();
	theCanvas->SaveAs(imagePath);

	delete theCanvas;

	for (int i = 0; i < numSubHists; i++)
	{
		delete subHists[i];

	}
	delete[] subHists;
	delete[] residualTotal;
	delete[] residualError;

}

void setHistNameTitle(TH1* inpHist, TString nameTitleString)
{
	inpHist->SetName(nameTitleString);
	inpHist->SetTitle(nameTitleString);
}

TH1* makeEmptyTH1Clone(TH1* inpHist, TString newName, TString newTitle)
{
	TH1* theClone = (TH1*) inpHist->Clone(newName);

	theClone->SetTitle(newTitle);
	theClone->Reset();
	return theClone;
}

int checkArrayForInt(int toCheckFor, const int numInArray, const int* theArray)
{
	for (int i = 0; i < numInArray; i++)
	{
		if(theArray[i] == toCheckFor) return i;
	}
	return -1;
}

void plotResidualsAndComparisonToImage(int numHistsA, TH1** histsA, int numHistsB, TH1** histsB, TString imageTitleString, TString imageFileName, double lowRange, double highRange, bool normalizedResidualsFlag)
{
	plotResidualsToImage(numHistsA, histsA, numHistsB, histsB, imageTitleString, addBeforeExtension(imageFileName, "_Residuals"), lowRange, highRange, normalizedResidualsFlag);
	stringstream theSStream;
	theSStream.precision(2);
	TString titleString;
	double integral, error;

	for (int i = 0; i < numHistsA; i++)
	{
		integral = histsA[i]->IntegralAndError(1, histsA[i]->GetNbinsX(), error);
		theSStream << histsA[i]->GetTitle() << ": " << scientific << integral << " #pm " << scientific << error;
		titleString = theSStream.str();
		histsA[i]->SetTitle(titleString);
		theSStream.str("");
	}

	for (int i = 0; i < numHistsB; i++)
	{
		integral = histsB[i]->IntegralAndError(1, histsB[i]->GetNbinsX(), error);
		theSStream << histsB[i]->GetTitle() << ": " << scientific << integral << " #pm " << scientific << error;
		titleString = theSStream.str();
		histsB[i]->SetTitle(titleString);
		theSStream.str("");
	}

	plotExpVersusMCToImage(numHistsA, histsA, numHistsB, histsB, imageTitleString, imageFileName);
}

TString getBranchVarType(const TBranch* inpBranch)
{
	TString branchWVar = inpBranch->GetTitle();
	Ssiz_t pos = branchWVar.Last('/');
	if(pos == kNPOS) return "";
	return branchWVar(pos + 1, branchWVar.Length() - pos - 1);
}

void changeGainOfHist(TH1D* inpHist, double gainFactor)
{
	TAxis* theAxis = inpHist->GetXaxis();
	theAxis->Set(theAxis->GetNbins(), theAxis->GetXmin(), theAxis->GetXmax() * gainFactor);
}

bool checkTolerance(double target, double variable, double tolerance)
{
	return abs((target - variable) / target) < tolerance;
}

int fitHistWithGaus(TH1D* inpHist, TF1* gaussFuncToFit)
{
	int fitStatus;
	double guessMean = gaussFuncToFit->GetParameter(1);
	double guessXmin = gaussFuncToFit->GetXmin();
	fitStatus = inpHist->Fit(gaussFuncToFit, "IENMR");        //Fit once on imperfect range

	gaussFuncToFit->SetRange(gaussFuncToFit->GetParameter(1) - (guessMean - guessXmin), gaussFuncToFit->GetParameter(1) + (guessMean - guessXmin));
	fitStatus = inpHist->Fit(gaussFuncToFit, "IENMR");        //Fit on a better range

	if((fitStatus != 0 && fitStatus != 4000))
	{
		cout << "Error in fit for peak " << gaussFuncToFit->GetTitle() << " in hist: " << inpHist->GetTitle() << endl;
	}

	return fitStatus;

}

double getBasicErrorProgation(double x, double y, double xErr, double yErr, BASIC_MATH_OPERATION inpOperation)
{
	if(inpOperation == ADDITION || inpOperation == SUBTRACTION)
	{
		return sqrt(xErr * xErr + yErr * yErr);
	}
	else if(inpOperation == MULTIPLICATION)
	{
		return (x * y) * sqrt(pow(xErr / x, 2) + pow(yErr / y, 2));
	}
	else if(inpOperation == DIVISION)
	{
		return (x / y) * sqrt(pow(xErr / x, 2) + pow(yErr / y, 2));
	}
	else
	{
		return 0;
	}
}

TH1* calcAvgHistWithErrorFromHists(const int numHists, TH1** inpHists, TString histName)
{
	TH1* outHist = (TH1*) inpHists[0]->Clone();
	outHist->SetName(histName);
	outHist->Reset();
	outHist->Sumw2();
	for (int i = 1; i <= inpHists[0]->GetNbinsX(); i++) //Loop over bins
	{
		double binCounts[numHists];
		for (int j = 0; j < numHists; j++) //loop over hists
		{
			binCounts[j] = inpHists[j]->GetBinContent(i);

		}
		double mean = TMath::Mean(numHists, binCounts);
		double error = TMath::RMS(numHists, binCounts);

		outHist->SetBinContent(i, mean);
		outHist->SetBinError(i, error);
	}
	return outHist;
}

TString getCurrentDateString()
{
	time_t t = time(0);   // get time now
	struct tm * now = localtime(&t);
	TString year = TString::Format("%d", now->tm_year + 1900 - 2000);
	TString month = TString::Format("%d", now->tm_mon + 1);
	TString day = TString::Format("%d", now->tm_mday);
	if(month.Length() == 1) month = "0" + month;
	if(day.Length() == 1) day = "0" + day;

	TString str = year + month + day;
	return str;
}

void createBlankThesisPlot(TString imageName, TString imageTitleString, TString legendList)
{
	gROOT->cd();

	int numLegend = numItemsInStringList(legendList);

	TH1D blankHist("blankHist", imageTitleString, 100, -3, 3);
	blankHist.FillRandom("gaus", 10000);
	TCanvas* c2 = new TCanvas(imageTitleString, "blankThesisPlot", 1200, 960);
	c2->SetFillColor(kGray);
	c2->cd();
	c2->SetLeftMargin(.09);
	c2->SetRightMargin(.05);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetFillColor(kWhite);

	TLegend *theLegend;
	double legendSize = .1 * numLegend;
	if(legendSize < .15) legendSize = .15;
	theLegend = new TLegend(0.55, .9 - legendSize, 0.8705314, 0.8381924);

	for (int i = 0; i < numLegend; i++)
	{
		TString legendString;
		stringFromList(legendList, i, legendString);
		theLegend->AddEntry(&blankHist, legendString);
	}

	blankHist.Draw("P E1 X0");

	theLegend->Draw();
	theLegend->SetBorderSize(1);
	theLegend->SetFillColor(kWhite);

	TPaveText* textBox = new TPaveText(.14, .945, .8705314, .9825073, "NDC");
	textBox->SetTextSize(.035);
	textBox->AddText(blankHist.GetTitle());
	textBox->SetFillColor(kWhite);
	textBox->SetLineColor(kWhite);

	textBox->SetShadowColor(kWhite);
	// textBox->SetFillStyle(4000);//transparent
	// textBox->SetLineStyle(4000);//transparent
	textBox->SetLineWidth(0);
	textBox->Draw();
	blankHist.SetTitle("");

	c2->Update();

	TString imagePath = "/home/mjbales/school/rdk/thesis/matt/pics/";
	imagePath += imageName;
	c2->SaveAs(imagePath);
	delete c2;
}

double GetIdealMaxForPlotting(TH1* inpHist, double sigmaLimit)
{
	double max = numeric_limits<double>::min();
	double maxError = 0;
	double meanError;
	double mean = getCountMeanMinusOutliers(inpHist, sigmaLimit, meanError);
	for (int i = 1; i <= inpHist->GetNbinsX(); i++)
	{
		if(inpHist->GetBinContent(i) > max && inpHist->GetBinContent(i) < mean + sigmaLimit * meanError)
		{
			max = inpHist->GetBinContent(i);
			maxError = inpHist->GetBinError(i);
		}
	}
	return max + maxError;
}

double GetIdealMinForPlotting(TH1* inpHist, double sigmaLimit)
{
	double min = numeric_limits<double>::max();
	double minError = 0;
	double meanError;
	double mean = getCountMeanMinusOutliers(inpHist, sigmaLimit, meanError);
	for (int i = 1; i <= inpHist->GetNbinsX(); i++)
	{
		if(inpHist->GetBinContent(i) < min && inpHist->GetBinContent(i) > mean - sigmaLimit * meanError)
		{
			min = inpHist->GetBinContent(i);
			minError = inpHist->GetBinError(i);
		}
	}
	return min - minError;
}

double getCountMeanMinusOutliers(TH1* inpHist, double sigmaLimit, double& meanError)
{
	const int n = inpHist->GetNbinsX();
	bool outliers[n];
	for (int i = 1; i <= n; i++)
	{
		double minusIthMean = 0;
		double minusIthMeanError = 0;
		for (int j = 1; j <= n; j++) //Calc mean minus the ith bin
		{
			if(j != i)
			{
				minusIthMean += inpHist->GetBinContent(j);
			}
		}
		minusIthMean /= n;
		for (int j = 1; j <= n; j++) //Calc mean minus the ith bin
		{
			if(j != i)
			{
				minusIthMeanError += pow(inpHist->GetBinContent(j) - minusIthMean, 2);
			}
		}
		minusIthMeanError = sqrt(minusIthMeanError / n);

		double theContent = inpHist->GetBinContent(i);

		if(abs(theContent) < minusIthMean + minusIthMeanError * sigmaLimit)
		{
			outliers[i - 1] = false;
		}
		else
		{
			outliers[i - 1] = true;
		}
	}

	double mean = 0;
	meanError = 0;
	for (int i = 1; i <= n; i++)
	{
		if(!outliers[i - 1])
		{
			mean += inpHist->GetBinContent(i);
		}

	}
	mean /= n;

	for (int i = 1; i <= n; i++)
	{
		if(!outliers[i - 1])
		{
			meanError += pow(inpHist->GetBinContent(i) - mean, 2);
		}

	}
	meanError = sqrt(meanError / n);

	return mean;
}

TCanvas* getThesisPlotCanvas(TString name)
{
	TCanvas* theCanvas = new TCanvas(name, "thesisPlot", 10, 10, 1200, 700);
	theCanvas->SetLeftMargin(.09);
	theCanvas->SetRightMargin(.05);
	theCanvas->SetTopMargin(.05);
	theCanvas->SetBottomMargin(.12);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetFillColor(kWhite);
	return theCanvas;
}

TCanvas* getHalfPresentationPlot(TString name)
{
	TCanvas* theCanvas = new TCanvas(name, "thesisPlot2", 10, 10, 600, 700);
	theCanvas->SetLeftMargin(.15);
	theCanvas->SetRightMargin(.05);
	theCanvas->SetTopMargin(.05);
	theCanvas->SetBottomMargin(.12);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetFillColor(kWhite);
	return theCanvas;
}

TCanvas* getHalfPresentationPlot2(TString name)
{
	TCanvas* theCanvas = new TCanvas(name, "thesisPlot3", 10, 10, 1200, 500);
	theCanvas->SetLeftMargin(.09);
	theCanvas->SetRightMargin(.05);
	theCanvas->SetTopMargin(.05);
	theCanvas->SetBottomMargin(.12);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetFillColor(kWhite);
	return theCanvas;
}

TH1* trimHist(TH1* inpHist, TString newName, int start, int end)
{
	int oBins = inpHist->GetNbinsX();

	int numTrim = start + end;

	if(oBins <= numTrim)
	{
		return nullptr;
	}

	//determine new edges

	const int numEdges = oBins - numTrim + 1;

	double edges[numEdges];

	for (int i = 1; i < numEdges; i++)
	{
		edges[i - 1] = inpHist->GetBinLowEdge(i + start);
	}

	edges[numEdges - 1] = inpHist->GetBinLowEdge(oBins - end) + inpHist->GetBinWidth(oBins - end);

	return inpHist->Rebin(numEdges - 1, newName, edges);
}

TH1D* makeTH1DFromDim(TString name, TString title, HistDim inpHistDim)
{
	return new TH1D(name, title, inpHistDim.numBins, inpHistDim.binLowEdge, inpHistDim.binHighEdge);
}

TH2D* makeTH2DFromDim(TString name, TString title, HistDim inpXHistDim, HistDim inpYHistDim)
{
	return new TH2D(name, title, inpXHistDim.numBins, inpXHistDim.binLowEdge, inpXHistDim.binHighEdge, inpYHistDim.numBins, inpYHistDim.binLowEdge, inpYHistDim.binHighEdge);
}

TH3D* makeTH3DFromDim(TString name, TString title, HistDim inpXHistDim, HistDim inpYHistDim, HistDim inpZHistDim)
{
	return new TH3D(name, title, inpXHistDim.numBins, inpXHistDim.binLowEdge, inpXHistDim.binHighEdge, inpYHistDim.numBins, inpYHistDim.binLowEdge, inpYHistDim.binHighEdge, inpZHistDim.numBins, inpZHistDim.binLowEdge, inpZHistDim.binHighEdge);
}

HistDim getDimFromHist(TH1* inpHist)
{
	HistDim outDim;
	outDim.numBins = inpHist->GetNbinsX();
	outDim.binLowEdge = inpHist->GetBinLowEdge(1);
	outDim.binHighEdge = inpHist->GetBinWidth(outDim.numBins) + inpHist->GetBinLowEdge(outDim.numBins);
	return outDim;
}

int getVarArrayFromHist(TH1* inpHist, double** outArray)
{
	const int numBins = inpHist->GetNbinsX();
	*outArray = new double[numBins + 1];

	for (int i = 0; i < numBins; i++)
	{
		(*outArray)[i] = inpHist->GetBinLowEdge(i + 1);

	}
	(*outArray)[numBins] = inpHist->GetBinLowEdge(numBins) + inpHist->GetBinWidth(numBins);

	return numBins;
}

//Creates a new histogram with the same data and bin count but with a new start and end
TH1D* adjustXScaleHist(TH1D* inpHist, double newStart, double newEnd, TString adjustedHistName)
{
	if(adjustedHistName == "")
	{
		adjustedHistName = TString(inpHist->GetName()) + "_adjusted";
	}
	int binCount = inpHist->GetNbinsX();
	TH1D* outHist = new TH1D(inpHist->GetName() + TString("_adjusted"), inpHist->GetTitle(), binCount, newStart, newEnd);
	for (int i = 0; i < binCount; ++i)
	{
		outHist->SetBinContent(i + 1, inpHist->GetBinContent(i + 1));
	}
	return outHist;
}

//Creates a new histogram with the same data and bin count but recenters the histogram by a number of bins amount. Assumes constant bin width
//numBinsShift can be postive or negative
TH1D* recenterPeriodicHist(TH1D* inpHist, int numBinsShift, TString adjustedHistName)
{
	if(adjustedHistName == "")
	{
		adjustedHistName = TString(inpHist->GetName()) + "_adjusted";
	}
	int binCount = inpHist->GetNbinsX();
	int correctedBinShift = numBinsShift - binCount * (numBinsShift / binCount);  //Note integer math
	double binWidth = inpHist->GetBinWidth(1); //Assumes constant bin width
	double shift = binWidth * correctedBinShift;
	double newStart = getHistExtentLow(inpHist) + shift;
	double newEnd = getHistExtentHigh(inpHist) + shift;
	TH1D* outHist = new TH1D(inpHist->GetName() + TString("_adjusted"), inpHist->GetTitle(), binCount, newStart, newEnd);
	outHist->GetXaxis()->SetTitle(inpHist->GetXaxis()->GetTitle());
	outHist->GetYaxis()->SetTitle(inpHist->GetYaxis()->GetTitle());
	for (int i = 0; i < binCount; ++i)
	{ //Loop goes through the inpHist and sets the appropriate loc in out hist
		int j = i - correctedBinShift;
		if(j <= 0)
		{
			j += binCount;
		}
		else if(j >= binCount)
		{
			j -= binCount;
		}

		outHist->SetBinContent(j + 1, inpHist->GetBinContent(i + 1));
	}

	return outHist;

}

double reducePeriodicNumber(double inp, double start, double end)
{
	double period = end - start;
	double answer = inp;
	double fractPart, intPart;

	fractPart = modf((inp - start) / period, &intPart);

	if(inp > end)
	{
		answer = fractPart * period + start;
	}
	else if(inp < start)
	{
		answer = fractPart * period + end;
	}
	return answer;
}

TVector3 getRandomDirection()
{

	double phi = gRandom->Rndm() * 2. * TMath::Pi();
	double cosTheta = gRandom->Rndm() * 2. - 1.;
	double sinTheta = sqrt(1. - cosTheta * cosTheta);
	return TVector3(sinTheta * cos(phi), sinTheta * sin(phi), cosTheta);
}

TVector2 getRandomDirection2D()
{

	double phi = gRandom->Rndm() * 2. * TMath::Pi();
	return TVector2(cos(phi), sin(phi));
}

TVector3 getRandomPointInCylinder(double radius, double height)
{
	double r = sqrt(gRandom->Rndm()) * radius;
	double theta = gRandom->Rndm() * 2. * TMath::Pi();
	double z = (gRandom->Rndm() - 0.5) * height;
	return TVector3(r * cos(theta), r * sin(theta), z);
}

TVector2 getRandomPointInCircle(double radius)
{
	double r = sqrt(gRandom->Rndm()) * radius;
	double theta = gRandom->Rndm() * 2. * TMath::Pi();
	return TVector2(r * cos(theta), r * sin(theta));
}

double meanArray(int numData, double* theData)
{
	double mean = 0;
	for (int i = 0; i < numData; ++i)
	{
		mean += theData[i];
	}
	mean /= numData;
	return mean;
}

double stDevArray(int numData, double* theData, bool useBesselCorrection)
{
	double mean = meanArray(numData, theData);
	double stdev = 0;

	for (int i = 0; i < numData; ++i)
	{
		stdev += pow(theData[i] - mean, 2);
	}

	if(useBesselCorrection)
	{
		stdev /= numData - 1;
	}
	else
	{
		stdev /= numData;
	}

	stdev = sqrt(stdev);

	return stdev;

}

double absMinArray(int numData, double* theData)
{
	double min = DBL_MAX;
	for (int i = 0; i < numData; ++i)
	{
		if(theData[i] < min)
		{
			min = theData[i];
		}
	}

	return min;
}

double absMaxArray(int numData, double* theData)
{
	double max = 0;
	for (int i = 0; i < numData; ++i)
	{
		if(theData[i] > max)
		{
			max = theData[i];
		}
	}

	return max;
}

//Default y axis unles forX=true;
double maxFromTGraphs(vector<TGraph*> theGraphs, bool includeError, bool forX)
{
	double theMax = -DBL_MAX;
	for (auto aGraph : theGraphs)
	{
		double* vals;
		double* errors;
		if(forX)
		{
			vals = aGraph->GetX();
			errors = aGraph->GetEX();
		}
		else
		{
			vals = aGraph->GetY();
			errors = aGraph->GetEY();
		}
		for (int i = 0; i < aGraph->GetN(); i++)
		{
			if(includeError)
				theMax = max(vals[i] + errors[i], theMax);
			else
				theMax = max(vals[i], theMax);
		}
	}
	return theMax;
}

double minFromTGraphs(std::vector<TGraph*> theGraphs, bool includeError, bool forX)
{
	double theMin = DBL_MAX;
	for (auto aGraph : theGraphs)
	{
		double* vals;
		double* errors;
		if(forX)
		{
			vals = aGraph->GetX();
			errors = aGraph->GetEX();
		}
		else
		{
			vals = aGraph->GetY();
			errors = aGraph->GetEY();
		}
		for (int i = 0; i < aGraph->GetN(); i++)
		{
			if(includeError)
				theMin = min(vals[i] - errors[i], theMin);
			else
				theMin = min(vals[i], theMin);
		}
	}
	return theMin;
}

TH1I* histFromHist(TH1I* inpHist)
{
	int max = inpHist->GetMaximum();
	int min = inpHist->GetMinimum();
	TH1I* newHist = new TH1I(TString("histOf_") + inpHist->GetName(), TString("Hist of ") + inpHist->GetTitle(), 1000, min, max);
	for (int i = 0; i < inpHist->GetNbinsX(); ++i)
	{
		newHist->Fill(inpHist->GetBinContent(i + 1));
	}
	return newHist;
}

void scaleTGraphErrors(TGraphErrors* inpGraph, double scaleFactor)
{
	double* y = inpGraph->GetY();
	double* eY = inpGraph->GetEY();
	for (int i = 0; i < inpGraph->GetN(); ++i)
	{
		y[i] *= scaleFactor;
		eY[i] *= scaleFactor;
	}
}

double* makeEvenSampleOverLog(const int numPoints, const double begin, const double end)
{
	double logBegin = log10(begin);
	double logEnd = log10(end);
	double logIncrement = (logEnd - logBegin) / (numPoints - 1);

	double* theArray = new double[numPoints];

	for (int i = 0; i < numPoints; ++i)
	{
		theArray[i] = pow(10, logBegin + logIncrement * i);
	}
	return theArray;
}

double meanVector(const vector<double>& theData)
{
	double mean = 0;
	for (double x : theData)
	{
		mean += x;
	}
	mean /= theData.size();
	return mean;
}

//Worried about precision
double carefulMeanVector(const vector<double>& theData)
{
	//First calc mean to use in subtraction
	double tempMean = meanVector(theData);

	double mean = 0;
	for (double x : theData)
	{
		mean += (x - tempMean);
	}
	mean /= theData.size();
	mean += tempMean;
	return mean;
}

double stDevVector(const vector<double>& theData, const bool useBesselCorrection)
{
	double mean = meanVector(theData);
	double stdev = 0;

	for (double x : theData)
	{
		stdev += pow(x - mean, 2);
	}

	if(useBesselCorrection)
	{
		stdev /= theData.size() - 1;
	}
	else
	{
		stdev /= theData.size();
	}

	stdev = sqrt(stdev);

	return stdev;

}

double carefullStDevVector(const vector<double>& theData, const bool useBesselCorrection)
{
	double mean = carefulMeanVector(theData);
	double stdev = 0;

	for (double x : theData)
	{
		stdev += pow(x - mean, 2);
	}

	if(useBesselCorrection)
	{
		stdev /= theData.size() - 1;
	}
	else
	{
		stdev /= theData.size();
	}

	stdev = sqrt(stdev);

	return stdev;

}

double minVector(const vector<double>& theData)
{
	return *(min_element(theData.begin(), theData.end()));
}

double maxVector(const vector<double>& theData)
{
	return *(max_element(theData.begin(), theData.end()));
}

double reducePeriodicToMeanInVector(vector<double>& theData)
{
	double tempMean = carefulMeanVector(theData);
	double mean = 0;
	for (double& x : theData)
	{
		x = reducePeriodicNumber(x, -TMath::Pi() + tempMean, TMath::Pi() + tempMean);
		mean += x;
	}
	mean /= theData.size();

	return mean;

}

TH1* makeVarHistScaledByWidth(TH1* inpHist, vector<double> binSpacing)
{
	TH1* outHist=inpHist->Rebin(binSpacing.size()-1,inpHist->GetName()+TString("VarHist"),binSpacing.data());

	for(int i = 1;i <= outHist->GetNbinsX();i++)
	{
		double w=outHist->GetBinWidth(i);
		outHist->SetBinContent(i,outHist->GetBinContent(i)/w);
		outHist->SetBinError(i,outHist->GetBinError(i)/w);
	}

	return outHist;
}

