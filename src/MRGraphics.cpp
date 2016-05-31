#include "MRGraphics.h"

//Standard Libraries
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cctype>
#include <algorithm>
#include <sys/stat.h>
#include <cmath>

//Root Libraries
#include "TSystem.h"
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TStyle.h"
#include "TPad.h"
#include "TROOT.h"
#include "TPaveText.h"
#include "TLegend.h"

#include "MRMisc.h"
#include "MRText.h"
#include "MRIO.h"

using namespace std;

TCanvas* plotTGraphs(std::vector<TGraph*> theGraphs, bool autoFormat, TString titleString, bool logX, bool logY)
{
	gROOT->cd();

	TCanvas* theCanvas = new TCanvas(titleString, titleString, 1600, 1200);
	theCanvas->cd();
	theCanvas->SetLeftMargin(.13);
	theCanvas->SetRightMargin(.06);
	theCanvas->SetGrid(1, 1);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetLogx(logX);
	gPad->SetLogy(logY);

	//Create frame based on the min and max from all the histograms;
	double minX = minFromTGraphs(theGraphs, true, true);
	double minY = minFromTGraphs(theGraphs, true, false);
	double maxX = maxFromTGraphs(theGraphs, true, true);
	double maxY = maxFromTGraphs(theGraphs, true, false);
	double lengthX = maxX - minX;
	double lengthY = maxY - minY;
	double widen = 0.05;
	if(!logX)
	{
		minX -= widen * lengthX;
		maxX += widen * lengthX;
	}
	else
	{
		minX *= widen;
		maxX /= widen;
	}

	if(!logY)
	{
		minY -= widen * lengthY;
		maxY += widen * lengthY;
	}
	else
	{
		minY *= widen;
		maxY /= widen;
	}
	TH1F* theFrameHist = theCanvas->DrawFrame(minX, minY, maxX, maxY);
	theFrameHist->SetTitleOffset(1.5, "y");
	theFrameHist->SetTitleOffset(1.3, "X");
	theFrameHist->SetTitle(titleString);

	//If it doesn't have a title string use the first histogram
	if(TString(theFrameHist->GetXaxis()->GetTitle()) == "")
	{
		theFrameHist->SetXTitle(theGraphs[0]->GetXaxis()->GetTitle());
		theFrameHist->SetYTitle(theGraphs[0]->GetYaxis()->GetTitle());
	}

	theCanvas->Update();
	theCanvas->Modified();

	TLegend* theLegend = nullptr;

	if(theGraphs.size() > 1)
	{
		double legendSize = .05 * theGraphs.size();
		if(legendSize < .15) legendSize = .15;

//		theLegend = new TLegend(0.55, .8 - legendSize, 0.87, 0.8);  //Top Right
//		theLegend = new TLegend(0.2, .8 - legendSize, 0.47, 0.8);  //Top Left
		theLegend = new TLegend(0.55, .4 - legendSize, 0.87, 0.4); //Bottom Right
	}

	for (unsigned int i = 0; i < theGraphs.size(); i++)
	{
		if(autoFormat)
		{
			theGraphs[i]->SetLineColor(MR_GRAPH_COLOR_LIST[i % MR_GRAPH_NUMCOLORS]);
			theGraphs[i]->SetMarkerColor(MR_GRAPH_COLOR_LIST[i % MR_GRAPH_NUMCOLORS]);
			theGraphs[i]->SetMarkerStyle(MR_MARKER_STYLE_LIST[i % MR_MARKER_NUMSTYLES]);
		}
		theGraphs[i]->Draw("P E1 same");
		if(theLegend != nullptr)
		{
			TString legendString = theGraphs[i]->GetTitle();
			legendString.Resize(legendString.First(';'));
			theLegend->AddEntry(theGraphs[i], legendString, "P E1");
		}
	}

	if(theLegend != nullptr) theLegend->Draw("same");

	return theCanvas;
}

void plotTGraphs(std::vector<TGraph*> theGraphs, bool autoFormat, TString titleString, TString imagePath, bool logX, bool logY)
{
	TCanvas* theCanvas = plotTGraphs(theGraphs, autoFormat, titleString, logX, logY);
	theCanvas->SaveAs(imagePath);
	delete theCanvas;

}

void plotTGraphs(std::vector<TString> filePaths, TString titleString, TString imagePath, bool logX, bool logY)
{
	vector<TGraph*> theGraphs;

	for (auto aPath : filePaths)
	{
		theGraphs.push_back((TGraph*) getTabSeperatedTGraphErrors(aPath));
	}

	plotTGraphs(theGraphs, true, titleString, imagePath, logX, logY);
}

void plotTGraph(TString filePath, TString titleString, TString imagePath, bool logX, bool logY)
{
	vector<TString> thePaths = { filePath };
	plotTGraphs(thePaths, titleString, imagePath, logX, logY);
}

TCanvas* plotExpVersusMC(int numExpHists, TH1** expHists, int numMCHists, TH1** mcHists, TString titleString)
{
	gROOT->cd();
	int totalNumHists = numExpHists + numMCHists;

	if(totalNumHists < 1 || (expHists == nullptr && mcHists == nullptr))
	{
		return nullptr;
	}

	const int numColors = 11;
	const int colorList[numColors] = { kBlack, kBlue, kRed, kGreen, kCyan, kMagenta, kPink, kOrange, kViolet, kAzure, kYellow };
	const int markerList[14] = { 20, 22, 23, 25, 24, 32, 33, 3, 5, 28, 26, 30, 29 };
	//gStyle->SetTitleColor(kWhite);
	TCanvas* c2 = new TCanvas(titleString, "plotManyHistCanvas", 1600, 1200);
	c2->SetFillColor(kGray);
	c2->cd();
	c2->SetLeftMargin(.13);
	c2->SetRightMargin(.06);
	c2->SetGrid(1, 1);
	gPad->SetTickx(1);
	gPad->SetTicky(1);
	gPad->SetFillColor(kWhite);

	TLegend *theLegend;

	if(totalNumHists >= 1)
	{
		double legendSize = .1 * totalNumHists;
		if(legendSize < .15) legendSize = .15;
		theLegend = new TLegend(0.55, .8 - legendSize, 0.87, 0.8);  //Top Right
//        theLegend = new TLegend(0.2,.8-legendSize,0.47,0.8);  //Top Left
//        theLegend = new TLegend(0.55,.4-legendSize,0.87,0.4); //Bottom Right
	}

	//Prep MCHists
	for (int i = 0; i < numMCHists; i++)
	{
		mcHists[i]->SetLineColor(colorList[(i) % numColors]);
		mcHists[i]->SetMarkerColor(colorList[(i) % numColors]);
		mcHists[i]->SetLineStyle((i - numExpHists) % 5 + 1);
		mcHists[i]->SetLineWidth(2);
		mcHists[i]->GetYaxis()->SetTitleOffset(1.5);
	}

	//Prep MCHists
	for (int i = 0; i < numExpHists; i++)
	{
		expHists[i]->SetLineColor(colorList[(i + numMCHists) % numColors]);
		expHists[i]->SetMarkerColor(colorList[(i + numMCHists) % numColors]);
		expHists[i]->SetMarkerStyle(markerList[(i) % 14]);
		expHists[i]->SetLineWidth(2);
		expHists[i]->SetMarkerSize(2);
		expHists[i]->GetYaxis()->SetTitleOffset(1.5);
	}

	bool drawExpFirst = true;

	for (int i = 0; i < totalNumHists; i++)
	{
		if(drawExpFirst && numExpHists > 0)
		{
			if(i < numExpHists)
			{
				theLegend->AddEntry(expHists[i], expHists[i]->GetTitle(), "PE");
				if(i == 0)
				{
					expHists[i]->SetTitle("");
					expHists[i]->Draw("P E1 X0");
				}
				else
				{
					expHists[i]->Draw("P E1 X0 SAME");
				}
			}
			else
			{
				theLegend->AddEntry(mcHists[i - numExpHists], mcHists[i - numExpHists]->GetTitle(), "L");
				mcHists[i - numExpHists]->Draw("HIST SAME");
			}
		}
		else
		{
			if(i < numMCHists)
			{
				theLegend->AddEntry(mcHists[i], mcHists[i]->GetTitle(), "L");
				if(i == 0)
				{
					mcHists[i]->SetTitle("");
					mcHists[i]->Draw("HIST");
				}
				else
				{
					mcHists[i]->Draw("HIST SAME");
				}
			}
			else
			{
				theLegend->AddEntry(expHists[i - numMCHists], expHists[i - numMCHists]->GetTitle(), "PE");
				expHists[i - numMCHists]->Draw("P E1 X0 SAME");
			}
		}
	}

	if(totalNumHists >= 1)
	{
		theLegend->Draw();
		theLegend->SetBorderSize(1);
		theLegend->SetFillColor(kWhite);
	}

	TPaveText* textBox = new TPaveText(.14, .945, .8705314, .9825073, "NDC");
	textBox->SetTextSize(.035);
	textBox->AddText(titleString);
	textBox->SetFillColor(kWhite);
	textBox->SetLineColor(kWhite);

	textBox->SetShadowColor(kWhite);
	// textBox->SetFillStyle(4000);//transparent
	// textBox->SetLineStyle(4000);//transparent
	textBox->SetLineWidth(0);
	textBox->Draw();

	c2->Update();

	return c2;
}

void plotExpVersusMCToImage(int numExpHists, TH1** expHists, int numMCHists, TH1** mcHists, TString titleString, TString imagePath, bool useLog)
{
	TString tempString;
	if(numMCHists > 0)
		tempString = mcHists[0]->GetTitle();
	else
		tempString = expHists[0]->GetTitle();
	TCanvas* theCanvas = plotExpVersusMC(numExpHists, expHists, numMCHists, mcHists, titleString);
	if(useLog) theCanvas->SetLogy();

	cout << "Saving image to " << imagePath << endl;
	theCanvas->SaveAs(imagePath);

	delete theCanvas;
	if(numMCHists > 0)
		mcHists[0]->SetTitle(tempString);
	else
		expHists[0]->SetTitle(tempString);

}

void plotExpVersusMCToImage(TString fileNameExp, TString fileNameMC, HistDim inpDim, TString titleString, TString imagePath, bool useLog)
{
	TH1* expHists[1];
	TH1* mcHists[1];

	expHists[0] = getTabSeperatedHist(fileNameExp, inpDim);
	expHists[0]->SetTitle("Exp");
	expHists[0]->SetName("Exp");
	mcHists[0] = getTabSeperatedHist(fileNameMC, inpDim);
	mcHists[0]->SetTitle("MC");
	expHists[0]->SetName("MC");

	expHists[0]->Rebin(5);
	mcHists[0]->Rebin(5);

	plotExpVersusMCToImage(1, expHists, 1, mcHists, titleString, imagePath, useLog);
	delete expHists[0];
	delete mcHists[0];

}
