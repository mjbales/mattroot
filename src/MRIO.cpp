#include "MRIO.h"
#include <algorithm>

//Standard Libraries
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits>
#include <vector>
#include <sstream>
#include <iomanip>
#include <sys/types.h>

//ROOT Libraries
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TRandom3.h"
#include "TList.h"
#include "TString.h"
#include "TGraphErrors.h"

//Matt ROOT Libraries
#include "MRText.h"
#include "MRPhys.h"
#include "MRMisc.h"
#include "MRHistDim.h"

using namespace std;

void convertTH1ToTXT(TH1* inpHist, TString outputFileName, bool printWidths)
{
	ofstream outFile;
	outFile.open(outputFileName);
	double binCenter, data, error, width;
	if(printWidths)
	{
		outFile << "#Bin Center\tBin Width\tBin Content\tBin Error" << endl;
	}
	else
	{
		outFile << "#Bin Center\tBin Content\tBin Error" << endl;
	}
	for (int i = 1; i <= inpHist->GetNbinsX(); i++)
	{
		outFile.precision(15);
		data = inpHist->GetBinContent(i);
		binCenter = inpHist->GetBinCenter(i);
		error = inpHist->GetBinError(i);
		width=inpHist->GetBinWidth(i);
		if(printWidths)
		{
			outFile << scientific << binCenter << "\t" << scientific << width << "\t" << scientific << data << "\t" << scientific << error;
		}
		else
		{
			outFile << scientific << binCenter << "\t" << scientific << data << "\t" << scientific << error;
		}
		if(i != inpHist->GetNbinsX()) outFile << endl;
	}
	outFile.close();
}

void convertXYDataWErrorToTXT(int numData, double* x, double* y, double* eX, double* eY, TString titleString, TString outputFileName)
{
	ofstream outFile;
	outFile.open(outputFileName);

	outFile << titleString << endl;
	for (int i = 0; i < numData; i++)
	{
		outFile.precision(15);
		outFile << scientific << x[i];
		if(eX != nullptr)
		{
			outFile << "\t" << scientific << eX[i];
		}
		else
		{
			outFile << "\t" << scientific << 0;
		}

		outFile << "\t" << scientific << y[i];
		if(eY != nullptr)
		{
			outFile << "\t" << scientific << eY[i];
		}
		else
		{
			outFile << "\t" << scientific << 0;
		}

		if(i != numData - 1)
		{
			outFile << endl;
		}

	}
	outFile.close();
}

void convertTGraphErrorsToTXT(TGraphErrors* inpGraph, TString outputFileName)
{
	ofstream outFile;
	outFile.open(outputFileName);

	double* x = inpGraph->GetX();
	double* y = inpGraph->GetY();
	double* eX = inpGraph->GetEX();
	double* eY = inpGraph->GetEY();
	TString titleString = TString("#") + inpGraph->GetTitle() + ";" + inpGraph->GetXaxis()->GetTitle() + ";" + inpGraph->GetYaxis()->GetTitle();
	convertXYDataWErrorToTXT(inpGraph->GetN(), x, y, eX, eY, titleString, outputFileName);

}

TH1D* getTabSeperatedHist(TString fileName, HistDim inpHistDim, char delim)
{
	return getTabSeperatedHist(fileName, inpHistDim.numBins, inpHistDim.binLowEdge, inpHistDim.binHighEdge, delim);
}

TH1D* getTabSeperatedHist(TString fileName, int numBins, double start, double end, char delim)
{
	TH1D* histPrototype = new TH1D("histPrototype", "histPrototype", numBins, start, end);
	TH1D* outHist = getTabSeperatedHist(fileName, histPrototype, delim);
	delete histPrototype;
	return outHist;
}

TH1D* getTabSeperatedHist(TString fileName, TH1D* histPrototype, char delim)
{
	if(!FileExists(fileName))
	{
		cout << "Hist: " << fileName << " does not exist." << endl;
		return nullptr;
	}
	ifstream inpFile;
	inpFile.open(fileName);

	TH1D* theHist = (TH1D*) histPrototype->Clone(fileName);
	theHist->Sumw2();
	theHist->SetTitle(fileName);

	TString temp, line;

	double binCenter, counts, unc;
	for (int i = 0; i < histPrototype->GetNbinsX(); i++)
	{

		getNonCommentLine(inpFile, line, delim);
		stringstream theSS(line.Data());
		theSS >> binCenter >> counts >> unc;

		theHist->SetBinContent(i + 1, counts);
		theHist->SetBinError(i + 1, unc);
	}

	return theHist;
}

//IF IT USES COUT!  NOT PRINTF
void convertTObjectPrint2StringStream(TObject* theObject, stringstream* theStream)
{
	if(theObject == nullptr || theStream == nullptr) return;

	streambuf *psbuf, *backup;

	backup = cout.rdbuf();     // back up cout's streambuf

	psbuf = theStream->rdbuf();   // get file's streambuf
	cout.rdbuf(psbuf);         // assign streambuf to cout
	cout << "Test Line 2" << endl;

	theObject->Print();

	cout.rdbuf(backup);        // restore cout's original streambuf
}

//Returns branch
void convertTTree2Text(TString rootFileName, TString treeName, TString outputFileName, int numEntriesToRun)
{
	TFile* rootFile = new TFile(rootFileName, "read");
	if(rootFile->IsZombie())
	{
		cout << "Error in opening " << rootFileName << endl;
		return;
	}

	TTree* theTree = (TTree*) rootFile->Get(treeName);

	int numToActuallyRun;
	if(numEntriesToRun > theTree->GetEntries())
		numToActuallyRun = theTree->GetEntries();
	else
		numToActuallyRun = numEntriesToRun;

	const int numBranches = theTree->GetNbranches();

	TBranch* branchArray[numBranches];
	TString branchNameArray[numBranches];
	Char_t branchCharArray[numBranches];

	int intArray[numBranches];
	double doubleArray[numBranches];
	Char_t charArray[numBranches];
	TString branchTitleArray[numBranches];

	ofstream outputFile;
	outputFile.open(outputFileName);
	outputFile << "$ Tree Name: " << theTree->GetName() << endl;
	outputFile << "$ Num Entries: " << numEntriesToRun << endl;
	outputFile << "$ Entry";

	for (int i = 0; i < numBranches; i++)
	{
		branchArray[i] = theTree->GetBranch(((theTree->GetListOfBranches())[0][i])->GetName());
		branchTitleArray[i] = branchArray[i]->GetTitle();
		branchNameArray[i] = branchArray[i]->GetName();
		branchCharArray[i] = branchTitleArray[i][branchTitleArray[i].Length() - 1];

		if(branchCharArray[i] == 'I') //Integer
		{
			theTree->SetBranchAddress(branchNameArray[i], &(intArray[i]));
		}
		else if(branchCharArray[i] == 'D') //Double
		{
			theTree->SetBranchAddress(branchNameArray[i], &(doubleArray[i]));
		}
		else if(branchCharArray[i] == 'B') //Char
		{
			theTree->SetBranchAddress(branchNameArray[i], &(charArray[i]));
		}

		outputFile << "\t" << branchNameArray[i];

	}

	outputFile.precision(10);

	int tempInt;

	for (int j = 0; j < numToActuallyRun; j++)
	{
		outputFile << endl;
		theTree->GetEntry(j);
		outputFile << j << "\t";
		for (int i = 0; i < numBranches; i++)
		{
			if(branchCharArray[i] == 'I') //Integer
			{
				outputFile << intArray[i];
			}
			else if(branchCharArray[i] == 'D') //Double
			{
				outputFile << doubleArray[i];
			}
			else if(branchCharArray[i] == 'B') //Char
			{
				tempInt = charArray[i];
				outputFile << tempInt;
			}
			if(i != (numBranches - 1)) outputFile << "\t";

		}

	}
	outputFile.close();

	rootFile->Close();
	delete rootFile;

}

TH1D* getCSVHist(TString fileName, int numBins, double start, double end, char endLineChar)
{
	if(!FileExists(fileName)) return nullptr;
	ifstream inpFile;
	inpFile.open(fileName);

	TH1D* theHist = new TH1D(fileName, fileName, numBins, start, end);
	theHist->Sumw2();

	string temp, line;

	double counts;
	double binError;

	getline(inpFile, line, endLineChar);

	string binCenterString, countString, errorString;

	for (int i = 0; i < numBins; i++)
	{
		getline(inpFile, binCenterString, ',');
		getline(inpFile, countString, ',');
		getline(inpFile, errorString, endLineChar);
		counts = str2double(countString);
		binError = str2double(errorString);

		theHist->SetBinContent(i, counts);
		theHist->SetBinError(i, binError);
	}

	return theHist;

}

TH1D* getCSVHist2(TString fileName, int numBins, double start, double end)
{
	ifstream inpFile;
	inpFile.open(fileName);

	TH1D* theHist = new TH1D(fileName, fileName, numBins, start, end);

	string temp, line;

	double counts;

	for (int i = 0; i < numBins; i++)
	{
		getline(inpFile, line, ',');
		getline(inpFile, line, ',');

		sscanf(line.data(), "%lf", &counts);
		theHist->SetBinContent(i, counts * 1e7);
		getline(inpFile, line, '\n');
	}
	inpFile.close();
	return theHist;

}

void convertTH2ToTXT(TH2* inpHist, TString outputFilePath)
{
	ofstream outFile;
	outFile.open(outputFilePath);
	double data, error;
	for (int k = 0; k < 2; k++)
	{
		if(k == 0)
		{
			outFile << "#Values:Rows = constant y, columns = constant x  Error histogram follows value histogram" << endl;
		}
		else
		{
			outFile << "#Error:Rows = constant y, columns = constant x  Value histogram precedes value histogram" << endl;
		}
		for (int i = 1; i <= inpHist->GetNbinsY(); i++)
		{
			for (int j = 1; j <= inpHist->GetNbinsX(); j++)
			{
				outFile.precision(15);
				if(k == 0)
				{
					data = inpHist->GetBinContent(j, i);
					outFile << scientific << data;
				}
				else
				{
					error = inpHist->GetBinError(j, i);
					outFile << scientific << error;
				}
				if(j != inpHist->GetNbinsX()) outFile << "\t";
			}
			if(i != inpHist->GetNbinsY() || k == 0) outFile << endl;

		}
	}
	outFile.close();
}

//Max four columns
TGraphErrors* getTabSeperatedTGraphErrors(TString filePath, char delim)
{
	std::vector<double> values[4];  //x,xE,y,yE depending on numColumns;

	if(!FileExists(filePath))
	{
		cout << "File, " << filePath << ", does not exist." << endl;
		return nullptr;
	}

	ifstream inpFile;
	inpFile.open(filePath);

	//If comment line is present presume that it is the title
	TString titleString = fileNameFromFullPath(filePath);
	if(inpFile.peek() == '#')
	{
		titleString = "";
		titleString.ReadToDelim(inpFile, delim);
		titleString.Remove(0, 1);
		titleString = titleString.Strip(TString::EStripType::kTrailing, '\t');
	}

	TString temp, line;
	double value;

	//Get first line to check for column count and then return to beginning;
	getNonCommentLine(inpFile, line, delim);
	const int numColumns = countTString(line, '\t') + 1; //counts how many delims there are to determine number of columns
	inpFile.seekg(0);

	if(numColumns < 2 || numColumns > 4)
	{
		cout << "Data expected in X,Y or X,Y,Yerror   or  X,Xerror,Y,Yerror format!" << endl;
		cout << "Number of columns set, " << numColumns << ", is incorrected." << endl;
		return nullptr;
	}

	//Read data
	for (int i = 0; i < 1000000; i++)  //Arbitrary 1 mill. data point max
	{

		if(getNonCommentLine(inpFile, line, delim))
		{
			break;  //Usually end of file
		}

		stringstream theSS(line.Data());

		theSS >> value;

		values[0].push_back(value); //X

		if(numColumns == 2)
		{
			theSS >> value;
			values[2].push_back(value); //Y;
			values[1].push_back(0); //X error is 0
			values[3].push_back(0); //Y error is 0
		}
		else if(numColumns == 3)
		{
			theSS >> value;
			values[2].push_back(value); //Y;
			theSS >> value;
			values[3].push_back(value); //Y error
			values[1].push_back(0); //X error is 0
		}
		else //4 columns
		{
			theSS >> value;
			values[1].push_back(value); //X error
			theSS >> value;
			values[2].push_back(value); //Y;
			theSS >> value;
			values[3].push_back(value); //Y error
		}
	}

	int numData = values[0].size();

	TGraphErrors* outGraph;

	if(numColumns == 2)
	{
		outGraph = new TGraphErrors(numData, &values[0][0], &values[2][0], nullptr, nullptr);
	}
	else if(numColumns == 3)
	{
		outGraph = new TGraphErrors(numData, &values[0][0], &values[2][0], nullptr, &values[3][0]);
	}
	else if(numColumns == 4)
	{
		outGraph = new TGraphErrors(numData, &values[0][0], &values[2][0], &values[1][0], &values[3][0]);
	}

	outGraph->SetName(fileNameFromFullPath(filePath));
	outGraph->SetTitle(titleString);

	cout << "File Path: " << filePath << endl;
	cout << "Title: " << outGraph->GetTitle() << endl;
	cout << "X Title: " << outGraph->GetXaxis()->GetTitle() << endl;
	cout << "Y Title: " << outGraph->GetYaxis()->GetTitle() << endl;

	return outGraph;
}
