#include "MRText.h"

#include <sstream>
#include <sys/stat.h>

#include "TH1.h"

using namespace std;

using std::stringstream;

int skipCommentLines(ifstream& inpFileStream)
{
	TString sBuffer;
	if(!inpFileStream.is_open() || inpFileStream.fail() || inpFileStream.eof()) return -1;
	while (inpFileStream.peek() == '#' || inpFileStream.peek() == '%')
	{
		sBuffer.ReadLine(inpFileStream);
	}

	if(!inpFileStream.is_open() || inpFileStream.fail() || inpFileStream.eof()) return -1;

	return 0;
}

int getNonCommentLine(ifstream& inpFileStream, TString& outString, char delim)
{
	outString = "";
	if(!inpFileStream.is_open() || inpFileStream.fail() || inpFileStream.eof()) return -1;

	TString sBuffer("#");

	while (sBuffer[0] == '#' || sBuffer[0] == '%')
	{
		sBuffer.ReadToDelim(inpFileStream, delim);
		if(inpFileStream.fail()) return -1;
	}

	outString = sBuffer;
	return 0;
}

int doubleFromList(TString theList, int pos, double& result)
{
	istringstream theStringStream(theList.Data());

	TString tempString;
	result = 0;

	for (int i = 0; i <= pos; i++)
	{
		theStringStream >> tempString;
		if(theStringStream.fail()) return -1;
	}
	istringstream resultStream(tempString.Data());
	resultStream >> result;
	return 0;
}

int intFromList(TString theList, int pos, int& result)
{
	double tempDouble;
	int errorCode;
	errorCode = doubleFromList(theList, pos, tempDouble);
	result = (int) tempDouble;
	return errorCode;
}

int boolFromList(TString theList, int pos, bool& result)
{
	double tempDouble;
	int errorCode;
	errorCode = doubleFromList(theList, pos, tempDouble);
	result = (bool) tempDouble;
	return errorCode;
}

int stringFromList(TString theList, int pos, TString& result)
{
	istringstream theStringStream(theList.Data());
	TString temp;

	for (int i = 0; i <= pos; i++)
	{
		if(theStringStream.fail() || theList == "") return -1;
		theStringStream >> temp;

	}
	result = temp;
	return 0;
}

TString stringFromList(TString theList, int pos)
{
	istringstream theStringStream(theList.Data());
	TString result;

	for (int i = 0; i <= pos; i++)
	{
		if(theStringStream.fail() || theList == "") return -1;
		theStringStream >> result;

	}
	return result;
}

int numItemsInStringList(TString theList)
{
	int numItems = 0;
	TString temp;

	while (!stringFromList(theList, numItems, temp)) //As long as it can find an item
	{
		numItems++;
	}
	return numItems;
}

bool FileExists(TString strFileName)
{
	struct stat stFileInfo;
	bool blnReturn;
	int intStat;
//#ifndef __CINT__
	// Attempt to get the file attributes
	intStat = stat(strFileName.Data(), &stFileInfo);
	if(intStat == 0)
	{
		// We were able to get the file attributes
		// so the file obviously exists.
		blnReturn = true;
	}
	else
	{
		// We were not able to get the file attributes.
		// This may mean that we don't have permission to
		// access the folder which contains this file. If you
		// need to do that level of checking, lookup the
		// return values of stat which will give you
		// more details on why stat failed.
		blnReturn = false;
	}

	return (blnReturn);
}
//#else //__CINT__
//    Long_t *id,*size,*flags,*mt;
//    return !(gSystem->GetPathInfo(strFileName,id,size,flags,mt));
//}
//#endif //__CINT__

TString addBeforeExtension(TString baseString, TString addString)
{
	TString outString = baseString;

	Ssiz_t found = outString.Last('.');
	if(found == kNPOS)
	{
		return outString;
	}
	outString.Insert(found, addString);
	return outString;
}

TString int2str(int i)
{
	stringstream out;
	out << i;
	return out.str();
}

int getMaxNumberByDigits(int digits)
{
	int maxNum = 1;
	for (int k = 0; k < digits; k++)
		maxNum *= 10;
	maxNum -= 1;
	return maxNum;
}

TString int2str(int i, int digits, bool strict)
{
	stringstream out;
	int maxNum = getMaxNumberByDigits(digits);

	if(strict && (i > maxNum || i < 0)) return out.str();

	for (int j = digits; j >= 1; j--)
	{
		if(i < getMaxNumberByDigits(j - 1) + 1)
			out << "0";
		else
		{
			out << i;
			break;
		}
	}
	return out.str();
}

TString d2str(double d, int inpPrecision, int opt)
{
	stringstream out;
	if(inpPrecision > 0)
	{
		out.precision(inpPrecision);
		out.unsetf(ios::floatfield);
	}

	if(opt == 0)
	{

	}
	else if(opt == 1)
	{
		out << std::fixed;
	}
	else if(opt == 2)
	{
		out << std::scientific;
	}
	else
	{
		cout << "Error!" << endl;
	}
	out << d;
	return out.str();
}

double str2double(TString s)
{
	return atof(s);
}

int str2int(TString s)
{
	return atoi(s);
}
unsigned int str2uint(TString s)
{
	return atoi(s);
}

int convertITXtoTH1I(TString itxFileName, TFile* theRootFile)
{
	theRootFile->cd();

	int numBins;
	double binStart, binEnd, binSize;
	TString s, waveName;
	stringstream ss;
	int data;
	TString endString = "\tEND";

	ifstream itxFile(itxFileName);

	if(!itxFile.is_open() || itxFile.fail() || itxFile.eof()) return -1;

	s.ReadToDelim(itxFile, '\r'); //IGOR
	s.ReadToDelim(itxFile, '\r'); //WAVE <WAVENAME>
	ss << s;
	ss >> s >> waveName;
	ss.clear();
	s.ReadToDelim(itxFile, '\r'); //Begin
	//First going to find the bin information
	for (numBins = 0; s.ReadToDelim(itxFile, '\r') && s != endString && numBins < 1000000 && !itxFile.eof(); numBins++)
	{

	}
	numBins--;

	sscanf(s, "%*s %*s %*s %lf,%lf", &binStart, &binSize);

	binEnd = binStart + binSize * numBins;

	cout << waveName << " in " << itxFileName << " will be converted to histogram in " << theRootFile->GetPath() << endl;
	cout << "Num Bins: " << numBins << "  Bin Start: " << binStart << "  Bin End: " << binEnd << endl;
	itxFile.close();

	TH1I theHist(waveName, waveName, numBins, binStart, binEnd);
	theHist.Sumw2();
	itxFile.open(itxFileName);

	s.ReadToDelim(itxFile, '\r'); //IGOR
	s.ReadToDelim(itxFile, '\r'); //WAVE <WAVENAME>
	s.ReadToDelim(itxFile, '\r'); //Begin
	for (int i = 0; i < numBins; i++)
	{
		s.ReadToDelim(itxFile, '\r'); //Data line
		sscanf(s, "%d", &data);
		for (int j = 0; j < data; j++)
			theHist.Fill((i + .5) * binSize);
	}

	theHist.Write();

	itxFile.close();

	return 0;
}

void typeAnythingToContinue(TString message)
{
	TString tempString;
	cout << message << " Type Anything to continue: ";
	cin >> tempString;
	return;
}

TString addExeDir(TString inpString)
{
	TString returnString = "/media/mjbexternal/exe/" + inpString;
	return returnString;
}

TString fileNameFromFullPath(TString fullPath)  //return fullPath if not found
{
	Ssiz_t found;
	found = fullPath.Last('/');
	if(found == kNPOS) return fullPath;
	return fullPath.Remove(0, found + 1);
}

TString filePathFromFullPath(TString fullPath)
{
	Ssiz_t found;
	found = fullPath.Last('/');
	if(found == kNPOS) return "";
	return fullPath.Remove(found + 1, fullPath.Length() - found);
}

TString fileExtensionFromFullPath(TString fullPath)
{
	Ssiz_t found;
	found = fullPath.Last('.');
	if(found == kNPOS) return "";
	return fullPath.Remove(0, found + 1);
}

int countTString(TString inpString, char inpChar)
{
	int numCharFound = 0;
	for (int i = 0; i < inpString.Sizeof(); i++)
	{
		if(inpString[i] == inpChar)
		{
			numCharFound++;
		}
	}
	return numCharFound;
}
