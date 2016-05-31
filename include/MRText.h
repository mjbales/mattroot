#ifndef CMRKTEXT_H_INCLUDED
#define CMRKTEXT_H_INCLUDED

//////////////////////////////////////////////////
//MRText.h - by Matthew Bales
//misc string/text related functions
//
///////////////////////////////////////////////////

#include <iostream>
#include <fstream>

#include "TString.h"
#include "TFile.h"

int getNonCommentLine(std::ifstream& inpFileStream, TString& outString, char delim = '\n');  /// given and input file stream, extracts lines until a non comment line is reached.  A comment line is designed by starting with a #
int doubleFromList(TString theList, int pos, double& result);         /// Given a sting based list of doubles seperated by spaces/tabs, it will extact the pos'th one
int intFromList(TString theList, int pos, int& result);               /// Same only for ints
int boolFromList(TString theList, int pos, bool& result);             /// Same only for bools
int stringFromList(TString theList, int pos, TString& result);            /// Same only for strings
TString stringFromList(TString theList, int pos);
int numItemsInStringList(TString theList);                              /// using stringFromList, determine how many strings are in list
int skipCommentLines(std::ifstream&);                                    /// Given and input stream it skips lines until the next line will be a noncomment line
bool FileExists(TString strFilename);                                /// Checks to see whether file exists.  Has IFDEF for working with CINT
TString addBeforeExtension(TString baseString, TString addString);     /// Adds addTString before the extension in baseTString (for filenames)

void typeAnythingToContinue(TString message);                        /// Simple cin input to type anything and press enter before continuing...useful for debugging

TString int2str(int i);              /// Converts an int to a str (SHOULD USE std::to_string)
TString int2str(int i, int digits, bool strict = true);   /// Converts an int to a str ensuring it has the appropriate digits (for simple alphabatetical order)
TString d2str(double d, int inpPrecision = -1, int opt = 0);        /// Converts a double to a string (SHOULD USE std::to_string for simple situations)

double str2double(TString s);        /// Converts a TString to a double (SHOULD USE std::stod)
int str2int(TString s);              /// Converts a TString to an int (SHOULD USE std::stoi)
unsigned int str2uint(TString s);              /// Converts a TString to an int (SHOULD use std::stoul)

int convertITXtoTH1I(TString itxFileName, TFile* theRootFile);  /// Conversion of ITX to TH1I. Not very robust, must be customized conversion of an IGOR ITX file to a ROOT histogram

TString fileNameFromFullPath(TString fullPath);  /// Get file name from the full path
TString filePathFromFullPath(TString fullPath); /// Get file path (no file name) from the full path
TString fileExtensionFromFullPath(TString fullPath); /// Get the file extension from the full path
int countTString(TString inpString, char inpChar);  /// Count the number of occurances of the character in the string
inline void getline(std::istream& inpFileStream, TString inpString) //Creating getline function similar to std::string version
{
	inpString.ReadLine(inpFileStream);
}

#endif // CMRKTEXT_H_INCLUDED
