//////////////////////////////////////////////////
//MRPHYS - by Matthew Bales
//misc math and physics functions that are useful in multiple applications
//
///////////////////////////////////////////////////
#ifndef MATTPHYS_H_INCLUDED
#define MATTPHYS_H_INCLUDED

//Standard Libraries
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

//Matt Libraries

//ROOT Libraries
#include "TRandom3.h"
#include "TMath.h"
#include "TCut.h"
#include "TH1.h"
#include "TVector3.h"

TVector3 convertMomentumToVelocity(TVector3 pIn, double massEnergy);  //Converts momentum in keV to m/s
double convertMomentumToVelocity(double inpP, double massEnergy);        //Converts momentum in keV to velocity in m/s
double convertKEToMomentum(double inpKE, double massEnergy);            //Converts kinetic energy in keV to momentum in keV/c
double convertMomentumToKE(double inpMom, double massEnergy);           //Converts momentum in keV to kinetic energy in keV
double convertMomentumToKE(TVector3 inpMom, double massEnergy);        //Converts momentum in keV to kinetic energy in keV
double convertVelocityToKE(TVector3 vIn, double massEnergy);           //Converts velocity in m/s to kinetic energy in keV
double convertVelocityToKE(double vx, double vy, double vz, double massEnergy);   //Converts veloccity in m/s to kinetic energy in keV
TVector3 convertVelocityToMomentum(TVector3 vIn, double massEnergy);   //Converts velocity in m/s to momentum in keV
double convertVelocityToMomentum(double vIn, double massEnergy);         //Converts velocity in m/s to momentum in keV

double relGamma(double vMagIn);     //calculates special relativity's gamma constant assuming v is in m/s
double relGamma(TVector3 vIn);     //calculates special relativity's gamma constant assuming v is in m/s

#endif // MATTPHYS_H_INCLUDED
