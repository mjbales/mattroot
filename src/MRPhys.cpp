#include "MRPhys.h"

using namespace std;

TVector3 convertMomentumToVelocity(TVector3 inpP, double massEnergy)
{
	return inpP * (TMath::C() / sqrt(inpP.X() * inpP.X() + inpP.Y() * inpP.Y() + inpP.Z() * inpP.Z() + massEnergy * massEnergy));
}

double convertMomentumToVelocity(double inpP, double massEnergy)
{
	return TMath::C() * inpP / sqrt(inpP * inpP + pow(massEnergy, 2));
}

double convertKEToMomentum(double inpKE, double massEnergy)
{
	return sqrt(pow(inpKE + massEnergy, 2) - pow(massEnergy, 2));
}

double convertMomentumToKE(double inpMom, double massEnergy)
{
	return sqrt(inpMom * inpMom + massEnergy * massEnergy) - massEnergy;
}

double convertMomentumToKE(TVector3 inpMom, double massEnergy)
{
	return sqrt(inpMom.X() * inpMom.X() + inpMom.Y() * inpMom.Y() + inpMom.Z() * inpMom.Z() + massEnergy * massEnergy) - massEnergy;

}

TVector3 convertVelocityToMomentum(TVector3 vIn, double massEnergy)
{
	vIn.SetMag(relGamma(vIn) * massEnergy / TMath::C()); //Gamma*mass*velocity
	return vIn;
}

double convertVelocityToMomentum(double vIn, double massEnergy)
{
	return relGamma(vIn) * massEnergy * vIn / TMath::C();
}

double convertVelocityToKE(TVector3 vIn, double massEnergy)
{

	return convertMomentumToKE(convertVelocityToMomentum(vIn, massEnergy), massEnergy);
}

double convertVelocityToKE(double vx, double vy, double vz, double massEnergy)
{
	TVector3 vIn(vx, vy, vz);
	return convertVelocityToKE(vIn, massEnergy);
}

double relGamma(double vMagIn)
{
	return 1. / sqrt(1. - vMagIn * vMagIn / (TMath::C() * TMath::C()));
}

double relGamma(TVector3 vIn)
{
	return 1. / sqrt(1. - (vIn.X() * vIn.X() + vIn.Y() * vIn.Y() + vIn.Z() * vIn.Z()) / (TMath::C() * TMath::C()));
}

