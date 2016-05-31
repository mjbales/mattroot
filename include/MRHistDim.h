#ifndef MATTROOT_INCLUDE_MRHISTDIM_H_
#define MATTROOT_INCLUDE_MRHISTDIM_H_

//////////////////////////////////////////////////
//MRHistDim.h - by Matthew Bales
//Simple histogram dimensions
//
///////////////////////////////////////////////////

struct HistDim
{
	int numBins;
	double binLowEdge;
	double binHighEdge;

	inline double getBinWidth() const {return (binHighEdge-binLowEdge)/(double) numBins;}
};

#endif /* MATTROOT_INCLUDE_MRHISTDIM_H_ */
