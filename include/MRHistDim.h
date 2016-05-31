/*
 * MRHistDim.h
 *
 *  Created on: Oct 7, 2014
 *      Author: mjbales
 */

#ifndef MATTROOT_INCLUDE_MRHISTDIM_H_
#define MATTROOT_INCLUDE_MRHISTDIM_H_

struct HistDim
{
	int numBins;
	double binLowEdge;
	double binHighEdge;

	inline double getBinWidth(){return (binHighEdge-binLowEdge)/(double) numBins;}
};

#endif /* MATTROOT_INCLUDE_MRHISTDIM_H_ */
