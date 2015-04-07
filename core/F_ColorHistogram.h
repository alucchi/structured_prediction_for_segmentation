/////////////////////////////////////////////////////////////////////////
// This program is free software; you can redistribute it and/or       //
// modify it under the terms of the GNU General Public License         //
// version 2 as published by the Free Software Foundation.             //
//                                                                     //
// This program is distributed in the hope that it will be useful, but //
// WITHOUT ANY WARRANTY; without even the implied warranty of          //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   //
// General Public License for more details.                            //
//                                                                     //
// Written and (C) by Aurelien Lucchi                                  //
// Contact <aurelien.lucchi@gmail.com> for comments & bug reports      //
/////////////////////////////////////////////////////////////////////////

#ifndef F_COLHISTO_H
#define F_COLHISTO_H

#include "F_Histogram.h"
#include "Histogram.h"
#include "Slice.h"
#include "Slice3d.h"
#include "Feature.h"

//-------------------------------------------------------------------------TYPES

//-------------------------------------------------------------------------CLASS

class F_ColorHistogram : public Feature
{
 public:	

  /**
   * @param img an image can be passed to set the correct number of channels.
   * If 0, look at the config file.
   */
  F_ColorHistogram(int _nb_bins,
                   int _max_pixel_value,
                   eHistogramType _histoType,
                   bool _useColorImage=true,
                   IplImage* img = 0);

  int getSizeFeatureVectorForOneSupernode();

  /**
   * Extract a feature vector for a given supernode in a 2d slice
   */
  bool getFeatureVectorForOneSupernode(osvm_node *x,
                                       Slice* slice,
                                       const int supernodeId);

  /**
   * Extract a feature vector for a given supernode in a 3d volume
   */
  bool getFeatureVectorForOneSupernode(osvm_node *x,
                                       Slice3d* slice3d,
                                       const int supernodeId);

  int nChannels;
  int nBinsPerSupernode;
  int nBins;
  int nTotalBins; //total number of bins for supernode and its neighbors
  int max_pixel_value;
  eHistogramType histoType;
  bool useColorImage;

 private:
  int offsetNeighbors;
  int nBinsPerSupernode2; // nBinsPerSupernode squared
  float valToBin;
  int nLocations;
};

#endif // F_COLHISTO_H
