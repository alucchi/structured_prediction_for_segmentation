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

#ifndef F_HISTO_H
#define F_HISTO_H

#include "Histogram.h"
#include "Slice_P.h"
#include "Slice.h"
#include "Slice3d.h"
#include "Feature.h"

//-------------------------------------------------------------------------TYPES

enum eHistogramType
{
  NO_NEIGHBORS = 0,
  INCLUDE_NEIGHBORS,
  INCLUDE_NEIGHBORS_IN_SEPARATE_BINS,
  INCLUDE_NEIGHBORS_PLUS_LOCATION
};

//-------------------------------------------------------------------------CLASS

class F_Histogram : public Feature
{
 public:	

  /**
   * @param img an image can be passed to set the correct number of channels.
   * If 0, look at the config file.
   */
  F_Histogram(int _nb_bins,
              int _max_pixel_value,
              eHistogramType _histoType,
              bool _useColorImage = false,
              IplImage* img = 0);

 protected:
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

 public:
  int nChannels;
  int nBinsPerSupernode;
  int nBins; //total number of bins for supernode and its neighbors
  int max_pixel_value;
  eHistogramType histoType;
  bool useColorImage;
  bool normalize_l1_norm;
  bool normalize_l2_norm;
};

#endif // F_HISTO_H
