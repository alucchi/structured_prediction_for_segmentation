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

#ifndef F_Filter_H
#define F_Filter_H

#include "Feature.h"
#include "oSVM.h"
#include "Slice.h"
#include "Slice3d.h"
#include "Slice_P.h"

//-------------------------------------------------------------------------TYPES

//-------------------------------------------------------------------------CLASS

class F_Filter : public Feature
{
 public:	

  F_Filter(Slice_P& slice);

  ~F_Filter();

  void createSupernodeBasedFeatures(Slice_P& slice, uchar* node_features, int featIdx);

  int getSizeFeatureVectorForOneSupernode();

  bool getFeatureVector(osvm_node *n,
                        const int x,
                        const int y);

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

  bool getFeatureVector(osvm_node *x,
                        Slice3d* slice3d,
                        const int gx,
                        const int gy,
                        const int gz);

  void precomputeFeatures(Slice_P& slice);

private:
  uchar** features;
  int sizeFV;
};

#endif // F_Filter_H
