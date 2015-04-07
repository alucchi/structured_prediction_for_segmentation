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
// Written and (C) by Aurelien Lucchi and Kevin Smith                  //
// Contact <aurelien.lucchi@gmail.com> for comments & bug reports      //
// Contact <kevin.smith@epfl.ch> for comments & bug reports            //
/////////////////////////////////////////////////////////////////////////

#ifndef F_PRECOMPUTED_H
#define F_PRECOMPUTED_H

#include "Feature.h"

#include <string>
#include <vector>

//------------------------------------------------------------------------------

class F_Precomputed : public Feature
{
 public:	

  F_Precomputed(map<sidType, osvm_node*>* _features, int _feature_size);

  ~F_Precomputed();

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

  void setFeatureSize(int _feature_size) { feature_size = _feature_size; }

 private:

  int feature_size;

  // precomputed quantities for nodes
  map<sidType, osvm_node*>* features;
};

#endif // F_PRECOMPUTED_H
