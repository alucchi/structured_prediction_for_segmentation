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

#ifndef F_COMBO_H
#define F_COMBO_H

// SliceMe
#include "Feature.h"
#include "globalsE.h"
#include "Slice.h"
#include "Slice3d.h"

#include <vector>

using namespace std;

/**
 * Combo class combines different features
 */
class F_Combo : public Feature
{
 public:	

  /*
  // AL28
  F_Combo(vector<eFeatureType>& feature_types,
          const char* image_name);
  */

  F_Combo(vector<eFeatureType>& types,
          Slice* slice);

  F_Combo(vector<eFeatureType>& types,
          Slice3d* slice3d);

  ~F_Combo();

  inline int getSizeFeatureVectorForOneSupernode();

  const vector<Feature*>& getFeatures() { return features; }

  bool getFeatureVectorForOneSupernode(osvm_node *n,
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

  bool getFeatureVectorForOneSupernode(osvm_node *n,
                                       Slice3d* slice3d,
                                       const int x,
                                       const int y,
                                       const int z);

  void init();

 private:
  vector<Feature*> features;
  int normalize_features;
  int sizeFV;

  // use to store features
  osvm_node** feature_buffer;
};

#endif // F_COMBO_H
