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

#include "F_Precomputed.h"
#include "utils.h"
#include "Config.h"

using namespace std;

//-----------------------------------------------------------------------------

//--------------------------------------------------------------------- METHODS

F_Precomputed::F_Precomputed(map<sidType, osvm_node*>* _features, int _feature_size)
{
   features = _features;
   feature_size = _feature_size;
}

F_Precomputed::~F_Precomputed()
{
}

int F_Precomputed::getSizeFeatureVectorForOneSupernode()
{
  return feature_size;
}

bool F_Precomputed::getFeatureVectorForOneSupernode(osvm_node *x, Slice* slice, int supernodeId)
{
  // not implemented yet
  return false;
}

bool F_Precomputed::getFeatureVectorForOneSupernode(osvm_node *x,
                                                    Slice3d* slice3d,
                                                    const int supernodeId)
{
  for(int i = 0; i < feature_size; i++) {
    x[i].value = (*features)[supernodeId][i].value;
  }
  return true;
}
