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

#include "F_Lbp.h"
#include "Config.h"

#define LBP_N_BINS 256

F_Lbp::~F_Lbp()
{
}

F_Lbp::F_Lbp()
{
}

int F_Lbp::getSizeFeatureVector()
{
  return LBP_N_BINS;
}

bool F_Lbp::getFeatureVector(osvm_node *x, Slice* slice, int supernodeId)
{
  supernode* s = slice->getSupernode(supernodeId);
  double value, n_value;
  node n;

  const int w = slice->getWidth();
  const int h = slice->getHeight();
  const int z = 0;

  int sizeFV = getSizeFeatureVector();
  const int neighborhoodSize = 8;
  int pattern;

  // initialize feature vector
  for(int i = 0; i < sizeFV; ++i) {
    x[i].value = 0;
  }

  nodeIterator ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd()) {
    ni.get(n);
    ni.next();

    value = slice->getIntensity(n.x, n.y, n.z);
    pattern = 0;

    // check if neighbors are part of the supernode.
    int min_x = max(0, n.x - 1);
    int max_x = min(w, n.x + 1);
    int min_y = max(0, n.y - 1);
    int max_y = min(h, n.y + 1);
    int count = 0;
    for(int _x = min_x; _x <= max_x; ++_x) {
      for(int _y = min_y; _y <= max_y; ++_y) {
        sidType nsid = slice->getSid(_x, _y, _z);
        if(supernodeId == nsid) {
          n_value = slice->getIntensity(_x, _y, _z);
          pattern |= (value > n_value) << count;
          ++count;
        }
      }
    }
    if(count == neighborhoodSize) {
      // valid
      ++x[pattern].value;
    }
  }

  return true;
}

bool F_Lbp::getFeatureVector(osvm_node *x, Slice3d* slice3d, int supernodeId)
{
  assert(0);
}
