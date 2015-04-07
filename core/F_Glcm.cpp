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

#include "F_Glcm.h"
#include "Config.h"

#define GLCM_DEFAULT_N_LEVELS 8

F_Glcm::~F_Glcm()
{
}

F_Glcm::F_Glcm()
{
  maxIntensity = 255;
  nItensityLevels = GLCM_DEFAULT_N_LEVELS;
  valToIdx = nItensityLevels/(double)maxIntensity;
}

int F_Glcm::getSizeFeatureVectorForOneSupernode()
{
  return nItensityLevels*nItensityLevels;
}

bool F_Glcm::getFeatureVectorForOneSupernode(osvm_node *x, Slice_P* slice, int supernodeId)
{
  supernode* s = slice->getSupernode(supernodeId);
  double value, n_value;
  int idx1, idx2, idx;
  node n;

  int w = slice->getWidth();
  int h = slice->getHeight();
  int d = slice->getDepth();

  // allocate memory to store data (done inside this function for multi-thread code)
  int sizeFV = getSizeFeatureVectorForOneSupernode();
  double* glcm_data = new double[sizeFV];
  memset(glcm_data, 0, sizeof(double)*sizeFV);

  nodeIterator ni = s->getIterator();
  ni.goToBegin();

  while(!ni.isAtEnd()) {
    ni.get(n);
    ni.next();

    value = slice->getIntensity(n.x, n.y, n.z);
    //idx1 = floor(value * valToIdx);
    idx1 = (int)(value * valToIdx);
    //assert(idx1 <= nItensityLevels);
    if(idx1 == nItensityLevels) {
      idx1 = nItensityLevels - 1;
    }

    // check if neighbors are part of the supernode.
    int min_x = max(0, n.x - 1);
    int max_x = min(w - 1, n.x + 1);
    int min_y = max(0, n.y - 1);
    int max_y = min(h - 1, n.y + 1);
    int min_z = max(0, n.z - 1);
    int max_z = min(d - 1, n.z + 1);
    for(int _x = min_x; _x <= max_x; ++_x) {
      for(int _y = min_y; _y <= max_y; ++_y) {
        for(int _z = min_z; _z <= max_z; ++_z) {
          sidType nsid = slice->getSid(_x, _y, _z);
	  if(_x != n.x || _y != n.y || _z != n.z) {
	    // make sure both pixels belong to the same supernode
	    if(supernodeId == nsid) {
	      n_value = slice->getIntensity(_x, _y, _z);
	      //idx2 = floor(n_value * valToIdx);
	      idx2 = (int)(n_value * valToIdx);
	      //assert(idx2 <= nItensityLevels);
	      if(idx2 == nItensityLevels) {
		idx2 = nItensityLevels - 1;
	      }
	      idx = (idx1*nItensityLevels + idx2);
	      ++glcm_data[idx];
	    }
	  }
	}
      }
    }
  }

  // copy feature vector
  for(int i = 0; i < sizeFV; ++i) {
    x[i].value = glcm_data[i];
  }

  delete[] glcm_data;

  return true;
}

bool F_Glcm::getFeatureVectorForOneSupernode(osvm_node *x, Slice* slice, int supernodeId)
{
  Slice_P* slice_p = static_cast<Slice_P*>(slice);
  return getFeatureVectorForOneSupernode(x, slice_p, supernodeId);
}

bool F_Glcm::getFeatureVectorForOneSupernode(osvm_node *x, Slice3d* slice3d, int supernodeId)
{
  Slice_P* slice_p = static_cast<Slice_P*>(slice3d);
  return getFeatureVectorForOneSupernode(x, slice_p, supernodeId);
}
