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

#include "F_Combo.h"
#include "oSVM.h"

//--------------------------------------------------------------------- METHODS

F_Combo::F_Combo(vector<eFeatureType>& feature_types,
                 Slice* slice)
{
  init();
  Feature* _feature = 0;
  for(vector<eFeatureType>::iterator it = feature_types.begin();
      it != feature_types.end(); it++) {
    _feature = Feature::getFeature(slice, *it);
    features.push_back(_feature);
  }

  // Initialize feature vectors
  printf("[F_Combo] Allocating memory to store %ld different features\n", features.size());
  feature_buffer = new osvm_node*[features.size()];

  uint fidx = 0;
  for(vector<Feature*>::iterator iFeature = features.begin();
      iFeature != features.end(); iFeature++) {
    int max_index = (*iFeature)->getSizeFeatureVectorForOneSupernode()+1;
    feature_buffer[fidx] = new osvm_node[max_index];

    int i = 0;
    for(i = 0;i < max_index-1; i++)
      feature_buffer[fidx][i].index = i+1;
    feature_buffer[fidx][i].index = -1;

    ++fidx;
  }

  sizeFV = 0;
  for(vector<Feature*>::iterator iFeature = features.begin();
      iFeature != features.end(); iFeature++) {
    sizeFV += (*iFeature)->getSizeFeatureVectorForOneSupernode();
  }
}

F_Combo::F_Combo(vector<eFeatureType>& feature_types,
                 Slice3d* slice)
{
  init();
  Feature* _feature = 0;
  for(vector<eFeatureType>::iterator it = feature_types.begin();
      it != feature_types.end(); it++) {
    _feature = Feature::getFeature(slice, *it);
    features.push_back(_feature);
  }

  // Initialize feature vectors
  printf("[F_Combo] Allocating memory to store %ld different features\n", features.size());
  feature_buffer = new osvm_node*[features.size()];

  uint fidx = 0;
  for(vector<Feature*>::iterator iFeature = features.begin();
      iFeature != features.end(); iFeature++) {
    int max_index = (*iFeature)->getSizeFeatureVectorForOneSupernode()+1;
    feature_buffer[fidx] = new osvm_node[max_index];

    int i = 0;
    for(i = 0;i < max_index-1; i++)
      feature_buffer[fidx][i].index = i+1;
    feature_buffer[fidx][i].index = -1;

    ++fidx;
  }

  sizeFV = 0;
  for(vector<Feature*>::iterator iFeature = features.begin();
      iFeature != features.end(); iFeature++) {
    sizeFV += (*iFeature)->getSizeFeatureVectorForOneSupernode();
  }
}


void F_Combo::init()
{
  normalize_features = false;
  string config_tmp;
  if(Config::Instance()->getParameter("normalize_combo_features", config_tmp)) {
    normalize_features = atoi(config_tmp.c_str());
    printf("[F_Combo] normalize_features = %d\n", (int)normalize_features);
  }
}

F_Combo::~F_Combo()
{
  uint fidx = 0;
  for(vector<Feature*>::iterator iFeature = features.begin();
      iFeature != features.end(); iFeature++) {
    delete *iFeature;
    delete[] feature_buffer[fidx];
    ++fidx;
  }
}

int F_Combo::getSizeFeatureVectorForOneSupernode()
{
  return sizeFV;
}

bool F_Combo::getFeatureVectorForOneSupernode(osvm_node *x, Slice* slice, int supernodeId)
{
  uint idx = 0;
  uint fidx = 0;
  for(vector<Feature*>::iterator iFeature = features.begin();
      iFeature != features.end(); iFeature++) {

    osvm_node *sx = feature_buffer[fidx];
    ++fidx;
    (*iFeature)->getFeatureVectorForOneSupernode(sx, slice, supernodeId);

    if(normalize_features > 0) {
      double norm = oSVM::norm(sx, normalize_features);
      if(fabs(norm) < 1e-20) {
        norm = 1.0;
      }
      for(int i = 0;sx[i].index != -1; i++) {
        x[idx].value = sx[i].value/norm;
        idx++;
      }
    } else {
      for(int i = 0;sx[i].index != -1; i++) {
        x[idx].value = sx[i].value;
        idx++;
      }
    }
  }

  return true;
}

bool F_Combo::getFeatureVectorForOneSupernode(osvm_node *n,
                                              const int x,
                                              const int y)
{
  uint idx = 0;
  uint fidx = 0;
  for(vector<Feature*>::iterator iFeature = features.begin();
      iFeature != features.end(); iFeature++) {
    osvm_node *sx = feature_buffer[fidx];
    ++fidx;
    (*iFeature)->getFeatureVectorForOneSupernode(sx,x,y);

    if(normalize_features > 0) {
      double norm = oSVM::norm(sx, normalize_features);
      for(int i = 0;sx[i].index != -1; i++) {
        n[idx].value = sx[i].value/norm;
        idx++;
      }
    } else {
      for(int i = 0;sx[i].index != -1; i++) {
        n[idx].value = sx[i].value;
        idx++;
      }
    }
  }

  return true;
}

bool F_Combo::getFeatureVectorForOneSupernode(osvm_node *x, Slice3d* slice3d, int supernodeId)
{
  uint idx = 0;
  uint fidx = 0;
  for(vector<Feature*>::iterator iFeature = features.begin();
      iFeature != features.end(); iFeature++) {
    osvm_node *sx = feature_buffer[fidx];
    ++fidx;
    (*iFeature)->getFeatureVectorForOneSupernode(sx, slice3d, supernodeId);

    if(normalize_features > 0) {
      double norm = oSVM::norm(sx, normalize_features);
      for(int i = 0;sx[i].index != -1; i++) {
        x[idx].value = sx[i].value/norm;
        idx++;
      }
    } else {
      for(int i = 0;sx[i].index != -1; i++) {
        x[idx].value = sx[i].value;
        idx++;
      }
    }
  }

  return true;
}

bool F_Combo::getFeatureVectorForOneSupernode(osvm_node *n, Slice3d* slice3d,
                                              const int x, const int y, const int z)
{
  uint idx = 0;
  uint fidx = 0;
  for(vector<Feature*>::iterator iFeature = features.begin();
      iFeature != features.end(); iFeature++) {
    osvm_node *sx = feature_buffer[fidx];
    ++fidx;
    (*iFeature)->getFeatureVectorForOneSupernode(sx,slice3d,x,y,z);

    if(normalize_features > 0) {
      double norm = oSVM::norm(sx, normalize_features);
      for(int i = 0;sx[i].index != -1; i++) {
        n[idx].value = sx[i].value/norm;
        idx++;
      }
    } else {
      for(int i = 0;sx[i].index != -1; i++) {
        n[idx].value = sx[i].value;
        idx++;
      }
    }
  }

  return true;
}
