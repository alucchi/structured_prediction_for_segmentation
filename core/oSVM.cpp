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

// standard libraries
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <float.h>
#include <math.h>
#include <fstream>
#include <algorithm>

// SliceMe
#include "Config.h"
#include "oSVM.h"
#include "utils.h"
#include "globalsE.h"

using namespace std;

//------------------------------------------------------------------------------

void oSVM::initSVMNode(osvm_node*& x, int d)
{
  x = new osvm_node[d+1];
  int i = 0;
  for(i = 0; i < d; i++)
    x[i].index = i+1;
  x[i].index = -1;
}

void oSVM::initSVMNodes(osvm_node**& xs, uint nNodes, int d)
{
  xs = new osvm_node*[nNodes];
  for(uint k=0; k < nNodes; k++)
    {
      xs[k] = new osvm_node[d+1];
      int i = 0;
      for(i = 0;i < d; i++)
        xs[k][i].index = i+1;
      xs[k][i].index = -1;
    }
}

void oSVM::print(osvm_node *x, const char* title)
{
  if(title)
    printf("%s: \n", title);
  for(int i = 0;x[i].index != -1; i++) {
    printf("%d:%g ",x[i].index,x[i].value);
  }
  printf("\n");
}

void oSVM::printNonZeros(osvm_node *x, const char* title)
{
  if(title)
    printf("%s: \n", title);
  for(int i = 0;x[i].index != -1; i++) {
    if (x[i].value != 0) {
      printf("%d:%g ",x[i].index,x[i].value);
    }
  }
  printf("\n");
}

double oSVM::norm(osvm_node* x, int norm_type)
{
  double _norm = 0;
  if(norm_type == L1_NORM) {
    for(int i = 0;x[i].index != -1; i++) {
      _norm += x[i].value;
    }
  } else {
    for(int i = 0;x[i].index != -1; i++) {
      _norm += x[i].value*x[i].value;
    }
    _norm = sqrt(_norm);
  }
  return _norm;
}

void oSVM::compute_range(Slice_P* slice, Feature* feature, map<int, double>& feature_min, map<int, double>& feature_max)
{
  int fvSize = feature->getSizeFeatureVector();
  int max_index = fvSize + 1;

  // first index=1
  int offset_index = 1;

  feature_min.clear();
  feature_max.clear();
  for(int idx = 0;idx < fvSize; idx++) {
    feature_max[offset_index+idx]=-DBL_MAX;
    feature_min[offset_index+idx]=DBL_MAX;
  }

  osvm_node* x = new osvm_node[max_index];
  int i = 0;
  for(i = 0;i < max_index-1; i++)
    x[i].index = i+1;
  x[i].index = -1;

  // copy all the features in an array of std::vectors
  vector<float>* features = new vector<float>[fvSize];
  for(int sid = 0; sid < (int)slice->getNbSupernodes(); sid++) {
    feature->getFeatureVector(x, slice, sid);
    for(int i=0;i<fvSize;i++) {
      features[i].push_back(x[i].value);
    }
  }

  ulong idx_min = 0.01*slice->getNbSupernodes();
  ulong idx_max = 0.99*slice->getNbSupernodes();
  for(int i = 0;i < fvSize;i++) {
    nth_element(features[i].begin(), features[i].begin()+idx_min, features[i].end());
    feature_min[offset_index+i] = features[i][idx_min];

    nth_element(features[i].begin(), features[i].begin()+idx_max, features[i].end());
    feature_max[offset_index+i] = features[i][idx_max];
  }

  delete[] x;
  delete[] features;
}

void oSVM::rescale(osvm_node *x, map<int, double>& feature_min, map<int, double>& feature_max)
{
  const double lower_value = 0; // lower value used when rescaling feature vectors
  const double upper_value = 1; // upper value used when rescaling feature vectors

  int n = 0;
  int index = x[n].index;
  while(index != -1)
    {
      if((feature_max[index] - feature_min[index])>0.0001)
        {
          if(x[n].value == feature_min[index])
            x[n].value = lower_value;
          else if(x[n].value == feature_max[index])
            x[n].value = upper_value;
          else
            x[n].value = lower_value + (upper_value-lower_value) * 
              (x[n].value-feature_min[index])/
              (feature_max[index]-feature_min[index]);
        }
      n++;
      index = x[n].index;
    }
}
