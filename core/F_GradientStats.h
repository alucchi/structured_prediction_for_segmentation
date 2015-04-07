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

#ifndef F_GRADIENTSTATS_H
#define F_GRADIENTSTATS_H

#include "Feature.h"

#include <string>
#include <vector>

enum eGradientType
{
  GRADIENT_X = 0,
  GRADIENT_Y,
  GRADIENT_Z,
  GRADIENT_MAG
};

//------------------------------------------------------------------------------

// Use 7 scales
static const int numScales = 7;
static const double scales[numScales] = {0.3, 0.7, 1.0, 1.6, 3.5, 5.0, 10.0};  // see ilastik

//------------------------------------------------------------------------------

class F_GradientStats : public Feature
{
 public:	

  F_GradientStats(Slice3d& slice3d,
                  eGradientType gradientType);

  ~F_GradientStats();

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

 private:

  void computeGradientCube(Slice3d& slice3d,
                           eGradientType gradientType);

  float* minGradientValue;
  float* maxGradientValue;

  int nBinsPerScale;

  /**
   * Cube containing the gradient
   */
  float** fGradientCube;
};

#endif // F_GRADIENTSTATS_H
