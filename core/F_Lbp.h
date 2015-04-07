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

#ifndef F_LBP_H
#define F_LBP_H

#include "Slice.h"
#include "Slice3d.h"
#include "Feature.h"
#include "oSVM.h"

//-------------------------------------------------------------------------TYPES

//-------------------------------------------------------------------------CLASS

class F_Lbp : public Feature
{
 public:	

  F_Lbp();
  ~F_Lbp();

  int getSizeFeatureVector();

  bool getFeatureVector(osvm_node *x,
                        Slice_P* slice,
                        int supernodeId);

  /**
   * Extract a feature vector for a given supernode in a 2d slice
   */
  bool getFeatureVector(osvm_node *x,
                        Slice* slice,
                        const int supernodeId);

  /**
   * Extract a feature vector for a given supernode in a 3d volume
   */
  bool getFeatureVector(osvm_node *x,
                        Slice3d* slice3d,
                        const int supernodeId);

 private:
  int nItensityLevels;
  int maxIntensity;
  double* glcm_data;
};

#endif // F_LBP_H
