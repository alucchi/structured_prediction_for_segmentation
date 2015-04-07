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

#ifndef F_SIFT_H
#define F_SIFT_H

#include "cv.h"
#include "sift.hpp"

// SliceMe
#include "Feature.h"
#include "oSVM.h"
#include "Slice.h"

class F_Sift : public Feature
{
 public:	
  VL::Sift* sift;
  float sigma0;

  /**
   * @param aoctaves is the number of octaves in the pyramid
   * @param levels : default value should be 3
   * @param aomin is the first octave. Note that setting it to -1, -2, ...
   * will cause the base of the pyramid to be
   * two, three, ... times larger than the input image.
   */
  F_Sift(const char* image_name,
         int aoctaves,
         int alevels,
         int aomin);
  ~F_Sift();

  void addScale(float _scale);

  bool getFeatureVectorForOneSupernode(osvm_node *n, Slice* slice, int supernodeId);

  bool getFeatureVectorForOneSupernode(osvm_node *n,
                                       Slice* slice,
                                       const int x,
                                       const int y);

  int getSizeFeatureVectorForOneSupernode();

  void setScales(vector<float>& _scales);

 private:
  VL::pixel_t* im_pt; // image data
  VL::float_t* descr_pt; //descriptor data

  int octaves;
  int levels;
  int omin;
  vector<float> scales;

  static const int desc_size = 128;

};

#endif // F_SIFT_H
