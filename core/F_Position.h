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

#ifndef F_Position_H
#define F_Position_H

#include "Feature.h"
#include "Slice.h"
#include "Slice3d.h"

//-------------------------------------------------------------------------TYPES

//-------------------------------------------------------------------------CLASS

class F_Position : public Feature
{
 public:	

  F_Position(const char* image_name=0);

  ~F_Position();

  int getSizeFeatureVector();

  bool getFeatureVector(osvm_node *n,
                        const int x,
                        const int y);
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

  bool getFeatureVector(osvm_node *x,
                        Slice3d* slice3d,
                        const int gx,
                        const int gy,
                        const int gz);

private:
  IplImage* img;

};

#endif // F_Position_H
