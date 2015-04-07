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

#ifndef F_ORIENTEDHISTOGRAM_H
#define F_ORIENTEDHISTOGRAM_H

#include "Slice.h"
#include "Slice3d.h"

//-------------------------------------------------------------------------TYPES

//-------------------------------------------------------------------------CLASS

class F_OrientedHistogram
{
 public:	

  F_OrientedHistogram(int _nOrientation);

  int getSizeFeatureVector();

  /**
   * Extract a feature vector for a given supernode in a 2d slice
   */
  bool getFeatureVector(osvm_node *x,
                        Slice* slice,
                        int sid, int nsid);

  /**
   * Extract a feature vector for a given supernode in a 3d volume
   */
  bool getFeatureVector(osvm_node *x,
                        Slice3d* slice3d,
                        int sid, int nsid);

 private:
  int nOrientation;
};

#endif // F_ORIENTEDHISTOGRAM_H
