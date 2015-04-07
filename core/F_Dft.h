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

#ifndef F_Dft_H
#define F_Dft_H

#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

// SliceMe
#include "Feature.h"
#include "oSVM.h"
#include "Slice.h"

class F_Dft : public Feature
{
 public:	

  /**
   * Constructor
   */
  F_Dft(const char* image_name,
        int _n_rows, int _n_cols);

  ~F_Dft();

  bool getFeatureVectorForOneSupernode(osvm_node *n, Slice* slice, int supernodeId);

  int getSizeFeatureVectorForOneSupernode();

 private:
  cv::Mat mag;
  int n_cols;
  int n_rows;

};

#endif // F_Dft_H
