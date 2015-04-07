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

#ifndef HISTO_H
#define HISTO_H

#include "StatModel.h"

#include <cv.h>
#include <highgui.h>


class Histogram : public StatModel
{
 public:

  /**
  * Constructor. The init has to be called later on.
  */
  Histogram() {};

  Histogram(int _nBins);

  ~Histogram();

  void draw(IplImage* histImage);

  double getNbEntries();

  double getValue(int bin);

  void init(int _nBins);

  void load(const char* filename);

  void normalize();

  void print();

  void print(int nBinsPerLine);

  void reset();

  void save(const char* filename);

  // private:
  double * histData;
  int nBins;
};

#endif // HISTO_H
