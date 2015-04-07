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

#ifndef COLORMAP_H
#define COLORMAP_H

#include "Supernode.h"

#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <string>

//------------------------------------------------------------------------------

#define COLORMAP_SIZE 64

/* Colormap from Matlab :
autumn=colormap(autumn)*255;
save('autumn.txt','-ascii','-tabs','autumn');
bone=colormap(bone)*255;
save('bone.txt','-ascii','-tabs','bone');
jet=colormap(jet)*255;
save('jet.txt','-ascii','-tabs','jet');
*/

enum eColorMapType
  {
    COLORMAP_PROBS,
    COLORMAP_AUTUMN,
    COLORMAP_BONE,
    COLORMAP_JET
  };

extern float colormap_autumn[COLORMAP_SIZE*3];

extern float colormap_bone[COLORMAP_SIZE*3];

extern float colormap_jet[COLORMAP_SIZE*3];

class Colormap
{
 public:

  static Colormap* pInstance;

  static Colormap* Instance()
  {    
    if (pInstance == 0)  // is it the first call?
      {
        //printf("[Colormap] Error : you have to load a configuration file\n");
        //exit(-1);
        pInstance = new Colormap(); // create unique instance
      }
    return pInstance; // address of unique instance
  }

  std::map<ulong, labelType>& get() { return classIdxToLabel; }

  void set(const char* colormapFilename)
  {
    printf("[Colormap] Setting colormap = %s\n", colormapFilename);
    getClassToLabelMap(colormapFilename, classIdxToLabel);
  }

 private:
  std::map<ulong, labelType> classIdxToLabel;

  // FIXME : functions below were duplicated from utils.cpp

  // FIXME : change basis, should use 256 instead of 255...
  //classIdx = b*1 + g*255 + r*255*255
  void classIdxToRGB(ulong classIdx, uchar& r, uchar& g, uchar& b)
  {
    int ir = classIdx / (int)pow(255.0f,2);
    r = (ir > 255)?255:(uchar)ir;
    classIdx -= r*pow(255.0f,2);
    int ig = classIdx / 255;
    g = (ig > 255)?255:(uchar)ig;
    classIdx -= g*255;
    int ib = classIdx;
    b = (ib > 255)?255:(uchar)ib;
  }


  void getClassToLabelMap(const char* colormapFilename, map<ulong, labelType>& classIdxToLabel)
  {
    // Load colormap information
    std::string line;
    int label = 0;
    ulong classIdx = 0;

    std::ifstream ifsCol(colormapFilename);
    if(ifsCol.fail()) {
      printf("[colormap] Error while loading %s\n", colormapFilename);
      exit(-1);
    }

    while(std::getline(ifsCol, line)) {
      sscanf(line.c_str(),"%d %lu", &label, &classIdx);
      uchar g,r,b;
      classIdxToRGB(classIdx,r,g,b);
      classIdxToLabel[classIdx]= (labelType)label;
    }
    ifsCol.close();
  }
};

#endif
