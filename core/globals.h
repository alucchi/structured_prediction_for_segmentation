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

#ifndef GLOBALS_H
#define GLOBALS_H

#include <string>

bool verbose = false;

int BACKGROUND = 0;
int FOREGROUND = 1;
int BOUNDARY = 2;
int BACKGROUND_HYBRID = 3; // do not sample this type
int FOREGROUND_HYBRID = 4; // do not sample this type
int UNKNOWN_LABEL = 5;
int OTHER_LABEL = 6;

int MAX_INTENSITY_GRADIENT = 115;

int SUPERPIXEL_DEFAULT_STEP_SIZE = 7;
int SUPERPIXEL_DEFAULT_M = 10;
int DEFAULT_VOXEL_STEP = 7;
int SUPERVOXEL_DEFAULT_CUBENESS = 40.0;

// minimum percent pixels require to assign a label to a supernode
float MIN_PERCENT_TO_ASSIGN_LABEL = 0.25;

int DEFAULT_FEATURE_DISTANCE = 2;

int DEFAULT_LONG_RANGE_EDGES_DISTANCE = 2;

typedef unsigned long ulong;

#endif //GLOBALS_H
