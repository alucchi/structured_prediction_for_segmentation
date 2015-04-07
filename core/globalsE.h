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

#ifndef GLOBALSE_H
#define GLOBALSE_H

#include <vector>

using namespace std;

//-----------------------------------------------------------------------------
// Experimenting...

#define USE_LONG_RANGE_EDGES 0
#define MAX_SQ_DISTANCE_LONG_RANGE_EDGES 500

#define USE_SPARSE_VECTORS 0

//-----------------------------------------------------------------------------

enum eFeatureType
{
  F_UNKNOWN = 0,
  F_HISTOGRAM = 1,
  F_COLOR_HISTOGRAM = 2,
  F_LOADFROMFILE = 4,
  F_GLCM = 8,
  F_POSITION = 16,
  F_FILTER = 32,
  F_GAUSSIAN = 64,
  F_GRADIENTSTATS = 128,
  F_BIAS = 256,
  F_DFT = 512,
  F_SIFT = 1024,
  F_END_FEATURETYPE = 2048
};

// Background, foreground and boundary have to be assigned to the first 3 labels.
enum eType
{
  T_BACKGROUND = 0,
  T_FOREGROUND,
  T_BOUNDARY,
  T_BACKGROUND_HYBRID, // do not sample this type
  T_FOREGROUND_HYBRID, // do not sample this type
  T_UNKNOWN_LABEL,
  T_OTHER_LABEL,
  NUMBER_TYPE
};

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned long ulong;

struct label_index
{
  //int backgroundIndex;
  vector<int> backgroundIndex;
  int boundaryIndex;
  int mitochondriaIndex; // or any object that has to be segmented

  label_index()
  {
    //backgroundIndex=-1;
   boundaryIndex=-1;
   mitochondriaIndex=-1; 
  }
};

struct point_3df
{
  float x;
  float y;
  float z;
};

struct s_RGB
{
  uchar r;
  uchar g;
  uchar b;
};

//typedef ulong sizeSliceType;
typedef int sizeSliceType;

//-----------------------------------------------------------------------------

extern int BACKGROUND;
extern int BOUNDARY;
extern int FOREGROUND;
extern int BACKGROUND_HYBRID; // do not sample this type
extern int FOREGROUND_HYBRID; // do not sample this type
extern int UNKNOWN_LABEL;
extern int OTHER_LABEL;

extern int MAX_INTENSITY_GRADIENT;

extern int SUPERPIXEL_DEFAULT_STEP_SIZE;
extern int SUPERPIXEL_DEFAULT_M;
extern int DEFAULT_VOXEL_STEP;
extern int SUPERVOXEL_DEFAULT_CUBENESS;

extern float MIN_PERCENT_TO_ASSIGN_LABEL;

// define the extent (i.e. number of supernodes) to which the features are computed
extern int DEFAULT_FEATURE_DISTANCE;

extern int DEFAULT_LONG_RANGE_EDGES_DISTANCE;

//-----------------------------------------------------------------------------

#define DEFAULT_FEATURE_TYPE (F_HISTOGRAM | F_FILTER | F_BIAS)
//#define DEFAULT_FEATURE_TYPE (F_HISTOGRAM | F_FILTER)

#define MAX_PIXEL_VALUE 256
#define NB_BINS 10

#define MAX_INTENSITY 255
//#define MAX_INTENSITY_GRADIENT 115 // gradient for imagesLeft_35-199_top
//#define DELTA_PB 1e-37 // probability smoothing factor
#define DELTA_PB 1e-10 // probability smoothing factor

#define INVALID_LABEL 255 //uchar

// number of pixels sampled in each supervoxels (around 5%)
//#define NB_SAMPLED_PIXELS 50
#define NB_SAMPLED_PIXELS 10

// white pixels indicate foreground region
#define FOREGROUND_MASKVALUE 255
//#define FOREGROUND_ADVANCED_MASKVALUE 90
#define FOREGROUND_ADVANCED_MASKVALUE 79
#define BACKGROUND_ADVANCED_MASKVALUE 13
#define BACKGROUND_MASKVALUE 0

#define SUPERPIXEL_IMAGE -2
#define SUPERPIXEL_CUBE -1

#define SVM_MAX_LINE_LENGTH 500

#define ORIENTED_HIST_NORIENTATION 8

#define SVM_DEFAULT_EPS 1e-3

#define PI 3.141592653589793

#define FEATURE_FIELD_SEPARATOR ':'

// assume features files are 1-based index (which is the case for files created
// by libsvmwrite)
#define FEATURE_FIRST_INDEX 1

extern bool verbose;

#define PRINT_MESSAGE(format, ...) if(verbose) printf (format, ## __VA_ARGS__)

#endif //GLOBALSE_H
