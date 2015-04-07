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

#ifndef INFERENCE_GLOBALS_H
#define INFERENCE_GLOBALS_H

#include "globalsE.h"

#ifdef _WIN32
#define isnan(x) _isnan(x)
#define isinf(x) (!_finite(x))
#endif

//------------------------------------------------------------------------------

typedef float nodeCoeffType;
typedef float edgeCoeffType;

//------------------------------------------------------------------------------

// set this to 1 to lear an offset for each class
// comment line to deactivate offset (do not set to 0)
#define W_OFFSET 0

#define INCLUDE_LOCAL_EDGES 1
//#define INCLUDE_LOCAL_EDGES 0

#define INCLUDE_GLOBAL_CLASSIFIER_IN_UNARY_TERM 1

#define SAMPLING_MUL_COEFF 10.0
//#define SAMPLING_MUL_COEFF 2.0

#define LOSS_VOC 2
#define LOSS_NODE_BASED 3

#define MAX_DISTANCE_DT 7

//------------------------------------------------------------------------------
// DO NOT EDIT

// macros used to index w vector
#if USE_LONG_RANGE_EDGES
#define SVM_FEAT_INDEX0(param) (((param)->nDistances*(param)->nGradientLevels*(param)->nOrientations*(param)->nClasses*(param)->nClasses) + \
                                (param)->nUnaryWeights + ((param)->nGradientLevels==0 && (param)->includeLocalEdges))
#else
#define SVM_FEAT_INDEX0(param) (((param)->nGradientLevels*(param)->nOrientations*(param)->nClasses*(param)->nClasses) + \
                                (param)->nUnaryWeights + ((param)->nGradientLevels==0 && (param)->includeLocalEdges))
#endif

#define SVM_FEAT_NUM_CLASSES(param) ((param)->nUnaryWeights)
#define SVM_FEAT_INDEX(param,c,f) (SVM_FEAT_INDEX0(param) + ((f)*SVM_FEAT_NUM_CLASSES(param))+(c))

#define ORIENTATION_SENSITIVE 1

#define COMPUTE_GLOBAL_LEVELS 1

// HACK for MSRC database in order to remove the global nodes for the last 3 classes
#define MSRC_NB_CLASSES_TO_REMOVE 3

#define MVC_MAX_ITER 200
#define INFERENCE_MAX_ITER 100

#define ONE_PARAM_PER_CLASS 0

// Metric
#define METRIC_SUPERNODE_BASED_01 0
#define METRIC_NODE_BASED_01 1

#define INFERENCE_VERBOSE 1
#define INFERENCE_PRINT(format, ...) if(INFERENCE_VERBOSE) printf (format, ## __VA_ARGS__)

//---------------------------------------------------------------------FUNCTIONS

#endif //INFERENCE_GLOBALS_H
