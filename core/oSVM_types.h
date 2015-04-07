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

#ifndef OSVM_TYPES_H
#define OSVM_TYPES_H

#include "svm.h"

//------------------------------------------------------------------------------

#ifdef USE_LIBLINEAR
typedef struct feature_node osvm_node;
typedef struct parameter osvm_parameter;
typedef struct problem osvm_problem;
typedef struct model osvm_model;
typedef int osvm_label_type;

#define osvm_destroy_param destroy_param
#define osvm_load_model load_model
#define osvm_save_model save_model
#define osvm_free_and_destroy_model free_and_destroy_model
#define osvm_destroy_param destroy_param
#define osvm_predict predict
#define osvm_predict_probability predict_probability
#define osvm_cross_validation cross_validation
#define osvm_get_labels get_labels

//#define osvm_train(p, m, s, o) train(p, m)
struct model* osvm_train(const struct problem *prob, const struct parameter *param,
                         bool** lSV = 0, float* obj = 0);

#else
typedef struct svm_node osvm_node;
typedef struct svm_parameter osvm_parameter;
typedef struct svm_problem osvm_problem;
typedef struct svm_model osvm_model;
typedef double osvm_label_type;

#define osvm_destroy_param svm_destroy_param
#define osvm_load_model svm_load_model
#define osvm_save_model svm_save_model
#define osvm_free_and_destroy_model svm_free_and_destroy_model
#define osvm_destroy_param svm_destroy_param
#define osvm_predict svm_predict
#define osvm_predict_probability svm_predict_probability
#define osvm_cross_validation svm_cross_validation
#define osvm_get_labels svm_get_labels
#define osvm_train svm_train
#endif

#endif //OSVM_TYPES_H
