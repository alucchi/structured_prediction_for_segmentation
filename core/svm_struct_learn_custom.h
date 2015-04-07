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

#ifndef SVM_STRUCT_LEARN_CUSTOM_H
#define SVM_STRUCT_LEARN_CUSTOM_H

#include "svm_struct_api.h"
#include "svm_light/svm_common.h"
#include "svm_struct/svm_struct_common.h"
#include "svm_struct/svm_struct_learn.h"

#include "constraint_set.h"

enum eInitWeightType
  {
    INIT_WEIGHT_0,
    INIT_WEIGHT_1,
    INIT_WEIGHT_RANDOM
  };

enum eGradientType
  {
    GRADIENT_GT = 0,
    GRADIENT_DIRECT_ADD, // add the loss : "away-from-worse"
    GRADIENT_DIRECT_SUBTRACT // subtract the loss : "toward-better"
  };

enum eLossType
  {
    LINEAR_LOSS = 0,
    LOG_LOSS,
    HINGE_LOSS,
    SQUARE_LOG_LOSS,
    SQUARE_HINGE_LOSS
  };

enum eUpdateType
  {
    UPDATE_CONSTANT = 0,
    UPDATE_DECREASING,
    UPDATE_MOMENTUM,
    UPDATE_MOMENTUM_DECREASING,
    UPDATE_PASSIVE_AGGRESSIVE
  };

typedef struct struct_gradient_parm {
  double learning_rate_0;
  double learning_rate;
  double learning_rate_exponent;
  bool enforce_submodularity;
  eUpdateType update_type;
  eLossType loss_type;
  float momentum_weight;
  float regularization_weight;
  eGradientType gradient_type;
  bool use_history;
  eSortingType constraint_set_type;
  bool ignore_loss;
  EXAMPLE* examples_all;
  int n_total_examples;
  int n_batch_examples;
  double* momentum;
  double max_norm_w;
  bool use_random_weights;
} GRADIENT_PARM;


//---------------------------------------------------------------------FUNCTIONS

/**
 * Compute gradient using history.
 */
double compute_gradient_with_history(STRUCT_LEARN_PARM *sparm, STRUCTMODEL *sm,
                                     EXAMPLE* ex,
                                     GRADIENT_PARM* gparm, SWORD* fy_to,
                                     double *dfy, double* loss);

double compute_gradient_with_history(STRUCT_LEARN_PARM *sparm, STRUCTMODEL *sm,
                                     EXAMPLE* ex, LABEL* y_bar, LABEL* y_direct,
                                     GRADIENT_PARM* gparm, SWORD* fy_to, SWORD* fy_away,
                                     double *dfy, double* loss);

/**
 * Compute m = delta_score + loss.
 */
double compute_m(STRUCT_LEARN_PARM *sparm,
                 STRUCTMODEL *sm, EXAMPLE *ex, long nExamples,
                 GRADIENT_PARM* gparm, LABEL* ybar, LABEL* y_direct,
                 SWORD* fy, SWORD* fybar, double *dfy);

double do_gradient_step(STRUCT_LEARN_PARM *sparm,
                        STRUCTMODEL *sm, EXAMPLE *ex, long nExamples,
                        GRADIENT_PARM* gparm,
                        double* momentum, double& dscore,
                        LABEL* y_bar);

double do_gradient_step(STRUCT_LEARN_PARM *sparm,
                        STRUCTMODEL *sm, EXAMPLE *ex, long nExamples,
                        GRADIENT_PARM* gparm,
                        double* momentum,
                        SWORD* fy, SWORD* fybar, double *dfy, double& dscore,
                        LABEL* y_bar);

void exportLabels(STRUCT_LEARN_PARM *sparm, EXAMPLE* ex,
                  LABEL* y, const char* dir_name);

void init_sampling(STRUCT_LEARN_PARM *sparm, STRUCTMODEL *sm, EXAMPLE* examples, long nExamples);

void init_w(STRUCT_LEARN_PARM *sparm, STRUCTMODEL *sm, GRADIENT_PARM* gparm, EXAMPLE* ex, int init_type);

void init_gradient_param(GRADIENT_PARM& gparm, Config* config,
                         ConstraintSet* constraint_set);

void enforce_submodularity(STRUCT_LEARN_PARM *sparm, double* smw);

void project_w(STRUCTMODEL *sm, GRADIENT_PARM* gparm);

void write_vector(const char* filename, double* v, int size_v);

void write_vector(const char* filename, SWORD* v);

//void write_scalar(const char* filename, double* v, int size_v);

#endif //SVM_STRUCT_LEARN_CUSTOM_H
