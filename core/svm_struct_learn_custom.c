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

// This code is based on the template provided by Thorsten Joachims.

#include <iomanip>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#ifdef _WIN32
#include <io.h>
#include <direct.h>
#define mkdir(x,y) _mkdir(x)
#include "gettimeofday.h"
#else
#include <unistd.h>
#include <sys/time.h>
#endif

#include "constraint_set.h"
#include "label_cache.h"
#include "svm_struct_learn_custom.h"
#include "svm_struct_api.h"
#include "svm_light/svm_common.h"
#include "svm_struct/svm_struct_common.h"
#include "svm_struct/svm_struct_learn.h"
#include "svm_struct_globals.h"

#include "Config.h"

#include "highgui.h"

#include "energyParam.h"
#include "graphInference.h"
#include "inference.h"

//------------------------------------------------------------------------MACROS

#define BUFFER_SIZE 250

// If greater than 1, output dscore, norm(dfy), loss
// If greater than 2, output dfy
#define CUSTOM_VERBOSITY 3

#define CUSTOM_VERBOSITY_F(X, Y) if(CUSTOM_VERBOSITY > X) { Y }

//---------------------------------------------------------------------FUNCTIONS

void write_vector(const char* filename, double* v, int size_v)
{
  ofstream ofs(filename, ios::app);
  for(int i = 0; i < size_v; ++i) {
    ofs << v[i] << " ";
  }
  ofs << endl;
  ofs.close();
}

/**
 * Write vector to a file (don't overwrite but append new line).
 */
void write_vector(const char* filename, SWORD* v)
{
  ofstream ofs(filename, ios::app);
  SWORD* w = v;
  while (w->wnum) {
    ofs << w->weight << " ";
    ++w;
  }
  ofs << endl;
  ofs.close();
}

/**
 * Write set of scalar values to a file (don't overwrite but append new lines).
 */
void write_scalars(const char* filename, double* v, int size_v)
{
  ofstream ofs(filename, ios::app);
  for(int i = 0; i < size_v; ++i) {
    ofs << v[i] << endl;
  }
  ofs.close();
}

/**
 * Write scalar value to a file (don't overwrite but append new line).
 */
void write_scalar(const char* filename, double v)
{
  ofstream ofs(filename, ios::app);
  ofs << v << endl;
  ofs.close();
}

/**
 * Returns squared norm
 */
double get_sq_norm(double* v, int _sizePsi)
{
  double sq_norm_v = 0;
  for(int i = 1; i < _sizePsi; ++i) {
    sq_norm_v += v[i]*v[i];
  }
  return sq_norm_v;
}

double get_norm(double* v, int _sizePsi)
{
  double norm_v = 0;
  for(int i = 1; i < _sizePsi; ++i) {
    norm_v += v[i]*v[i];
  }
  return sqrt(norm_v);
}

/**
 * Compute the average norm of psi over the training data
 */
double get_norm_psi_gt(STRUCT_LEARN_PARM *sparm, STRUCTMODEL *sm, EXAMPLE *examples, long nExamples)
{
  int _sizePsi = sm->sizePsi + 1;
  SWORD* fy_to = new SWORD[_sizePsi];
  double avg_norm = 0;

  for(long i = 0; i < nExamples; ++i) {
    computePsi(fy_to, examples[i].x, examples[i].y, sm, sparm);
    double norm_wy_to = 0;
    SWORD* wy_to = fy_to;
    while (wy_to->wnum) {
      norm_wy_to += wy_to->weight*wy_to->weight;
      ++wy_to;
    }
    norm_wy_to = sqrt(norm_wy_to);
    avg_norm += norm_wy_to;
  }
  avg_norm /= nExamples;

  delete[] fy_to;
  return avg_norm;
}

/**
 * accumulate gradient in dfy
 */
void compute_gradient_accumulate(STRUCTMODEL *sm, GRADIENT_PARM* gparm,
                                 SWORD* fy_to, SWORD* fy_away, double *dfy,
                                 const double loss, const double dfy_weight)
{
#if CUSTOM_VERBOSITY > 2
  double score_y = 0;
  double score_y_away = 0;
#endif

  SWORD* wy_to = fy_to;
  SWORD* wy_away = fy_away;
  switch(gparm->loss_type)
    {
    case LOG_LOSS:
      {
        // L(w) = log(1+e(m(x)))
        // where m(x) = (loss(y,y_bar) + score(x,y_bar)) - score(x,y)
        // and score(x,y) = w^T*psi(x,y)
        // dL(w)/dw = ( m'(x) e(m(x)) ) / ( 1 + e(m(x)))
        // m'(x) = psi(x,y_bar) - psi(x,y)
        double m = 0;
        double dm;
        while (wy_to->wnum) {
          while(wy_away->wnum && (wy_away->wnum < wy_to->wnum)) {
            ++wy_away;
          }

          if(wy_to->wnum == wy_away->wnum) {
            dm = wy_away->weight - wy_to->weight;
          } else {
            dm = - wy_to->weight;
          }
          m += (sm->w[wy_to->wnum]*dm);
          ++wy_to;
        }
        m += loss;
        double e_m = 0;
        if(m < 100) {
          e_m = exp(m);
        }

        wy_to = fy_to;
        wy_away = fy_away;
        while (wy_to->wnum) {
          while(wy_away->wnum && (wy_away->wnum < wy_to->wnum)) {
            ++wy_away;
          }
          if(wy_to->wnum == wy_away->wnum) {
            dm = wy_away->weight - wy_to->weight;
          } else {
            dm = - wy_to->weight;
          }

          if(m >= 100) {
            dfy[wy_to->wnum] += dfy_weight * dm;
          } else {
            dfy[wy_to->wnum] += dfy_weight * (dm*e_m / (e_m + 1));
          }
#if CUSTOM_VERBOSITY > 2
          score_y += sm->w[wy_to->wnum]*wy_to->weight;
          score_y_away += sm->w[wy_to->wnum]*wy_away->weight;
#endif
          ++wy_to;
        }
      }
      break;
    case HINGE_LOSS:
      {
        // L(w) = (loss(y,y_bar) + score(x,y_bar)) - score(x,y)
        // where score(x,y) = w^T*psi(x,y)
        // dL(w)/dw = psi(x,y_bar) - psi(x,y)
        double dm;
        while (wy_to->wnum) {
          while(wy_away->wnum && (wy_away->wnum < wy_to->wnum)) {
            ++wy_away;
          }

          if(wy_to->wnum == wy_away->wnum) {
            dm = wy_away->weight - wy_to->weight;
          } else {
            dm = - wy_to->weight;
          }

          dfy[wy_to->wnum] += dfy_weight * dm;
#if CUSTOM_VERBOSITY > 2
          score_y += sm->w[wy_to->wnum]*wy_to->weight;
          score_y_away += sm->w[wy_to->wnum]*wy_away->weight;
#endif
          ++wy_to;
        }
      }
      break;
    case SQUARE_HINGE_LOSS:
      {
        // L(w) = log(1+e(m(x)))
        // where m(x) = (loss(y,y_bar) + score(x,y_bar)) - score(x,y)
        // and score(x,y) = w^T*psi(x,y)
        // dL(w)/dw = ( m'(x) e(m(x)) ) / ( 1 + e(m(x)))
        // m'(x) = psi(x,y_bar) - psi(x,y)
        double m = 0;
        double dm;
        while (wy_to->wnum) {
          while(wy_away->wnum && (wy_away->wnum < wy_to->wnum)) {
            ++wy_away;
          }

          if(wy_to->wnum == wy_away->wnum) {
            dm = wy_away->weight - wy_to->weight;
          } else {
            dm = - wy_to->weight;
          }
          m += (sm->w[wy_to->wnum]*dm);
          ++wy_to;
        }
        m += loss;

        wy_to = fy_to;
        wy_away = fy_away;
        while (wy_to->wnum) {
          while(wy_away->wnum && (wy_away->wnum < wy_to->wnum)) {
            ++wy_away;
          }
          if(wy_to->wnum == wy_away->wnum) {
            dm = wy_away->weight - wy_to->weight;
          } else {
            dm = - wy_to->weight;
          }

          dfy[wy_to->wnum] += 1e-30 * dfy_weight * dm * m;
#if CUSTOM_VERBOSITY > 2
          score_y += sm->w[wy_to->wnum]*wy_to->weight;
          score_y_away += sm->w[wy_to->wnum]*wy_away->weight;
#endif
          ++wy_to;
        }
      }
      break;
    default:
      printf("[svm_struct_custom] Unknown loss type %d\n", gparm->loss_type);
      exit(-1);
      break;
    }

#if CUSTOM_VERBOSITY > 2
  ofstream ofs_score_y("score_y.txt", ios::app);
  ofs_score_y << score_y << endl;
  ofs_score_y.close();

  ofstream ofs_score_y_away("score_y_away.txt", ios::app);
  ofs_score_y_away << score_y_away << endl;
  ofs_score_y_away.close();
#endif
}

void compute_psi(STRUCT_LEARN_PARM *sparm, STRUCTMODEL *sm,
                   EXAMPLE* ex, LABEL* y_bar, LABEL* y_direct,
                   GRADIENT_PARM* gparm, SWORD* fy_to, SWORD* fy_away,
                   double* loss)
{
  labelType* y_to = 0;
  labelType* y_away = 0;
  switch(gparm->gradient_type) {
  case GRADIENT_GT:
    // moves toward ground truth, away from larger loss
    y_to = ex->y.nodeLabels;
    y_away = y_bar->nodeLabels;
    computePsi(fy_to, ex->x, ex->y, sm, sparm);
    computePsi(fy_away, ex->x, *y_bar, sm, sparm);
    break;
  case GRADIENT_DIRECT_ADD:
    // moves away from larger loss
    y_to = y_direct->nodeLabels;
    y_away = y_bar->nodeLabels;
    computePsi(fy_to, ex->x, *y_direct, sm, sparm);
    computePsi(fy_away, ex->x, *y_bar, sm, sparm);
    break;
  case GRADIENT_DIRECT_SUBTRACT:
    // moves toward better label
    y_to = y_direct->nodeLabels;
    y_away = y_bar->nodeLabels;
    computePsi(fy_to, ex->x, *y_direct, sm, sparm);
    computePsi(fy_away, ex->x, *y_bar, sm, sparm);
    break;
  default:
    printf("[svm_struct_custom] Unknown gradient type\n");
    exit(-1);
    break;
  }

  if(!gparm->ignore_loss) {
    int nDiff;
    double _loss;
    computeLoss(y_to, y_away, ex->y.nNodes, sparm, _loss, nDiff);
    if(loss) {
      *loss = _loss;
    }
  }
}

void compute_psi_to(STRUCT_LEARN_PARM *sparm, STRUCTMODEL *sm,
                    EXAMPLE* ex, GRADIENT_PARM* gparm, SWORD* fy_to)
{
  switch(gparm->gradient_type) {
  case GRADIENT_GT:
    // moves toward ground truth, away from larger loss
    computePsi(fy_to, ex->x, ex->y, sm, sparm);
    break;
    /*
  case GRADIENT_DIRECT_ADD:
    // moves away from larger loss
    computePsi(fy_to, ex->x, *y_direct, sm, sparm);
    break;
  case GRADIENT_DIRECT_SUBTRACT:
    // moves toward better label
    computePsi(fy_to, ex->x, *y_direct, sm, sparm);
    break;
    */
  default:
    printf("[svm_struct_custom] Unknown gradient type\n");
    exit(-1);
    break;
  }
}

double compute_gradient_accumulate(STRUCT_LEARN_PARM *sparm, STRUCTMODEL *sm,
                                   EXAMPLE* ex, LABEL* y_bar, LABEL* y_direct,
                                   GRADIENT_PARM* gparm, SWORD* fy_to, SWORD* fy_away,
                                   double *dfy, double* loss, const double dfy_weight)
{
  int _sizePsi = sm->sizePsi + 1;
  double _loss;
  compute_psi(sparm, sm, ex, y_bar, y_direct, gparm, fy_to, fy_away, &_loss);
  if(loss) {
    *loss = _loss;
  }

  compute_gradient_accumulate(sm, gparm, fy_to, fy_away, dfy, _loss, dfy_weight);

#if CUSTOM_VERBOSITY > 3
  write_vector("dfy.txt", dfy, _sizePsi);
#endif

  double dscore = 0;
  // do not add +1 here as dfy also has an additional dummy entry at index 0.
  double* smw = sm->w;
  for(int i = 0; i < _sizePsi; ++i) {
    dscore += smw[i]*dfy[i];
  }
  return dscore;
}

double compute_gradient(STRUCT_LEARN_PARM *sparm, STRUCTMODEL *sm,
                        EXAMPLE* ex, LABEL* y_bar, LABEL* y_direct,
                        GRADIENT_PARM* gparm, SWORD* fy_to, SWORD* fy_away,
                        double *dfy, double* loss, const double dfy_weight)
{
  // initialize dfy to 0
  int _sizePsi = sm->sizePsi + 1;
  for(int i = 0; i < _sizePsi; ++i) {
    dfy[i] = 0;
  }

  return compute_gradient_accumulate(sparm, sm, ex, y_bar, y_direct, gparm, fy_to, fy_away, dfy, loss, dfy_weight);
}

double compute_gradient(STRUCTMODEL *sm, GRADIENT_PARM* gparm,
                        SWORD* fy_to, SWORD* fy_away, double *dfy,
                        const double loss, const double dfy_weight)
{
  // initialize dfy to 0
  int _sizePsi = sm->sizePsi + 1;
  for(int i = 0; i < _sizePsi; ++i) {
    dfy[i] = 0;
  }

  compute_gradient_accumulate(sm, gparm, fy_to, fy_away, dfy, loss, dfy_weight);

  double dscore = 0;
  // do not add +1 here as dfy also has an additional dummy entry at index 0.
  double* smw = sm->w;
  for(int i = 0; i < _sizePsi; ++i) {
    dscore += smw[i]*dfy[i];
  }
  return dscore;
}

void exportLabels(STRUCT_LEARN_PARM *sparm, EXAMPLE* ex,
                  LABEL* y, const char* dir_name)
{
  string paramSlice3d;
  Config::Instance()->getParameter("slice3d", paramSlice3d);
  bool useSlice3d = paramSlice3d.c_str()[0] == '1';
  string paramVOC;
  Config::Instance()->getParameter("voc", paramVOC);
  bool useVOC = paramVOC.c_str()[0] == '1';

  stringstream ss_dir;
  ss_dir << dir_name;
  mkdir(ss_dir.str().c_str(), 0777);
  if(!useSlice3d) {
    //TODO : Remove !useVOC
    if(useVOC) {
      ss_dir << "x" << sparm->iterationId;
    }
    else {
      ss_dir << "x" << sparm->iterationId << "/";
    }
  }
  mkdir(ss_dir.str().c_str(), 0777);

  stringstream soutColoredImage;
  soutColoredImage << ss_dir.str();
  if(useSlice3d) {
    soutColoredImage << getNameFromPathWithoutExtension(ex->x.slice->getName());
    soutColoredImage << "_";
    soutColoredImage << sparm->iterationId;
  } else {
    soutColoredImage << ex->x.slice->getName();
  }

  ex->x.slice->exportSupernodeLabels(soutColoredImage.str().c_str(),
                                     sparm->nClasses,
                                     y->nodeLabels,
                                     y->nNodes,
                                     &(sparm->labelToClassIdx));

  if(useSlice3d) {
    zipAndDeleteCube(soutColoredImage.str().c_str());
  }
}

double do_gradient_step(STRUCT_LEARN_PARM *sparm,
                        STRUCTMODEL *sm, EXAMPLE *ex, long nExamples,
                        GRADIENT_PARM* gparm,
                        double* momentum, double& dscore, LABEL* y_bar)
{
  int _sizePsi = sm->sizePsi + 1;
  SWORD* fy_to = new SWORD[_sizePsi];
  SWORD* fy_away = new SWORD[_sizePsi];
  double* dfy = new double[_sizePsi];
  memset((void*)dfy, 0, sizeof(double)*(_sizePsi));

  double m = do_gradient_step(sparm, sm, ex, nExamples, gparm,
                              momentum, fy_to, fy_away, dfy, dscore, y_bar);
  delete[] fy_to;
  delete[] fy_away;
  delete[] dfy;
  return m;
}

double compute_gradient_with_history(STRUCT_LEARN_PARM *sparm, STRUCTMODEL *sm,
                                     EXAMPLE* ex,
                                     GRADIENT_PARM* gparm, SWORD* fy_to,
                                     double *dfy, double* loss)
{
  ConstraintSet* cs = ConstraintSet::Instance();
  const vector< constraint >* constraints = cs->getConstraints(ex->x.id);
  assert(constraints != 0);
  int n_cs = constraints->size();

  double* dfy_weights = new double[n_cs];
  if(gparm->use_random_weights) {
    double total_weights = 0;
    for(int c = 0; c < n_cs; ++c) {
      dfy_weights[c] = rand() * ((double)n_cs/(double)RAND_MAX);
      total_weights += dfy_weights[c];
    }
    for(int c = 0; c < n_cs; ++c) {
      dfy_weights[c] /= total_weights;
    }
  } else {
    double total_weights = 0;
    for(int c = 0; c < n_cs; ++c) {
      dfy_weights[c] = 1.0/(double)(n_cs);
      total_weights += dfy_weights[c];
    }
    for(int c = 0; c < n_cs; ++c) {
      dfy_weights[c] /= total_weights;
    }
  }

  int _sizePsi = sm->sizePsi + 1;
  // initialize dfy to 0
  for(int i = 0; i < _sizePsi; ++i) {
    dfy[i] = 0;
  }

  if(loss) {
    *loss = 0;
  }

  // add gradient for history of constraints
  if(gparm->loss_type != HINGE_LOSS && gparm->loss_type != SQUARE_HINGE_LOSS) {
    // use all the constraints in the set
    int c = 0;
    for(vector<constraint>::const_iterator it = constraints->begin();
        it != constraints->end(); ++it) {
      compute_gradient_accumulate(sm, gparm, fy_to,
                                  it->first->w, dfy, it->first->loss, dfy_weights[c]);
      if(loss) {
        *loss += it->first->loss;
      }      
      ++c;
    }
  } else {
    // only use violated constraints

    double score_gt = computeScore(sm, fy_to);
    int c = 0;
    for(vector<constraint>::const_iterator it = constraints->begin();
        it != constraints->end(); ++it) {
      // check if constraint is violated
      double score_cs = computeScore(sm, it->first->w);
      bool positive_margin = (score_cs - score_gt + it->first->loss) > 0;
      //printf("Margin constraint %d: score_cs = %g, score_gt = %g, loss = %g, margin = %g\n",
      //       c, score_cs, score_gt, it->first->loss, score_cs - score_gt + it->first->loss);

      if(positive_margin) {
        compute_gradient_accumulate(sm, gparm, fy_to,
                                    it->first->w, dfy, it->first->loss, dfy_weights[c]);
        if(loss) {
          *loss += it->first->loss;
        }
      }     
      ++c;
    }
  }

  double total_dscore = 0;
  // do not add +1 here as dfy also has an additional dummy entry at index 0.
  double* smw = sm->w;
  for(int i = 0; i < _sizePsi; ++i) {
    total_dscore += smw[i]*dfy[i];
  }

  delete[] dfy_weights;

  return total_dscore;
}

double compute_gradient_with_history(STRUCT_LEARN_PARM *sparm, STRUCTMODEL *sm,
                                     EXAMPLE* ex, LABEL* y_bar, LABEL* y_direct,
                                     GRADIENT_PARM* gparm, SWORD* fy_to, SWORD* fy_away,
                                     double *dfy, double* loss)
{
  double dfy_weight = 1.0;
  int _sizePsi = sm->sizePsi + 1;
  // initialize dfy to 0
  for(int i = 0; i < _sizePsi; ++i) {
    dfy[i] = 0;
  }

  double _loss;
  double _dscore = compute_gradient(sparm, sm, ex, y_bar, y_direct, gparm, fy_to,
                                    fy_away, dfy, &_loss, dfy_weight);
  if(loss) {
    *loss += _loss;
  }

  // add gradient for history of constraints
  ConstraintSet* cs = ConstraintSet::Instance();
  const vector< constraint >* constraints = cs->getConstraints(ex->x.id);
  if(constraints) {
    dfy_weight = 1.0/(double)(constraints->size()+1.0);
    for(vector<constraint>::const_iterator it = constraints->begin();
        it != constraints->end(); ++it) {
      compute_gradient_accumulate(sm, gparm, fy_to,
                                  it->first->w, dfy, it->first->loss, dfy_weight);
      if(loss) {
        *loss += it->first->loss;
      }
    }
  }

  double total_dscore = 0;
  // do not add +1 here as dfy also has an additional dummy entry at index 0.
  double* smw = sm->w;
  for(int i = 0; i < _sizePsi; ++i) {
    total_dscore += smw[i]*dfy[i];
  }
  total_dscore += _dscore;
  return total_dscore;
}

void update_w(STRUCT_LEARN_PARM *sparm, STRUCTMODEL *sm, GRADIENT_PARM* gparm,
              double* momentum, double *dfy)
{
  int _sizePsi = sm->sizePsi + 1;

  // do not add +1 here as dfy also has an additional dummy entry at index 0.
  double* smw = sm->w;
  if(momentum) {
    // update momentum
    for(int i = 1; i < _sizePsi; ++i) {
      momentum[i] = (gparm->learning_rate*(dfy[i] + (gparm->regularization_weight*smw[i])) + gparm->momentum_weight*momentum[i]);
    }
    for(int i = 1; i < _sizePsi; ++i) {
      smw[i] -= momentum[i];
    }
  } else {
    for(int i = 1; i < _sizePsi; ++i) {
      smw[i] -= gparm->learning_rate*(dfy[i]+(gparm->regularization_weight*smw[i]));
    }
  }
}

double do_gradient_step(STRUCT_LEARN_PARM *sparm,
                        STRUCTMODEL *sm, EXAMPLE *ex, long nExamples,
                        GRADIENT_PARM* gparm,
                        double* momentum,
                        SWORD* fy_to, SWORD* fy_away, double *dfy,
                        double& dscore,
                        LABEL* y_bar)
{
  int _sizePsi = sm->sizePsi + 1;
  LABEL* y_direct = 0;

  double* _lossPerLabel = sparm->lossPerLabel;
  if(gparm->ignore_loss) {
    sparm->lossPerLabel = 0;
  }

  // setting this to 1 will make the example loop below single thread so that
  // several threads can be run for different temperature while running the
  // samplign code.
#define USE_SAMPLING 0

#if USE_SAMPLING
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
#endif

  /*** precomputation step ***/
  for(int i = 0; i < nExamples; i++) {

#if USE_SAMPLING
#ifdef USE_OPENMP
    int threadId = omp_get_thread_num();
    printf("[svm_struct_custom] Thread %d/%d\n", threadId,omp_get_num_threads());
#endif
#endif

    if(sparm->loss_type == SLACK_RESCALING) {
      y_bar[i] = find_most_violated_constraint_slackrescaling(ex[i].x, ex[i].y,
                                                             sm, sparm);
    } else {
      y_bar[i] = find_most_violated_constraint_marginrescaling(ex[i].x, ex[i].y,
                                                              sm, sparm);
    }
  }

  if(gparm->ignore_loss) {
    sparm->lossPerLabel = _lossPerLabel;
  }

  if(gparm->gradient_type == GRADIENT_DIRECT_ADD ||
     gparm->gradient_type == GRADIENT_DIRECT_SUBTRACT) {

    // allocate memory
    y_direct = new LABEL[nExamples];

    // temporarily remove loss
    double* _lossPerLabel = sparm->lossPerLabel;
    sparm->lossPerLabel = 0;

    for(int il = 0; il < nExamples; il++) {

#ifdef USE_OPENMP
      int threadId = omp_get_thread_num();
      printf("[svm_struct_custom] Thread %d/%d\n", threadId,omp_get_num_threads());
#else
      int threadId = 0;
#endif

      // check if labels are stored in the cache
      int cacheId = nExamples + ex[il].x.id;
      bool labelFound = LabelCache::Instance()->getLabel(cacheId, *y_direct);
      if(!labelFound) {
        // allocate memory
        y_direct->nNodes = ex[il].y.nNodes;
        y_direct->nodeLabels = new labelType[y_direct->nNodes];
        for(int n = 0; n < ex[il].y.nNodes; n++) {
          y_direct->nodeLabels[n] = ex[il].y.nodeLabels[n];
        }
        y_direct->cachedNodeLabels = false;
        labelFound = true;
      }

      runInference(ex[il].x, ex[il].y, sm, sparm, y_direct[il], threadId, labelFound, cacheId);
      //exportLabels(sparm, &ex[il], y_bar, "direct/");

    }
    sparm->lossPerLabel = _lossPerLabel;
  }

#if CUSTOM_VERBOSITY > 2
  ofstream ofs_cs_dscore("constraint_set_dscore.txt", ios::app);
#endif
  int n_satisfied = 0;
  int n_not_satisfied = 0;

  const double dfy_weight = 1.0;
  double total_dscore = 0;
  double total_dloss = 0;

  if(gparm->constraint_set_type == CS_USE_MVC) {

    ConstraintSet* cs = ConstraintSet::Instance();

    for(int il = 0; il < nExamples; il++) { /*** example loop ***/

      double _loss = 0;
      compute_gradient(sparm, sm, &ex[il], &y_bar[il], &y_direct[il], gparm,
                       fy_to, fy_away, dfy, &_loss, dfy_weight);

      // add the current constraint first
      if(gparm->constraint_set_type == CS_MARGIN || gparm->constraint_set_type == CS_MARGIN_DISTANCE) {
        double margin = total_dscore + total_dloss;
        double sorting_value = (fabs(margin) < 1e-38)?0 : 1.0/margin;
        cs->add(ex[il].x.id, fy_away, _loss, _sizePsi, sorting_value);
      } else {
        cs->add(ex[il].x.id, fy_away, _loss, _sizePsi);
      }

      const constraint* c = cs->getMostViolatedConstraint(ex[il].x.id, sm->w);
      double dscore_cs = compute_gradient(sm, gparm, fy_to,
                                          c->first->w, dfy, c->first->loss, dfy_weight);
      bool positive_margin = (dscore_cs + c->first->loss) > 0;

      if( (gparm->loss_type != HINGE_LOSS && gparm->loss_type != SQUARE_HINGE_LOSS) || positive_margin) {
        update_w(sparm, sm, gparm, momentum, dfy);
      }

#if CUSTOM_VERBOSITY > 2
      ofs_cs_dscore << dscore_cs << "," << c->first->loss << endl;
#endif
    }

  } else {
    // use all the constraints in the working set instead of the most violated one

    for(int il = 0; il < nExamples; il++) { /*** example loop ***/

      // compute gradient for last generated constraint
      double _loss;
      double _dscore = compute_gradient(sparm, sm, &ex[il], &y_bar[il], &y_direct[il], gparm,
                                        fy_to, fy_away, dfy, &_loss, dfy_weight);

      if(gparm->use_history) {

        // add last generated constraint to the working set
        ConstraintSet* cs = ConstraintSet::Instance();
        if(gparm->constraint_set_type == CS_MARGIN || gparm->constraint_set_type == CS_MARGIN_DISTANCE) {
          double margin = _dscore + _loss;
          double sorting_value = (fabs(margin) < 1e-38)?0 : 1.0/margin;
          cs->add(ex[il].x.id, fy_away, _loss, _sizePsi, sorting_value);
        } else {
          cs->add(ex[il].x.id, fy_away, _loss, _sizePsi);
        }

        const vector< constraint >* constraints = cs->getConstraints(ex[il].x.id);

        if(constraints) {
          for(vector<constraint>::const_iterator it = constraints->begin();
              it != constraints->end(); ++it) {
            double dscore_cs = compute_gradient(sm, gparm, fy_to,
                                                it->first->w, dfy, it->first->loss, dfy_weight);
            total_dloss += it->first->loss;
            bool positive_margin = (dscore_cs + it->first->loss) > 0;

            if( (gparm->loss_type != HINGE_LOSS && gparm->loss_type != SQUARE_HINGE_LOSS) || positive_margin) {
              update_w(sparm, sm, gparm, momentum, dfy);
              total_dscore += dscore_cs;
            }

            if(positive_margin) {
              ++n_not_satisfied;
            } else {
              ++n_satisfied;
            }
#if CUSTOM_VERBOSITY > 2
            ofs_cs_dscore << dscore_cs << "," << it->first->loss << " ";
#endif
          }
        }

#if CUSTOM_VERBOSITY > 2
        ofs_cs_dscore << " , " << _dscore << endl;
#endif

        if(gparm->constraint_set_type == CS_MARGIN || gparm->constraint_set_type == CS_MARGIN_DISTANCE) {
          // compute margin
          double margin = total_dscore + total_dloss;
          double sorting_value = (fabs(margin) < 1e-38)?0 : 1.0/margin;
          cs->add(ex[il].x.id, fy_away, _loss, _sizePsi, sorting_value);
        } else {
          cs->add(ex[il].x.id, fy_away, _loss, _sizePsi);
        }
      } else {
        bool positive_margin = (_dscore + _loss) > 0;
        if( (gparm->loss_type != HINGE_LOSS && gparm->loss_type != SQUARE_HINGE_LOSS) || positive_margin) {
          update_w(sparm, sm, gparm, momentum, dfy);
          total_dscore += _dscore;
        }
      }
    }
  }

#if CUSTOM_VERBOSITY > 2
  ofs_cs_dscore.close();
#endif

#if CUSTOM_VERBOSITY > 1
  ofstream ofs_cs_card("constraint_set_card.txt", ios::app);
  ofs_cs_card << n_satisfied << " " << n_not_satisfied << " " << n_satisfied+n_not_satisfied << endl;
  ofs_cs_card.close();
#endif

#if CUSTOM_VERBOSITY > 1
  double sq_norm_dfy = 0;
  for(int i = 1; i < _sizePsi; ++i) {
    sq_norm_dfy += dfy[i]*dfy[i];
  }
  ofstream ofs_norm_dfy("norm_dfy.txt", ios::app);
  ofs_norm_dfy << sqrt(sq_norm_dfy) << endl;
  ofs_norm_dfy.close();

  ofstream ofs_dscore("dscore.txt", ios::app);
  ofs_dscore << total_dscore << endl;
  ofs_dscore.close();

  if(sparm->giType == T_GI_SAMPLING) {
    ofstream ofs_temp("temperature.txt", ios::app);
    ofs_temp << sparm->sampling_temperature_0 << endl;
    ofs_temp.close();
  }

#endif // CUSTOM_VERBOSITY

  dscore = total_dscore;

  double m = compute_m(sparm, sm, ex, nExamples, gparm, y_bar, y_direct, fy_to, fy_away, dfy);

  if(y_direct) {
    delete[] y_direct;
  }

  return m;
}

double compute_m(STRUCT_LEARN_PARM *sparm,
                 STRUCTMODEL *sm, EXAMPLE *ex, long nExamples,
                 GRADIENT_PARM* gparm, LABEL* y_bar, LABEL* y_direct,
                 SWORD* fy_to, SWORD* fy_away, double *dfy)
{
  const double dfy_weight = 1.0;
  double total_loss = 0; // cumulative loss for all examples
  double total_dscore = 0;

  if(gparm->use_history) {

    // use history of constraints
    ConstraintSet* cs = ConstraintSet::Instance();

    for(int il = 0; il < nExamples; il++) { /*** example loop ***/
      const vector< constraint >* constraints = cs->getConstraints(ex[il].x.id);
      if(constraints) {
        for(vector<constraint>::const_iterator it = constraints->begin();
            it != constraints->end(); ++it) {
          double dscore_cs = compute_gradient(sm, gparm, fy_to,
                                              it->first->w, dfy, it->first->loss, dfy_weight);

          // do not add if negative to avoid adding and subtracting values.
          // this score is just logged, not used in any computation.
          if(dscore_cs > 0) {
            total_dscore += dscore_cs;
            total_loss += it->first->loss;
          }
        }
      }
    }
  } else {
    for(int il = 0; il < nExamples; il++) { /*** example loop ***/
      double _loss  = 0;
      total_dscore += compute_gradient(sparm, sm, &ex[il], &y_bar[il],
                                       &y_direct[il], gparm, fy_to, fy_away,
                                       dfy, &_loss, dfy_weight);
      total_loss += _loss;
    }
  }

#if CUSTOM_VERBOSITY > 1
  // loss computed after updating weight vector
  ofstream ofs("loss.txt", ios::app);
  ofs << total_loss << endl;
  ofs.close();

  // dscore computed after updating weight vector
  ofstream ofs_dscore("a_dscore.txt", ios::app);
  ofs_dscore << total_dscore << endl;
  ofs_dscore.close();
#endif

  double m = total_dscore + total_loss;
  return m;
}

/**
 * compute an estimate of the score
 * return max among a subset of sampled superpixels and only compute score for class 0
 */
double compute_score_estimate(STRUCT_LEARN_PARM *sparm, STRUCTMODEL *sm, GRADIENT_PARM& gparm, EXAMPLE* example)
{
  const int c = 0; // class 0
  labelType* groundTruthLabels = example->y.nodeLabels;
  EnergyParam param;
  sparmToEnergyParam(*sparm, &param);
  Slice_P* slice = example->x.slice;
  GraphInference gi(slice, (const EnergyParam*)&param, sm->w+1, example->x.feature, 0, 0);

  double max_potential = 0;

  string config_tmp;
  double sampling_rate = 0.3;
  if(Config::Instance()->getParameter("sampling_rate", config_tmp)) {
    sampling_rate = atof(config_tmp.c_str()); 
    printf("[svm_struct] sampling_rate = %g\n", sampling_rate);
  }
  bool draw_samples = sampling_rate != 1.0;
  ulong nSupernodes = slice->getNbSupernodes();
  int sid = 0;
  for(int i = 0; i < nSupernodes*sampling_rate; ++i) {

    if(draw_samples) {
      // Select a pixel at random
      sid = rand() * ((double)nSupernodes/(double)RAND_MAX);
    } else {
      sid = i;
    }

    supernode* s = slice->getSupernode(sid);
    vector < supernode* >* lNeighbors = &(s->neighbors);

    double potential = gi.computeUnaryPotential(slice, sid, c);
    for(vector < supernode* >::iterator itN = lNeighbors->begin();
        itN != lNeighbors->end(); itN++) {
      const int c2 = rand()*sparm->nClasses / (double)RAND_MAX;
      double pairwisePotential = gi.computePairwisePotential(slice, s, (*itN),
                                                             c, c2);
      potential += pairwisePotential;
    }

    if(!gparm.ignore_loss) {
      if(c != groundTruthLabels[sid]) {
        // add loss of the ground truth label
        potential += sparm->lossPerLabel[groundTruthLabels[sid]];
      }
    }

    if(potential > max_potential) {
      max_potential = potential;
    }
  }

  return max_potential;
}

void set_sampling_temperature(STRUCT_LEARN_PARM *sparm,
                              STRUCTMODEL *sm,
                              GRADIENT_PARM& gparm,
                              EXAMPLE* ex,
                              long nExamples,
                              double* dscores,
                              double* ddfy)
{

  int schedulingType = 2;
  string config_tmp;
  if(Config::Instance()->getParameter("schedulingType", config_tmp)) {
    schedulingType = atoi(config_tmp.c_str()); 
  }
  printf("[svm_struct_custom] schedulingType=%d\n", schedulingType);

  int nItems = sparm->stepForOutputFiles;

  switch(schedulingType) {
  case 0:
    {
    // try to lower the temperature as much as we can

    double score = compute_score_estimate(sparm, sm, gparm, &ex[0]);

    // compute average dscore
    double avg_dscore = 0;
    for(int i = 0; i < nItems; ++i) {
      avg_dscore += dscores[i];
    }

#if CUSTOM_VERBOSITY > 2
    {
      ofstream ofs_temp("avg_dscore.txt", ios::app);
      ofs_temp << sparm->iterationId << " " << avg_dscore << endl;
      ofs_temp.close();
    }
#endif

    bool exp_is_inf = false;
    while(!exp_is_inf) {
      exp_is_inf = isinf(exp(score/(sparm->sampling_temperature_0/SAMPLING_MUL_COEFF)));
      if(!exp_is_inf) {
        // decrease randomness
        sparm->sampling_temperature_0 /= SAMPLING_MUL_COEFF;
      }
    }
    }
    break;
  case 1:
    {
    // decrease temperature if dscore < 0
    // increase temperature if not enough changes

    double score = compute_score_estimate(sparm, sm, gparm, &ex[0]);
        
    // compute average dscore
    double avg_dscore = 0;
    for(int i = 0; i < nItems; ++i) {
      avg_dscore += dscores[i];
    }

#if CUSTOM_VERBOSITY > 2
    {
      ofstream ofs_temp("avg_dscore.txt", ios::app);
      ofs_temp << sparm->iterationId << " " << avg_dscore << endl;
      ofs_temp.close();
    }
#endif

    if(avg_dscore < 0) {
      bool exp_is_inf = isinf(exp(score/(sparm->sampling_temperature_0/SAMPLING_MUL_COEFF)));
      if(!exp_is_inf && sparm->sampling_temperature_0 > score*1e-10) {
        // decrease randomness
        sparm->sampling_temperature_0 /= SAMPLING_MUL_COEFF;

        // reset labels in the cache to ground-truth
        for(int i = 0; i < nExamples; i++) {
          int cacheId = ex[i].x.id;
          LabelCache::Instance()->setLabel(cacheId, ex[i].y);
        }

#if CUSTOM_VERBOSITY > 2
        ofstream ofs_temp("temperature_change.txt", ios::app);
        ofs_temp << sparm->iterationId << " " << score << " ";
        ofs_temp << score/sparm->sampling_temperature_0 << " " << exp(score/(sparm->sampling_temperature_0));
        ofs_temp << " " << sparm->sampling_temperature_0 << endl;
        ofs_temp.close();
#endif
      } else {

        // decrease sampling rate
        sparm->sampling_rate = max(0.1, sparm->sampling_rate-0.1);

#if CUSTOM_VERBOSITY > 2
        ofstream ofs_temp("sampling_rate.txt", ios::app);
        ofs_temp << sparm->sampling_rate << endl;
        ofs_temp.close();
#endif

      }
    }

    }
    break;
  case 2:
    {
      // decrease temperature
      if((sparm->iterationId % 20) == 0) {
      //if((sparm->iterationId % sparm->stepForOutputFiles) == 0) {
        double score = compute_score_estimate(sparm, sm, gparm, &ex[0]);
        bool exp_is_inf = isinf(exp(score/(sparm->sampling_temperature_0/SAMPLING_MUL_COEFF)));
        printf("[svm_learn] Set temp %g %d\n", score, (int)exp_is_inf);
        if(!exp_is_inf) {
          // decrease randomness
          sparm->sampling_temperature_0 /= SAMPLING_MUL_COEFF;
        }
      }
    }
    break;
  default:
    // do not change the temperature
    break;
  }
}

// use -w 9 to call this function
void svm_learn_struct_joint_custom(SAMPLE sample, STRUCT_LEARN_PARM *sparm,
				   LEARN_PARM *lparm, KERNEL_PARM *kparm, 
				   STRUCTMODEL *sm)
     /* Input: sample (training examples)
	       sparm (structural learning parameters)
               lparm (svm learning parameters)
               kparm (kernel parameters)
	       Output: sm (learned model) */
{
  long        nTotalExamples = sample.n;
  EXAMPLE     *examples = sample.examples;

  CONSTSET    cset;
  int         cached_constraint = 0;
  int         numIt = 0;
  double      *alpha = NULL;
  double ceps;
  bool finalized = false;
  string config_tmp;

  init_struct_model(sample, sm, sparm, lparm, kparm); 

  Config* config = Config::Instance();
  GRADIENT_PARM gparm;
  init_gradient_param(gparm, config, ConstraintSet::Instance());
  gparm.examples_all = examples;
  gparm.n_total_examples = nTotalExamples;

  // 0 = initialize parameter vector to 0
  // 1 = initialize parameter vector using random values
  int init_type = INIT_WEIGHT_0;
  if(gparm.ignore_loss == 1) {
    init_type = INIT_WEIGHT_RANDOM;
  }
  if(config->getParameter("sgd_init_type", config_tmp)) {
    init_type = atoi(config_tmp.c_str());
  }
  printf("[SVM_struct_custom] init_type = %d\n", init_type);

  sm->svm_model = 0;
  //sm->w = new double[sm->sizePsi+1];
  // use C style to be compatible with svm-light
  // and to make sure there is no error in free_struct_model
  sm->w = (double *)my_malloc(sizeof(double)*(sm->sizePsi+1));
  init_w(sparm, sm, &gparm, examples, init_type);

  if(sparm->giType == T_GI_SAMPLING) {
    init_sampling(sparm, sm, examples, nTotalExamples);
  }

  int n_iterations_update_learning_rate = 1000;
  if(config->getParameter("n_iterations_update_learning_rate", config_tmp)) {
    n_iterations_update_learning_rate = atoi(config_tmp.c_str());
  }
  printf("[SVM_struct_custom] n_iterations_update_learning_rate = %d\n", n_iterations_update_learning_rate);

  if(gparm.ignore_loss) {
    sparm->lossPerLabel = 0;
  }

  double* momentum = 0;
  if( (gparm.update_type == UPDATE_MOMENTUM) || (gparm.update_type == UPDATE_MOMENTUM_DECREASING)) {
    int _sizePsi = sm->sizePsi + 1;
    momentum = new double[_sizePsi];
    for(int i = 0; i < _sizePsi; ++i) {
      momentum[i] = 0;
    }
  }

  double last_obj = 0;
  int _sizePsi = sm->sizePsi + 1;
  SWORD* fy_to = new SWORD[_sizePsi];
  SWORD* fy_away = new SWORD[_sizePsi];
  double* dfy = new double[_sizePsi];
  memset((void*)dfy, 0, sizeof(double)*(_sizePsi));
  LABEL* y_bar = new LABEL[nTotalExamples];

  int nMaxIterations = (sparm->nMaxIterations!=-1)?sparm->nMaxIterations:1e7;
  int nItems = sparm->stepForOutputFiles;
  double* learning_rates = new double[nItems];
  double* norm_ws = new double[nItems];
  double* ms = new double[nItems];
  double* objs = new double[nItems];
  double* dscores = new double[nItems];
  int idx = 0;
  int example_id = 0;

  // sampling: keep previous dfy and check how much change
  // occurs to decide if we should raise the temperature
  double* dfy_p = 0;
  double* ddfy = 0;
  if(sparm->giType == T_GI_SAMPLING) {
    dfy_p = new double[_sizePsi];
    memset((void*)dfy_p, 0, sizeof(double)*(_sizePsi));
    ddfy = new double[nItems];
  }

  do {

    EXAMPLE* _ex = examples;

    int _nBatchExamples = nTotalExamples;

    if(gparm.n_batch_examples != -1) {
      _ex += example_id;

      // make sure next batch doesn't go over size of training set
      if(example_id + gparm.n_batch_examples > nTotalExamples) {
        _nBatchExamples = nTotalExamples - example_id;
        printf("[SVM_struct_custom] Batch size for iteration %d is too large. Reducing to %d\n", sparm->iterationId, _nBatchExamples);
      } else {
        _nBatchExamples = gparm.n_batch_examples;
      }
    }

    double m = do_gradient_step(sparm, sm, _ex, _nBatchExamples,
                                &gparm, momentum, fy_to, fy_away, dfy,
                                dscores[idx], y_bar);

    // projection
    if(gparm.enforce_submodularity) {
      double* smw = sm->w + 1;
      enforce_submodularity(sparm, smw);
    }
    if(gparm.max_norm_w > 0) {
      project_w(sm, &gparm);
    }

    double obj = 0;
    switch(gparm.loss_type) {
    case LINEAR_LOSS:
      obj = m;
      break;
    case LOG_LOSS:
      if(m > 100) {
        obj = log(1 + std::exp(100.0));
      } else {
        obj = log(1 + std::exp(m));
      }
      break;
    case HINGE_LOSS:
      obj = max(0.0, m);
      break;
    case SQUARE_HINGE_LOSS:
      obj = max(0.0, m*m);
      break;
    default:
      printf("[svm_struct_custom] Unknown loss type %d\n", gparm.loss_type);
      exit(-1);
      break;
    }

    if(numIt == 0) {
      ceps = fabs(obj);
    } else {
      ceps = fabs(last_obj - obj);
    }
    last_obj = obj;

    finalized = finalize_iteration(ceps,cached_constraint,sample,sm,cset,alpha,sparm);

    switch(gparm.update_type)
      {
      case UPDATE_DECREASING:
      case UPDATE_MOMENTUM_DECREASING:
        {
          if(numIt > n_iterations_update_learning_rate) {
            gparm.learning_rate = gparm.learning_rate_0/(double)pow(numIt-n_iterations_update_learning_rate,
                                                                    gparm.learning_rate_exponent);
          } else {
            if(numIt > n_iterations_update_learning_rate) {
              gparm.learning_rate = gparm.learning_rate_0/(double)pow(numIt-n_iterations_update_learning_rate,
                                                                      gparm.learning_rate_exponent);
            }
          }
        }
        break;
      default:
        break;
      }

    int _sizePsi = sm->sizePsi + 1;
    double norm_w = 0;
    for(int i = 1; i < _sizePsi; ++i) {
      norm_w += sm->w[i]*sm->w[i];
    }

    norm_ws[idx] = sqrt(norm_w);
    ms[idx] = m;
    objs[idx] = obj;
    learning_rates[idx] = gparm.learning_rate;

    if(ddfy) {
      // compute difference between successive dfy vectors
      ddfy[idx] = 0;
      for(int i = 0; i < _sizePsi; ++i) {
        ddfy[idx] += fabs(dfy_p[i] - dfy[i]);
        dfy_p[i] = dfy[i];
      }
    }

    // dump stats
    if(idx == nItems-1) {
      idx = 0;
      write_vector("learning_rate.txt", learning_rates, nItems);
      write_vector("m.txt", ms, nItems);
      write_vector("norm_w.txt", norm_ws, nItems);
      write_vector("obj.txt", objs, nItems);

      if(momentum) {
        stringstream sout;
        sout << "step_" << nItems << ".txt";
        string s_it = sout.str();
        string momentum_filename = "momentum_" + s_it;
        write_vector(momentum_filename.c_str(), momentum, _sizePsi);
      }

    } else {
      ++idx;
    }

    // dynamically change temperature for sampling
    if(sparm->giType == T_GI_SAMPLING) {
      set_sampling_temperature(sparm, sm, gparm, examples, gparm.n_total_examples, dscores, ddfy);
    }

    fflush(stdout);

    example_id += gparm.n_batch_examples;
    if(example_id >= gparm.n_total_examples) {
      example_id = 0;
    }

    ++numIt;

  } while( (numIt < nMaxIterations) &&
           (!finalized)
	 );

  ConstraintSet::Instance()->save("constraint_set.txt");

  if(numIt >= nMaxIterations) {
    printf("[svm_struct_custom] Reached max number of iterations %d\n",
           nMaxIterations);
  }
  if (ceps <= sparm->epsilon) {
    printf("[svm_struct_custom] Reached epsilon value %g/%g\n",
           ceps, sparm->epsilon);
  }

  if(momentum) {
    delete[] momentum;
    momentum = 0;
  }

  delete[] fy_to;
  delete[] fy_away;
  delete[] dfy;
  delete[] objs;
  delete[] ms;
  delete[] learning_rates;
  delete[] norm_ws;
  delete[] dscores;

  for(int s = 0; s < nTotalExamples; ++s) {
    free_label(y_bar[s]);
  }
  delete[] y_bar;

  if(dfy_p) {
    delete[] dfy_p;
  }
  if(ddfy) {
    delete[] ddfy;
  }
}

void init_sampling(STRUCT_LEARN_PARM *sparm, STRUCTMODEL *sm, EXAMPLE* examples, long nExamples)
{
  string config_tmp;
  // set temperature for sampling
  double norm_wy_to = get_norm_psi_gt(sparm, sm, examples, nExamples);

  sparm->sampling_temperature_0 = norm_wy_to*1e4;
  printf("[SVM_struct_custom] sampling_temperature_0 = %g\n", sparm->sampling_temperature_0);

  sparm->sampling_rate = 0.3;
  if(Config::Instance()->getParameter("sampling_rate", config_tmp)) {
    sparm->sampling_rate = atof(config_tmp.c_str()); 
  }
  printf("[SVM_struct_custom] sampling_rate = %g\n", sparm->sampling_rate);

  //sparm->lossScale = norm_wy_to;
  sparm->lossScale = 1;
  if(Config::Instance()->getParameter("sgd_loss_scale", config_tmp)) {
    double lossScale = atof(config_tmp.c_str());
    if(lossScale > 0) {
      sparm->lossScale = lossScale;
    }
  }
  printf("[SVM_struct_custom] sparm->lossScale = %g\n", sparm->lossScale);

  for(int c = 0; c < sparm->nClasses; ++c) {
    sparm->lossPerLabel[c] *= sparm->lossScale;
  }
}

void init_w(STRUCT_LEARN_PARM *sparm, STRUCTMODEL *sm, GRADIENT_PARM* gparm, EXAMPLE* examples, int init_type)
{
  vector<string> paramFiles;
  getFilesInDir("parameter_vector0/", paramFiles, "txt", true);

  if(paramFiles.size() > 0) {
    string initial_parameter_vector = paramFiles[paramFiles.size()-1].c_str();
    printf("[SVM_struct_custom] Loading parameter file %s\n",
           initial_parameter_vector.c_str());
    sm->w[0] = 0;
    EnergyParam _param(initial_parameter_vector.c_str());
    for(int i = 1; i < sm->sizePsi+1; ++i) {
      sm->w[i] = _param.weights[i-1];
    }
    rename(initial_parameter_vector.c_str(), "initial_parameter_vector.txt");
  } else {

    if(init_type == 0) {
      // easier to debug cause we know exactly what to expect for the first
      // most violated constraint.
      printf("[SVM_struct_custom] Initializing weight vector to 0\n");
      for(int i = 0; i < sm->sizePsi+1; ++i) {
        sm->w[i] = 0;
      }
    } else {
      printf("[SVM_struct_custom] Initializing weight vector randomly\n");
      srand(time(NULL));
      double norm_w = 0;
      sm->w[0] = 0;
      for(int i = 1; i < sm->sizePsi+1; ++i) {
        sm->w[i] = rand() / (double)RAND_MAX;
        norm_w += sm->w[i];
      }

      int _sizePsi = sm->sizePsi + 1;
      // todo ; loop over all examples ?
      SWORD* fy_to = new SWORD[_sizePsi];
      computePsi(fy_to, examples[0].x,examples[0].y,sm,sparm);
      double norm_wy_to = 0;
      SWORD* wy_to = fy_to;
      while (wy_to->wnum) {
        norm_wy_to += wy_to->weight*wy_to->weight;
        ++wy_to;
      }
      delete[] fy_to;

      if(gparm->max_norm_w > 0 && norm_w > 0) {
        double scale = gparm->max_norm_w/norm_w;
        printf("[SVM_struct_custom] Rescaling w by %g\n", scale);
        for(int i = 1; i < sm->sizePsi+1; ++i) {
          sm->w[i] *= scale;
        }
      }
    }

    stringstream sout;
    sout << "parameter_vector0";
    sout << "/iteration_0.txt";
    save_parameters(sout.str().c_str(), sparm, sm);
  }
}

void project_w(STRUCTMODEL *sm, GRADIENT_PARM* gparm)
{
  double norm_w = 0;
  for(int i = 1; i < sm->sizePsi+1; ++i) {
    norm_w += sm->w[i];
  }

  if(gparm->max_norm_w > 0 && norm_w > gparm->max_norm_w) {
    double scale = gparm->max_norm_w/norm_w;
    printf("[SVM_struct_custom] Rescaling w by %g/%g = %g\n",
           gparm->max_norm_w, norm_w, scale);
    for(int i = 1; i < sm->sizePsi+1; ++i) {
      sm->w[i] *= scale;
    }
  }
}

void enforce_submodularity(STRUCT_LEARN_PARM *sparm, double* smw)
{
  double* pw = smw + sparm->nUnaryWeights;

  int nloops = 1;
  int offset = 0;
  if(sparm->nClasses > 2) {
    nloops = 2;
    offset = 1;
  }

  // For 2 classes, use graph-cuts if pairwise potential is attractive/submodular.
  // Let E be the energy and S the score (E=-S)
  // Recall that potential is submodular if :
  // E(0,0) + E(1,1) =< E(0,1) + E(1,0)
  // S(0,0) + S(1,1) >= S(0,1) + S(1,0) or d >= a

  for(int l = 0 ; l < nloops; ++l) {

    double d = 0; // diagonal
    double a = 0; //anti-diagonal
    double d2 = 0; // diagonal
    double a2 = 0; //anti-diagonal

#if USE_LONG_RANGE_EDGES
    for(int p = 0; submodularEnergy && (p < sparm->nDistances); ++p)
      {
        a = 0;
        d = 0;
#else
        {
#endif

          for(int g = 0; g < sparm->nGradientLevels; ++g) {
            d += pw[0] + pw[3+offset];
            a += pw[1] + pw[2+offset];
            if(a > d) {
              // compute accumulated terms up to gradient - 1
              a2 = 0;
              d2 = 0;
              double* pw2 = smw + sparm->nUnaryWeights;
              if(l > 0) {
                pw2 += 4;
              }
              for(int g2 = 0; g2 < g; ++g2) {
                d2 += pw2[0] + pw2[3+offset];
                a2 += pw2[1] + pw2[2+offset];
                pw2 += sparm->nClasses*sparm->nClasses;
              }

              // make sure that a < d (remember that submodular means that the sum of the terms on the diagonal is higher)
              pw[1] = (d-a2)/2.0 - 1e-2;
              pw[2+offset] = pw[1];

              // add last term
              d2 += pw2[0] + pw2[3+offset];
              a2 += pw2[1] + pw2[2+offset];
          
              a = a2;
              d = d2;
            } else {
              // a <= d
              // make sure matrix is symmetric
              pw[1] = (pw[1] + pw[2+offset])/2.0;
              pw[2+offset] = pw[1];
            }

            pw += sparm->nClasses*sparm->nClasses;
          }
        }

        pw = smw + sparm->nUnaryWeights;
        pw += 4;
      }
}

void init_gradient_param(GRADIENT_PARM& gparm, Config* config,
                         ConstraintSet* constraint_set)
{
  string config_tmp;

  double learning_rate = 1e-9;
  if(config->getParameter("learning_rate", config_tmp)) {
    learning_rate = atof(config_tmp.c_str());
  }
  printf("[SVM_struct_custom] learning_rate = %g\n", learning_rate);

  eUpdateType sgd_update_type = UPDATE_MOMENTUM_DECREASING;
  if(config->getParameter("sgd_update_type", config_tmp)) {
    sgd_update_type = (eUpdateType)atoi(config_tmp.c_str());
  }
  printf("[SVM_struct_custom] sgd_update_type = %d\n", (int)sgd_update_type);

  float momentum_weight = 0;
  if(config->getParameter("sgd_momentum_weight", config_tmp)) {
    momentum_weight = atof(config_tmp.c_str());
  }
  printf("[SVM_struct_custom] momentum_weight = %g\n", momentum_weight);
  assert(momentum_weight <= 1.0);

  float regularization_weight = 0;
  if(config->getParameter("sgd_regularization_weight", config_tmp)) {
    regularization_weight = atof(config_tmp.c_str());
  }
  printf("[SVM_struct_custom] regularization_weight = %g\n", regularization_weight);

  double learning_rate_exponent = 0.5;
  if(config->getParameter("learning_rate_exponent", config_tmp)) {
    learning_rate_exponent = atof(config_tmp.c_str());
  }
  printf("[SVM_struct_custom] learning_rate_exponent = %g\n", learning_rate_exponent);

  bool update_loss_function = false;
  if(config->getParameter("update_loss_function", config_tmp)) {
    update_loss_function = config_tmp.c_str()[0] == '1';
  }
  printf("[SVM_struct_custom] update_loss_function=%d\n", (int)update_loss_function);

  eLossType sgd_loss_type = HINGE_LOSS;
  if(config->getParameter("sgd_loss_type", config_tmp)) {
    sgd_loss_type = (eLossType)atoi(config_tmp.c_str());
  }
  printf("[SVM_struct_custom] sgd_loss_type = %d\n", (int)sgd_loss_type);
  assert(sgd_loss_type != LINEAR_LOSS);

  eGradientType sgd_gradient_type = GRADIENT_GT;
  if(config->getParameter("sgd_gradient_type", config_tmp)) {
    sgd_gradient_type = (eGradientType)atoi(config_tmp.c_str());
  }
  printf("[SVM_struct_custom] sgd_gradient_type = %d\n", (int)sgd_gradient_type);

  int sgd_n_batch_examples = -1; // use all examples at once
  if(config->getParameter("sgd_n_batch_examples", config_tmp)) {
    sgd_n_batch_examples = atoi(config_tmp.c_str());
  }
  printf("[SVM_struct_custom] sgd_n_batch_examples = %d\n", sgd_n_batch_examples);

  int max_number_constraints = 10000;
  bool sgd_use_history = true;
  if(config->getParameter("cs_max_number_constraints", config_tmp)) {
    max_number_constraints = atoi(config_tmp.c_str());
    sgd_use_history = max_number_constraints > 0;
  }
  printf("[SVM_struct_custom] sgd_use_history = %d\n", (int)sgd_use_history);

  eSortingType constraint_set_type = CS_DISTANCE;
  if(config->getParameter("constraint_set_type", config_tmp)) {
    constraint_set_type = (eSortingType)atoi(config_tmp.c_str());
  }
  printf("[SVM_struct_custom] constraint_set_type = %d\n", (int)constraint_set_type);
  constraint_set->setSortingAlgorithm(constraint_set_type);

  bool sgd_ignore_loss = false;
  if(config->getParameter("sgd_ignore_loss", config_tmp)) {
    sgd_ignore_loss = config_tmp.c_str()[0] == '1';
  }
  printf("[SVM_struct_custom] sgd_ignore_loss = %d\n", (int)sgd_ignore_loss);

  double max_norm_w = -1;
  if(Config::Instance()->getParameter("max_norm_w", config_tmp)) {
    max_norm_w = atof(config_tmp.c_str());
  }
  printf("[SVM_struct_custom] max_norm_w = %g\n", max_norm_w);

  bool enforce_submodularity = true;
  if(Config::Instance()->getParameter("project_w", config_tmp)) {
    enforce_submodularity = config_tmp.c_str()[0] == '1';
  }
  if(Config::Instance()->getParameter("enforce_submodularity", config_tmp)) {
    enforce_submodularity = config_tmp.c_str()[0] == '1';
  }
  printf("[SVM_struct_custom] enforce_submodularity = %d\n", (int)enforce_submodularity);

  gparm.learning_rate = learning_rate;
  gparm.learning_rate_0 = learning_rate; // initial learning rate
  gparm.learning_rate_exponent = learning_rate_exponent;
  gparm.max_norm_w = max_norm_w;
  gparm.enforce_submodularity = enforce_submodularity;
  gparm.update_type = sgd_update_type;
  gparm.loss_type = sgd_loss_type;
  gparm.momentum_weight = momentum_weight;
  gparm.regularization_weight = regularization_weight;
  gparm.gradient_type = sgd_gradient_type;
  gparm.use_history = sgd_use_history;
  gparm.constraint_set_type = constraint_set_type;
  gparm.ignore_loss = sgd_ignore_loss;
  gparm.n_batch_examples = sgd_n_batch_examples;
}
