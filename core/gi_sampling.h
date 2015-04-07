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

#ifndef GI_SAMPLING_H
#define GI_SAMPLING_H

// SliceMe
#include "Slice.h"

#include "graphInference.h"
#include "energyParam.h"

#include <map>
#include <vector>

//------------------------------------------------------------------------------

template<template <typename> class P = std::less >
struct compare_pair_second {
    template<class T1, class T2> bool operator()(const std::pair<T1, T2>& left, const std::pair<T1, T2>& right) {
        return P<T2>()(left.second, right.second);
    }
};

typedef std::pair<int, double> prob_pair;

//------------------------------------------------------------------------------

class GI_sampling : public GraphInference
{
 public:
  /**
   * Constructor for SSVM framework
   */
  GI_sampling(Slice_P* _slice,
              const EnergyParam* _param,
              double* _smw,
              labelType* _groundTruthLabels,
              double* _lossPerLabel,
              Feature* _feature,
              std::map<sidType, nodeCoeffType>* _nodeCoeffs,
              double _sampling_rate);

  /**
   * This function can start running multiple chains in parallel
   * You have to initialize the seed for the time function before calling the run function
   * For example :
   * #if USE_GSL_DEBUG
   *  gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
   * #else
   *  srand(time(NULL));
   * #endif
   */
  double run(labelType* inferredLabels,
             int id,
             size_t maxiter,
             labelType* nodeLabelsGroundTruth = 0,
             bool computeEnergyAtEachIteration = false,
             double* _loss = 0);

  double run_VOC(labelType* inferredLabels,
                 int id,
                 size_t maxiter,
                 labelType* nodeLabelsGroundTruth,
                 ulong* TPs, ulong* FPs, ulong* FNs);

  /**
   * Run a single chain
   */
  double runOnce(labelType* inferredLabels,
                 size_t maxiter,
                 labelType* nodeLabelsGroundTruth,
                 double* _loss,
                 double temperature);

  /**
   * Run a single chain
   * Optimize VOC score
   */
  double runOnce_VOC(labelType* inferredLabels,
                     size_t maxiter,
                     labelType* nodeLabelsGroundTruth,
                     double temperature,
                     ulong* TPs, ulong* FPs, ulong* FNs);

  void setInitializedLabels(bool value) { initializedLabels = value; }

 private:
  bool initializedLabels;
  double sampling_rate;

  bool replaceVoidMSRC;
  labelType voidLabel;
  labelType moutainLabel;
  labelType horseLabel;

};

#endif // GI_SAMPLING_H
