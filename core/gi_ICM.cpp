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

#include "gi_ICM.h"

// SliceMe
#include "Config.h"
#include "utils.h"

#include "inference_globals.h"

#if USE_GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif

//------------------------------------------------------------------------------

GI_ICM::GI_ICM(Slice_P* _slice, 
               const EnergyParam* _param,
               double* _smw,
               labelType* _groundTruthLabels,
               double* _lossPerLabel,
               Feature* _feature,
               map<sidType, nodeCoeffType>* _nodeCoeffs)
{
  GraphInference::init();
  slice = _slice;
  param = _param;
  smw = _smw;
  lossPerLabel = _lossPerLabel;
  groundTruthLabels = _groundTruthLabels;
  lossPerLabel = _lossPerLabel;
  feature = _feature;
  nodeCoeffs = _nodeCoeffs;
}


double GI_ICM::run(labelType* inferredLabels,
                   int id,
                   size_t maxiter,
                   labelType* nodeLabelsGroundTruth,
                   bool computeEnergyAtEachIteration,
                   double* _loss)
{
  bool useLossFunction = lossPerLabel!=0;
  string paramMSRC;
  Config::Instance()->getParameter("msrc", paramMSRC);
  bool useMSRC = paramMSRC.c_str()[0] == '1';
  bool replaceVoidMSRC = false;
  labelType voidLabel = 0;
  labelType moutainLabel = 0;
  labelType horseLabel = 0;
  if(!useLossFunction && useMSRC) {
    Config::Instance()->getParameter("msrc_replace_void", paramMSRC);
    replaceVoidMSRC = paramMSRC.c_str()[0] == '1';
    voidLabel = classIdxToLabel[0];
    moutainLabel = classIdxToLabel[4161600];
    horseLabel = classIdxToLabel[8323328];
    printf("[GI_ICM] MSRC void=%d, moutain=%d, horse=%d\n",
           (int)voidLabel, (int)moutainLabel, (int)horseLabel);
  } else {
    printf("[GI_ICM] Do not replace void labels\n");
  }

  double *buf = new double[param->nClasses];
  int sid = 0;
  double maxScore = 0;
  double totalScore_old = 0;
  double totalScore = 10;
  supernode *s;
  node cs;

  // allocate memory to store features
  int fvSize = feature->getSizeFeatureVector();
  int max_index = fvSize + 1;
  osvm_node* n = new osvm_node[max_index];
  int i = 0;
  for(i = 0;i < max_index-1; i++)
    n[i].index = i+1;
  n[i].index = -1;

  const map<int, supernode* >& _supernodes = slice->getSupernodes();

  // random selection
  int maxIter = 10;
  for(int iter = 0; iter < maxIter && (totalScore - totalScore_old) > 1.0; ++iter) {
    
    printf("[GI_ICM] Iteration %d/%d\n", iter, maxIter);

    totalScore_old = totalScore;
    totalScore = 0;

    // Select a pixel at random
    //int sid = gsl_rng_uniform_int(r, nSupernodes);

    int nLabelsChanged = 0;

    for(map<int, supernode* >::const_iterator its = _supernodes.begin();
        its != _supernodes.end(); its++) {
      sid = its->first;
      s = its->second;

      if(param->nOrientations > 1) {
        s->getCenter(cs);
      }

      if(param->nClasses != 2) {

        for(int c = 0; c < (int)param->nClasses; c++) {
          double unaryPotential = computeUnaryPotential(slice, sid, i);
          buf[c] = unaryPotential;

          // add pairwise potential
          if(param->includeLocalEdges) {
            vector < supernode* >* lNeighbors = &(s->neighbors);
            for(vector < supernode* >::iterator itN = lNeighbors->begin();
                itN != lNeighbors->end(); itN++) {
              double pairwisePotential = computePairwisePotential(slice, s, (*itN),
                                                               c, inferredLabels[(*itN)->id]);
              buf[c] += pairwisePotential;
            }
          }

        }
      } else {
        buf[T_FOREGROUND] = 0;

        const int c = T_BACKGROUND;
        double unaryPotential = computeUnaryPotential(slice, sid, c);
        buf[c] = unaryPotential;

        // add pairwise potential
        if(param->includeLocalEdges) {
          vector < supernode* >* lNeighbors = &(s->neighbors);
          for(vector < supernode* >::iterator itN = lNeighbors->begin();
              itN != lNeighbors->end(); itN++) {
            double pairwisePotential = computePairwisePotential(slice, s, (*itN),
                                                             c, inferredLabels[(*itN)->id]);

            buf[c] += pairwisePotential;
          }
        }

      }

      if(useLossFunction) {
        // loss function
        for(int s = 0; s < (int)param->nClasses; s++) {
          if(s != groundTruthLabels[sid]) {
            // add loss of the ground truth label
            buf[s] = buf[s] + lossPerLabel[groundTruthLabels[sid]];
          }
        }

        // loss function : -1 if correct label (same as adding +1 to all incorrect labels)
        //buf[groundTruthLabels[sid]] = std::exp(std::log(buf[groundTruthLabels[sid]])-1);
      }

      // pick max
      //bool initialized = false;
      maxScore = buf[inferredLabels[sid]];
      for(int s = 0; s < (int)param->nClasses; s++) {
        if(!replaceVoidMSRC || (s != voidLabel && s != moutainLabel && s != horseLabel)) {
          //if(!initialized || maxScore < buf[s]) {
          if(maxScore < buf[s]) {
            printf("[GI_ICM] Changing sid %d label %d - > %d score %g -> %g\n", sid, inferredLabels[sid], s, maxScore, buf[s]);
            computeEnergy(inferredLabels);
            maxScore = buf[s];
            inferredLabels[sid] = s;
            printf("[GI_ICM] Iteration %d/%d %g\n", iter, maxIter, computeEnergy(inferredLabels));
            ++nLabelsChanged;
            //initialized = true;
          }
        }
      }
      totalScore += maxScore;
    }
    printf("[GI_ICM] Iteration %d/%d. Total score = %g. nLabelsChanged = %d\n", iter, maxIter,
           totalScore, nLabelsChanged);
  }

  // cleaning  
  delete[] n;
  delete[] buf;
  return computeEnergy(inferredLabels);
}
