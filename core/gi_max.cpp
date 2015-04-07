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

#include "gi_max.h"

// SliceMe
#include "Config.h"
#include "utils.h"

#include "inference_globals.h"

//------------------------------------------------------------------------------

GI_max::GI_max(Slice_P* _slice, 
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
  feature = _feature;
  nodeCoeffs = _nodeCoeffs;

  if(param->includeLocalEdges == 1) {
    printf("[GI_max] Warning : param->includeLocalEdges=%d\n", param->includeLocalEdges);
  }
}


double GI_max::run(labelType* inferredLabels,
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
    printf("[gi_max] MSRC void=%d, moutain=%d, horse=%d\n",
           (int)voidLabel, (int)moutainLabel, (int)horseLabel);
  } else {
    printf("[gi_max] Do not replace void labels\n");
  }

  double *buf = new double[param->nClasses];
  int sid = 0;
  double maxScore = 0;

  // allocate memory to store features
  int fvSize = feature->getSizeFeatureVector();
  int max_index = fvSize + 1;
  osvm_node* n = new osvm_node[max_index];
  int i = 0;
  for(i = 0;i < max_index-1; i++)
    n[i].index = i+1;
  n[i].index = -1;

  const map<int, supernode* >& _supernodes = slice->getSupernodes();
  for(map<int, supernode* >::const_iterator its = _supernodes.begin();
      its != _supernodes.end(); its++) {
    sid = its->first;

    if(param->nClasses != 2) {
      for(int i = 0; i < (int)param->nClasses; i++) {

        double unaryPotential = computeUnaryPotential(slice, sid, i);

        if(nodeCoeffs) {
          unaryPotential *= (*nodeCoeffs)[sid];
        }

        buf[i] = unaryPotential;
      }
    } else {
      buf[T_FOREGROUND] = 0;

      const int i = T_BACKGROUND;
      double unaryPotential = computeUnaryPotential(slice, sid, i);

      if(nodeCoeffs) {
        unaryPotential *= (*nodeCoeffs)[sid];
      }

      buf[i] = unaryPotential;
    }

    if(useLossFunction) {
      // loss function
      for(int s = 0; s < (int)param->nClasses; s++) {
        if(s != groundTruthLabels[sid]) {
          // add loss of the ground truth label
          if(nodeCoeffs) {
            buf[s] += (*nodeCoeffs)[sid]*lossPerLabel[groundTruthLabels[sid]];
          } else {
            buf[s] += lossPerLabel[groundTruthLabels[sid]];
          }
        }
      }

      // loss function : -1 if correct label (same as adding +1 to all incorrect labels)
      //buf[groundTruthLabels[sid]] = std::exp(std::log(buf[groundTruthLabels[sid]])-1);
    }

    // pick max
    bool initialized = false;
    maxScore = buf[0];
    inferredLabels[sid] = 0;
    for(int s = 0; s < (int)param->nClasses; s++) {
      if(!replaceVoidMSRC || (s != voidLabel && s != moutainLabel && s != horseLabel)) {
        if(!initialized || maxScore < buf[s]) {
          maxScore = buf[s];
          inferredLabels[sid] = s;
          initialized = true;
        }
      }
    }
  }
  
  delete[] n;
  delete[] buf;
  return computeEnergy(inferredLabels);
}
