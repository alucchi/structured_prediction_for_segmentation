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

#include "graphInference.h"

#include "globalsE.h"
#include "Config.h"

//------------------------------------------------------------------------------

map<ulong, labelType> GraphInference::classIdxToLabel;

//------------------------------------------------------------------------------

GraphInference::GraphInference(Slice_P* _slice,
                               const EnergyParam* _param,
                               double* _smw,
                               Feature* _feature,
                               map<sidType, nodeCoeffType>* _nodeCoeffs,
                               map<sidType, edgeCoeffType>* _edgeCoeffs
                               )
{
  init();

  slice = _slice;
  param = _param;
  smw = _smw;
  feature = _feature;
  nodeCoeffs = _nodeCoeffs;
  edgeCoeffs = _edgeCoeffs;
}

GraphInference::~GraphInference()
{
}

void GraphInference::init()
{
  nodeCoeffs = 0;
  edgeCoeffs = 0;
  lossPerLabel = 0;
}

/**
 * Compute energy for a given configuration nodeLabels
 */
double GraphInference::computeEnergy(labelType* nodeLabels)
{
  // Compute energy
  double energyU = 0.0;
  double energyP = 0.0;
  double loss = 0.0;
  int nDiff = 0;

  int fvSize = feature->getSizeFeatureVector();
  int label = 0;
  int nSupernodes = slice->getNbSupernodes();
  double energySupernode = 0;
  double energyEdge = 0;
  for(int sid = 0; sid < nSupernodes; sid++)
    {
      label = nodeLabels[sid];

      if(lossPerLabel) {
        if(label != groundTruthLabels[sid]) {
          if(nodeCoeffs) {
            loss += (*nodeCoeffs)[sid]*lossPerLabel[groundTruthLabels[sid]];
          } else {
            loss += lossPerLabel[groundTruthLabels[sid]];
          }
          ++nDiff;
        }
      }

      if (param->nUnaryWeights == 1 && label == T_FOREGROUND) {
	// Only accumulates weights for BACKGROUND class.
	continue;
      }

      osvm_node *n = slice->getFeature(sid);
      energySupernode = 0;
      for(int s = 0; s < fvSize; s++) {
        energySupernode -= smw[SVM_FEAT_INDEX(param, label,s)]*n[s].value;
      }

#ifdef W_OFFSET
      energySupernode -= smw[label];
#endif
      
      if(nodeCoeffs) {
        energySupernode *= (*nodeCoeffs)[sid];
      }

      energyU += energySupernode;
    }

  if(param->includeLocalEdges)
    {
      // add energy for pairwize term
      if(param->nGradientLevels == 0) {
        ulong edgeId = 0;
        const map<int, supernode* >& _supernodes = slice->getSupernodes();
        for(map<sidType, supernode* >::const_iterator itNode = _supernodes.begin();
            itNode != _supernodes.end(); itNode++) {
          for(vector<supernode*>::iterator itNode2 = itNode->second->neighbors.begin();
              itNode2 != itNode->second->neighbors.end(); itNode2++) {
            // set edges once
            if(itNode->first < (*itNode2)->id) {
              continue;
            }
            
            if(nodeLabels[itNode->first] == nodeLabels[(*itNode2)->id]) {
              energyEdge = smw[param->nUnaryWeights];
              if(edgeCoeffs) {
                energyEdge *= (*edgeCoeffs)[edgeId];
              }
              energyP -= energyEdge;
            }
            ++edgeId;
          }
        }
      } else {
#if USE_LONG_RANGE_EDGES
        int distanceIdx;
#endif
        int gradientIdx;
        int orientationIdx;
        int w_edgeIdx;
        ulong edgeId = 0;
        const map<int, supernode* >& _supernodes = slice->getSupernodes();
        for(map<sidType, supernode* >::const_iterator itNode = _supernodes.begin();
            itNode != _supernodes.end(); itNode++) {
          for(vector<supernode*>::iterator itNode2 = itNode->second->neighbors.begin();
              itNode2 != itNode->second->neighbors.end(); itNode2++) {
            // set edges once
            if(itNode->first < (*itNode2)->id) {
              continue;
            }
          
            gradientIdx = slice->getGradientIdx(itNode->first,(*itNode2)->id);
            orientationIdx = slice->getOrientationIdx(itNode->first, (*itNode2)->id);

            int offset = (orientationIdx*param->nClasses*param->nClasses);

#if USE_LONG_RANGE_EDGES
            distanceIdx = slice->getDistanceIdx(itNode->first,
                                                (*itNode2)->id);
            offset += distanceIdx*param->nGradientLevels*param->nClasses*param->nClasses*param->nOrientations;
#endif

            energyEdge = 0;
            for(int i =0; i <= gradientIdx; i++) {
              w_edgeIdx = (i*param->nClasses*param->nClasses*param->nOrientations) + offset + nodeLabels[itNode->first]*param->nClasses + nodeLabels[(*itNode2)->id];
              energyEdge -= smw[w_edgeIdx + param->nUnaryWeights];
            }

            if(edgeCoeffs) {
              energyEdge *= (*edgeCoeffs)[edgeId];
            }

            energyP += energyEdge;
            ++edgeId;
          }
        }
      }
    }

  double energy = energyU + energyP;

  INFERENCE_PRINT("[graphInference] Energy unary=%g, pairwise=%g, nDiff=%d/%d, loss=%g, totalE=%g, totalE-loss=%g\n",
                  energyU, energyP, nDiff, nSupernodes, loss, energy, energy-loss);

  return energy - loss;
}

void GraphInference::computeNodePotentials(double**& unaryPotentials, double& maxPotential)
{
  // allocate memory to store features
  int fvSize = feature->getSizeFeatureVector();
  int max_index = fvSize + 1;
  osvm_node* n = new osvm_node[max_index];
  int i = 0;
  for(i = 0;i < max_index-1; i++)
    n[i].index = i+1;
  n[i].index = -1;

  bool useLossFunction = lossPerLabel!=0;
  int sid = 0;
  double p = 0;
  const map<int, supernode* >& _supernodes = slice->getSupernodes();
  maxPotential = 0;
  for(map<int, supernode* >::const_iterator its = _supernodes.begin();
      its != _supernodes.end(); its++) {
    sid = its->first;
    double* buf = unaryPotentials[sid];

    if(param->nClasses != 2) {
      for(int i = 0; i < (int)param->nClasses; i++) {
        //INFERENCE_PRINT("1.%d:%e:%e:%e ", i, buf[i],std::log(buf[i]),std::log(buf[i])*smw[SVM_PROB_INDEX]);

        feature->getFeatureVector(n, slice, sid);
        p = 0;
        for(int s=0; s < fvSize; s++) {
          p += n[s].value*smw[SVM_FEAT_INDEX(param,i,s)];
        }

#ifdef W_OFFSET
        p += smw[i];
#endif
        if(nodeCoeffs) {
          p *= (*nodeCoeffs)[sid];
        }
        buf[i] = p;

        if (fabs(buf[i]) > maxPotential) {
          maxPotential = fabs(buf[i]);
        }
      }
    } else {
      // Only 2 classes.
      buf[T_FOREGROUND] = 0;

      int i = T_BACKGROUND;
      feature->getFeatureVector(n, slice, sid);
      p = 0;
      for(int s=0; s < fvSize; s++) {
        assert(s+1 == n[s].index);
        p += n[s].value*smw[SVM_FEAT_INDEX(param,i,s)];
      }
#ifdef W_OFFSET
      p += smw[i];
#endif

      if(nodeCoeffs) {
        p *= (*nodeCoeffs)[sid];
      }

      buf[i] = p;

      if (fabs(buf[i]) > maxPotential) {
        maxPotential = fabs(buf[i]);
      }
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

          if (fabs(buf[s]) > maxPotential) {
            maxPotential = fabs(buf[s]);
          }
        }
      }

      // loss function : -1 if correct label (same as adding +1 to all incorrect labels)
      //buf[groundTruthLabels[sid]] = buf[groundTruthLabels[sid]]-1;
    }

    for(int c = 0; c < (int)param->nClasses; c++) {
      unaryPotentials[sid][c] = buf[c];
    }
  }
  delete[] n;
}
