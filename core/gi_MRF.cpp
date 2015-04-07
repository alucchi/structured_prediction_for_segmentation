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

#include "gi_MRF.h"

#include <cmath>
#include <limits.h>

// SliceMe
#include "Config.h"
#include "Slice_P.h"
#include "utils.h"

#define INF_COST INT_MAX/100000
#define MAX_COST INT_MAX/10000000

#include "inference_globals.h"

//#define SHIFT_AND_RESCALE_W 1
#define DEBUG 1

//-----------------------------------------------------------------------GLOBALS

//------------------------------------------------------------------------------


MRF::CostVal dCost(int pix, int i, void* ptr)
{
  GI_MRF* gi = (GI_MRF*) ptr;
  return gi->unaryCosts[pix][i];
}

MRF::CostVal fnCost(int pix1, int pix2, int i, int j, void* ptr)
{
  GI_MRF* gi = (GI_MRF*) ptr;

  if(!gi->param->includeLocalEdges)
    return 1;

  // compute index associated to edge (pix1, pix2)
  int edgeIdx = pix1*gi->nNodes + pix2;
  //INFERENCE_PRINT("edge %d %d %d %d\n",pix1,pix2,edgeIdx,gi->edgeIdxToEdgeType[edgeIdx]);
  // compute index associated to the pair of labels (i,j)
  int labelIdx;

  MRF::CostVal cost;
  if(gi->param->nGradientLevels == 0) {
    if(i == j) // different labels
      cost = gi->pairwiseCosts_Local[0][0];
    else
      cost = 0; // labels are the same
  } else {
    labelIdx = i*gi->param->nClasses + j;
    cost = gi->pairwiseCosts_Local[gi->edgeIdxToGradientIdx[edgeIdx]][labelIdx];
  }

  return cost;
}

//------------------------------------------------------------------------------

GI_MRF::GI_MRF(Slice_P* _slice, 
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

  mrf = 0;
  useQPBO = false;
  nNodes = slice->getNbSupernodes(); // local nodes
  nEdges = slice->getNbEdges(); // edges between local nodes

  if(param->includeLocalEdges) {
    if(param->nGradientLevels == 0) {
      nPairwiseCosts_Local = 1;
      nPairwiseStates = 1;
    } else {
      nPairwiseCosts_Local = (param->nOrientations*param->nGradientLevels);
      nPairwiseStates = param->nClasses*param->nClasses;
    }
  } else {
    nPairwiseCosts_Local = 0;
    nPairwiseStates = 0;
  }

  // add +1 causes indices start at 1
  int sizePsi = (param->nOrientations*param->nGradientLevels*param->nClasses*param->nClasses) + (param->nGradientLevels==0 && param->includeLocalEdges) + param->nUnaryWeights + param->nScalingCoefficients;
  INFERENCE_PRINT("[GI_MRF] param->nOrientations=%d, param->nGradientLevels=%d, param->nClasses=%d, sizePsi = %d\n",
         param->nOrientations, param->nGradientLevels, param->nClasses, sizePsi);

  double w_min = _smw[0];
  double w_max = _smw[0];
  for(int i=1; i< sizePsi; i++) {
    if(w_min > _smw[i])
      w_min = _smw[i];
    if(w_max < _smw[i])
      w_max = _smw[i];
    }
  
  // negate : cost = -score
  double c_min = -w_max;
  double c_max = -w_min;
#ifdef SHIFT_AND_RESCALE_W
  double c_diff = c_max - c_min;
  double t;
#endif

  cost = new energyType[sizePsi];
  for(int i=0; i< sizePsi; i++) {

#ifdef SHIFT_AND_RESCALE_W
      t = -_smw[i]; // negate
      t -= c_min;  //shift so that min value is 0
      if(c_diff != 0)
        t /= c_diff;  //rescale
      cost[i] = t*MAX_COST; //rescale

      // REMOVE THIS WITH DOUBLES !
      assert(cost[i]>=0);
#else
      cost[i] = -_smw[i]*MAX_COST;
#endif
    }

  INFERENCE_PRINT("[gi_MRF] c_min=%g, c_max=%g\n", c_min,c_max);

#ifdef DEBUG
  INFERENCE_PRINT("[gi_MRF] w\n");
  for(int i=0; i< sizePsi; i++)
    if(_smw[i] != 0)
      INFERENCE_PRINT("%d:%g ", i, _smw[i]);
  INFERENCE_PRINT("\n");
  INFERENCE_PRINT("[gi_MRF] cost\n");
  for(int i=0; i< sizePsi; i++)
    if(_smw[i] != 0)
      INFERENCE_PRINT("%d:%g ", i, (double)cost[i]);
  INFERENCE_PRINT("\n");
#endif

  if(lossPerLabel)
    {
      lossPerLabelRescaled = new energyType[param->nClasses];
      for(int c=0; c < (int)param->nClasses; c++)
        {

#ifdef SHIFT_AND_RESCALE_W
          lossPerLabelRescaled[c] = lossPerLabel[c]*MAX_COST;

          if(c_diff != 0)
            lossPerLabelRescaled[c]/=c_diff;
#else
          lossPerLabelRescaled[c] = lossPerLabel[c]*MAX_COST;
#endif

        }

#ifdef DEBUG
      INFERENCE_PRINT("[gi_MRF] lossPerLabel\n");
      for(int c=0; c < param->nClasses; c++)
        INFERENCE_PRINT("%d:%g ", c, (double)lossPerLabel[c]);
      INFERENCE_PRINT("\n");
      INFERENCE_PRINT("[gi_MRF] lossPerLabelRescaled\n");
      for(int c=0; c < param->nClasses; c++)
        INFERENCE_PRINT("%d:%g ", c, (double)lossPerLabelRescaled[c]);
      INFERENCE_PRINT("\n");
#endif

    }
  else
    {
      INFERENCE_PRINT("[gi_MRF] No loss specified\n");
      lossPerLabelRescaled = 0;
    }

  createGraph();
}

GI_MRF::~GI_MRF()
{
  for(int i=0; i < nNodes; i++)
    delete[] unaryCosts[i];
  delete[] unaryCosts;

  if(mrf)
    delete mrf;
  if(cost)
    delete[] cost;
  if(lossPerLabelRescaled)
    delete[] lossPerLabelRescaled;

  if(pairwiseCosts_Local) {
    for(int i=0; i < nPairwiseCosts_Local; i++)
      delete[] pairwiseCosts_Local[i];
    delete[] pairwiseCosts_Local;
  } 

  if(data)
    delete data;
  if(smooth)
    delete smooth;
  if(energy)
    delete energy;
}

void GI_MRF::createGraph()
{
  data = new DataCost(dCost);  
  smooth = new SmoothnessCost(fnCost); 
  energy = new EnergyFunction(data,smooth);

  // nNodes*nClasses table that contains unary costs for global and local nodes
  unaryCosts = new energyType*[nNodes];
  for(int i=0; i < nNodes; i++)
    unaryCosts[i] = new energyType[param->nClasses];

  if(param->includeLocalEdges) {
    pairwiseCosts_Local = new energyType*[nPairwiseCosts_Local];
    for(int i=0; i < nPairwiseCosts_Local; i++)
      pairwiseCosts_Local[i] = new energyType[nPairwiseStates];
  } else {
    pairwiseCosts_Local = 0;
  }

  minUnaryCost = INT_MAX;
  minPairwiseCost = INT_MAX;

  addUnaryNodes();

  if(param->includeLocalEdges) {
    precomputeEdgeCosts();
  }

  energyType minCost = min(minUnaryCost, minPairwiseCost);
  energyType offset = -minCost;
  if(minCost < 0) {
    INFERENCE_PRINT("[gi_MRF] min unary cost is negative (%g). Shifting all the values...\n", (double)minCost);
    // shift costs so that all values are positives !
    for(int i=0; i < nNodes; i++)
      for(int c=0; c < param->nClasses; c++) {
        unaryCosts[i][c] += offset;
      }
  }
  else {
    INFERENCE_PRINT("[gi_MRF] min unary cost is positive (%g)\n",(double)minCost);
  }

  if(minCost < 0) {
    INFERENCE_PRINT("[gi_MRF] min cost for local edges is negative (%g). Shifting all the values...\n", (double)minCost);
    // shift costs so that all values are positives !
    for(int i=0; i < nPairwiseCosts_Local; i++)
      for(int c=0; c < nPairwiseStates; c++) {
        pairwiseCosts_Local[i][c] += offset;
      }
  }
  else {
    INFERENCE_PRINT("[gi_MRF] min cost for local edges is positive (%g)\n",(double)minCost);
  }

#ifdef DEBUG2
  INFERENCE_PRINT("[gi_MRF] Unary costs (%dx%d):\n",nNodes,param->nClasses);
  for(int i=0; i < nNodes; i++)
    {
      INFERENCE_PRINT("%d:", i);
      for(int s=0; s < param->nClasses; s++)
        {
          INFERENCE_PRINT(" %d:%e",s,(float)unaryCosts[i][s]);
        }
      INFERENCE_PRINT("\n");
    }

  INFERENCE_PRINT("[gi_MRF] Pairwise costs for local edges (%dx%d):\n",nPairwiseCosts_Local,nPairwiseStates);
  for(int i=0; i < nPairwiseCosts_Local; i++)
    {
      INFERENCE_PRINT("%d:", i);
      for(int s=0; s < nPairwiseStates; s++)
        {
          INFERENCE_PRINT(" %d:%e",s,(float)pairwiseCosts_Local[i][s]);
        }
      INFERENCE_PRINT("\n");
    }
#endif

}

void GI_MRF::addUnaryNodes()
{
  // allocate memory to store features
  int fvSize = feature->getSizeFeatureVector();
  int max_index = fvSize + 1;
  osvm_node* n = new osvm_node[max_index];
  int i = 0;
  for(i = 0;i < max_index-1; i++)
    n[i].index = i+1;
  n[i].index = -1;

  int pix;
  double p = 0;
  // local nodes
  const map<int, supernode* >& _supernodes = slice->getSupernodes();
  for(map<int, supernode* >::const_iterator its = _supernodes.begin();
      its != _supernodes.end(); its++) {
    pix = its->first;

    if(param->nClasses != 2) {
      for(int i = 0; i < (int)param->nClasses; i++) {
        feature->getFeatureVector(n, slice, pix);
        p = 0;
        for(int s=0; s < fvSize; s++) {
          p += n[s].value*cost[SVM_FEAT_INDEX(param,i,s)];
        }

#ifdef W_OFFSET
        p += cost[i];
#endif
        /*
        if(nodeCoeffs) {
          p *= (*nodeCoeffs)[sid];
        }
        */
        unaryCosts[pix][i] = p;
      }
    } else {
      unaryCosts[pix][T_FOREGROUND] = 0;

      int i = T_BACKGROUND;
      feature->getFeatureVector(n, slice, pix);
      p = 0;
      for(int s=0; s < fvSize; s++) {
        p += n[s].value*cost[SVM_FEAT_INDEX(param,i,s)];
      }
#ifdef W_OFFSET
      p += cost[i];
#endif
      unaryCosts[pix][i] = p;
    }

    if(lossPerLabelRescaled) {
      for(int s = 0; s < (int)param->nClasses; s++) {
        if(s != groundTruthLabels[pix]) {
          // add loss of the ground truth label
          unaryCosts[pix][s] -= lossPerLabelRescaled[groundTruthLabels[pix]];
        }
      }      
    }
  }

  delete[] n;

  // compute min cost value
  for(int i=0; i < nNodes; i++) {
    for(int c=0; c < param->nClasses; c++) {
      if(unaryCosts[i][c] < minUnaryCost) {
        minUnaryCost = unaryCosts[i][c];
      }
    }
  }
}

void GI_MRF::addLocalEdges()
{
  //MRFGraph::captype weight = 1;
  MRF::CostVal weight = 1;

  // edges between local nodes
  const map<int, supernode* >& _supernodes = slice->getSupernodes();
  for(map<int, supernode* >::const_iterator its = _supernodes.begin();
      its != _supernodes.end(); its++) {
    vector < supernode* >* lNeighbors = &(its->second->neighbors);
    for(vector < supernode* >::iterator itN = lNeighbors->begin();
        itN != lNeighbors->end(); itN++) {
      // set edges once
      if(its->first < (*itN)->id)
        continue;

      // only call this function once for each edge !
      mrf->setNeighbors(its->first, (*itN)->id, weight);
    }
  }
}

void GI_MRF::precomputeEdgeCosts()
{
  uint edgeIdx;
  uint gradientIdx;
  int pix1, pix2; // pixel ids
  int orientationIdx = 0;

  // edges between local nodes

  // Go over all the edges and compute the gradient
  const map<int, supernode* >& _supernodes = slice->getSupernodes();
  for(map<int, supernode* >::const_iterator its = _supernodes.begin();
      its != _supernodes.end(); its++) {
    vector < supernode* >* lNeighbors = &(its->second->neighbors);
    pix1 = its->first;

    for(vector < supernode* >::iterator itN = lNeighbors->begin();
        itN != lNeighbors->end(); itN++) {
      // set edges once
      if(its->first < (*itN)->id) {
        continue;
      }
        
      pix2 = (*itN)->id;

      gradientIdx = slice->getGradientIdx(pix1,pix2);
      orientationIdx = slice->getOrientationIdx(pix1, pix2);
        
      edgeIdx = pix1*nNodes + pix2;
      edgeIdxToGradientIdx[edgeIdx] = gradientIdx;
    }
  }

  if(param->nGradientLevels == 0) {
    if(param->includeLocalEdges) {
      pairwiseCosts_Local[0][0] = cost[param->nUnaryWeights]; //param->nUnaryWeights = offset unary terms
      INFERENCE_PRINT("pairwiseCosts_Local[0][0] %g\n", pairwiseCosts_Local[0][0]);
    }
    else {
      pairwiseCosts_Local[0][0] = 0;
    }
  } else {
    // Compute cost (param->nClasses*param->nClasses table) for each gradient
    for(int g=0; g < param->nGradientLevels; g++) {
      for(int i=0; i < param->nClasses; i++) {
        for(int j=0; j < param->nClasses; j++) {
          uint labelIdx = i*param->nClasses+j;

          // w[0..nClasses-1] contains the unary weights
          double w_sum = 0;
          int _idx;
          for(int dg = 0; dg <= g; dg++) {
            _idx = (dg*nPairwiseStates*param->nOrientations) + (orientationIdx*nPairwiseStates) + labelIdx;
            w_sum += cost[_idx+param->nUnaryWeights]; //param->nUnaryWeights = offset unary terms
          }
          pairwiseCosts_Local[g][labelIdx] = w_sum;
          //printf("pairwiseCosts_Local[%d][%d] %g\n", g, labelIdx, pairwiseCosts_Local[g][labelIdx]);
        }
      }
    }
  }

  // Shift pairwiseCosts_Local
  for(int i=0; i < nPairwiseCosts_Local; i++)
    for(int c=0; c < nPairwiseStates; c++) {
      if(pairwiseCosts_Local[i][c] < minPairwiseCost)
        minPairwiseCost = pairwiseCosts_Local[i][c];
    }
}

double GI_MRF::run(labelType* inferredLabels,
                   int id,
                   size_t maxiter,
                   labelType* nodeLabelsGroundTruth,
                   bool computeEnergyAtEachIteration,
                   double* _loss)
{
  if(mrf) {
    INFERENCE_PRINT("[GI_MRF] MRF already instantiated\n");
    exit(-1);
  } else {
    // Assumes all the nodes have the same number of labels
    // Since the global nodes have less labels than the local nodes,
    // we set the cost of the extra labels to infinity.
    //mrf = new GCoptimization(nNodes, param->nClasses,energy, this);
    mrf = new Expansion(nNodes, param->nClasses,energy, this);
    //mrf = new MaxProdBP(nNodes, param->nClasses,energy, this);
    
    if(param->includeLocalEdges) {
      addLocalEdges();
    }
  }

  mrf->initialize();
  mrf->clearAnswer();

  if(useQPBO) {
    INFERENCE_PRINT("[GI_MRF] setUseQPBO to true, %d edges\n", nEdges);
    ((Expansion*)mrf)->setUseQPBO(true,2*nEdges);
    INFERENCE_PRINT("[GI_MRF] setProbing to true\n");
    ((Expansion*)mrf)->setUseProbing(true);
  }

  // Use non random order
  ((Expansion*)mrf)->setLabelOrder(false);
  //((Expansion*)mrf)->setLabelOrder(true);
    
  MRF::EnergyVal E = mrf->totalEnergy();
  INFERENCE_PRINT("[GI_MRF] Energy at the Start= %g (%g,%g)\n", (float)E,
         (float)mrf->smoothnessEnergy(), (float)mrf->dataEnergy());

#ifdef COUNT_TRUNCATIONS
  truncCnt = totalCnt = 0;
#endif
  float t;
  /*
  float tot_t = 0;
  for (int iter=0; iter<6; iter++) {
    mrf->optimize(1, t);

    E = mrf->totalEnergy();
    tot_t = tot_t + t ;
    INFERENCE_PRINT("[GI_MRF] energy = %g (%f secs)\n", (float)E, tot_t);
  }
  */
  mrf->optimize(100,t);  // run for 5 iterations, store time t it took 
#ifdef COUNT_TRUNCATIONS
  if (truncCnt > 0)
    INFERENCE_PRINT("[GI_MRF] ***WARNING: %d terms (%.2f%%) were truncated to ensure regularity\n", 
           truncCnt, (float)(100.0 * truncCnt / totalCnt));
#endif

  int n = slice->getNbSupernodes();
  INFERENCE_PRINT("[GI_MRF] Copying %d labels\n", n);
  for (int pix = 0; pix < n; pix++) {
    inferredLabels[pix] = mrf->getLabel(pix);
    //assert(inferredLabels[pix]<2);
    //INFERENCE_PRINT("inferredLabels[pix] %d %d\n", pix, inferredLabels[pix]);
  }

  INFERENCE_PRINT("[GI_MRF] Done\n");

  //return mrf->totalEnergy();
  return GraphInference::computeEnergy(inferredLabels);
}
