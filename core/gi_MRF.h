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

// SliceMe
#include "Slice_P.h"

#include "graphInference.h"
#include "energyParam.h"

#include "mrf.h"
#include "GCoptimization.h"
//#include "MaxProdBP.h"

#include <map>
#include <vector>

typedef int giLabelType;

//typedef int energyType;
typedef double energyType;

//------------------------------------------------------------------------------

class GI_MRF : public GraphInference
{
 public:

  GI_MRF(Slice_P* _slice, 
         const EnergyParam* _param,
         double* _smw,
         labelType* _groundTruthLabels,
         double* _lossPerLabel,
         Feature* _feature,
         std::map<sidType, nodeCoeffType>* _nodeCoeffs);

  ~GI_MRF();

  void addLocalEdges();

  void createGraph();

  void addUnaryNodes();
  void precomputeEdgeCosts();

  double run(labelType* inferredLabels,
             int id,
             size_t maxiter,
             labelType* nodeLabelsGroundTruth,
             bool computeEnergyAtEachIteration,
             double* _loss = 0);

  energyType* lossPerLabelRescaled;

  int nNodes; // total number of nodes in the graph (local and global nodes)
  int nEdges;

  // Costs for graph nodes and edges are precomputed
  // and store in those 2 arrays for efficiency reasons
  energyType** unaryCosts;
  energyType** pairwiseCosts_Local; // cost associated to each gradient : nPairwiseCosts*nPairwiseStates
  int nPairwiseCosts_Local;
  int nPairwiseStates;

  map<uint,uint> edgeIdxToGradientIdx;

  void setUseQPBO(bool _val) { useQPBO = _val; }

 private:
  MRF* mrf;
  EnergyFunction *energy;
  DataCost *data;
  SmoothnessCost *smooth;

  /**
   * Store the cost computed from the w vector
   */
  energyType* cost;

  energyType minUnaryCost;
  energyType minPairwiseCost;

  bool useQPBO;
};
