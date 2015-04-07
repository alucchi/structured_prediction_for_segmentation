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

#ifndef GI_MAXFLOW_H
#define GI_MAXFLOW_H

// Include graph.h first to avoid compilation errors
// due to the re-definition of symbols
#include "kgraph.h"
#include "maxflow.h"

// SliceMe
#include "Feature.h"
#include "Slice_P.h"

#include "graphInference.h"
#include "energyParam.h"

// standard libraries
#include <map>
#include <vector>

typedef double maxflow_cap_type;
typedef maxflow::Graph<float,float,float> GraphType;


//------------------------------------------------------------------------------

class GI_maxflow : public GraphInference
{
 public:

  /**
   * Constructor for SSVM framework
   */
  GI_maxflow(Slice_P* _slice, 
             const EnergyParam* _param,
             double* _smw,
             labelType* _groundTruthLabels,
             double* _lossPerLabel,
             Feature* _feature,
             std::map<sidType, nodeCoeffType>* _nodeCoeffs,
             std::map<sidType, edgeCoeffType>* _edgeCoeffs
             );

  ~GI_maxflow();

  void addLocalEdges();

  void addUnaryNodes();

  void createGraph();

  void precomputeEdgePotentials();

  void precomputeUnaryPotentials();

  double run(labelType* inferredLabels,
             int id,
             size_t maxiter,
             labelType* nodeLabelsGroundTruth = 0,
             bool computeEnergyAtEachIteration = false,
             double* _loss = 0);

 private:
  GraphType* g;

  maxflow_cap_type** unaryPotentials;
  ulong nUnaryPotentials;
  maxflow_cap_type* edgePotentials;
  ulong nEdgePotentials;
  maxflow_cap_type minPotential;

};


class GI_multiobject : public GraphInference
{
 public:

  /**
   * Constructor for SSVM framework
   */
  GI_multiobject(Slice_P* _slice, 
             const EnergyParam* _param,
             double* _smw,
             labelType* _groundTruthLabels,
             double* _lossPerLabel,
             Feature* _feature,
             std::map<sidType, nodeCoeffType>* _nodeCoeffs,
             std::map<sidType, edgeCoeffType>* _edgeCoeffs
             );

  ~GI_multiobject();

  void addInterLayerEdges();

  void addLocalEdges();

  void addUnaryNodes();

  void createGraph();

  void precomputeEdgePotentials();

  void precomputeDataTerm();

  double run(labelType* inferredLabels,
             int id,
             size_t maxiter,
             labelType* nodeLabelsGroundTruth = 0,
             bool computeEnergyAtEachIteration = false,
             double* _loss = 0);

 private:
  GraphType* g;

  maxflow_cap_type** unaryPotentials;
  ulong nUnaryPotentials;
  maxflow_cap_type* edgePotentials;
  ulong nEdgePotentials;
  maxflow_cap_type* interLayer_edgePotentials;
  ulong interLayer_nEdgePotentials;

  ulong nNodes;
  ulong nEdges;

  std::vector<long> nextLayerEdgeId;

};

#endif //GI_MAXFLOW_H
