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

#include "gi_multiobject.h"

// SliceMe
#include "Config.h"
#include "utils.h"

#include "inference_globals.h"

//------------------------------------------------------------------------------

// this code was written for 2 layers only. Do not change this variable!
#define N_LAYERS 2

using namespace std;

#define OUTSIDE_LABEL 0
#define INSIDE_LABEL 1

//------------------------------------------------------------------------------

GI_multiobject::GI_multiobject(Slice_P* _slice,
                       const EnergyParam* _param,
                       double* _smw,
                       labelType* _groundTruthLabels,
                       double* _lossPerLabel,
                       Feature* _feature,
                       map<sidType, nodeCoeffType>* _nodeCoeffs,
                       map<sidType, edgeCoeffType>* _edgeCoeffs)
{
  GraphInference::init();
  slice = _slice;
  param = _param;
  smw = _smw;
  lossPerLabel = _lossPerLabel;
  groundTruthLabels = _groundTruthLabels;
  feature = _feature;
  nodeCoeffs = _nodeCoeffs;
  edgeCoeffs = _edgeCoeffs;
  g = 0;
  edgePotentials = 0;
  nEdgePotentials = 0;
  nNodes = 0;
  nEdges = 0;

  createGraph();

  assert(param->nGradientLevels > 0);
  assert(nodeCoeffs == 0);
  assert(edgeCoeffs == 0);
}

GI_multiobject::~GI_multiobject() {
  if(g) {
    delete g;
  }
  if(unaryPotentials) {
    for(uint p = 0; p < nUnaryPotentials; p++) {
      delete[] unaryPotentials[p];
    }
    delete[] unaryPotentials;
  }
  if(edgePotentials) {
    delete[] edgePotentials;
  }
}

void GI_multiobject::createGraph()
{
  if(g) {
    delete g;
  }

  nNodes = slice->getNbSupernodes();
  // compute estimate number of edges
  // this ONLY enforces a 1-supernode thickness boundary for containment constraint!
  ulong nInterLayerEdges = slice->getNbUndirectedEdges();
  nEdges = (nInterLayerEdges*(N_LAYERS+2*(N_LAYERS-1))) + nNodes*(N_LAYERS-1);
  INFERENCE_PRINT("[GI_multiobject] nNodes=%ld nEdges=%ld nTotalEdges=%ld\n", nNodes, nInterLayerEdges, nEdges);
  g = new GraphType(nNodes*N_LAYERS, nEdges);

  // allocate memory for unary and pairwise potentials
  nUnaryPotentials = nNodes*N_LAYERS;
  unaryPotentials = new maxflow_cap_type*[nUnaryPotentials];
  for(uint p = 0; p < nUnaryPotentials; p++) {
    unaryPotentials[p] = new maxflow_cap_type[param->nClasses];

    for(int c = 0; c < param->nClasses; ++c) {
      unaryPotentials[p][c] = 0;
    }
  }

  nEdgePotentials = nEdges;
  edgePotentials = new maxflow_cap_type[nEdgePotentials];
  for(uint e = 0; e < nEdgePotentials; ++e) {
    edgePotentials[e] = 0;
  }

  // pre-allocate array
  nextLayerEdgeId.resize(nNodes, -1);

  if(param->includeLocalEdges) {
    precomputeEdgePotentials();
  }
  precomputeDataTerm();

  addUnaryNodes();
  if(param->includeLocalEdges) {
    addLocalEdges();
  }
}

/**
 * Run maxflow on a given factor graph
 * @param inferredLabels is the output
 */
double GI_multiobject::run(labelType* inferredLabels,
                       int id,
                       size_t maxiter,
                       labelType* nodeLabelsGroundTruth,
                       bool computeEnergyAtEachIteration,
                       double* _loss)
{
  double flow = g->maxflow();
  INFERENCE_PRINT("[GI_multiobject] flow=%g\n", flow);
  INFERENCE_PRINT("[GI_multiobject] nNodes=%ld nEdges=%ld\n", nNodes, nEdges);

  // check binary labels of each supernode in both layers to infer what the label is
  const map<sidType, supernode* >& _supernodes = slice->getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    // check label of node in first layer
    if(g->what_segment(it->first) == GraphType::SOURCE) {
      // check label of node in second layer
      if(g->what_segment(it->first + nNodes) == GraphType::SOURCE) {
        inferredLabels[it->first] = BACKGROUND;
      } else {
        // should not happen as we assigned an infinite cost to this configuration
        assert(0);
      }
    } else {
      // check label of node in second layer
      if(g->what_segment(it->first + nNodes) == GraphType::SINK) {
        inferredLabels[it->first] = FOREGROUND;
      } else {
        inferredLabels[it->first] = BOUNDARY;
      }
    }
  }

  // sanity check
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    supernode* s = it->second;
    for(vector < supernode* >::iterator itN = s->neighbors.begin();
        itN != s->neighbors.end();itN++) {

      // set edges once
      if(it->first < (*itN)->id) {
        continue;
      }

      if( (inferredLabels[it->first] == BACKGROUND && inferredLabels[(*itN)->id] == FOREGROUND)
          || (inferredLabels[it->first] == FOREGROUND && inferredLabels[(*itN)->id] == BACKGROUND) ) {
        printf("[gi_multiobject] WARNING %d %d\n", it->first, (*itN)->id);
      }
    }
  }

  double energy = GraphInference::computeEnergy(inferredLabels);
  return energy;
}

void GI_multiobject::addUnaryNodes()
{
  const map<sidType, supernode* >& _supernodes = slice->getSupernodes();
  g->add_node(_supernodes.size()*N_LAYERS);

  ulong sid = 0;
  for(int l = 0; l < N_LAYERS; ++l) {
    for(uint i = 0; i < _supernodes.size(); ++i) {
      // source capacity is first

      maxflow_cap_type minPotential = min(unaryPotentials[sid][INSIDE_LABEL], unaryPotentials[sid][OUTSIDE_LABEL]);
      g->add_tweights(sid, unaryPotentials[sid][OUTSIDE_LABEL] - minPotential,
                      unaryPotentials[sid][INSIDE_LABEL] - minPotential);
      ++sid;
    }
  }
}

// todo: pass unary term
// 0 0: Background
// 0 1: K (does not matter)
// 1 0: Boundary
// 1 1: Foreground
void GI_multiobject::precomputeDataTerm()
{
  const map<sidType, supernode* >& _supernodes = slice->getSupernodes();
  bool useLossFunction = lossPerLabel!=0;

  double weightBackground = 0;
  double weightForeground = 0;
  double weightBoundary = 0;
  sidType sid = 0;
  double score[4];
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    sid = it->first;

    weightBackground = computeUnaryPotential(slice, sid, BACKGROUND);
    weightBoundary = computeUnaryPotential(slice, sid, BOUNDARY);
    weightForeground = computeUnaryPotential(slice, sid, FOREGROUND);

    if(useLossFunction) {
      // add loss of the ground truth label
      if(FOREGROUND == groundTruthLabels[sid]) {
        weightForeground -= lossPerLabel[groundTruthLabels[sid]];
      } else {
        if(BOUNDARY == groundTruthLabels[sid]) {
          weightBoundary -= lossPerLabel[groundTruthLabels[sid]];
        } else {
          weightBackground -= lossPerLabel[groundTruthLabels[sid]];
        }
      }
    }

    score[0] = weightBackground;
    score[1] = -1e20;
    score[2] = weightBoundary;
    score[3] = weightForeground;
    double D = score[0] + score[3] - score[1] - score[2];
    assert(D>=0); //submodularity condition

    unaryPotentials[sid][OUTSIDE_LABEL] += score[0]-score[2]; // A-C
    unaryPotentials[nNodes + sid][INSIDE_LABEL] += score[3] - score[2]; // D-C
    edgePotentials[nextLayerEdgeId[sid]] += D;
  }
}

void GI_multiobject::precomputeEdgePotentials()
{
  int multiObj_nClasses = 2;
  int nPairwiseStates = multiObj_nClasses*multiObj_nClasses;
  double *score = new double[nPairwiseStates];
  int gradientIdx;
  int idx;
  int oidx = 0;
  supernode *s;
  ulong edgeId = 0;
  double s33[3][3];

  for(int n = 0; n < N_LAYERS; ++n) {

    const map<int, supernode* >& _supernodes = slice->getSupernodes();
    for(map<int, supernode* >::const_iterator its = _supernodes.begin();
        its != _supernodes.end(); its++) {
      s = its->second;

      vector < supernode* >* lNeighbors = &(s->neighbors);
      for(vector < supernode* >::iterator itN = lNeighbors->begin();
          itN != lNeighbors->end(); itN++) {

        // set edges once
        if(its->first < (*itN)->id) {
          continue;
        }

        // get gradient index. Do not use any orientation index with maxflow
        // as it's only used for the EM dataset
        gradientIdx = slice->getGradientIdx(its->first,(*itN)->id);

#if USE_LONG_RANGE_EDGES
        int distanceIdx = slice->getDistanceIdx(s->id, (*itN)->id);
        oidx = distanceIdx*param->nGradientLevels*param->nClasses*param->nClasses*param->nOrientations;
#endif

        for(int r = 0; r < 3; ++r) {
          for(int c = 0; c < 3; ++c) {
            double w_sum = 0;
            for(int i = 0; i <= gradientIdx; i++) {
              int p = 3*r + c;
              idx = (i*param->nClasses*param->nClasses) + oidx + p;
              w_sum += smw[idx+param->nUnaryWeights]; // param->nUnaryWeights is the offset due to unary terms
            }
            s33[r][c] = w_sum;
          }
        }

        if(n == 0) {
          // layer 0 is inside + boundary
          score[0] = s33[BACKGROUND][BACKGROUND];
          score[1] = s33[BACKGROUND][BOUNDARY];
          score[2] = s33[BOUNDARY][BACKGROUND];
          score[3] = s33[BOUNDARY][BOUNDARY];
        } else {
          // layer 1 is inside
          score[0] = s33[BOUNDARY][BOUNDARY];
          score[1] = s33[BOUNDARY][FOREGROUND];
          score[2] = s33[FOREGROUND][BOUNDARY];
          score[3] = s33[FOREGROUND][FOREGROUND];
        }

        //double D = edgePotential[1] + edgePotential[2] - edgePotential[0] - edgePotential[3];
        double D = score[0] + score[3] - score[1] - score[2];
        if(D < 0) {
          printf("Score for layer %d, gradient %d = %g\na=%f, d=%f\n%f %f\n%f %f\n",
                 n, gradientIdx, D, score[0] + score[3], score[1] + score[2], score[0], score[1], score[2], score[3]);

          idx = 0;
          printf("pw\n");
          double* _pw = smw + param->nUnaryWeights;
          for(int g = 0; g < param->nGradientLevels; ++g) {
            for(int i = 0; i < param->nClasses; ++i) {
              for(int j = 0; j < param->nClasses; ++j) {
                printf("%g ", _pw[idx]);
                ++idx;
              }
              printf("\n");
            }
            //_pw += param->nClasses*param->nClasses;
          }
          printf("\n");

          for(int r = 0; r < 3; ++r) {
            for(int c = 0; c < 3; ++c) {
              double w_sum = 0;
              for(int i = 0; i <= gradientIdx; i++) {
                int p = 3*r + c;
                idx = (i*param->nClasses*param->nClasses) + oidx + p;
                w_sum += smw[idx+param->nUnaryWeights]; // param->nUnaryWeights is the offset due to unary terms
                printf("rc %d %d %d %g\n", r, c, idx, smw[idx+param->nUnaryWeights]);
              }
              s33[r][c] = w_sum;
            }
          }

          printf("s33\n");
          for(int r = 0; r < 3; ++r) {
            for(int c = 0; c < 3; ++c) {
              printf("%f ", s33[r][c]);
            }
            printf("\n");
          }
          printf("\n");

        }
        assert(D>=0); //submodularity condition

        unaryPotentials[its->first + (n*nNodes)][OUTSIDE_LABEL] += (score[0] - score[2]); // A-C
        unaryPotentials[(*itN)->id + (n*nNodes)][INSIDE_LABEL] += (score[3] - score[2]); // D-C

        edgePotentials[edgeId] += D;

        ++edgeId;
      }
    }
  }

  // layer 0 corresponds to foreground + boundary
  // layer 1 corresponds to foreground only
  const map<int, supernode* >& _supernodes = slice->getSupernodes();
  for(map<int, supernode* >::const_iterator its = _supernodes.begin();
      its != _supernodes.end(); its++) {
    s = its->second;

    score[0] = 0;
    score[1] = -1e30; // -infinity
    score[2] = 0;
    score[3] = 0;
    double D = score[0] + score[3] - score[1] - score[2];
    //printf("score = %g\n%g %g\n%g %g\n",D,score[0], score[1], score[2], score[3]);
    assert(D>=0); //submodularity condition

    // s-t links capacities are 0
    // so we don't update their t-weights (stored in "unaryPotentials") here

    // link between a node and itself in the next layer
    edgePotentials[edgeId] += D;
    nextLayerEdgeId[its->first] = edgeId;
    ++edgeId;

    vector < supernode* >* lNeighbors = &(s->neighbors);
    for(vector < supernode* >::iterator itN = lNeighbors->begin();
        itN != lNeighbors->end(); itN++) {
      edgePotentials[edgeId] += D;
      ++edgeId;
    }
  }

  delete[] score;

  INFERENCE_PRINT("[gi_multiobj] EDGE %ld %ld\n", edgeId, nEdgePotentials);
}

void GI_multiobject::addLocalEdges()
{
  const map<int, supernode* >& _supernodes = slice->getSupernodes();
  ulong edgeId = 0;
  supernode *s;

  // Must iterate in the exact same order as precomputeEdgePotentials!!

  for(int n = 0; n < N_LAYERS; ++n) {
    for(map<int, supernode* >::const_iterator its = _supernodes.begin();
        its != _supernodes.end(); its++) {
      s = its->second;
      vector < supernode* >* lNeighbors = &(s->neighbors);
      for(vector < supernode* >::iterator itN = lNeighbors->begin();
          itN != lNeighbors->end(); itN++) {

        // set edges once
        if(its->first < (*itN)->id) {
          continue;
        }

        // add edge in all layers
        g->add_edge(its->first + (n*nNodes), (*itN)->id + (n*nNodes), edgePotentials[edgeId], 0);

        ++edgeId;
      }
    }
  }

  // the following is valid for 2 layers only
  for(map<int, supernode* >::const_iterator its = _supernodes.begin();
      its != _supernodes.end(); its++) {
    s = its->second;

    // s-t links capacities are 0
    // so we don't update their t-weights (stored in "unaryPotentials") here

    // link between a node and itself in the next layer
    g->add_edge(its->first, its->first + nNodes, edgePotentials[edgeId], 0);
    ++edgeId;

    vector < supernode* >* lNeighbors = &(s->neighbors);
    for(vector < supernode* >::iterator itN = lNeighbors->begin();
        itN != lNeighbors->end(); itN++) {
      g->add_edge((*itN)->id, its->first + nNodes, edgePotentials[edgeId], 0);
      ++edgeId;

    }
  }
}
