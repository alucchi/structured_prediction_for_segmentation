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

#include "gi_maxflow.h"

// SliceMe
#include "Config.h"
#include "utils.h"

#include "inference_globals.h"

//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

GI_maxflow::GI_maxflow(Slice_P* _slice,
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
  unaryPotentials = 0;
  edgePotentials = 0;
  nUnaryPotentials = 0;
  nEdgePotentials = 0;
  createGraph();
}

GI_maxflow::~GI_maxflow() {
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

void GI_maxflow::createGraph()
{
  if(g) {
    delete g;
  }

  g = new GraphType(slice->getNbSupernodes(), slice->getNbUndirectedEdges());
  precomputeUnaryPotentials();
  if(param->includeLocalEdges) {
    precomputeEdgePotentials();
  }

  // Shift unary potentials
  // Note that there is no need to rescale the edge terms because of submodularity condition
  minPotential = unaryPotentials[0][0];
  maxflow_cap_type maxPotential = unaryPotentials[0][0];
  for(uint p = 0; p < nUnaryPotentials; p++) {
    for(int c = 0; c < param->nClasses; c++) {
      if(unaryPotentials[p][c] < minPotential) {
        minPotential = unaryPotentials[p][c];
      }
      if(fabs(unaryPotentials[p][c]) > maxPotential) {
        maxPotential = fabs(unaryPotentials[p][c]);
      }
    }
  }
  INFERENCE_PRINT("[GI_maxflow] minPotential for nodes = %g\n", minPotential);
  INFERENCE_PRINT("[GI_maxflow] maxPotential for nodes = %g\n", maxPotential);

  if(minPotential < 0) {
    for(uint p = 0; p < nUnaryPotentials; p++) {
      for(int c = 0; c < param->nClasses; c++) {
        unaryPotentials[p][c] -= minPotential;
      }
    }
  }

  addUnaryNodes();
  if(param->includeLocalEdges) {
    addLocalEdges();
  }
}

/**
 * Run BP on a given factor graph
 * @param nodeLabelsBP node inferred by BP
 */
double GI_maxflow::run(labelType* inferredLabels,
                       int id,
                       size_t maxiter,
                       labelType* nodeLabelsGroundTruth,
                       bool computeEnergyAtEachIteration,
                       double* _loss)
{
  double flow = g->maxflow();
  INFERENCE_PRINT("[GI_maxflow] flow=%g\n", flow);

  const map<sidType, supernode* >& _supernodes = slice->getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    if(g->what_segment(it->first) == GraphType::SOURCE) {
      inferredLabels[it->first] = BACKGROUND;
    } else {
      inferredLabels[it->first] = FOREGROUND;
    }
  }

  double energy = GraphInference::computeEnergy(inferredLabels);
  return energy;
}

void GI_maxflow::addUnaryNodes()
{
  const map<sidType, supernode* >& _supernodes = slice->getSupernodes();
  g->add_node(_supernodes.size());

  // a node connected to source is BACKGROUND
  // a node connected to sink is FOREGROUND
  int sid = 0;
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    //assert(sid==it->first); //sanity check
    sid = it->first;
    // source capacity is first
    g->add_tweights(sid, unaryPotentials[sid][T_BACKGROUND],
                    unaryPotentials[sid][T_FOREGROUND]);
  }
}

void GI_maxflow::precomputeUnaryPotentials()
{
  const map<sidType, supernode* >& _supernodes = slice->getSupernodes();
  nUnaryPotentials = _supernodes.size();
  unaryPotentials = new maxflow_cap_type*[nUnaryPotentials];
  for(uint p = 0; p < nUnaryPotentials; p++) {
    unaryPotentials[p] = new maxflow_cap_type[param->nClasses];
  }

  bool useLossFunction = lossPerLabel!=0;

  string config_tmp;
  int loss_function = 0;
  if(Config::Instance()->getParameter("loss_function", config_tmp)) {
    loss_function = atoi(config_tmp.c_str());
  }

  double weightToSource = 0;
  double weightToSink = 0;
  sidType sid = 0;
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    sid = it->first;

    // Source = background
    weightToSource = computeUnaryPotential(slice, sid, T_BACKGROUND);

    // Sink = foreground
    weightToSink = 0;

    if(useLossFunction) {
      // add loss of the ground truth label
      if(T_FOREGROUND != groundTruthLabels[sid]) {
        if(loss_function == LOSS_NODE_BASED) {
          weightToSink += lossPerLabel[sid];
        } else {
          weightToSink += lossPerLabel[groundTruthLabels[sid]];
        }
      } else {
        if(loss_function == LOSS_NODE_BASED) {
          weightToSource += lossPerLabel[sid];
        } else {
          weightToSource += lossPerLabel[groundTruthLabels[sid]];
        }
      }
    }

    if(nodeCoeffs) {
      weightToSource *= (*nodeCoeffs)[sid];
      weightToSink *= (*nodeCoeffs)[sid];
    }

    unaryPotentials[sid][T_BACKGROUND] = weightToSource;
    unaryPotentials[sid][T_FOREGROUND] = weightToSink;
    //printf("sid %d %g %g\n", sid, weightToSource, weightToSink);
  }
}

void GI_maxflow::precomputeEdgePotentials()
{
  nEdgePotentials = slice->getNbUndirectedEdges();
  edgePotentials = new maxflow_cap_type[nEdgePotentials];

  // Pairwise factors
  if(param->nGradientLevels > 0) {
    int nPairwiseStates = param->nClasses*param->nClasses;
    double *score = new double[nPairwiseStates];
    int gradientIdx;
    int idx;
    int oidx = 0;
    supernode *s;
    const map<int, supernode* >& _supernodes = slice->getSupernodes();
    ulong edgeId = 0;
    for(map<int, supernode* >::const_iterator its = _supernodes.begin();
        its != _supernodes.end(); its++) {
      s= its->second;

      vector < supernode* >* lNeighbors = &(s->neighbors);
      for(vector < supernode* >::iterator itN = lNeighbors->begin();
          itN != lNeighbors->end(); itN++) {
        //INFERENCE_PRINT("GI_maxflow::addLocalEdges %d %d\n", its->first,(*itN)->id);

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
        //printf("distanceIdx %d/%d\n",distanceIdx,param->nDistances);
#endif

        for(int p = 0; p < nPairwiseStates; p++ ) {
          // w[0..nClasses-1] contains the unary weights
          double w_sum = 0;
          for(int i = 0; i <= gradientIdx; i++) {
            idx = (i*param->nClasses*param->nClasses) + oidx + p;
            w_sum += smw[idx+param->nUnaryWeights]; // param->nUnaryWeights is the offset due to unary terms
          }

          if(edgeCoeffs) {
            w_sum *= (*edgeCoeffs)[edgeId];
          }

          score[p] = w_sum;
          
        }

        double D = score[0] + score[3] - score[1] - score[2];
        assert(D>=0); //submodularity condition

        unaryPotentials[its->first][T_BACKGROUND] += (score[0] - score[2]); // A-C
        unaryPotentials[(*itN)->id][T_FOREGROUND] += (score[3] - score[2]); // D-C

        edgePotentials[edgeId] = D;

        ++edgeId;
      }
    }
    delete[] score;
    assert(edgeId == nEdgePotentials);
  } else {
    assert(0);
    // potts model
    //smw[param->nUnaryWeights]);
  }
}

void GI_maxflow::addLocalEdges()
{
  ulong edgeId = 0;
  supernode *s;

  const map<int, supernode* >& _supernodes = slice->getSupernodes();
  for(map<int, supernode* >::const_iterator its = _supernodes.begin();
      its != _supernodes.end(); its++) {
    s= its->second;
    vector < supernode* >* lNeighbors = &(s->neighbors);
    for(vector < supernode* >::iterator itN = lNeighbors->begin();
        itN != lNeighbors->end(); itN++) {

      // set edges once
      if(its->first < (*itN)->id) {
        continue;
      }

      g->add_edge(its->first, (*itN)->id, edgePotentials[edgeId], 0);
      ++edgeId;
    }
  }
}
