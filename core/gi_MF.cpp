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

#include "gi_MF.h"

// SliceMe
#include "Config.h"
#include "globalsE.h"
#include "utils.h"

#include "inference_globals.h"

//------------------------------------------------------------------------------

#define MAX_POTENTIAL 1.0

#define EXP_DOMAIN 0

//------------------------------------------------------------------------------

GI_MF::GI_MF(Slice_P* _slice, 
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
  believes = 0;
  ownBelievesBuffer = true;
}

GI_MF::~GI_MF()
{
  if(ownBelievesBuffer && believes) {
    for(int c = 0; c < param->nClasses; ++c) {
      delete[] believes[c];
    }
    delete[] believes;
  }
}

double GI_MF::run(labelType* inferredLabels,
                   int id,
                   size_t maxiter,
                   labelType* nodeLabelsGroundTruth,
                   bool computenergyAtEachIteration,
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
    printf("[GI_MF] MSRC void=%d, moutain=%d, horse=%d\n",
           (int)voidLabel, (int)moutainLabel, (int)horseLabel);
  } else {
    printf("[GI_MF] Do not replace void labels\n");
  }

  int sid = 0;
  double maxScore = 0;
  double totalScore_old = 0;
  double totalScore = 10;
  supernode *s;

  const map<int, supernode* >& _supernodes = slice->getSupernodes();
  ulong nSupernodes = slice->getNbSupernodes();  

  // check if memory was already allocated for believes
  if(!believes) {
    believes = new double*[nSupernodes];
    for (uint k = 0; k < nSupernodes; ++k) {
      believes[k] = new double[param->nClasses];
    }
  }

  double maxPotential;
  computeNodePotentials(believes, maxPotential);
    
  double scale = 1.0;
  if (fabs(maxPotential) > 1e-30) {
    scale = MAX_POTENTIAL/maxPotential;
  }
  INFERENCE_PRINT("[gi_MF] maxPotential=%g, scale=%g\n", maxPotential, scale);

#if EXP_DOMAIN
  INFERENCE_PRINT("[gi_MF] Computing in exp domain\n");
  // go to exponential domain
  for (uint sid = 0; sid < nSupernodes; ++sid) {
    for(int c = 0; c < param->nClasses; ++c) {
      believes[sid][c] = std::exp(believes[sid][c]*scale);
    }
  }
#else
  INFERENCE_PRINT("[gi_MF] Computing in log domain\n");
  // stay in log domain
  for (uint sid = 0; sid < nSupernodes; ++sid) {
    for(int c = 0; c < param->nClasses; ++c) {
      believes[sid][c] *= scale;
    }
  }
#endif

  //exportBelieves("believes0");

  for(map<int, supernode* >::const_iterator its = _supernodes.begin();
      its != _supernodes.end(); its++) {
    sid = its->first;
    inferredLabels[sid] = 0;
  }

  // random selection
  for(uint iter = 0; iter < maxiter && (totalScore - totalScore_old) > 1.0; ++iter) {
    
    printf("[GI_MF] Iteration %d/%ld\n", iter, maxiter);

    totalScore_old = totalScore;
    totalScore = 0;

    // Select a pixel at random
    //int sid = gsl_rng_uniform_int(r, nSupernodes);

    int nLabelsChanged = 0;
    double* bs = 0;
    double maxBelief = 0;

    for(map<int, supernode* >::const_iterator its = _supernodes.begin();
        its != _supernodes.end(); its++) {
      sid = its->first;
      s = its->second;
      bs = believes[sid];

      if(param->nClasses != 2) {

        for(int c = 0; c < (int)param->nClasses; c++) {
          double unaryPotential = computeUnaryPotential(slice, sid, c) * scale;

          double pairwiseBelief = 0;
          if(param->includeLocalEdges) {
            vector < supernode* >* lNeighbors = &(s->neighbors);
            for(vector < supernode* >::iterator itN = lNeighbors->begin();
                itN != lNeighbors->end(); itN++) {

              // set edges once
              if(its->first < (*itN)->id) {
                continue;
              }

#if USE_LONG_RANGE_EDGES
               double pairwisePotential = computePairwisePotential_distance(slice, s, (*itN),
                                                                            c, inferredLabels[(*itN)->id]);
#else
               double pairwisePotential = computePairwisePotential(slice, s, (*itN),
                                                                   c, inferredLabels[(*itN)->id]);
#endif
               pairwisePotential *= scale;
               pairwiseBelief += believes[(*itN)->id][c] * pairwisePotential;
            }
          }
#if EXP_DOMAIN
        bs[c] = exp(unaryPotential)*exp(pairwiseBelief);
#else
        bs[c] = unaryPotential + pairwiseBelief;
#endif
          if(maxBelief < bs[c]) {
            maxBelief = bs[c];
          }
        }
      } else {
        int c = T_FOREGROUND;
#if EXP_DOMAIN
        bs[T_FOREGROUND] = 1;
#else
        bs[T_FOREGROUND] = 0;
#endif

        c = T_BACKGROUND;
        double unaryPotential = computeUnaryPotential(slice, sid, c) * scale;
        //printf("unaryPotential %d %g\n", sid, unaryPotential);

        double pairwiseBelief = 0;
        if(param->includeLocalEdges) {
          vector < supernode* >* lNeighbors = &(s->neighbors);
          for(vector < supernode* >::iterator itN = lNeighbors->begin();
              itN != lNeighbors->end(); itN++) {

            // set edges once
            if(its->first < (*itN)->id) {
              continue;
            }

#if USE_LONG_RANGE_EDGES
            double pairwisePotential = computePairwisePotential_distance(slice, s, (*itN),
                                                                         c, inferredLabels[(*itN)->id]);
#else
            double pairwisePotential = computePairwisePotential(slice, s, (*itN),
                                                                c, inferredLabels[(*itN)->id]);
#endif
            pairwisePotential *= scale;
            pairwiseBelief += believes[(*itN)->id][c] * pairwisePotential;
            //printf("pairwisePotential %d %d %d %g %g %g\n", sid, (*itN)->id, c, pairwisePotential, believes[(*itN)->id][c], pairwiseBelief);
          }
        }
#if EXP_DOMAIN
        bs[c] = exp(unaryPotential)*exp(pairwiseBelief);
#else
        bs[c] = unaryPotential + pairwiseBelief;
#endif
        //printf("bs[%d] %g\n",c,bs[c]);
        if(maxBelief < bs[c]) {
          maxBelief = bs[c];
        }
      }

      if(useLossFunction) {
        // loss function
        for(int c = 0; c < (int)param->nClasses; c++) {
          if(c != groundTruthLabels[sid]) {
            double _loss = lossPerLabel[groundTruthLabels[sid]] * scale;
            // add loss of the ground truth label
#if EXP_DOMAIN
            bs[c] *= exp(_loss);
#else
            bs[c] += _loss;
#endif

            if(maxBelief < bs[c]) {
              maxBelief = bs[c];
            }
          }
        }
      }

      // pick max
      maxScore = bs[inferredLabels[sid]];
      for(int c = 0; c < (int)param->nClasses; c++) {
        if(!replaceVoidMSRC || (c != voidLabel && c != moutainLabel && c != horseLabel)) {
          if(maxScore < bs[c]) {
            maxScore = bs[c];
            inferredLabels[sid] = c;
            ++nLabelsChanged;
          }
        }
      }
      totalScore += maxScore;
    }
    printf("[GI_MF] Iteration %d/%ld. Total score = %g. nLabelsChanged = %d\n", iter, maxiter,
           totalScore, nLabelsChanged);

    if(maxBelief > MAX_POTENTIAL) {
      // normalize believes
      printf("[GI_MF] Normalizing believes. maxBelief = %g\n", maxBelief);
      for (uint sid = 0; sid < nSupernodes; ++sid) {
        for(int c = 0; c < param->nClasses; ++c) {
          believes[sid][c] /= maxBelief;
        }
      }
    }
    //exportBelieves("believes_last");

  }

  // cleaning
  return computeEnergy(inferredLabels);
}

void GI_MF::exportBelieves(const char* filename)
{
  ofstream ofs(filename);
  ulong nSupernodes = slice->getNbSupernodes();
  for (uint sid = 0; sid < nSupernodes; ++sid) {
    ofs << sid;
    for(int c = 0; c < param->nClasses; ++c) {
      ofs << " " << believes[sid][c];
    }
    ofs << endl;
  }
  ofs.close();
}
