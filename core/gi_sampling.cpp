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

#include "gi_sampling.h"

// SliceMe
#include "Config.h"
#include "utils.h"

#include "inference_globals.h"

#if USE_GSL_DEBUG
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#else
#ifdef _WIN32
#include "gettimeofday.h"
#else
#include <sys/time.h>
#include <unistd.h>
#endif
#endif

//------------------------------------------------------------------------------

#if USE_GSL_DEBUG
#define GI_SAMPLING_GET_RAND_UNIFORM gsl_rng_uniform(rng);
#define GI_SAMPLING_GET_RAND(n) gsl_rng_uniform_int(rng, n);
#else
#define GI_SAMPLING_GET_RAND_UNIFORM rand() / (double)RAND_MAX;
#define GI_SAMPLING_GET_RAND(n) rand() * ((double)n/(double)RAND_MAX);
#endif


//------------------------------------------------------------------------------

GI_sampling::GI_sampling(Slice_P* _slice, 
                         const EnergyParam* _param,
                         double* _smw,
                         labelType* _groundTruthLabels,
                         double* _lossPerLabel,
                         Feature* _feature,
                         map<sidType, nodeCoeffType>* _nodeCoeffs,
                         double _sampling_rate)
{
  GraphInference::init();
  slice = _slice;
  param = _param;
  smw = _smw;
  lossPerLabel = _lossPerLabel;
  groundTruthLabels = _groundTruthLabels;
  feature = _feature;
  nodeCoeffs = _nodeCoeffs;
  sampling_rate = _sampling_rate;

  initializedLabels = false;

  bool useLossFunction = lossPerLabel!=0;
  string paramMSRC;
  Config::Instance()->getParameter("msrc", paramMSRC);
  bool useMSRC = paramMSRC.c_str()[0] == '1';
  replaceVoidMSRC = false;
  voidLabel = 0;
  moutainLabel = 0;
  horseLabel = 0;
  if(!useLossFunction && useMSRC) {
    Config::Instance()->getParameter("msrc_replace_void", paramMSRC);
    replaceVoidMSRC = paramMSRC.c_str()[0] == '1';
    voidLabel = classIdxToLabel[0];
    moutainLabel = classIdxToLabel[4161600];
    horseLabel = classIdxToLabel[8323328];
    printf("[gi_sampling] MSRC void=%d, moutain=%d, horse=%d\n",
           (int)voidLabel, (int)moutainLabel, (int)horseLabel);
  }
}

double GI_sampling::run(labelType* inferredLabels,
                        int id,
                        size_t maxiter,
                        labelType* nodeLabelsGroundTruth,
                        bool computeEnergyAtEachIteration,
                        double* _loss)
{
  // higher temperature means more randomness
  double temperature = 0.001;
  string config_tmp;
  if(Config::Instance()->getParameter("sampling_initial_temperature", config_tmp)) {
    temperature = atof(config_tmp.c_str()); 
    printf("[gi_sampling] sampling_initial_temperature = %g\n", temperature);
  }

  return runOnce(inferredLabels, maxiter, nodeLabelsGroundTruth, _loss,
                 temperature);
}

double GI_sampling::run_VOC(labelType* inferredLabels,
                            int id,
                            size_t maxiter,
                            labelType* nodeLabelsGroundTruth,
                            ulong* TPs, ulong* FPs, ulong* FNs)
{
  // higher temperature means more randomness
  double temperature = 0.001;
  string config_tmp;
  if(Config::Instance()->getParameter("sampling_initial_temperature", config_tmp)) {
    temperature = atof(config_tmp.c_str()); 
  }

  return runOnce_VOC(inferredLabels, maxiter, nodeLabelsGroundTruth,
                     temperature, TPs, FPs, FNs);
}

double GI_sampling::runOnce(labelType* inferredLabels,
                            size_t maxiter,
                            labelType* nodeLabelsGroundTruth,
                            double* _loss,
                            double temperature)
{
  double *buf = new double[param->nClasses];
  double *bufUnary = new double[param->nClasses];
  double *bufPairwise = new double[param->nClasses];
  double *bufLoss = new double[param->nClasses];
  int sid = 0;
  double totalScore_old = 0;
  double totalScore = 10; // different from 0 for first loop
  double totalUnaryScore = 0;
  double totalPairwiseScore = 0;
  double totalLoss = 0;
  supernode *s;
  bool useLossFunction = lossPerLabel!=0;

  // allocate memory to store features
  int fvSize = feature->getSizeFeatureVector();
  int max_index = fvSize + 1;
  osvm_node* n = new osvm_node[max_index];
  int i = 0;
  for(i = 0;i < max_index-1; i++)
    n[i].index = i+1;
  n[i].index = -1;

  const map<int, supernode* >& _supernodes = slice->getSupernodes();

  if(!initializedLabels) {
    if(nodeLabelsGroundTruth) {
      printf("[gi_sampling] Initializing labels to groundtruth\n");
      for(map<int, supernode* >::const_iterator its = _supernodes.begin();
          its != _supernodes.end(); its++) {
        sid = its->first;
        inferredLabels[sid] = nodeLabelsGroundTruth[sid];
      }
    } else {
      printf("[gi_sampling] Initializing labels to 0\n");
      for(map<int, supernode* >::const_iterator its = _supernodes.begin();
          its != _supernodes.end(); its++) {
        sid = its->first;
        inferredLabels[sid] = 0;
      }
    }
  }

  double *probs = new double[param->nClasses];
  double *cumulated_probs = new double[param->nClasses];

#if USE_GSL_DEBUG
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
#else
  //srand(time(NULL));

#ifndef NOTIME
  struct timeval _t;
  gettimeofday(&_t, NULL);
  srand(_t.tv_usec);
#endif

#endif


  int maxIter = 1;
  ulong nSupernodes = slice->getNbSupernodes();
  for(int iter = 0; iter < maxIter && (totalScore - totalScore_old) > 1.0; ++iter) {
    
    totalScore_old = totalScore;
    totalScore = 0;
    totalUnaryScore = 0;
    totalPairwiseScore = 0;
    totalLoss = 0;

    bool draw_samples = sampling_rate != 1.0;
    for(int i = 0; i < nSupernodes*sampling_rate; ++i) {

      if(draw_samples) {
        // Select a pixel at random
        sid = GI_SAMPLING_GET_RAND(nSupernodes);
      } else {
        sid = i;
      }
      s = slice->getSupernode(sid);

      vector < supernode* >* lNeighbors = &(s->neighbors);

      if(param->nClasses != 2) {

        for(int c = 0; c < (int)param->nClasses; c++) {
          bufUnary[c] = computeUnaryPotential(slice, sid, c);

          // add pairwise potential
          bufPairwise[c] = 0;
          for(vector < supernode* >::iterator itN = lNeighbors->begin();
              itN != lNeighbors->end(); itN++) {

            // set edges once
            if(s->id < (*itN)->id) {
              continue;
            }

#if USE_LONG_RANGE_EDGES
            double pairwisePotential = computePairwisePotential_distance(slice, s, (*itN),
                                                                         c, inferredLabels[(*itN)->id]);

#else
            double pairwisePotential = computePairwisePotential(slice, s, (*itN),
                                                                c, inferredLabels[(*itN)->id]);
#endif

            bufPairwise[c] += pairwisePotential;
          }
        }
      } else {
        bufUnary[T_FOREGROUND] = 0;
        bufPairwise[T_FOREGROUND] = 0;

        int c = T_BACKGROUND;
        bufUnary[c] = computeUnaryPotential(slice, sid, c);
        //printf("unary %d %d %g\n", sid, c, buf[c]);

        // add pairwise potential
        bufPairwise[c] = 0;
        for(vector < supernode* >::iterator itN = lNeighbors->begin();
            itN != lNeighbors->end(); itN++) {

          // set edges once
          if(s->id < (*itN)->id) {
            continue;
          }

#if USE_LONG_RANGE_EDGES
          double pairwisePotential = computePairwisePotential_distance(slice, s, (*itN),
                                                                       c, inferredLabels[(*itN)->id]);
#else
          double pairwisePotential = computePairwisePotential(slice, s, (*itN),
                                                              c, inferredLabels[(*itN)->id]);
#endif
          bufPairwise[c] += pairwisePotential;
          //printf("pairwise %d, %d %d %g\n", sid, (*itN)->id, c, pairwisePotential);
        }
      }

      if(useLossFunction) {
        for(int c = 0; c < (int)param->nClasses; c++) {
          if(c != groundTruthLabels[sid]) {
            // add loss of the ground truth label
            //buf[c] += lossPerLabel[groundTruthLabels[sid]];
            bufLoss[c] = lossPerLabel[groundTruthLabels[sid]];
          } else {
            bufLoss[c] = 0;
          }
        }
      }

      // compute posterior probability for each class
      double Z = 0;
      for(int c = 0; c < param->nClasses; ++c) {
        buf[c] = bufUnary[c] + bufPairwise[c] + bufLoss[c];
        probs[c] = std::exp(buf[c]/temperature);
        Z += probs[c];
      }

      // normalize probabilities
      if(fabs(Z) > 1e-30 && !isinf(Z)) {
        for(int c = 0; c < param->nClasses; ++c) {
          probs[c] /= Z;
        }
      }

      // sort probabilities and cumulate them
      double cumulated_prob = 0;
      for(int c = 0; c < param->nClasses; ++c) {
        cumulated_probs[c] = cumulated_prob;
        cumulated_prob += probs[c];
        //printf("cumulated_prob %g %g\n",probs[c], cumulated_prob);
      }

      bool labelSet = false;
      do {
        double rand_n = GI_SAMPLING_GET_RAND_UNIFORM;

        //printf("rand 0:(%g,%g) 1:(%g,%g) %g %g\n",buf[0],p0,buf[1],p1,prob,rand_n);
        int label = 0;
        int counter = param->nClasses - 1;
        while(counter >= 0) {
          //printf("it %d %g %g %d %g\n", it->first, it->second, rand_n, counter, cumulated_probs[counter]);
          if(rand_n > cumulated_probs[counter]) {
            label = counter;
            totalScore += buf[counter]; // include loss
            totalUnaryScore += bufUnary[counter];
            totalPairwiseScore += bufPairwise[counter];
            totalLoss += bufLoss[counter];
            break;
          }
          --counter;
        }

        if(!replaceVoidMSRC || (label != voidLabel && label != moutainLabel && label != horseLabel)) {
          inferredLabels[sid] = label;
          labelSet = true;
        }
      } while(!labelSet);

    }

    // Decrease temperature?
    // temperature = temperature*0.1;

    printf("[gi_sampling] Iteration %d/%d. Scores for sampled nodes: unaryScore = %g pairwiseScore = %g loss = %g score = %g\n", iter, maxIter,
           totalUnaryScore, totalPairwiseScore, totalLoss, totalScore);
  }

#if USE_GSL_DEBUG
    gsl_rng_free(rng);
#endif
  
  delete[] n;

  delete[] buf;
  delete[] bufUnary;
  delete[] bufPairwise;
  delete[] bufLoss;

  delete[] probs;
  delete[] cumulated_probs; 

  return computeEnergy(inferredLabels);
}

double GI_sampling::runOnce_VOC(labelType* inferredLabels,
                                size_t maxiter,
                                labelType* nodeLabelsGroundTruth,
                                double temperature,
                                ulong* TPs, ulong* FPs, ulong* FNs)
{
  double *buf = new double[param->nClasses];
  int sid = 0;
  double totalScore_old = 0;
  double totalScore = 10;
  supernode *s;
  bool useLossFunction = lossPerLabel!=0;

  // allocate memory to store features
  int fvSize = feature->getSizeFeatureVector();
  int max_index = fvSize + 1;
  osvm_node* n = new osvm_node[max_index];
  int i = 0;
  for(i = 0;i < max_index-1; i++)
    n[i].index = i+1;
  n[i].index = -1;

  const map<int, supernode* >& _supernodes = slice->getSupernodes();

  if(!initializedLabels) {
    if(nodeLabelsGroundTruth) {
      printf("[gi_sampling] Initializing labels using the ground truth\n");
      for(map<int, supernode* >::const_iterator its = _supernodes.begin();
          its != _supernodes.end(); its++) {
        sid = its->first;
        inferredLabels[sid] = nodeLabelsGroundTruth[sid];
      }
    } else {
      printf("[gi_sampling] Initializing labels to 0\n");
      for(map<int, supernode* >::const_iterator its = _supernodes.begin();
          its != _supernodes.end(); its++) {
        sid = its->first;
        inferredLabels[sid] = 0;
      }
    }
  }

  string config_tmp;
  double sampling_rate = 1.0;
  if(Config::Instance()->getParameter("sampling_rate", config_tmp)) {
    sampling_rate = atof(config_tmp.c_str()); 
    printf("[gi_sampling] sampling_rate = %g\n", sampling_rate);
  }

  double *probs = new double[param->nClasses];
  double *cumulated_probs = new double[param->nClasses];

#if USE_GSL_DEBUG
  gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
#else
  //srand(time(NULL));

#ifndef NOTIME
  struct timeval _t;
  gettimeofday(&_t, NULL);
  srand(_t.tv_usec);
#endif
#endif

  int maxIter = 1;
  ulong nSupernodes = slice->getNbSupernodes();
  for(int iter = 0; iter < maxIter && (totalScore - totalScore_old) > 1.0; ++iter) {
    
    printf("[gi_sampling] Iteration %d/%d\n", iter, maxIter);

    totalScore_old = totalScore;
    totalScore = 0;

    bool draw_samples = sampling_rate != 1.0;
    for(int i = 0; i < nSupernodes*sampling_rate; ++i) {

      if(draw_samples) {
        // Select a pixel at random
        sid = GI_SAMPLING_GET_RAND(nSupernodes);
      } else {
        sid = i;
      }
      s = slice->getSupernode(sid);

      vector < supernode* >* lNeighbors = &(s->neighbors);

      if(param->nClasses != 2) {

        for(int c = 0; c < (int)param->nClasses; c++) {
          buf[c] = computeUnaryPotential(slice, sid, c);

          if(iter >= 0) {
            // add pairwise potential
            for(vector < supernode* >::iterator itN = lNeighbors->begin();
                itN != lNeighbors->end(); itN++) {
#if USE_LONG_RANGE_EDGES
              double pairwisePotential = computePairwisePotential_distance(slice, s, (*itN),
                                                                  c, inferredLabels[(*itN)->id]);
#else
              double pairwisePotential = computePairwisePotential(slice, s, (*itN),
                                                                           c, inferredLabels[(*itN)->id]);
#endif

              buf[c] += pairwisePotential;
            }
          }
        }
      } else {
        buf[T_FOREGROUND] = 0;

        int c = T_BACKGROUND;
        buf[c] = computeUnaryPotential(slice, sid, c);

        if(iter >= 0) {
          // add pairwise potential
          for(vector < supernode* >::iterator itN = lNeighbors->begin();
              itN != lNeighbors->end(); itN++) {
#if USE_LONG_RANGE_EDGES
            double pairwisePotential = computePairwisePotential_distance(slice, s, (*itN),
                                                                c, inferredLabels[(*itN)->id]);
#else
            double pairwisePotential = computePairwisePotential(slice, s, (*itN),
                                                                c, inferredLabels[(*itN)->id]);
#endif
            buf[c] += pairwisePotential;
          }
        }
      }

      // debugging
      //printf("sid %d %d %d\n",sid,inferredLabels[sid],param->nClasses);
      assert(inferredLabels[sid]<param->nClasses);

      if(useLossFunction) {
        int g = nodeLabelsGroundTruth[sid];
        for(int c = 0; c < (int)param->nClasses; c++) {

          // update results of previous prediction
          int p = (c > 0)?c-1:inferredLabels[sid];

          if(p == g) {
            --TPs[p];
          } else {
            // -1 false positive for class p
            --FPs[p];
            // -1 false negative for class g
            --FNs[g];
          }

          if(c == g) {
            ++TPs[c];
          } else {
            // +1 false positive for class c
            ++FPs[c];
            // +1 false negative for class g
            ++FNs[g];
          }

          double score = 0;
          for(int _c = 0; _c < (int)param->nClasses; _c++) {
            double d = TPs[_c] + FPs[_c] + FNs[_c];
            if(fabs(d) > 1e-10) {
              score += ((double)TPs[_c])/d;
              //printf("score += %ld/%g -> %g\n", TPs[_c], d, score);
            }
          }
          // add loss
          //printf("[gi_sampling] sid %d BEFORE class %d -> %ld %ld %ld %g %g\n",
          //       sid, c, TPs[c], FPs[c], FNs[c], score, buf[c]);
          buf[c] += 1.0 - (score/param->nClasses);
          //printf("[gi_sampling] sid %d AFTER class %d -> %g\n", sid, c, buf[c]);
        }
      }

      // compute posterior probability for each class
      double Z = 0;
      for(int c = 0; c < param->nClasses; ++c) {
        probs[c] = std::exp(buf[c]/temperature);
        Z += probs[c];
      }

      // normalize probabilities
      if(fabs(Z) > 1e-30 && !isinf(Z)) {
        for(int c = 0; c < param->nClasses; ++c) {
          probs[c] /= Z;
        }
      }

      // sort probabilities and cumulate them
      double cumulated_prob = 0;
      for(int c = 0; c < param->nClasses; ++c) {
        cumulated_probs[c] = cumulated_prob;
        cumulated_prob += probs[c];
        //printf("cumulated_prob %g %g\n",probs[c].second, cumulated_prob);
      }

      bool labelSet = false;
      int p = param->nClasses - 1;
      do {
        double rand_n = GI_SAMPLING_GET_RAND_UNIFORM;

        //printf("rand 0:(%g,%g) 1:(%g,%g) %g %g\n",buf[0],p0,buf[1],p1,prob,rand_n);
        int label = 0;
        int counter = param->nClasses - 1;
        while(counter >= 0) {
          //printf("it %d %g %g %d %g\n", it->first, it->second, rand_n, counter, cumulated_probs[counter]);
          if(rand_n > cumulated_probs[counter]) {
            label = counter;
            totalScore += buf[counter];
            break;
          }
          --counter;
        }

        if(!replaceVoidMSRC || (label != voidLabel && label != moutainLabel && label != horseLabel)) {
          inferredLabels[sid] = label;
          labelSet = true;
        }
      } while(!labelSet);
      
      if(useLossFunction) {
        // update stats for VOC score
        int g = nodeLabelsGroundTruth[sid];

        // update results of previous prediction
        if(p == g) {
          --TPs[p];
        } else {
          // -1 false positive for class p
          --FPs[p];
          // -1 false negative for class g
          --FNs[g];
        }

        if(inferredLabels[sid] == g) {
          ++TPs[inferredLabels[sid]];
        } else {
          // +1 false positive for class inferredLabels[sid]
          ++FPs[inferredLabels[sid]];
          // +1 false negative for class g
          ++FNs[g];
        }
      }

    }

    printf("[gi_sampling] Iteration %d/%d. Total score = %g\n", iter, maxIter,
           totalScore);
    }

#if USE_GSL_DEBUG
    gsl_rng_free(rng);
#endif
  
  delete[] n;
  delete[] probs;
  delete[] cumulated_probs;

  return computeEnergy(inferredLabels);
}
