////////////////////////////////////////////////////////////////////////
// This program is free software; you can redistribute it and/or       //
// modify it under the terms of the GNU General Public License         //
// version 2 as published by the Free Software Foundation.             //
//                                                                     //
// This program is distributed in the hope that it will be useful, but //
// WITHOUT ANY WARRANTY; without even the implied warranty of          //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   //
// General Public License for more details.                            //
//                                                                     //
// Written and (C) by Aurelien Lucchi and Yunpeng Li                   //
// Contact <aurelien.lucchi@gmail.com> for comments & bug reports      //
/////////////////////////////////////////////////////////////////////////

/***********************************************************************/
/*                                                                     */
/*   svm_struct_api.c                                                  */
/*                                                                     */
/*   Definition of API for attaching implementing SVM learning of      */
/*   structures (e.g. parsing, multi-label classification, HMM)        */ 
/*                                                                     */
/*   Author: Thorsten Joachims                                         */
/*   Date: 03.07.04                                                    */
/*                                                                     */
/*   Copyright (c) 2004  Thorsten Joachims - All rights reserved       */
/*                                                                     */
/*   This software is available for non-commercial use only. It must   */
/*   not be modified and distributed without prior permission of the   */
/*   author. The author is not responsible for implications from the   */
/*   use of this software.                                             */
/*                                                                     */
/***********************************************************************/

#include "svm_struct/svm_struct_common.h"
#include "svm_struct_api.h"
#include "svm_struct_globals.h"
#include "label_cache.h"

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <map>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>
#ifdef _WIN32
#include <io.h>
#include <direct.h>
#define mkdir(x,y) _mkdir(x)
#else
#include <unistd.h>
#endif
#include <omp.h>

// inference
#include "inference.h"
#include "inference_globals.h"
#include "graphInference.h"
#include "gi_ICM.h"
#include "gi_max.h"
#include "gi_MF.h"
#include "gi_sampling.h"

#if USE_LIBDAI
#include "gi_libDAI.h"
using namespace dai;
#endif
#if USE_MAXFLOW
#include "gi_maxflow.h"
#endif
#if USE_MULTIOBJ
#include "gi_multiobject.h"
#endif

// SliceMe
#include "Config.h"
#include "Feature.h"
#include "F_Combo.h"
#include "F_LoadFromFile.h"
#include "F_Precomputed.h"
#include "Slice.h"
#include "utils.h"
#include "utils_ITK.h"

using namespace std;

#define fileExtension_Inference "png"
#define SEPARATOR '\t'

struct search_item {
  supernode* s;
  int distance;
};

//-----------------------------------------------------------------------GLOBALS

// global variable used to store the total loss at every iteration
double totalLoss = 0;
double totalLossMVC = 0;
double totalLhsXw = 0;
double totalWPsi = 0;
double totalWPsiGT = 0;

// Maximum number of threads you will be running.
// TODO : Make sure this is a correct number !!!
const int maxBuffers = 100;

// Memory buffer used to run inference after each iteration.
// maxBuffers buffers will be allocated
labelType** tempNodeLabels = 0;

ulong** TPs = 0;
ulong** FPs = 0;
ulong** FNs = 0;

// Memory buffer used to store the potentials
// Array of size maxBuffers*nNodes*nClasses 
double*** tempPotentials = 0;

// Optional : add label names in this vector if you want to see them printed in the log file
vector<string> labelNames;

// Directories containing images for training/test and validation set
string trainingDir;      // directory containing training images
string testDir;       // directory containing test images
string validationDir; // directory containing validation images

// Directory containg superpixel labels (dat files)
string labelsDir;
string maskTrainingDir;
string maskTestDir;
string maskValidationDir;

// save C value and compare it with the new one after each iteration.
// If the old and the new values are different, images are saved in a
// directory and the score is computed and appended to a file finalScore.txt 
double oldC;

// hack to remove last 3 classes for MSRC database
// set in read_struct_examples if useMSRC si set to 1
int NB_CLASSES_TO_REMOVE = 0;

EXAMPLE  *test_examples;
EXAMPLE  *validation_examples;

long nTestExamples;
uint maxNbTestNodes;
long nValidationExamples;
uint maxNbValidationNodes;

clock_t time_0;

bool useVOC = false;
bool useMSRC = false;
bool useSlice3d = true;
bool use01Loss = false;
bool generateFirstConstraint = false;
bool useGCForSubModularEnergy = true;
bool predictTrainingImages = true;
int nParallelChains = 12;

// output directories
string mostViolatedConstraintDir = "mostViolatedConstraint0";
string inferenceDir_training = "inference_training0";
string inferenceDir_test = "inference_test0";
string inferenceDir_validation = "inference_validation0";
string parameterDir = "parameter_vector0/";
string groundtruthDir = "groundTruth/";
string scoreDir = "scores0/";

//---------------------------------------------------------------------FUNCTIONS

/**
 * Compute loss function defined on nodes only.
 * You can also define the loss as a function of edges
 * (to get a smoother segmentation results for example)
 * @param y = ground truth
 * @param ybar = prediction
 */
void computeLoss(LABEL& y, labelType* ybar, const STRUCT_LEARN_PARM *sparm,
                 double& loss, int& nDiff)
{
  if(sparm->lossPerLabel == 0) {
    loss  = 0;
    nDiff = 0;
    return;
  }

  if(sparm->loss_function == LOSS_VOC) {
    computeVOCLoss(y, ybar, sparm->nClasses, loss, nDiff);
  } else {
    // use per node loss
    loss  = 0;
    nDiff = 0;

    // weight errors inversely proportional to the frequency of a class
    if(y.nodeCoeffs) {
      for(int n = 0; n < y.nNodes; n++) {
        if(y.nodeLabels[n] != ybar[n]) {
          if(sparm->loss_function == LOSS_NODE_BASED) {
          assert(sparm->lossPerLabel[n]>=0);
            loss += (*y.nodeCoeffs)[n]*sparm->lossPerLabel[n];
          } else {
            assert(sparm->lossPerLabel[y.nodeLabels[n]]>=0);
            loss += (*y.nodeCoeffs)[n]*sparm->lossPerLabel[y.nodeLabels[n]];
          }
          nDiff++;
        }
      }
    } else {
      for(int n = 0; n < y.nNodes; n++) {
        if(y.nodeLabels[n] != ybar[n]) {
          if(sparm->loss_function == LOSS_NODE_BASED) {
            assert(sparm->lossPerLabel[n]>=0);
            loss += sparm->lossPerLabel[n];
          } else {
            assert(sparm->lossPerLabel[y.nodeLabels[n]]>=0);
            loss += sparm->lossPerLabel[y.nodeLabels[n]];
          }
          nDiff++;
        }
      }
    }
  }
}

void computeLoss(labelType* y, labelType* ybar, int nNodes, const STRUCT_LEARN_PARM *sparm,
                 double& loss, int& nDiff)
{
  if(sparm->lossPerLabel == 0) {
    loss  = 0;
    nDiff = 0;
    return;
  }

  if(sparm->loss_function == LOSS_VOC) {
    computeVOCLoss(y, ybar, nNodes, sparm->nClasses, loss, nDiff, TPs[0], FPs[0], FNs[0]);
  } else {
    // use per node loss
    loss  = 0;
    nDiff = 0;

    for(int n = 0; n < nNodes; n++) {
      if(y[n] != ybar[n]) {
        if(sparm->loss_function == LOSS_NODE_BASED) {
          assert(sparm->lossPerLabel[n]>=0);
          loss += sparm->lossPerLabel[n];
        } else {
          assert(sparm->lossPerLabel[y[n]]>=0);
          loss += sparm->lossPerLabel[y[n]];
        }
        nDiff++;
      }
    }
  }
}

/**
 * Compute VOC score.
 * Note that SSVM can not minimize the VOC score.
 * This is only used for log purposes.
 * @param y = ground truth
 * @param ybar = prediction
 */
void computeVOCLoss(LABEL& y, labelType* ybar,
                    int nClasses,
                    double& loss, int& nDiff, SPATTERN x) {
  computeVOCLoss(y, ybar, nClasses, loss, nDiff, x.TPs, x.FPs, x.FNs);
}

void computeVOCLoss(LABEL& y, labelType* ybar,
                    int nClasses,
                    double& loss, int& nDiff) {
  computeVOCLoss(y, ybar, nClasses, loss, nDiff, TPs[0], FPs[0], FNs[0]);
}

void computeVOCLoss(LABEL& y, labelType* ybar,
                    int nClasses,
                    double& loss, int& nDiff,
                    ulong* TPs, ulong* FPs, ulong* FNs) {
  computeVOCLoss(y.nodeLabels, ybar, y.nNodes, nClasses, loss, nDiff, TPs, FPs,
                 FNs);
}

double computeScore(const STRUCTMODEL *sm, SWORD* fy)
{
  double score = 0;
  double* smw = sm->w + 1;

  SWORD* wfy = fy;
  for(int i = 0; i < sm->sizePsi; ++i) {
    score += smw[i]*wfy->weight;
    ++wfy;
  }

  return score;
}

/**
 * Print w vector
 */
void print_w(double* w, const STRUCT_LEARN_PARM *sparm,
             vector<string>& labelNames)
{
  SSVM_PRINT("[SVM_struct] w=[");

  //unary terms
  for(int c=0; c < sparm->nUnaryWeights; c++) {
    SSVM_PRINT("%d:%.2g ", c,w[c]);
  }

  for(int s=0; s < sparm->nScalingCoefficients; s++) {
    SSVM_PRINT("%d:%.2g ",SVM_FEAT_INDEX0(sparm) + s, w[SVM_FEAT_INDEX0(sparm)+s]);
  }
  SSVM_PRINT("\n");

  // pairwise terms
  int idx = sparm->nUnaryWeights;
  if(sparm->includeLocalEdges) {
    if(sparm->nGradientLevels == 0) {
      SSVM_PRINT("%d:%.2g ", idx,w[idx]);
      idx++;
    } else {
      SSVM_PRINT("\t");
      for(int c=0; c < sparm->nClasses; c++) {
	SSVM_PRINT("%d\t", c);
      }
      SSVM_PRINT("\n");
      
      for(int g=0; g < sparm->nGradientLevels; g++) {
	for(int o=0; o < sparm->nOrientations; o++) {
	  for(int c=0; c < sparm->nClasses; c++) {
	    SSVM_PRINT("%d", c);
	    for(int c2=0; c2 < sparm->nClasses; c2++) {
	      SSVM_PRINT("\t%.2g", w[idx]);
	      idx++;
	    }
	    SSVM_PRINT("\n");
	  }
	}
      }
    }
  }
}

/**
 * Compute score for a given image
 */
double computeVOCScore(SPATTERN x,
                       const string& maskDir,
                       map<ulong, labelType>& classIdxToLabel,           
                       int nClasses) {
  string groundtruthName =  maskDir + x.slice->getName();
  if(!fileExists(groundtruthName.c_str())) {
    // check if ground truth file for MSRC exists ?
    string baseName = getNameFromPathWithoutExtension(x.slice->getName());
    groundtruthName =  maskDir + baseName + "_GT.bmp";
    if(!fileExists(groundtruthName.c_str())) {
      printf("[SVM_struct] Could not find %s\n", groundtruthName.c_str());
      return -1;
    }
  } else {
    printf("[SVM_struct] Could not find %s\n", groundtruthName.c_str());
    return -1;
  }

  uchar* ptrM;
  uchar* ptrI;
  IplImage* img = (static_cast<Slice*>(x.slice))->img;
  IplImage* imgAnnotation = cvLoadImage(groundtruthName.c_str());
  if(!imgAnnotation) {
    printf("[SVM_struct] Error : input mask %s was not found\n", groundtruthName.c_str());
    exit(-1);
  }

  // Compute TP,FP and FN used to compute VOC score
  for(int c = 0; c < nClasses; c++) {
    x.TPs[c] = 0;
    x.FPs[c] = 0;
    x.FNs[c] = 0;
  }
  set<uint> lClasses;
  for(int u=0; u< img->width; u++) {
    for(int v=0; v< img->height; v++) {
      ptrM = &((uchar*)(img->imageData + img->widthStep*v))[u*img->nChannels];
      ptrI = &((uchar*)(imgAnnotation->imageData + imgAnnotation->widthStep*v))[u*imgAnnotation->nChannels];

      ulong classIdx_M = ptrM[0]*65025 + ptrM[1]*255 + ptrM[0];
      ulong classIdx_I = ptrI[0]*65025 + ptrI[1]*255 + ptrI[0];
      int labelM = classIdxToLabel[classIdx_M];
      int labelI = classIdxToLabel[classIdx_I];

      assert(labelI<nClasses);
      assert(labelM<nClasses);

      if(labelM == labelI) {
        // +1 true positive
        x.TPs[labelM]++;
      } else {
        // +1 false positive for predicted class
        x.FPs[labelI]++;
        // +1 false negative for groundtruth class
        x.FNs[labelM]++;
      }
      
      if(lClasses.count(labelI) == 0)
        lClasses.insert(labelI);
      if(lClasses.count(labelM) == 0)
        lClasses.insert(labelM);
    }
  }

  double score = 0;
  for(set<uint>::iterator it = lClasses.begin(); it != lClasses.end(); it++) {
    uint c = *it;
    float d = x.TPs[c]+x.FPs[c]+x.FNs[c];
    if(d>0) {
      score += x.TPs[c]/d;
    }
  }

  score /= lClasses.size();
  cvReleaseImage(&imgAnnotation);
  return score;
}

/**
 * Compute psi function. Caller is responsible for freeing memory.
 * psi contains elements for both nodes and edges
 */
SWORD* computePsi(SPATTERN x, LABEL y, const STRUCTMODEL *sm,
                 const STRUCT_LEARN_PARM *sparm,
                 double* _score = 0)
{
  SWORD* words = new SWORD[sm->sizePsi + 1];
  return computePsi(words, x, y, sm, sparm, _score);
}

SWORD* computePsi(SWORD* words, SPATTERN x, LABEL y, const STRUCTMODEL *sm,
                 const STRUCT_LEARN_PARM *sparm,
                 double* _score)
{
  double* smw = sm->w + 1;
  double* feats = new double[sm->sizePsi]; // 0-indexed

  // initialize feature vector
  for(int i = 0; i < sm->sizePsi; i++) {
    feats[i] = 0; 
  }

  int fvSize = x.feature->getSizeFeatureVector();

  // local nodes
  int label;
  int featIdx = 0;
  sidType sid = 0;
  const map<int, supernode* >& _supernodes = x.slice->getSupernodes();
  for(map<int, supernode* >::const_iterator itNode = _supernodes.begin();
      itNode != _supernodes.end(); itNode++) {
    // +1 for the appropriate label
    sid = itNode->first;
    label = y.nodeLabels[sid]; // label of the node
    
    if (sparm->nUnaryWeights == 1 && label == T_FOREGROUND) {
      // Only accumulates weights for BACKGROUND class.
      continue;
    }
    
#ifdef W_OFFSET
    if(x.nodeCoeffs) {
      feats[label] += (*x.nodeCoeffs)[sid];
    } else {
      feats[label]++;
    }
#endif
    
    osvm_node *n = x.slice->getFeature(sid);    

#if USE_SPARSE_VECTORS
    assert(0); // TODO...
#endif

    for(int s = 0; s < fvSize; s++) {
      featIdx = SVM_FEAT_INDEX(sparm, label, s);
      if(featIdx >= sm->sizePsi) {
        printf("[SVM_struct] featIdx>=sm->sizePsi %d %d %d %d %d %d %ld\n",featIdx,label,T_FOREGROUND,sparm->nUnaryWeights,s,fvSize,sm->sizePsi);
        exit(-1);
      }
      assert(s+1 == n[s].index);
      if(x.nodeCoeffs) {
        feats[featIdx] += (*x.nodeCoeffs)[sid]*n[s].value;
      } else {
        feats[featIdx] += n[s].value;
      }
    }

  }
  
  // edges
  if(sparm->includeLocalEdges) {
    if(sparm->nGradientLevels == 0) {
      edgeCoeffType edgeCoeff = 1.0;
      ulong edgeId = 0;
      // Only learn diagonal element.
      const map<int, supernode* >& _supernodes = x.slice->getSupernodes();
      for(map<int, supernode* >::const_iterator itNode = _supernodes.begin();
          itNode != _supernodes.end(); itNode++) {
        for(vector<supernode*>::iterator itNode2 = itNode->second->neighbors.begin();
            itNode2 != itNode->second->neighbors.end(); itNode2++) {
          // set edges once
          if(itNode->first < (*itNode2)->id) {
            continue;
          }

          if(x.edgeCoeffs) {
            edgeCoeff = (*x.edgeCoeffs)[edgeId];
          }

          if(y.nodeLabels[itNode->first] == y.nodeLabels[(*itNode2)->id]) {
            //sparm->nUnaryWeights is the offset due to unary terms
            feats[sparm->nUnaryWeights] += edgeCoeff;
          }

          ++edgeId;
        }
      }
    } else {
      // full pairwise model

#if USE_LONG_RANGE_EDGES
      int distanceIdx;
#endif
      int gradientIdx;
      int orientationIdx;
      int featIdx;
      ulong edgeId = 0;
      edgeCoeffType edgeCoeff = 1.0;
      const map<int, supernode* >& _supernodes = x.slice->getSupernodes();
      for(map<int, supernode* >::const_iterator itNode = _supernodes.begin();
          itNode != _supernodes.end(); itNode++) {
        for(vector<supernode*>::iterator itNode2 = itNode->second->neighbors.begin();
            itNode2 != itNode->second->neighbors.end(); itNode2++) {
          // set edges once
          if(itNode->first < (*itNode2)->id) {
            continue;
          }

          if(x.edgeCoeffs) {
            edgeCoeff = (*x.edgeCoeffs)[edgeId];
          }

          gradientIdx = x.slice->getGradientIdx(itNode->first,
                                                (*itNode2)->id);
          orientationIdx = x.slice->getOrientationIdx(itNode->first,
                                                      (*itNode2)->id);

          int offset = (orientationIdx*sparm->nClasses*sparm->nClasses);

#if USE_LONG_RANGE_EDGES
          distanceIdx = x.slice->getDistanceIdx(itNode->first,
                                                (*itNode2)->id);
          offset += distanceIdx*sparm->nGradientLevels*sparm->nClasses*sparm->nClasses*sparm->nOrientations;
#endif

          if(sparm->nUnaryWeights < 3) {
            // symmetric case : add +0.5 to both indices
            // sparm->nUnaryWeights is the offset due to unary terms
            for(int i = 0; i <= gradientIdx; i++)  {
              featIdx = (i*sparm->nClasses*sparm->nClasses*sparm->nOrientations) + offset + y.nodeLabels[itNode->first]*sparm->nClasses + y.nodeLabels[(*itNode2)->id];
              feats[featIdx+sparm->nUnaryWeights] += edgeCoeff/2.0;

              featIdx = (i*sparm->nClasses*sparm->nClasses*sparm->nOrientations) + offset + y.nodeLabels[itNode->first] + y.nodeLabels[(*itNode2)->id]*sparm->nClasses;
              feats[featIdx+sparm->nUnaryWeights] += edgeCoeff/2.0;
            }
          } else {
            for(int i = 0; i <= gradientIdx; i++)  {
              featIdx = (i*sparm->nClasses*sparm->nClasses*sparm->nOrientations) + offset + y.nodeLabels[itNode->first]*sparm->nClasses + y.nodeLabels[(*itNode2)->id];
              //sparm->nUnaryWeights is the offset due to unary terms
              feats[featIdx+sparm->nUnaryWeights] += edgeCoeff;
            }
          }

          ++edgeId;
        }
      }
    }
  }

  // make 'words' from the feat vector
  for(int i = 0; i < sm->sizePsi; i++) {
    words[i].wnum = i + 1;
    words[i].weight = feats[i];
  }
  words[sm->sizePsi].wnum = 0;  // termination symbol
  words[sm->sizePsi].weight = 0;

  // unary
  double scoreU = 0;
  for(int i = 0; i < sparm->nUnaryWeights; i++) {
    scoreU += smw[i]*feats[i];
  }

  for(int s = 0; s < sparm->nScalingCoefficients; s++) {
    scoreU += smw[SVM_FEAT_INDEX0(sparm) + s]*feats[SVM_FEAT_INDEX0(sparm) + s];
  }

  // pairwise
  double scoreP = 0;
  int sizePsi_pairwise = sparm->nClasses*sparm->nClasses;

  if(!sparm->includeLocalEdges) {
    sizePsi_pairwise = 0;
  } else {
    if(sparm->nGradientLevels == 0) {
        sizePsi_pairwise = 1;
    } else {
      sizePsi_pairwise *= sparm->nGradientLevels;
      sizePsi_pairwise *= sparm->nOrientations;
#if USE_LONG_RANGE_EDGES
      sizePsi_pairwise *= sparm->nDistances;
#endif
    }
  }

  for(int i=sparm->nUnaryWeights; i<sparm->nUnaryWeights+sizePsi_pairwise; i++) {
    scoreP += smw[i]*feats[i];
  }

  SSVM_PRINT("[SVM_struct]::computePsi ScoreU=%g, ScoreP=%g, Score=%g (This should be equal to -Energy)\n",
         scoreU, scoreP, scoreU + scoreP);

  if(_score) {
    *_score = 0;
    for(int i=0; i<sm->sizePsi; i++) {
      *_score += smw[i]*feats[i];
    }    
  }

  delete[] feats;
  return words;
}

/**
 * Initialize loss function
 * The loss for each label is stored in a global variable 'lossPerLabel'
 * This function has to be called after loading all the images
 */
void initLossFunction(EXAMPLE  *examples, const long nExamples,
                      double*& lossPerLabel,
                      int nClasses,
                      int loss_function,
                      STRUCT_LEARN_PARM *sparm)
{
  if(sparm->loss_function == LOSS_NODE_BASED) {
    SSVM_PRINT("[SVM_struct] Initializing node-based loss function\n");
    initLossFunction_nodeBased(examples, nExamples, lossPerLabel, nClasses,
                                loss_function, sparm);
  } else {
    SSVM_PRINT("[SVM_struct] Initializing class-based loss function\n");
    initLossFunction_classBased(examples, nExamples, lossPerLabel, nClasses,
                                loss_function, sparm);
  }
}

void computeDistance(Slice_P* slice_p, supernode* s, double* distanceMap, int starting_label, int max_distance)
{
  // add supernode s to a queue of supernodes to be visited
  deque<search_item> pending_list;
  set<sidType> visited;
  search_item si;
  si.s = s;
  si.distance = 0;
  pending_list.push_front(si);

  // iterate over all supernodes in the queue
  while(!pending_list.empty()) {

    // retrieve supernode and associated distance from s
    search_item sqi = pending_list.front();
    supernode* sq = sqi.s;
    int dist = sqi.distance;

    // remove sq from the queue
    pending_list.pop_front();

    if(sq->id != s->id) {
      if(distanceMap[sq->id] > dist) {
        distanceMap[sq->id] = dist;
      }
    } else {
      distanceMap[sq->id] = 0;
    }

    if(dist < max_distance) {
      // add neighbors to the queue
      // iterate over neighbors of sq and add them to the queue if there are not already in it
      for(vector<supernode*>::iterator itNq = sq->neighbors.begin();
          itNq != sq->neighbors.end(); ++itNq) {
        supernode* snq = *itNq;

          set<sidType>::iterator lookup_visited = visited.find(snq->id);
          if(lookup_visited == visited.end()) {
            int dist_snq = dist;
            ++dist_snq;
            search_item snqi;
            snqi.s = snq;
            snqi.distance = dist_snq;
            pending_list.push_back(snqi);
            visited.insert(sq->id);
          }
      }
    }
  }
}

void initLossFunction_nodeBased(EXAMPLE  *examples, const long nExamples,
                                 double*& lossPerLabel,
                                 int nClasses,
                                 int loss_function,
                                 STRUCT_LEARN_PARM *sparm)
{
  assert(nExamples == 1);

  sparm->lossScale = 1.0;

  Slice_P* slice = examples[0].x.slice;
  labelType* nodeLabels = examples[0].y.nodeLabels;

  // create distance map
  const int starting_label = BOUNDARY;
  ulong nSupernodes = slice->getNbSupernodes();
  lossPerLabel = new double[nSupernodes];
  for(ulong n = 0; n < nSupernodes; ++n) {
    lossPerLabel[n] = (double)MAX_DISTANCE_DT;
  }

  const map<sidType, supernode* >& _supernodes = slice->getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    supernode* s = it->second;
        
    if(nodeLabels[it->first] == starting_label) {
      computeDistance(slice, s, lossPerLabel, starting_label, MAX_DISTANCE_DT);
    }
  }

  // loss = distance + 1
  for(ulong i = 0; i < nSupernodes; ++i) {
    ++lossPerLabel[i];
  }
}

void initLossFunction_classBased(EXAMPLE  *examples, const long nExamples,
                                 double*& lossPerLabel,
                                 int nClasses,
                                 int loss_function,
                                 STRUCT_LEARN_PARM *sparm)
{
  // FIXME : delete is missing
  lossPerLabel = (double *)my_malloc(sizeof(double)*nClasses);
  ulong *countLabels = (ulong *)my_malloc(sizeof(ulong)*nClasses);
  for(int i = 0; i < nClasses; i++) {
    countLabels[i] = 0;
  }

  sparm->lossScale = 1.0;

  if(loss_function == 0 && !use01Loss) {
    // count labels in all images
    int label;
    for(long idx = 0; idx < nExamples; idx++) {
      const map<int, supernode* >& _supernodes = examples[idx].x.slice->getSupernodes();
      for(map<int, supernode* >::const_iterator itNode = _supernodes.begin();
	  itNode != _supernodes.end(); itNode++) {
	label = examples[idx].y.nodeLabels[itNode->first];
	countLabels[label]++;
      }
    }

    // normalize per class
    double total = 0;
    for(int i = 0; i < nClasses; i++) {
      if(countLabels[i] != 0) {
	lossPerLabel[i] = 1.0/countLabels[i];
	total += lossPerLabel[i];
      } else {
	lossPerLabel[i] = 0; //1
      }
    }

    for(int i = 0; i < nClasses; i++) {
      lossPerLabel[i] /= total;
    }
  } else {
    if(loss_function == 2) {
      // use VOC loss
      for(int i = 0; i < nClasses; i++){
        lossPerLabel[i] = 0;
      }
    } else {
      // 0-1 loss function
      for(int i = 0; i < nClasses; i++){
        lossPerLabel[i] = 1;
      }
    }
  }

  assert(nClasses>NB_CLASSES_TO_REMOVE);

  // no penalty for the last 3 classes of the MSRC database
  for(int i = nClasses-NB_CLASSES_TO_REMOVE; i < nClasses; i++) {
    lossPerLabel[i]  = 0;
  }

  SSVM_PRINT("[SVM_struct] countLabels:");
  for(int i = 0; i < nClasses; i++) {
    SSVM_PRINT(" %ld", countLabels[i]);
  }
  SSVM_PRINT("\n");
  SSVM_PRINT("[SVM_struct] LossPerLabel:");
  for(int i = 0; i < nClasses; i++) {
    SSVM_PRINT(" %f", lossPerLabel[i]);
  }
  SSVM_PRINT("\n");

  ofstream ofsLossPerLabel("lossPerLabel.txt");
  for(int i = 0; i < nClasses; i++)
    ofsLossPerLabel << lossPerLabel[i] << endl;
  ofsLossPerLabel.close();

  // Comuting maximum possible loss for each image
  for(long idx = 0; idx < nExamples; idx++) {
    double loss  = 0;
    for(int n = 0; n < examples[idx].y.nNodes; n++) {    
      loss += lossPerLabel[examples[idx].y.nodeLabels[n]];
    }
    SSVM_PRINT("[SVM_struct] Maximum possible loss for image %ld = %g\n", idx, loss);
  }
  free(countLabels);
}

void sparmToEnergyParam(const STRUCT_LEARN_PARM& sparm, EnergyParam* param)
{
  param->nUnaryWeights = sparm.nUnaryWeights;
  param->nClasses = sparm.nClasses;
  param->nDistances = sparm.nDistances;
  param->nGradientLevels = sparm.nGradientLevels;
  param->nOrientations = sparm.nOrientations;
  param->nScalingCoefficients = sparm.nScalingCoefficients;
  param->includeLocalEdges = sparm.includeLocalEdges;
  param->nScales = sparm.nScales;

#if USE_LONG_RANGE_EDGES
  param->sizePsi = (sparm.nDistances*sparm.nGradientLevels*sparm.nOrientations*sparm.nClasses*sparm.nClasses) + (sparm.nGradientLevels==0 && sparm.includeLocalEdges)
    + sparm.nUnaryWeights + sparm.nScalingCoefficients;
#else
  param->sizePsi = (sparm.nGradientLevels*sparm.nOrientations*sparm.nClasses*sparm.nClasses) + (sparm.nGradientLevels==0 && sparm.includeLocalEdges)
    + sparm.nUnaryWeights + sparm.nScalingCoefficients;
#endif

  param->weights = 0; // not used in graphInference class
}

void energyParamToSparm(const EnergyParam& param, STRUCT_LEARN_PARM* sparm)
{
  sparm->nUnaryWeights = param.nUnaryWeights;
  sparm->nClasses = param.nClasses;
  sparm->nDistances = param.nDistances;
  sparm->nGradientLevels = param.nGradientLevels;
  sparm->nOrientations = param.nOrientations;
  sparm->nScalingCoefficients = param.nScalingCoefficients;
  sparm->includeLocalEdges = param.includeLocalEdges;
  sparm->nScales = param.nScales;
}

void runInferenceWithCurrentParameters(SPATTERN x, LABEL y, STRUCTMODEL *sm,
                                       STRUCT_LEARN_PARM *sparm, int threadId)
{
  SSVM_PRINT("------------------------------------------------INFERENCE-%d-%d-START\n",
             x.id, sparm->iterationId);

  EnergyParam param;
  sparmToEnergyParam(*sparm, &param);

  // Run inference with current parameters (ignore loss function)
  GraphInference* gi_Inference;
  switch(sparm->giType)
    {
#if USE_LIBDAI
    case T_GI_LIBDAI:
      gi_Inference = new GI_libDAI(x.slice,
                                   &param,
                                   sm->w+1,
                                   y.nodeLabels, // groundtruth labels used to compute loss  
                                   0, // no loss function
                                   x.feature,
                                   x.nodeCoeffs, 0);
      gi_Inference->run(tempNodeLabels[threadId], //store node computed by BP in a dummy vector
                        x.id,
                        INFERENCE_MAX_ITER);

      break;
#endif

    default:
      printf("[svm_struct] Unknown algorithm type\n");
      exit(-1);
      break;
    }

  delete gi_Inference;

#if VERBOSITY > 1
  stringstream soutColoredImage;

  // create directory
  stringstream iterDir;
  if(useSlice3d) {
    iterDir  << "inference_training/";
  } else {
    iterDir  << "inference_training_" << sparm->iterationId << "/";
  }
  mkdir(iterDir.str().c_str(), 0777);
  soutColoredImage << iterDir.str();

  if(useSlice3d) {
    soutColoredImage << getNameFromPathWithoutExtension(x.slice->getName());
    soutColoredImage << "_";
    soutColoredImage << sparm->iterationId;
  } else {
    soutColoredImage << x.slice->getName();
  }

  x.slice->exportSupernodeLabels(soutColoredImage.str().c_str(),
                                 sparm->nClasses,
                                 tempNodeLabels[threadId],
                                 y.nNodes,
                                 &(sparm->labelToClassIdx));

  if(useSlice3d) {
    zipAndDeleteCube(soutColoredImage.str().c_str());
  }

#endif

  // Compute loss for infered image
  // FIXME : Change name : GT=groundtruth
  double loss  = 0;
  int nDiff = 0;
  computeLoss(y,tempNodeLabels[threadId],sparm,loss,nDiff);

  stringstream sLoss;
  sLoss << "lossGT_" << x.id << ".txt";
  ofstream ofs(sLoss.str().c_str(), ios::app);
  ofs << loss << " " << nDiff << endl;
  ofs.close();

  if (useVOC) {
    double lossVOC  = 0;
    int nDiffVOC = 0;
    computeVOCLoss(y,tempNodeLabels[threadId],sparm->nClasses,lossVOC,nDiffVOC,x);

    stringstream sLossVOC;
    sLossVOC << "lossGT_VOC_" << x.id << ".txt";
    ofstream ofsVOC(sLossVOC.str().c_str(), ios::app);
    ofsVOC << lossVOC << endl;
    ofsVOC.close();

    SSVM_PRINT("[Inference] Loss = %g , LossVOC = %g (%d nodes have different labels)\n",
               loss, lossVOC, nDiff);
  } else {
    SSVM_PRINT("[Inference] Loss = %g (%d nodes have different labels)\n",
               loss, nDiff);
  }

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  totalLoss += loss;
  SSVM_PRINT("------------------------------------------------INFERENCE-%d-%d-END\n",
             x.id, sparm->iterationId);
}

void        svm_struct_learn_api_init(int argc, char* argv[])
{
  /* Called in learning part before anything else is done to allow
     any initializations that might be necessary. */
}

void        svm_struct_learn_api_exit()
{
  /* Called in learning part at the very end to allow any clean-up
     that might be necessary. */
}

void        svm_struct_classify_api_init(int argc, char* argv[])
{
  /* Called in prediction part before anything else is done to allow
     any initializations that might be necessary. */
}

void        svm_struct_classify_api_exit()
{
  /* Called in prediction part at the very end to allow any clean-up
     that might be necessary. */
}

void load_3d_dataset(string imageDir,
		     string maskDir,
		     STRUCT_LEARN_PARM *sparm,
                     EXAMPLE*& examples,
                     long* nExamples,
                     uint* maxNbNodes,
		     int* featureSize,
                     Config* config)
                     
{
  *nExamples = 1;
  examples = (EXAMPLE *)my_malloc(sizeof(EXAMPLE)*(*nExamples));

  Slice3d* slice3d = new Slice3d(imageDir.c_str());
  if(!slice3d) {
    printf("[SVM_struct] Error while loading %s\n", imageDir.c_str());
    exit(-1);
  }

  string config_tmp;
  bool rescale_raw_data = false;
  if(Config::Instance()->getParameter("rescale_raw_data", config_tmp)) {
    rescale_raw_data = config_tmp[0] == '1';
  }
  if(rescale_raw_data) {
    slice3d->rescaleRawData();

#if USE_ITK
    exportTIFCube(slice3d->raw_data,
                  "rescaled_data.tif",
                  slice3d->depth,
                  slice3d->height,
                  slice3d->width);
#endif

  }

  slice3d->loadSupervoxels(imageDir.c_str());

#if USE_LONG_RANGE_EDGES
  slice3d->addLongRangeEdges_supernodeBased(sparm->nDistances);
#endif

  // Load features
  vector<eFeatureType> feature_types;
  int paramFeatureTypes = DEFAULT_FEATURE_TYPE;
  if(config->getParameter("featureTypes", config_tmp)) {
    paramFeatureTypes = atoi(config_tmp.c_str());
  }
  getFeatureTypes(paramFeatureTypes, feature_types);

  stringstream sout_feature_filename;
  sout_feature_filename << slice3d->inputDir << "/features_";
  sout_feature_filename << slice3d->getSupernodeStep() << "_" << slice3d->getCubeness();
  sout_feature_filename << "_" << paramFeatureTypes;
  sout_feature_filename << "_" << DEFAULT_FEATURE_DISTANCE;

  printf("[SVM_struct] Checking %s\n", sout_feature_filename.str().c_str());
  Feature* feature = 0;
  bool featuresLoaded = false;
  if(fileExists(sout_feature_filename.str())) {
    printf("[SVM_struct] Loading features from %s\n", sout_feature_filename.str().c_str());
    *featureSize = -1;
    if(slice3d->loadFeatures(sout_feature_filename.str().c_str(), featureSize)) {
      featuresLoaded = true;
      feature = new F_Precomputed(slice3d->getPrecomputedFeatures(), *featureSize/DEFAULT_FEATURE_DISTANCE);
      printf("[SVM_struct] Features Loaded succesfully\n");
    } else {
      printf("[SVM_struct] Features not loaded succesfully\n");
    }
  } else {
    printf("[SVM_struct] File %s does not exist\n", sout_feature_filename.str().c_str());
  }

  if(!featuresLoaded) {
    feature = Feature::getFeature(slice3d, feature_types);

    *featureSize = feature->getSizeFeatureVector();
    SSVM_PRINT("[SVM_struct] Feature size = %d\n", *featureSize);
    slice3d->precomputeFeatures(feature);

#if VERBOSITY > 1
    // Dump features
    if(feature) {
      feature->save(*slice3d, sout_feature_filename.str().c_str());
    }
#endif

  }

  // precompute gradient indices to avoid race conditions
  slice3d->precomputeGradientIndices(sparm->nGradientLevels);
  slice3d->precomputeOrientationIndices(sparm->nOrientations);


  SSVM_PRINT("[SVM_struct] Slice3d instantiated. %ld nodes. %ld edges\n",
             slice3d->getNbSupernodes(),
             slice3d->getNbEdges());
  slice3d->printEdgeStats();
  slice3d->printDistanceIndicesCount(sparm->nDistances);

#if VERBOSITY > 1
  // Print one feature vector
  const int sid_to_print = 100;
  osvm_node* n = slice3d->getFeature(sid_to_print);

  SSVM_PRINT("[SVM_struct] Feature(%d):\n", sid_to_print);
  for(int s = 0; s < *featureSize; s++) {
    SSVM_PRINT("%d:%g ", n[s].index, n[s].value);
  }
  SSVM_PRINT("\n");
#endif

  bool update_loss_function = false;
  if(config->getParameter("update_loss_function", config_tmp)) {
    update_loss_function = config_tmp.c_str()[0] == '1';
  }
  printf("[SVM_struct] update_loss_function=%d\n", (int)update_loss_function);

  const int idx = 0;
  examples[idx].x.id = idx;
  examples[idx].x.slice = slice3d;
  examples[idx].x.feature = feature;
  if(update_loss_function) {
    int nNodes = slice3d->getNbSupernodes();
    map<sidType, nodeCoeffType>* _nodeCoeffs = new map<sidType, nodeCoeffType>;
    for(int n = 0; n < nNodes; ++n) {
      (*_nodeCoeffs)[n] = 1.0;
    }
    examples[idx].x.nodeCoeffs = _nodeCoeffs;
    examples[idx].y.nodeCoeffs = _nodeCoeffs;

    int nEdges = slice3d->getNbEdges();
    map<sidType, nodeCoeffType>* _edgeCoeffs = new map<sidType, edgeCoeffType>;
    for(int e = 0; e < nEdges; ++e) {
      (*_edgeCoeffs)[e] = 1.0;
    }
    examples[idx].x.edgeCoeffs = _edgeCoeffs;

  } else {
    examples[idx].x.nodeCoeffs = 0;
    examples[idx].y.nodeCoeffs = 0;
    examples[idx].x.edgeCoeffs = 0;
  }

  examples[idx].x.imgAnnotation = 0; // no ground-truth image

  switch(sparm->metric_type) {
    case METRIC_SUPERNODE_BASED_01:
      examples[idx].x.cubeAnnotation = 0;
      break;
    case METRIC_NODE_BASED_01:
      {
        // load all the mask images
        int _width;
        int _height;
        int _depth = -1;
        loadFromDir(maskDir.c_str(), examples[idx].x.cubeAnnotation,
                    _width, _height, &_depth);
      }
      break;
  default:
    break;
  }

  examples[idx].x.nEdges = slice3d->getNbEdges();

  // load ground truth
  bool includeBoundaryLabels = true;
  if(config->getParameter("includeBoundaryLabels", config_tmp)) {
    includeBoundaryLabels = config_tmp.c_str()[0] == '1';
  }
  SSVM_PRINT("[SVM_struct] includeBoundaryLabels=%d\n", (int)includeBoundaryLabels);
  bool includeUnknownType = false;
  slice3d->setIncludeOtherLabel(false);
  slice3d->generateSupernodeLabelFromMaskDirectory(maskDir.c_str(),
                                                   includeBoundaryLabels,
                                                   includeUnknownType);
  SSVM_PRINT("[SVM_struct] Ground truth volume has been loaded\n");

  if(slice3d->isSupernodeLabelsLoaded()) {

    examples[idx].y.nNodes = slice3d->getNbSupernodes();
    examples[idx].y.nodeLabels = new labelType[examples[idx].y.nNodes];
    examples[idx].y.cachedNodeLabels = false;

    if(*maxNbNodes < (uint)examples[idx].y.nNodes) {
      *maxNbNodes = (uint)examples[idx].y.nNodes;
    }

    SSVM_PRINT("[SVM_struct] Loading ground truth labels in memory\n");

    int label;
    int sid = 0;
    const map<int, supernode* >& _supernodes = slice3d->getSupernodes();
    for(map<int, supernode* >::const_iterator its = _supernodes.begin();
        its != _supernodes.end(); its++) {
      label = its->second->getLabel();
      examples[idx].y.nodeLabels[sid] = label;
      sid++;
    }

#if VERBOSITY > 2
    // Save ground truth labels
    stringstream soutColoredImageGT;
    soutColoredImageGT << groundtruthDir;
    soutColoredImageGT << slice3d->getName();
    SSVM_PRINT("[SVM_struct] Saving ground truth labels to %s\n",
               soutColoredImageGT.str().c_str());

    slice3d->exportSupernodeLabels(soutColoredImageGT.str().c_str(),
                                   sparm->nClasses,
                                   examples[idx].y.nodeLabels,
                                   examples[idx].y.nNodes,
                                   &(sparm->labelToClassIdx));

    stringstream soutOverlayImageGT;
    soutOverlayImageGT << groundtruthDir;
    soutOverlayImageGT << getNameFromPathWithoutExtension(slice3d->getName());
    soutOverlayImageGT << "_overlay.png";
    PRINT_MESSAGE("[SVM_struct] Saving ground truth overlay to %s\n", soutOverlayImageGT.str().c_str());
    slice3d->exportOverlay(soutOverlayImageGT.str().c_str(), examples[idx].y.nodeLabels);

#endif
  } else {
    examples[idx].y.nNodes = 0;
    examples[idx].y.nodeLabels = 0;
    examples[idx].y.cachedNodeLabels = false;
  }

  examples[idx].x.TPs = new ulong[sparm->nClasses];
  examples[idx].x.FPs = new ulong[sparm->nClasses];
  examples[idx].x.FNs = new ulong[sparm->nClasses];
  examples[idx].x.count = 0;

#if VERBOSITY > 1
  {
    // Print one feature vector
    const int sid_to_print = 100;
    osvm_node* n = slice3d->getFeature(sid_to_print);
      
    SSVM_PRINT("[SVM_struct] Feature(%d):\n", sid_to_print);
    for(int s = 0; s < *featureSize; s++) {
      SSVM_PRINT("%d:%g ", n[s].index, n[s].value);
    }
    SSVM_PRINT("\n");
  }
#endif
}


void load_2d_dataset(string imageDir,
                     int nImages,
                     string maskDir,
                     STRUCT_LEARN_PARM *sparm,
                     EXAMPLE*& examples,
                     long* nExamples,
                     uint* maxNbNodes,
		     int* featureSize,
                     Config* config)                     
{
  // load images
  vector<string> lFiles;
  int nFiles = 0;
  if( (strcmp(sparm->fileExtension, "png")==0) || (strcmp(sparm->fileExtension, "bmp")==0)
      || (strcmp(sparm->fileExtension, "jpg")==0) || (strcmp(sparm->fileExtension, "tif")==0)) {
    SSVM_PRINT("[SVM_struct] Listing files in %s\n", imageDir.c_str());
    getFilesInDir(imageDir.c_str(), lFiles, sparm->fileExtension);
    nFiles = lFiles.size();
    SSVM_PRINT("[SVM_struct] Found %d files in %s\n", nFiles, imageDir.c_str());
  } else {
    printf("[SVM_struct] Unknown file extension %s\n", sparm->fileExtension);
    exit(-1);
  }

  int maxNumberImages = 0;
  string config_tmp;
  if(config->getParameter("maxNumberImages", config_tmp)) {
    maxNumberImages = atoi(config_tmp.c_str());
  }
  SSVM_PRINT("[SVM_struct] maxNumberImages=%d\n", maxNumberImages);

  if(maxNumberImages != 0) {
    nFiles = min(maxNumberImages, nFiles);
  }

  if(nImages != -1) {
    nFiles = min(nImages, nFiles);
  }

  *nExamples = nFiles;
  examples = (EXAMPLE *)my_malloc(sizeof(EXAMPLE)*nFiles);

  int idx = 0;
  for(int i = 0; i < nFiles; i++) {
    string file = lFiles[i];
    string baseName = getNameFromPathWithoutExtension(file);

    string groundtruthName;
    bool found_GT = getGroundTruthName(groundtruthName, maskDir, file);
    if(!found_GT) {
      printf("[SVM_struct] Ground truth file associated to %s not found in directory %s\n",
             file.c_str(), maskDir.c_str());
    }

    // load image and superpixels
    stringstream imageName;
    imageName << imageDir << baseName << "." << sparm->fileExtension;
    stringstream superpixelLabels;

    superpixelLabels << labelsDir;
    superpixelLabels << baseName;
    superpixelLabels << ".dat";

    if(labelsDir != "" && !fileExists(superpixelLabels.str().c_str())) {
      SSVM_PRINT("[SVM_struct] Superpixel file %s not found\n", superpixelLabels.str().c_str());
      //continue;
    }

    SSVM_PRINT("[SVM_struct] Loading %d-th image file %s\n", idx, imageName.str().c_str());
    Slice* slice = new Slice(imageName.str().c_str(),
                             superpixelLabels.str().c_str());
    // Load features
    vector<eFeatureType> feature_types;
    int paramFeatureTypes = DEFAULT_FEATURE_TYPE;
    if(config->getParameter("featureTypes", config_tmp)) {
      paramFeatureTypes = atoi(config_tmp.c_str());
    }
    getFeatureTypes(paramFeatureTypes, feature_types);
    Feature* feature = Feature::getFeature(slice, feature_types);
    if(paramFeatureTypes & F_LOADFROMFILE) {
      F_LoadFromFile* fLoadFromFile = 0;
      if(paramFeatureTypes != F_LOADFROMFILE) {
        F_Combo* feat_combo = (F_Combo*)feature;
        const vector<Feature*>& feat_list = feat_combo->getFeatures();
        for (vector<Feature*>::const_iterator it = feat_list.begin();
             it != feat_list.end(); ++it) {
          if ((*it)->getFeatureType() == F_LOADFROMFILE) {
            fLoadFromFile = (F_LoadFromFile*)(*it);
          }
        }
      } else {
        fLoadFromFile = (F_LoadFromFile*)feature;
      }

      string featureFile;
      if(!config->getParameter("feature_file", featureFile)) {
        printf("[SVM_struct] Error : path for feature file was not set\n");
        exit(-1);
      }
      if(!fileExists(featureFile) || (isDirectory(featureFile))) {
        SSVM_PRINT("[SVM_struct] featureFile %s does not exist or is a directory...\n", featureFile.c_str());
        featureFile += getLastDirectoryFromPath(imageDir);
        featureFile += "/";
      }
      printf("[SVM_struct] Loading featureFile=%s\n", featureFile.c_str());
      fLoadFromFile->init(*slice, featureFile.c_str());
    }
    *featureSize = feature->getSizeFeatureVector();
    SSVM_PRINT("[SVM_struct] Feature size = %d\n", *featureSize);

    SSVM_PRINT("[SVM_struct] %ld node predictions loaded. %ld edges\n",
               slice->mSupernodes.size(),
               slice->getNbEdges());

    bool update_loss_function = false;
    if(config->getParameter("update_loss_function", config_tmp)) {
      update_loss_function = config_tmp.c_str()[0] == '1';
    }
    SSVM_PRINT("[SVM_struct] update_loss_function=%d\n", (int)update_loss_function);

    // precompute gradient indices to avoid race conditions
    slice->precomputeGradientIndices(sparm->nGradientLevels);
    slice->precomputeOrientationIndices(sparm->nOrientations);
    slice->precomputeFeatures(feature);
#if USE_LONG_RANGE_EDGES
    slice->addLongRangeEdges_supernodeBased(sparm->nDistances);
    //slice->precomputeDistanceIndices(sparm->nDistances);
#endif

    examples[idx].x.id = idx;
    examples[idx].x.slice = slice;
    examples[idx].x.feature = feature;
    examples[idx].x.cubeAnnotation = 0;
    if(update_loss_function) {
      int nNodes = slice->getNbSupernodes();
      map<sidType, nodeCoeffType>* _nodeCoeffs = new map<sidType, nodeCoeffType>;
      for(int n = 0; n < nNodes; ++n) {
        (*_nodeCoeffs)[n] = 1.0;
      }
      examples[idx].x.nodeCoeffs = _nodeCoeffs;
      examples[idx].y.nodeCoeffs = _nodeCoeffs;

      int nEdges = slice->getNbEdges();
      map<sidType, nodeCoeffType>* _edgeCoeffs = new map<sidType, edgeCoeffType>;
      for(int e = 0; e < nEdges; ++e) {
        (*_edgeCoeffs)[e] = 1.0;
      }
      examples[idx].x.edgeCoeffs = _edgeCoeffs;


    } else {
      examples[idx].x.nodeCoeffs = 0;
      examples[idx].y.nodeCoeffs = 0;
      examples[idx].x.edgeCoeffs = 0;
    }

    examples[idx].x.nEdges = slice->getNbEdges();

    if(found_GT) {
      // load labels from ground truth
      Slice* sliceGT = new Slice(imageName.str().c_str(),
                                 superpixelLabels.str().c_str());

      bool GT_TextFile = getExtension(groundtruthName) == "labels";
      if(GT_TextFile) {
        sliceGT->generateSupernodeLabelFromTextFile(groundtruthName.c_str(),
                                                    sparm->nClasses);
        examples[idx].x.imgAnnotation = 0; // no ground-truth image
      } else {
        bool includeBoundaryLabels = false;
        if(config->getParameter("includeBoundaryLabels", config_tmp)) {
          includeBoundaryLabels = config_tmp.c_str()[0] == '1';
        }
        SSVM_PRINT("[SVM_struct] includeBoundaryLabels=%d\n", (int)includeBoundaryLabels);

        if(sparm->classIdxToLabel.size() < 3 || includeBoundaryLabels) {
          SSVM_PRINT("[SVM_struct] Generating supernode labels from mask image\n");
          sliceGT->generateSupernodeLabelFromMaskImage(groundtruthName.c_str(), includeBoundaryLabels);
        } else {
          SSVM_PRINT("[SVM_struct] Generating supernode labels from multiclass image %s containing %ld labels\n",
                     groundtruthName.c_str(), sparm->classIdxToLabel.size());
          sliceGT->generateSupernodeLabelsFromMultiClassMaskImage(groundtruthName.c_str(),
                                                                  sparm->classIdxToLabel);
        }

        IplImage* imgAnnotation = cvLoadImage(groundtruthName.c_str());
        if(!imgAnnotation) {
          printf("[SVM_struct] Error : input mask %s was not found\n", groundtruthName.c_str());
          exit(-1);
        }
        examples[idx].x.imgAnnotation = imgAnnotation;
      }

      examples[idx].y.nNodes = sliceGT->getNbSupernodes();
      examples[idx].y.nodeLabels = new labelType[examples[idx].y.nNodes];
      examples[idx].y.cachedNodeLabels = false;

      if(*maxNbNodes < (uint)examples[idx].y.nNodes) {
        *maxNbNodes = examples[idx].y.nNodes;
      }

      int label;
      int sid = 0;
      const map<int, supernode* >& _supernodes = sliceGT->getSupernodes();
      for(map<int, supernode* >::const_iterator its = _supernodes.begin();
          its != _supernodes.end(); its++) {
        label = its->second->getLabel();
        examples[idx].y.nodeLabels[sid] = label;
        sid++;
      }

#if VERBOSITY > 2
      // Save ground truth labels
      string soutColoredImageGT = groundtruthDir + file;
      SSVM_PRINT("[Ground truth] Saving labels to %s\n", soutColoredImageGT.c_str());
      examples[idx].x.slice->exportSupernodeLabels(soutColoredImageGT.c_str(),
                                                   sparm->nClasses,
                                                   examples[idx].y.nodeLabels,
                                                   examples[idx].y.nNodes,
                                                   &(sparm->labelToClassIdx));
#endif

      delete sliceGT;
    } else {
      examples[idx].x.imgAnnotation = 0; // no ground-truth image
      examples[idx].y.nNodes = 0;
      examples[idx].y.nodeLabels = 0;
      examples[idx].y.cachedNodeLabels = false;
    }

    examples[idx].x.TPs = new ulong[sparm->nClasses];
    examples[idx].x.FPs = new ulong[sparm->nClasses];
    examples[idx].x.FNs = new ulong[sparm->nClasses];
    examples[idx].x.count = new ulong[sparm->nClasses];

    idx++;
  }
}

void load_learn_parm(STRUCT_LEARN_PARM *sparm, Config* config)
{
  string config_tmp;

  sparm->nClasses = 3;
  if(config->getParameter("nClasses", config_tmp)) {
    sparm->nClasses = atoi(config_tmp.c_str());
  }
  SSVM_PRINT("[SVM_struct] sparm->nClasses=%d\n", sparm->nClasses);

  if (sparm->nClasses == 2) {
    // Special case for which we only need to weight one class vs the other one
    sparm->nUnaryWeights = 1;
  } else {
    sparm->nUnaryWeights = sparm->nClasses;
  }
  SSVM_PRINT("[SVM_struct] sparm->nUnaryWeights=%d\n", sparm->nUnaryWeights);

  string fileExtension;
  config->getParameter("fileExtension", fileExtension);
  strncpy(sparm->fileExtension,fileExtension.c_str(),MAX_SIZE_EXT);

  // nGradientLevels
  if(config->getParameter("nGradientLevels", config_tmp)) {
    sparm->nGradientLevels = atoi(config_tmp.c_str());
  } else {
    sparm->nGradientLevels = 5;
  }
  SSVM_PRINT("[SVM_struct] nGradientLevels=%d\n", sparm->nGradientLevels);

  // nDistances
  if(config->getParameter("nDistances", config_tmp)) {
    sparm->nDistances = atoi(config_tmp.c_str());
  } else {
    sparm->nDistances = 0;
  }

#if USE_LONG_RANGE_EDGES == 0
  sparm->nDistances = 0;
#endif

  SSVM_PRINT("[SVM_struct] nDistances=%d\n", sparm->nDistances);

  // maxSqEdgeDistance
  if(config->getParameter("maxSqEdgeDistance", config_tmp)) {
    sparm->maxSqEdgeDistance = atoi(config_tmp.c_str());
  } else {
    sparm->maxSqEdgeDistance = 500;
  }
  SSVM_PRINT("[SVM_struct] maxSqEdgeDistance=%d\n", sparm->maxSqEdgeDistance);

  // nOrientations
  if(config->getParameter("nOrientations", config_tmp)) {
    sparm->nOrientations = atoi(config_tmp.c_str());

    if(sparm->nOrientations != 1 && sparm->nOrientations != 3) {
      printf("[SVM_struct] Error: sparm->nOrientations != 1 && sparm->nOrientations != 3\n");
      exit(-1);
    }
  } else {
    sparm->nOrientations = 1;
  }
  SSVM_PRINT("[SVM_struct] nOrientations=%d\n", sparm->nOrientations);

  // includeLocalEdges  
  sparm->includeLocalEdges = true;
  if(config->getParameter("includeLocalEdges", config_tmp)) {
    // 1 means include edges
    if(atoi(config_tmp.c_str()) == 0) {
      sparm->includeLocalEdges = false;
      SSVM_PRINT("[SVM_struct] No local edges will be included\n");
    }
  }
  SSVM_PRINT("[SVM_struct] sparm->includeLocalEdges=%d\n", (int)sparm->includeLocalEdges);

  if(config->getParameter("giType", config_tmp)) {
    sparm->giType = atoi(config_tmp.c_str());
    SSVM_PRINT("[SVM_struct] giType = %d\n", sparm->giType);
  }

  if(sparm->includeLocalEdges == 0) {
    sparm->giType = T_GI_MAX;
  }

  if(sparm->giType == T_GI_MULTIOBJ && sparm->nClasses <=2) {
    sparm->giType = T_GI_LIBDAI;
  }

  if(sparm->giType == T_GI_MAX && sparm->includeLocalEdges == 1) {
    printf("[SVM_struct] Error: you specified giType == T_GI_MAX that can only be used for the D-model (graph without edges)\n");
    exit(-1);
  }

  // nMaxIterations
  if(config->getParameter("nMaxIterations", config_tmp)) {
      sparm->nMaxIterations = atoi(config_tmp.c_str());
    }
  else {
    sparm->nMaxIterations = 800;
  }
  SSVM_PRINT("[SVM_struct] nMaxIterations=%d\n", sparm->nMaxIterations);

  // stepForOutputFiles
  sparm->stepForOutputFiles = 20;
  if(config->getParameter("stepForOutputFiles", config_tmp)) {
    sparm->stepForOutputFiles = atoi(config_tmp.c_str());
    if(sparm->stepForOutputFiles == 0) {
      printf("[SVM_struct] Error: stepForOutputFiles=0\n");
      exit(-1);
    }
  }
  SSVM_PRINT("[SVM_struct] stepForOutputFiles=%d\n", sparm->stepForOutputFiles);

  // stepForParameterFiles
  sparm->stepForParameterFiles = 20;
  if(config->getParameter("stepForParameterFiles", config_tmp)) {
    sparm->stepForParameterFiles = atoi(config_tmp.c_str());
    if(sparm->stepForParameterFiles == 0) {
      printf("[SVM_struct] Error: stepForParameterFiles=0\n");
      exit(-1);
    }
  }
  SSVM_PRINT("[SVM_struct] stepForParameterFiles=%d\n", sparm->stepForParameterFiles);

  sparm->nScales = 1;
  if(config->getParameter("nScales", config_tmp)) {
    sparm->nScales = atoi(config_tmp.c_str());
  }
  SSVM_PRINT("[SVM_struct] nScales=%d\n", sparm->nScales);

  string colormapFilename;
  config->getParameter("colormapFilename", colormapFilename);
  SSVM_PRINT("[SVM_struct] colormapFilename=%s\n", colormapFilename.c_str());

  // load colormap information with more than 3 classes
  if(colormapFilename != "") {
    // Load colormap information
    SSVM_PRINT("[SVM_struct] Loading colormap %s\n", colormapFilename.c_str());
    ifstream ifsCol(colormapFilename.c_str());
    string line;
    int label = 0;
    int r,g,b;
    ulong classIdx;
    char labelName[50];
    while(getline(ifsCol, line)) {
      sscanf(line.c_str(),"%d %lu %d %d %d %s", &label, &classIdx,&r,&g,&b,labelName);
      sparm->classIdxToLabel[classIdx] = label;
      sparm->labelToClassIdx[(labelType)label] = classIdx;

      SSVM_PRINT("label=%d, classIdx=%lu, rgb=(%d,%d,%d)\n",
             label, classIdx, (int)r, (int)g, (int)b);
    }
    ifsCol.close();
  } else {
    SSVM_PRINT("[SVM_struct] Generating random colormap\n");
    srand(time(NULL));

    for(int label = 0; label < sparm->nClasses; ++label) {
      uchar r = (uchar)(rand()*255.0 / (double)RAND_MAX);
      uchar g = (uchar)(rand()*255.0 / (double)RAND_MAX);
      uchar b = (uchar)(rand()*255.0 / (double)RAND_MAX);
      ulong classIdx = RGBToclassIdx(r, g, b);
      sparm->classIdxToLabel[classIdx] = label;
      sparm->labelToClassIdx[(labelType)label] = classIdx;
    }
  }

  sparm->ssvm_iteration = 0;
  sparm->iterationId = 1;

  // DEPRECATED
  sparm->nLocalScales = 1;
  sparm->featureSize = 1;
}

/**
 * Format of the input file is the following :
 * nClasses
 * image directory
 * Mask directory
 * Prediction directory
 * Label directory
 * File extension
 * Colormap filename
 * nGlobalStates
 * nGradientLevels
 * nOrientations
 * includeLocalEdges
 * Test directory
 * Validation directory
 * File extensions for prediction files
 * nMaxIterations
 * stepForOutputFiles
 */
SAMPLE      read_struct_examples(char *file, STRUCT_LEARN_PARM *sparm)
{
  /* Reads struct examples and returns them in sample. The number of
     examples must be written into sample.n */
  SAMPLE   sample;  /* sample */
  EXAMPLE  *examples;
  string config_tmp;

  verbose = true;

  SSVM_PRINT("[SVM_struct] Loading examples from %s\n", file);

  Config* config = new Config(file);
  Config::setInstance(config);

  set_default_parameters(config);

  string paramVOC;
  if(config->getParameter("voc", paramVOC)) {
    useVOC = paramVOC.c_str()[0] == '1';
  } else {
    useVOC = false;
  }

  string paramMSRC;
  if(config->getParameter("msrc", paramMSRC)) {
    useMSRC = paramMSRC.c_str()[0] == '1';
  } else {
    useMSRC = false;
  }

  string paramSlice3d;
  if(config->getParameter("slice3d", paramSlice3d)) {
    useSlice3d = paramSlice3d.c_str()[0] == '1';
  } else {
    useSlice3d = true;
  }

  string param01Loss;
  if(config->getParameter("01loss", param01Loss)) {
    use01Loss = param01Loss.c_str()[0] == '1';
  } else {
    use01Loss = false;
  }
  SSVM_PRINT("[SVM_struct] use01Loss=%d\n", (int)use01Loss);

  if(!useVOC && !useMSRC && !useSlice3d) {
    printf("[SVM_struct] Set voc=1 XOR msrc=1 XOR slice3d=1 in your config file\n");
    exit(-1);
  }

  SSVM_PRINT("[SVM_struct] useMSRC=%d, useVOC=%d useSlice3d=%d\n",
             useMSRC, useVOC, useSlice3d);

  if(useMSRC) {
    NB_CLASSES_TO_REMOVE = 3;
  } else {
    NB_CLASSES_TO_REMOVE = 0;
  }
  SSVM_PRINT("[SVM_struct] NB_CLASSES_TO_REMOVE=%d\n", NB_CLASSES_TO_REMOVE);

  sparm->alg_type = DEFAULT_ALG_TYPE;
  if(Config::Instance()->getParameter("alg_type", config_tmp)) {
    sparm->alg_type = atoi(config_tmp.c_str());
  }
  SSVM_PRINT("[SVM_struct] alg_type=%d\n", sparm->alg_type);

  if(Config::Instance()->getParameter("obj_c", config_tmp)) {
    sparm->C = atof(config_tmp.c_str());
  }
  SSVM_PRINT("[SVM_struct] obj_c=%f\n", sparm->C);

  load_learn_parm(sparm, config);

  Config::Instance()->getParameter("trainingDir", trainingDir);
  Config::Instance()->getParameter("maskTrainingDir", maskTrainingDir);
  Config::Instance()->getParameter("labelsDir", labelsDir);


  // test directory
  if(Config::Instance()->getParameter("testDir", testDir)) {
    SSVM_PRINT("[SVM_struct] testDir=%s\n", testDir.c_str());
  }

  // validation directory
  if(Config::Instance()->getParameter("validationDir", validationDir)) {
    SSVM_PRINT("[SVM_struct] validationDir=%s\n", validationDir.c_str());
  }

  // mask test directory
  if(Config::Instance()->getParameter("maskTestDir", maskTestDir)) {
    SSVM_PRINT("[SVM_struct] maskTestDir=%s\n", maskTestDir.c_str());
  }

  // mask validation directory
  if(Config::Instance()->getParameter("maskValidationDir", maskValidationDir)) {
    SSVM_PRINT("[SVM_struct] validationDir=%s\n", maskValidationDir.c_str());
  }

  if(Config::Instance()->getParameter("generateFirstConstraint", config_tmp)) {
    if(atoi(config_tmp.c_str())) {
      generateFirstConstraint = true;
    }
  }

  if(Config::Instance()->getParameter("useGCForSubModularEnergy", config_tmp)) {
    if(atoi(config_tmp.c_str())) {
      useGCForSubModularEnergy = true;
    }
  }
  SSVM_PRINT("[SVM_struct] useGCForSubModularEnergy=%d\n", (int)useGCForSubModularEnergy);

  if(Config::Instance()->getParameter("sampling_nParallelChains", config_tmp)) {
    nParallelChains = atoi(config_tmp.c_str());
  }
  SSVM_PRINT("[SVM_struct] nParallelChains=%d\n", nParallelChains);


  if(Config::Instance()->getParameter("predictTrainingImages", config_tmp)) {
    if(atoi(config_tmp.c_str())==1) {
      predictTrainingImages = true;
    } else {
      predictTrainingImages = false;
    }
  }
  SSVM_PRINT("[SVM_struct] predictTrainingImages=%d\n", (int)predictTrainingImages);

  SSVM_PRINT("[SVM_struct] SVM_FEAT_INDEX0=%d\n", SVM_FEAT_INDEX0(sparm));

  if(useMSRC) {
    GraphInference::setColormap(sparm->classIdxToLabel);
  }

  sparm->metric_type = METRIC_NODE_BASED_01;
  if(Config::Instance()->getParameter("metric_type", config_tmp)) {
    sparm->metric_type = atoi(config_tmp.c_str());
  }
  SSVM_PRINT("[SVM_struct] metric_type = %d\n", sparm->metric_type);

  // feature type
  int paramFeatureTypes = DEFAULT_FEATURE_TYPE;
  vector<eFeatureType> feature_types;
  if(Config::Instance()->getParameter("featureTypes", config_tmp)) {
    paramFeatureTypes = atoi(config_tmp.c_str());
  }
  getFeatureTypes(paramFeatureTypes, feature_types);
  
  if(Config::Instance()->getParameter("loss_function", config_tmp)) {
    sparm->loss_function = atoi(config_tmp.c_str());
  }
  SSVM_PRINT("[SVM_struct] loss_function=%d\n", sparm->loss_function);

  if(feature_types.size() == 0) {
    printf("[SVM_struct] Error : no features specified in the config file\n");
    exit(-1);
  }

#if VERBOSITY > 2
  mkdir(groundtruthDir.c_str(), 0777);
#endif

  if(sparm->nClasses == 3) {
    printf("[svm_struct] Set class labels\n");
    BACKGROUND = 0;
    BOUNDARY = 1;
    FOREGROUND = 2;
  }

  long nExamples = 0;
  uint maxNbNodes = 0;

  if (useSlice3d) {
    int featureSize = 0;
    load_3d_dataset(trainingDir, maskTrainingDir, sparm, examples, &nExamples, &maxNbNodes, &featureSize, config);
    sparm->featureSize = featureSize;
    //if(testDir != "0" && isDirectory(testDir)) {
    if(testDir != "0" && !testDir.empty()) {
      load_3d_dataset(testDir, maskTestDir, sparm, test_examples, &nTestExamples, &maxNbTestNodes, &featureSize, config);
    } else {
      printf("[svm_struct] Test directory %s is not valid\n", testDir.c_str());
    }
  } else {
    int featureSize = 0;
    int nImages = -1; // load all
    load_2d_dataset(trainingDir, nImages, maskTrainingDir, sparm, examples, &nExamples, &maxNbNodes, &featureSize, config);
    sparm->featureSize = featureSize;
    if(isDirectory(testDir)) {
      load_2d_dataset(testDir, nImages, maskTestDir, sparm, test_examples, &nTestExamples, &maxNbTestNodes, &featureSize, config);
    }
    if(isDirectory(validationDir)) {
      load_2d_dataset(validationDir, nImages, maskValidationDir, sparm, validation_examples, &nValidationExamples, &maxNbValidationNodes, &featureSize, config);
    }
  }

  SSVM_PRINT("[SVM_struct] Loaded %ld examples\n", nExamples);

  sparm->nScalingCoefficients = sparm->nScales*sparm->nUnaryWeights*sparm->featureSize;

  SSVM_PRINT("[SVM_struct] %d scales detected including %d local scales and %d scaling coefficients\n",
         sparm->nScales,sparm->nLocalScales,sparm->nScalingCoefficients);

  SSVM_PRINT("[SVM_struct] giType = %d\n", sparm->giType);

  // rescale features
  bool rescale_features = true;
  if(Config::Instance()->getParameter("rescale_features", config_tmp)) {
    rescale_features = config_tmp.c_str()[0] == '1';
  }
  SSVM_PRINT("[SVM_struct] rescale_features %d\n", (int)rescale_features);
  if(rescale_features) {
    SSVM_PRINT("[SVM_struct] Rescaling precomputed features for %ld samples\n", nExamples);

    assert(useSlice3d); // not implemented for 2d images yet

    const char* scale_filename = "scale.txt";

    // training images
    for(int i = 0; i < nExamples; ++i) {
      examples[i].x.slice->rescalePrecomputedFeatures(scale_filename);

#if VERBOSITY > 1
      const int sid_to_print = 100;
      osvm_node* n = examples[i].x.slice->getFeature(sid_to_print);

      SSVM_PRINT("[SVM_struct] Feature(%d):\n", sid_to_print);
      for(int s = 0; s < sparm->featureSize; s++) {
        SSVM_PRINT("%d:%g ", n[s].index, n[s].value);
      }
      SSVM_PRINT("\n");
#endif
    }    

    // testing images
    for(int i = 0; i < nTestExamples; ++i) {
      test_examples[i].x.slice->rescalePrecomputedFeatures(scale_filename);

#if VERBOSITY > 1
      const int sid_to_print = 100;
      osvm_node* n = test_examples[i].x.slice->getFeature(sid_to_print);

      SSVM_PRINT("[SVM_struct] Feature(%d):\n", sid_to_print);
      for(int s = 0; s < sparm->featureSize; s++) {
        SSVM_PRINT("%d:%g ", n[s].index, n[s].value);
      }
      SSVM_PRINT("\n");
#endif
    }

  }

  // Allocate memory for max number of nodes.
  SSVM_PRINT("[SVM_struct] Allocating temporary memory to run inference after each iteration. maxBuffers=%d, maxNbNodes=%d\n", maxBuffers, maxNbNodes);
  tempNodeLabels = new labelType*[maxBuffers];
  for(int il = 0; il < maxBuffers; il++) {
    tempNodeLabels[il] = new labelType[maxNbNodes];
  }

  TPs = new ulong*[maxBuffers];
  FPs = new ulong*[maxBuffers];
  FNs = new ulong*[maxBuffers];
  for(int il = 0; il < maxBuffers; il++) {
    TPs[il] = new ulong[sparm->nClasses];
    FPs[il] = new ulong[sparm->nClasses];
    FNs[il] = new ulong[sparm->nClasses];
  }

  printf("[svm_struct] Initializing loss function %d\n", sparm->loss_function);
  initLossFunction(examples, nExamples,
                   sparm->lossPerLabel,
                   sparm->nClasses,
                   sparm->loss_function,
                   sparm);

  if(sparm->giType == T_GI_MF) {
    // Allocate memory for potentials
    SSVM_PRINT("[SVM_struct] Allocating temporary memory to store potentials. maxBuffers=%d, maxNbNodes=%d\n", maxBuffers, maxNbNodes);
    tempPotentials = new double**[maxBuffers];
    for(int il = 0; il < maxBuffers; il++) {
      tempPotentials[il] = new double*[maxNbNodes];
      for(uint n = 0; n < maxNbNodes; n++) {
        tempPotentials[il][n] = new double[sparm->nClasses];
      }
    }
  }

  oldC = sparm->C;
  sparm->startC = sparm->C;
  sample.n = nExamples;
  sample.examples = examples;
  sample.nTest = nTestExamples;
  sample.test_examples = test_examples;
  sample.nValidation = nValidationExamples;
  sample.validation_examples = validation_examples;

#ifdef USE_OPENMP
  int nThreads = NTHREADS;
  if(config->getParameter("nThreads", config_tmp)) {
    nThreads = atoi(config_tmp.c_str());
  }
  printf("[svm_struct] Setting number of threads to %d\n", nThreads);
  omp_set_num_threads(nThreads);  
#endif

#if USE_LONG_RANGE_EDGES
  printf("[svm_struct] USE_LONG_RANGE_EDGES is ON\n");
#endif

  time_0 = clock();

  return(sample);
}

void        init_struct_model(SAMPLE sample, STRUCTMODEL *sm, 
			      STRUCT_LEARN_PARM *sparm, LEARN_PARM *lparm, 
			      KERNEL_PARM *kparm)
{
  /* Initialize structmodel sm. The weight vector w does not need to be
     initialized, but you need to provide the maximum size of the
     feature space in sizePsi. This is the maximum number of different
     weights that can be learned. Later, the weight vector w will
     contain the learned weights for the model. */

  // Check size of the feature vector
  //sparm->featureSize = sample.examples[0].x.feature->getSizeFeatureVector();
  SSVM_PRINT("[SVM_struct] featureSize = %d\n", sparm->featureSize);

  sparm->nScalingCoefficients = sparm->nScales*sparm->nUnaryWeights*sparm->featureSize;

  // K = sparm->nClasses
  // L = sparm->nGlobalStates
  // sizePsi = K*K + K + 1
  // + 1 for weight of the unary cost from the RBF SVM.

#if USE_LONG_RANGE_EDGES
  sm->sizePsi = (sparm->nDistances*sparm->nGradientLevels*sparm->nOrientations*sparm->nClasses*sparm->nClasses) + (sparm->nGradientLevels==0 && sparm->includeLocalEdges)
    + sparm->nUnaryWeights + sparm->nScalingCoefficients;
#else
  sm->sizePsi = (sparm->nGradientLevels*sparm->nOrientations*sparm->nClasses*sparm->nClasses) + (sparm->nGradientLevels==0 && sparm->includeLocalEdges)
    + sparm->nUnaryWeights + sparm->nScalingCoefficients;
#endif

  SSVM_PRINT("[SVM_struct] %d scales detected including %d local scales and %d scaling coefficients\n",
         sparm->nScales,sparm->nLocalScales,sparm->nScalingCoefficients);
  SSVM_PRINT("[SVM_struct] sizePsi = %ld\n", sm->sizePsi);

  update_output_dir(sparm->ssvm_iteration);
}

void update_output_dir(const int iteration)
{
  // Update directory names
  mostViolatedConstraintDir = StringPrintf("mostViolatedConstraint", iteration)  + "/";
  inferenceDir_training = StringPrintf("inference_training", iteration) + "/";
  inferenceDir_test = StringPrintf("inference_test", iteration) + "/";
  inferenceDir_validation = StringPrintf("inference_validation", iteration) + "/";
  parameterDir = StringPrintf("parameter_vector", iteration) + "/";
  scoreDir = StringPrintf("scores", iteration)  + "/";

  mkdir(scoreDir.c_str(), 0777);
  mkdir(parameterDir.c_str(), 0777);
  mkdir(inferenceDir_training.c_str(), 0777);
  mkdir(inferenceDir_test.c_str(), 0777);
  mkdir(inferenceDir_validation.c_str(), 0777);
  mkdir(mostViolatedConstraintDir.c_str(), 0777);
}

CONSTSET    init_struct_constraints(SAMPLE sample, STRUCTMODEL *sm, 
				    STRUCT_LEARN_PARM *sparm)
{
  /* Initializes the optimization problem. Typically, you do not need
     to change this function, since you want to start with an empty
     set of constraints. However, if for example you have constraints
     that certain weights need to be positive, you might put that in
     here. The constraints are represented as lhs[i]*w >= rhs[i]. lhs
     is an array of feature vectors, rhs is an array of doubles. m is
     the number of constraints. The function returns the initial
     set of constraints. */

  CONSTSET c;
  long     sizePsi=sm->sizePsi;
  //long     i;
  SWORD     words[2];

  if(1) { /* normal case: start with empty set of constraints */
    c.lhs=NULL;
    c.rhs=NULL;
    c.m=0;
  }
  else { /* add constraints so that all learned weights are
            positive. WARNING: Currently, they are positive only up to
            precision epsilon set by -e. */
    c.lhs=(DOC **)my_malloc(sizeof(DOC *)*sizePsi);
    c.rhs=(double*)my_malloc(sizeof(double)*sizePsi);
    //for(i=0; i<sizePsi; i++) {
    words[0].wnum=SVM_FEAT_INDEX0(sparm) + 1;
    words[0].weight=1e4;
    words[1].wnum=0;
    /* the following slackid is a hack. we will run into problems,
       if we have more than 1000000 slack sets (ie examples) */
    c.lhs[0]=create_example(0,0,1000000,1,create_svector(words,(char*)"",1.0));
    c.rhs[0]=1e4;
    //}
    c.m = 0;
  }

  return(c);
}

LABEL       classify_struct_example(SPATTERN x, STRUCTMODEL *sm, 
				    STRUCT_LEARN_PARM *sparm)
{
  /* Finds the label yhat for pattern x that scores the highest
     according to the linear evaluation function in sm, especially the
     weights sm.w. The returned label is taken as the prediction of sm
     for the pattern x. The weights correspond to the features defined
     by psi() and range from index 1 to index sm->sizePsi. If the
     function cannot find a label, it shall return an empty label as
     recognized by the function empty_label(y). */
  LABEL y;

  /* insert your code for computing the predicted label y here */
  assert(0);

  return(y);
}

LABEL       find_most_violated_constraint_slackrescaling(SPATTERN x, LABEL y, 
                                                         const STRUCTMODEL *sm, 
                                                         const STRUCT_LEARN_PARM *sparm)
{
  /* Finds the label ybar for pattern x that that is responsible for
     the most violated constraint for the slack rescaling
     formulation. For linear slack variables, this is that label ybar
     that maximizes

            argmax_{ybar} loss(y,ybar)*(1-psi(x,y)+psi(x,ybar)) 

     Note that ybar may be equal to y (i.e. the max is 0), which is
     different from the algorithms described in
     [Tschantaridis/05]. Note that this argmax has to take into
     account the scoring function in sm, especially the weights sm.w,
     as well as the loss function, and whether linear or quadratic
     slacks are used. The weights in sm.w correspond to the features
     defined by psi() and range from index 1 to index
     sm->sizePsi. Most simple is the case of the zero/one loss
     function. For the zero/one loss, this function should return the
     highest scoring label ybar (which may be equal to the correct
     label y), or the second highest scoring label ybar, if
     Psi(x,ybar)>Psi(x,y)-1. If the function cannot find a label, it
     shall return an empty label as recognized by the function
     empty_label(y). */
  /* insert your code for computing the label ybar here */
  return find_most_violated_constraint_marginrescaling(x,y,sm,sparm);
}

LABEL       find_most_violated_constraint_marginrescaling(SPATTERN x, LABEL y, 
                                                          const STRUCTMODEL *sm, 
                                                          const STRUCT_LEARN_PARM *sparm)
{
  /* Finds the label ybar for pattern x that that is responsible for
     the most violated constraint for the margin rescaling
     formulation. For linear slack variables, this is that label ybar
     that maximizes

            argmax_{ybar} loss(y,ybar)+sm.w*psi(x,ybar)

            argmin_{ybar} Energy(X,ybar) - loss(y,ybar)

     Note that ybar may be equal to y (i.e. the max is 0), which is
     different from the algorithms described in
     [Tschantaridis/05]. Note that this argmax has to take into
     account the scoring function in sm, especially the weights sm.w,
     as well as the loss function, and whether linear or quadratic
     slacks are used. The weights in sm.w correspond to the features
     defined by psi() and range from index 1 to index
     sm->sizePsi. Most simple is the case of the zero/one loss
     function. For the zero/one loss, this function should return the
     highest scoring label ybar (which may be equal to the correct
     label y), or the second highest scoring label ybar, if
     Psi(x,ybar)>Psi(x,y)-1. If the function cannot find a label, it
     shall return an empty label as recognized by the function
     empty_label(y). */
  LABEL ybar;

  clock_t t = clock() - time_0;
  SSVM_PRINT("[svm_struct] find_most_violated_constraint_marginrescaling::START %ld %f\n", t, t/(float)CLOCKS_PER_SEC);

  int threadId = omp_get_thread_num();
  SSVM_PRINT("[svm_struct] Thread %d/%d\n", threadId, omp_get_num_threads());
  if(threadId >= maxBuffers) {
    printf("[svm_struct] ThreadId=%d > maxBuffers=%d\n", threadId, maxBuffers);
    exit(-1);
  }

#if VERBOSITY > 3
  // sm->w[0] is a dummy variable
  double* smw = sm->w + 1;

  // Compute energy of the ground truth for current w
  EnergyParam param;
  sparmToEnergyParam(*sparm, &param);
  GraphInference gi(x.slice,
                    &param,
                    smw,
                    x.feature,
                    x.nodeCoeffs, 0);
  double energyGT = gi.computeEnergy(y.nodeLabels);

  stringstream sEnergy;
  sEnergy << "energy_x" << x.id << ".txt";
  ofstream ofs(sEnergy.str().c_str(), ios::app);
  ofs << energyGT << endl;
  ofs.close();

#ifdef USE_OPENMP
  #pragma omp atomic
#endif
  totalWPsiGT -= energyGT;

#endif

#if VERBOSITY > 3
  // output psi to a file
  SWORD* words = computePsi(x,y,sm,sparm);
  stringstream sPsi;
  sPsi << "psiGT_" << x.id << ".txt";
  ofstream ofsPsi(sPsi.str().c_str(), ios::app);
  for(int i = 0; i<sm->sizePsi; i++) {
    ofsPsi << words[i].weight;
    ofsPsi << " "; 
  }
  ofsPsi << endl;;
  ofsPsi.close();
  delete[] words;

#endif

  // check if labels are stored in the cache
  int cacheId = x.id;
  bool labelFound = LabelCache::Instance()->getLabel(cacheId, ybar);
  if(!labelFound) {
    // allocate memory
    ybar.nNodes = y.nNodes;
    ybar.nodeLabels = new labelType[ybar.nNodes];
    for(int n = 0; n < y.nNodes; n++) {
      ybar.nodeLabels[n] = y.nodeLabels[n];
    }
    ybar.cachedNodeLabels = false;
    labelFound = true;
  }

  if(sparm->loss_function == 2) {
    double lossVOC  = 0;
    int nDiffVOC = 0;
    computeVOCLoss(y,ybar.nodeLabels,sparm->nClasses,lossVOC,nDiffVOC,x);
  }

  SSVM_PRINT("------------------------------------------------MOSTVIOLATED-%d-%d-START\n",
             x.id, sparm->iterationId);

  if(sparm->iterationId == 1 && sparm->includeLocalEdges && generateFirstConstraint) {
    printf("[SVM_struct] Generating noisy pattern for first violated constraint\n");

    for(int n = 0; n < ybar.nNodes; ++n) {
      if(ybar.nodeLabels[n] <= 1) {
        ybar.nodeLabels[n] = 1-y.nodeLabels[n];
      } else {
        ybar.nodeLabels[n] = rand() % sparm->nClasses;
      }
    }
  } else {
    runInference(x, y, sm, sparm, ybar, threadId, labelFound, cacheId);
  }

#if VERBOSITY > 2

  if((sparm->iterationId%sparm->stepForOutputFiles)==0) {
    stringstream MVCDir;
    MVCDir << mostViolatedConstraintDir;
    if(!useSlice3d) {
      //TODO : Remove !useVOC
      if(useVOC) {
        MVCDir << "x" << sparm->iterationId;
      }
      else {
        MVCDir << "x" << sparm->iterationId << "/";
      }
    }
    mkdir(MVCDir.str().c_str(), 0777);

    stringstream soutColoredImage;
    soutColoredImage << MVCDir.str();
    if(useSlice3d) {
      soutColoredImage << getNameFromPathWithoutExtension(x.slice->getName());
      soutColoredImage << "_";
      soutColoredImage << sparm->iterationId;
    } else {
      soutColoredImage << getNameFromPathWithoutExtension(x.slice->getName());
      soutColoredImage << ".png";
    }

    x.slice->exportSupernodeLabels(soutColoredImage.str().c_str(),
                                   sparm->nClasses,
                                   ybar.nodeLabels,
                                   ybar.nNodes,
                                   &(sparm->labelToClassIdx));

    if(useSlice3d) {
      zipAndDeleteCube(soutColoredImage.str().c_str());
    }
  }

#endif

  SSVM_PRINT("------------------------------------------------MOSTVIOLATED-%d-%d-END\n",
             x.id, sparm->iterationId);

  t = clock() - time_0;
  SSVM_PRINT("[svm_struct] find_most_violated_constraint_marginrescaling::END %ld %f\n", t, t/(float)CLOCKS_PER_SEC);

  return(ybar);
}

bool isSubmodular(const STRUCT_LEARN_PARM *sparm, double* smw)
{
  bool submodularEnergy = false;
  double* pw = smw + sparm->nUnaryWeights;
  double d = 0; // diagonal
  double a = 0; //anti-diagonal

  if(useGCForSubModularEnergy) {
    // For 2 classes, use graph-cuts if pairwise potential is attractive/submodular.
    // Let E be the energy and S the score (E=-S)
    // Recall that potential is submodular if :
    // E(0,0) + E(1,1) =< E(0,1) + E(1,0)
    // S(0,0) + S(1,1) >= S(0,1) + S(1,0) or d >= a
    if(sparm->includeLocalEdges && sparm->nClasses == 2) {
      if(sparm->nGradientLevels == 0) {
        submodularEnergy = true;
      } else {
        submodularEnergy = true;
#if USE_LONG_RANGE_EDGES
        for(int p = 0; submodularEnergy && (p < sparm->nDistances); ++p)
          {
            a = 0;
            d = 0;
#else
          {
#endif
            for(int g = 0; g < sparm->nGradientLevels; ++g) {
              d += pw[0] + pw[3];
              a += pw[1] + pw[2];
              if(a > (d - 1e-8)) { // add a small margin because of numerical errors
              //if(a > d) {
                submodularEnergy = false;
                break;
              }
              pw += sparm->nClasses*sparm->nClasses;
            }
          }
      }
    }
  }
  return submodularEnergy;
}

void runInference(SPATTERN x, LABEL y, 
                  const STRUCTMODEL *sm, 
                  const STRUCT_LEARN_PARM *sparm,
                  LABEL& ybar, const int threadId, bool labelFound, int cacheId)
{
  GraphInference* gi_MVC;
  bool computeEnergyAtEachIteration = true;
  // sm->w[0] is a dummy variable
  double* smw = sm->w + 1;

  EnergyParam param;
  sparmToEnergyParam(*sparm, &param);

  switch(sparm->giType) {
  case T_GI_LIBDAI:
    {
      bool useGC = isSubmodular(sparm, smw);
      if(useGC) {
        double* pw = smw + sparm->nUnaryWeights;
        SSVM_PRINT("[SVM_struct] Use graph-cuts for attractive potential (%g<%g) at iteration %d\n",
                   pw[1]+pw[2], pw[0]+pw[3], sparm->iterationId);

#if USE_MAXFLOW
        gi_MVC = new GI_maxflow(x.slice,
                                &param,
                                smw,
                                y.nodeLabels, // groundtruth labels used to compute loss  
                                sparm->lossPerLabel,
                                x.feature,
                                x.nodeCoeffs,
                                x.edgeCoeffs
                                );

        double energy = gi_MVC->run(ybar.nodeLabels, // inferred labels
                                    x.id,
                                    MVC_MAX_ITER,
                                    y.nodeLabels, // ground truth
                                    computeEnergyAtEachIteration);

        SSVM_PRINT("[MostViolatedConstraint] maxflow energy=%g (This should be equal to -score)\n", energy);
#else
        printf("[svm_struct] USE_MAXFLOW is set to 0\n");
        exit(-1);
#endif
      } else {

#if USE_LIBDAI
        gi_MVC = new GI_libDAI(x.slice,
                               &param,
                               smw,
                               y.nodeLabels, // groundtruth labels used to compute loss  
                               sparm->lossPerLabel,
                               x.feature,
                               x.nodeCoeffs,
                               x.edgeCoeffs
                               );

        double energy = gi_MVC->run(ybar.nodeLabels, // inferred labels
                                    x.id,
                                    MVC_MAX_ITER,
                                    y.nodeLabels, // ground truth
                                    computeEnergyAtEachIteration);

        SSVM_PRINT("[MostViolatedConstraint] libDAI energy=%g (This should be equal to -score)\n", energy);

#else
        printf("[svm_struct] USE_LIBDAI is set to 0\n");
        exit(-1);
#endif
      }
    }
    break;

  case T_GI_LIBDAI_ICM:
    {
      bool useGC = isSubmodular(sparm, smw);
      double energy = 0;
      if(useGC) {
        double* pw = smw + sparm->nUnaryWeights;
        SSVM_PRINT("[SVM_struct] Use graph-cuts for attractive potential (%g<%g) at iteration %d\n",
                   pw[1]+pw[2], pw[0]+pw[3], sparm->iterationId);

#if USE_MAXFLOW
        gi_MVC = new GI_maxflow(x.slice,
                                &param,
                                smw,
                                y.nodeLabels, // groundtruth labels used to compute loss  
                                sparm->lossPerLabel,
                                x.feature,
                                x.nodeCoeffs,
                                x.edgeCoeffs
                                );

        energy = gi_MVC->run(ybar.nodeLabels, // inferred labels
                             x.id,
                             MVC_MAX_ITER,
                             y.nodeLabels, // ground truth
                             computeEnergyAtEachIteration);
#else
        printf("[svm_struct] USE_MAXFLOW is set to 0\n");
        exit(-1);
#endif

        SSVM_PRINT("[MostViolatedConstraint] maxflow energy=%g (This should be equal to -score)\n", energy);
      } else {

#if USE_LIBDAI
        gi_MVC = new GI_libDAI(x.slice,
                               &param,
                               smw,
                               y.nodeLabels, // groundtruth labels used to compute loss  
                               sparm->lossPerLabel,
                               x.feature,
                               x.nodeCoeffs,
                               x.edgeCoeffs
                               );

        energy = gi_MVC->run(ybar.nodeLabels, // inferred labels
                             x.id,
                             MVC_MAX_ITER,
                             y.nodeLabels, // ground truth
                             computeEnergyAtEachIteration);
#else
        printf("[svm_struct] USE_LIBDAI is set to 0\n");
        exit(-1);
#endif
        SSVM_PRINT("[MostViolatedConstraint] libDAI energy=%g (This should be equal to -score)\n", energy);
      }

      delete gi_MVC;

#if VERBOSITY > 3

      if((sparm->iterationId%sparm->stepForOutputFiles)==0) {
        stringstream MVCDir;
        MVCDir << mostViolatedConstraintDir;
        if(!useSlice3d) {
          //TODO : Remove !useVOC
          if(useVOC) {
            MVCDir << "xx" << sparm->iterationId;
          }
          else {
            MVCDir << "xx" << sparm->iterationId << "/";
          }
        }
        mkdir(MVCDir.str().c_str(), 0777);

        stringstream soutColoredImage;
        soutColoredImage << MVCDir.str();
        if(useSlice3d) {
          soutColoredImage << getNameFromPathWithoutExtension(x.slice->getName());
          soutColoredImage << "_";
          soutColoredImage << sparm->iterationId;
        } else {
          soutColoredImage << x.slice->getName();
        }

        // SSVM_PRINT("[SVM_struct] Saving labels for most violated constraint to %s\n",
        //           soutColoredImage.str().c_str());
        x.slice->exportSupernodeLabels(soutColoredImage.str().c_str(),
                                       sparm->nClasses,
                                       ybar.nodeLabels,
                                       ybar.nNodes,
                                       &(sparm->labelToClassIdx));

        if(useSlice3d) {
          zipAndDeleteCube(soutColoredImage.str().c_str());
        }
      }

#endif

      gi_MVC = new GI_ICM(x.slice,
                          &param,
                          smw,
                          y.nodeLabels, // groundtruth labels used to compute loss  
                          sparm->lossPerLabel,
                          x.feature,
                          x.nodeCoeffs
                          );

      double energy2 = gi_MVC->run(ybar.nodeLabels, // inferred labels
                                   x.id,
                                   MVC_MAX_ITER,
                                   y.nodeLabels, // ground truth
                                   computeEnergyAtEachIteration);

      SSVM_PRINT("[MostViolatedConstraint] libDAI+ICM energy=%g (This should be equal to -score)\n", energy2);
      //assert(energy2 < (energy+1e-3));
      if(energy2 > (energy+1e-3)) {
        printf("[svm_struct] Error Energy(ICM) = %g > Energy(BP) = %g\n", energy2, energy);
        exit(-1);
      }
    }
    break;

#if USE_MRF && USE_LIBDAI
  case T_GI_LIBDAI_ICM_QPBO:
    {
      // run QPBO first
      GI_MRF* _gi_mrf = new GI_MRF(x.slice,
                                   &param,
                                   smw,
                                   y.nodeLabels, // groundtruth labels used to compute loss  
                                   sparm->lossPerLabel,
                                   x.feature,
                                   x.nodeCoeffs
                                   );
      _gi_mrf->setUseQPBO(true);

      double energy_QPBO = _gi_mrf->run(//ybar.nodeLabels, // inferred labels
                                        tempNodeLabels[0],
                                        x.id,
                                        MVC_MAX_ITER,
                                        y.nodeLabels, // ground truth
                                        computeEnergyAtEachIteration);
      delete _gi_mrf;
      SSVM_PRINT("[MostViolatedConstraint] QPBO energy=%g (This should be equal to -score)\n", energy_QPBO);

      gi_MVC = new GI_libDAI(x.slice,
                             &param,
                             smw,
                             y.nodeLabels, // groundtruth labels used to compute loss  
                             sparm->lossPerLabel,
                             x.feature,
                             x.nodeCoeffs,
                             x.edgeCoeffs
                             );

      double energy = gi_MVC->run(ybar.nodeLabels, // inferred labels
                                  x.id,
                                  MVC_MAX_ITER,
                                  y.nodeLabels, // ground truth
                                  computeEnergyAtEachIteration);

      SSVM_PRINT("[MostViolatedConstraint] libDAI energy=%g (This should be equal to -score)\n", energy);

      delete gi_MVC;

      if(energy_QPBO < energy) {
        for(int n = 0; n < ybar.nNodes; ++n) {
          ybar.nodeLabels[n] = tempNodeLabels[0][n];
        }
      }

      // try to improve with ICM
      gi_MVC = new GI_ICM(x.slice,
                          &param,
                          smw,
                          y.nodeLabels, // groundtruth labels used to compute loss  
                          sparm->lossPerLabel,
                          x.feature,
                          x.nodeCoeffs
                          );

      double energy2 = gi_MVC->run(ybar.nodeLabels, // inferred labels
                                   x.id,
                                   MVC_MAX_ITER,
                                   y.nodeLabels, // ground truth
                                   computeEnergyAtEachIteration);

      SSVM_PRINT("[MostViolatedConstraint] libDAI+ICM energy=%g (This should be equal to -score)\n", energy2);
      //assert(energy2 < (energy+1e-3));
      if(energy2 > (energy+1e-3)) {
        printf("[svm_struct] Error Energy(ICM) = %g > Energy(BP) = %g\n", energy2, energy);
      }
    }
    break;
#else
    printf("[svm_struct] USE_LIBDAI and/or USE_MRF are set to 0\n");
    exit(-1);
#endif


#if USE_OPENGM
  case T_GI_OPENGM:
    {
      bool useGC = isSubmodular(sparm, smw);

      if(useGC) {
        double* pw = smw + sparm->nUnaryWeights;
        SSVM_PRINT("[SVM_struct] Use graph-cuts for attractive potential (%g<%g) at iteration %d\n",
                   pw[1]+pw[2], pw[0]+pw[3], sparm->iterationId);

#if USE_MAXFLOW
        gi_MVC = new GI_maxflow(x.slice,
                                &param,
                                smw,
                                y.nodeLabels, // groundtruth labels used to compute loss  
                                sparm->lossPerLabel,
                                x.feature,
                                x.nodeCoeffs,
                                x.edgeCoeffs
                                );

        double energy = gi_MVC->run(ybar.nodeLabels, // inferred labels
                                    x.id,
                                    MVC_MAX_ITER,
                                    y.nodeLabels, // ground truth
                                    computeEnergyAtEachIteration);

        SSVM_PRINT("[MostViolatedConstraint] maxflow energy=%g (This should be equal to -score)\n", energy);
#else
        printf("[svm_struct] USE_MAXFLOW is set to 0\n");
        exit(-1);
#endif
      } else {

        gi_MVC = new GI_opengm(x.slice,
                               &param,
                               smw,
                               y.nodeLabels, // groundtruth labels used to compute loss  
                               sparm->lossPerLabel,
                               x.feature
                               );
        double energy = gi_MVC->run(ybar.nodeLabels, // inferred labels
                                    //tempNodeLabels[0],
                                    x.id,
                                    MVC_MAX_ITER,
                                    y.nodeLabels, // ground truth
                                    computeEnergyAtEachIteration);

        SSVM_PRINT("[MostViolatedConstraint] opengm energy=%g (This should be equal to -score)\n", energy);
      }

      break;
    }
#endif

  case T_GI_MAX:
    {
      gi_MVC = new GI_max(x.slice,
                          &param,
                          smw,
                          y.nodeLabels, // groundtruth labels used to compute loss  
                          sparm->lossPerLabel,
                          x.feature,
                          x.nodeCoeffs);
      double energy = gi_MVC->run(ybar.nodeLabels, // inferred labels
                                  x.id,
                                  MVC_MAX_ITER,
                                  y.nodeLabels, // ground truth
                                  computeEnergyAtEachIteration);
      SSVM_PRINT("[MostViolatedConstraint] MAX energy=%g (This should equal to -score)\n", energy);
    }
    break;

  case T_GI_MF:
    {
      GI_MF* gi_MF = new GI_MF(x.slice,
                               &param,
                               smw,
                               y.nodeLabels, // groundtruth labels used to compute loss  
                               sparm->lossPerLabel,
                               x.feature,
                               x.nodeCoeffs);
      gi_MVC = gi_MF;
      gi_MF->setBelieves(tempPotentials[threadId]);
      double energy = gi_MVC->run(ybar.nodeLabels, // inferred labels
                                  x.id,
                                  MVC_MAX_ITER,
                                  y.nodeLabels, // ground truth
                                  computeEnergyAtEachIteration);
      SSVM_PRINT("[MostViolatedConstraint] MF energy=%g (This should equal to -score)\n", energy);
    }
    break;

  case T_GI_SAMPLING:
    {
      //if(nParallelChains == 1) {
      if(nParallelChains == 0) {
        // temperature is picked automatically

        GI_sampling* gi_sampling
          = new GI_sampling(x.slice,
                            &param,
                            smw,
                            y.nodeLabels, // groundtruth labels used to compute loss  
                            sparm->lossPerLabel,
                            x.feature,
                            x.nodeCoeffs,
                            sparm->sampling_rate);
        gi_MVC = gi_sampling;

        gi_sampling->setInitializedLabels(labelFound);

        double energy = 0;
        if(sparm->loss_function == 2) {
          energy = gi_sampling->run_VOC(ybar.nodeLabels, // inferred labels
                                        x.id,
                                        MVC_MAX_ITER,
                                        y.nodeLabels, // ground truth
                                        x.TPs, x.FPs, x.FNs);
        } else {
          energy = gi_MVC->run(ybar.nodeLabels, // inferred labels
                               x.id,
                               MVC_MAX_ITER,
                               y.nodeLabels, // ground truth
                               computeEnergyAtEachIteration);
        }

        SSVM_PRINT("[MostViolatedConstraint] SAMPLING energy=%g (This should equal to -score)\n", energy);

        if(cacheId != -1) {
          LabelCache::Instance()->setLabel(cacheId, ybar);
          ybar.cachedNodeLabels = true;
        }

      } else {
        // run several chains in parallel
        double *energies = (double *)my_malloc(sizeof(double)*nParallelChains);

        // temperature of the first chain
        const double temperature_chain0 = sparm->sampling_temperature_0*(pow(10.0,(nParallelChains/2)));
        SSVM_PRINT("[MostViolatedConstraint] temperature=%g temperature_chain0=%g\n",
                   sparm->sampling_temperature_0, temperature_chain0);

        // TODO : find a better way...
        if(sparm->nClasses > 3) {
          // no parallel for
          for(int iChain = 0; iChain < nParallelChains; ++iChain) {
            int bufferId = threadId*nParallelChains + iChain;
            assert(bufferId < maxBuffers);
            double temperature = temperature_chain0/(pow(10.0,iChain));
            GI_sampling* gi_sampling
              = new GI_sampling(x.slice, &param, smw, y.nodeLabels,
                                sparm->lossPerLabel, x.feature, x.nodeCoeffs,
                                sparm->sampling_rate);
            gi_sampling->setInitializedLabels(labelFound);
            if(labelFound) {
              for(int n = 0; n < ybar.nNodes; ++n) {
                tempNodeLabels[bufferId][n] = ybar.nodeLabels[n];
              }
            }

            if(sparm->loss_function == LOSS_VOC) {
              // copy TPS, FPs and FNs
              for(int c = 0; c < sparm->nClasses; ++c) {
                TPs[bufferId][c] = x.TPs[c];
                FPs[bufferId][c] = x.FPs[c];
                FNs[bufferId][c] = x.FNs[c];
              }

              energies[iChain] = gi_sampling->runOnce_VOC(tempNodeLabels[bufferId],
                                                          MVC_MAX_ITER, y.nodeLabels,
                                                          temperature,
                                                          TPs[bufferId], FPs[bufferId], FNs[bufferId]);
            } else {
              energies[iChain] = gi_sampling->runOnce(tempNodeLabels[bufferId],
                                                      MVC_MAX_ITER, y.nodeLabels,
                                                      0, temperature);
            }
            delete gi_sampling;
          }
        } else {
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
          for(int iChain = 0; iChain < nParallelChains; ++iChain) {
            int bufferId = threadId*nParallelChains + iChain;

            int _threadId = omp_get_thread_num();
            SSVM_PRINT("[svm_struct] Thread sampling %d/%d\n", _threadId, omp_get_num_threads());

            assert(bufferId < maxBuffers);
            double temperature = temperature_chain0/(pow(10.0,iChain));
            GI_sampling* gi_sampling
              = new GI_sampling(x.slice, &param, smw, y.nodeLabels,
                                sparm->lossPerLabel, x.feature, x.nodeCoeffs,
                                sparm->sampling_rate);
            
            /*
            // re-use previous sampling output
            gi_sampling->setInitializedLabels(labelFound);
            if(labelFound) {
              for(int n = 0; n < ybar.nNodes; ++n) {
                tempNodeLabels[bufferId][n] = ybar.nodeLabels[n];
              }
            }
            */

            if(sparm->loss_function == LOSS_VOC) {
              // copy TPS, FPs and FNs
              for(int c = 0; c < sparm->nClasses; ++c) {
                TPs[bufferId][c] = x.TPs[c];
                FPs[bufferId][c] = x.FPs[c];
                FNs[bufferId][c] = x.FNs[c];
              }

              energies[iChain] = gi_sampling->runOnce_VOC(tempNodeLabels[bufferId],
                                                          MVC_MAX_ITER, y.nodeLabels,
                                                          temperature,
                                                          TPs[bufferId], FPs[bufferId], FNs[bufferId]);
            } else {
              energies[iChain] = gi_sampling->runOnce(tempNodeLabels[bufferId],
                                                      MVC_MAX_ITER, y.nodeLabels,
                                                      0, temperature);
            }
            delete gi_sampling;
          }
        }

        gi_MVC = 0;

        // copy the best chain
        double minEnergy = energies[0];
        int bestChain = 0;
        for(int iChain = 0; iChain < nParallelChains; ++iChain) {
          double temperature = temperature_chain0/(pow(10.0,iChain));
          SSVM_PRINT("[MostViolatedConstraint] sampling (%d, %g) energy = %g\n",
                     iChain, temperature, energies[iChain]);
          if(energies[iChain] < minEnergy) {
            minEnergy = energies[iChain];
            bestChain = iChain;
          }
        }
        // check me
        int bufferId = threadId*nParallelChains + bestChain;
        assert(bufferId < maxBuffers);
        ybar.nodeLabels = tempNodeLabels[bufferId];

        //------------------------------------------------------------------------------
        // TODO: debugging
        // TODO: debugging
        int _sizePsi = sm->sizePsi + 1;
        SWORD* fy_to = new SWORD[_sizePsi];
        computePsi(fy_to, x, y, sm, sparm);
        double score_gt = computeScore(sm, fy_to);
        printf("[Sampling] temperature=%g scoreGT=%g\n",
                   sparm->sampling_temperature_0, score_gt);

        SWORD* fy_away = new SWORD[_sizePsi];
        computePsi(fy_away, x, ybar, sm, sparm);
        double score_c = computeScore(sm, fy_away);
        printf("[Sampling] temperature=%g scoreC=%g\n",
                   sparm->sampling_temperature_0, score_c);

        /*
        stringstream sPsi;
        sPsi << "psi_" << x.id << ".txt";
        ofstream ofsPsi(sPsi.str().c_str(), ios::app);
        for(int i = 0; i<sm->sizePsi; i++) {
          ofsPsi << fy_away[i].weight;
          ofsPsi << " "; 
        }
        ofsPsi << endl;;
        ofsPsi.close();
        */
        delete[] fy_away;
        //------------------------------------------------------------------------------

#if VERBOSITY > 2
        ofstream ofs_temp("best_temperature.txt", ios::app);
        ofs_temp << bestChain << endl;
        ofs_temp.close();
#endif

        double temperature = temperature_chain0/(pow(10.0,bestChain));
        SSVM_PRINT("[MostViolatedConstraint] best sampling (%d, %d, %g) energy = %g (This should equal to -score)\n",
                   bestChain, bufferId, temperature, minEnergy);

        if(cacheId != -1) {
          LabelCache::Instance()->setLabel(cacheId, ybar);
          ybar.cachedNodeLabels = true;
        }

        free(energies);
      }
    }
    break;

#if USE_MULTIOBJ
  case T_GI_MULTIOBJ:
    {

      gi_MVC = new GI_multiobject(x.slice,
                                  &param,
                                  smw,
                                  y.nodeLabels, // groundtruth labels used to compute loss  
                                  sparm->lossPerLabel,
                                  x.feature,
                                  x.nodeCoeffs, 0
                                  );

      double energy = gi_MVC->run(ybar.nodeLabels, // inferred labels
                                  x.id,
                                  MVC_MAX_ITER,
                                  y.nodeLabels, // ground truth
                                  computeEnergyAtEachIteration);

      SSVM_PRINT("[MostViolatedConstraint] multiobj energy=%g (This should be equal to -score)\n", energy);
      break;
    }
#endif

  default:
    printf("[svm_struct] Unknown type %d\n", sparm->giType);
    exit(-1);
    break;
  }

  if(gi_MVC) {
    delete gi_MVC;
  }
}

int         empty_label(LABEL y)
{
  /* Returns true, if y is an empty label. An empty label might be
     returned by find_most_violated_constraint_???(x, y, sm) if there
     is no incorrect label that can be found for x, or if it is unable
     to label x at all */
  return(0);
}

SVECTOR     *psi(SPATTERN x, LABEL y, const STRUCTMODEL *sm,
		 const STRUCT_LEARN_PARM *sparm)
{
  /* Returns a feature vector describing the match between pattern x
     and label y. The feature vector is returned as a list of
     SVECTOR's. Each SVECTOR is in a sparse representation of pairs
     <featurenumber:featurevalue>, where the last pair has
     featurenumber 0 as a terminator. Featurenumbers start with 1 and
     end with sizePsi. Featuresnumbers that are not specified default
     to value 0. As mentioned before, psi() actually returns a list of
     SVECTOR's. Each SVECTOR has a field 'factor' and 'next'. 'next'
     specifies the next element in the list, terminated by a NULL
     pointer. The list can be though of as a linear combination of
     vectors, where each vector is weighted by its 'factor'. This
     linear combination of feature vectors is multiplied with the
     learned (kernelized) weight vector to score label y for pattern
     x. Without kernels, there will be one weight in sm.w for each
     feature. Note that psi has to match
     find_most_violated_constraint_???(x, y, sm) and vice versa. In
     particular, find_most_violated_constraint_???(x, y, sm) finds
     that ybar!=y that maximizes psi(x,ybar,sm)*sm.w (where * is the
     inner vector product) and the appropriate function of the
     loss + margin/slack rescaling method. See that paper for details. */

  SVECTOR *fvec=NULL;

  clock_t t = clock() - time_0;
  SSVM_PRINT("[svm_struct] psi::START %ld %f\n", t, t/(float)CLOCKS_PER_SEC);

  /* insert code for computing the feature vector for x and y here */

  SSVM_PRINT("[MostViolatedConstraint] computePsi: ");
  double lhsXw = 0;
  SWORD* words = computePsi(x,y,sm,sparm,&lhsXw);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  totalWPsi += lhsXw;

#if VERBOSITY > 3
  stringstream sPsi;
  sPsi << "psiMVC_" << x.id << ".txt";
  ofstream ofsPsi(sPsi.str().c_str(), ios::app);
  ofsPsi << "Iteration " << sparm->iterationId << "\n";

#if VERBOSITY > 3
  SSVM_PRINT("[MostViolatedConstraint] psi:");
#endif

  for(int i = 0; i < sm->sizePsi; i++) {
#if VERBOSITY > 3
    SSVM_PRINT(" %d:%g", words[i].wnum, words[i].weight);
#endif  
    //ofsPsi << i << ":" << words[i].weight << " ";
    ofsPsi << words[i].weight << " ";
  }

#if VERBOSITY > 3
  SSVM_PRINT("\n");
#endif

  ofsPsi << endl;;
  ofsPsi.close();
#endif

  fvec = create_svector(words, (char*)"", 1.0);
  delete[] words;

  t = clock() - time_0;
  SSVM_PRINT("[svm_struct] psi::END %ld %f\n", t, t/(float)CLOCKS_PER_SEC);

  return fvec;
}

double      loss(LABEL y, LABEL ybar, const STRUCT_LEARN_PARM *sparm)
{
  /* loss for correct label y and predicted label ybar. The loss for
     y==ybar has to be zero. sparm->loss_function is set with the -l option. */
  double loss  = 0;
  int nDiff = 0;
  computeLoss(y,ybar.nodeLabels,sparm,loss,nDiff);

#if VERBOSITY > 3
  ofstream ofs("lossMVC.txt", ios::app);
  ofs << loss << endl;
  ofs.close();
#endif

  SSVM_PRINT("[MostViolatedConstraint] Loss = %g (%d/%d nodes have different labels)\n", loss, nDiff, y.nNodes);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  totalLossMVC += loss;

  return loss;
}

void save_parameters(const char* output_filename, STRUCT_LEARN_PARM *sparm, STRUCTMODEL *sm)
{
  ofstream ofs(output_filename);
  ofs << sm->sizePsi << endl;
  ofs << sparm->nClasses << endl;
  ofs << sparm->nGradientLevels << endl;
  ofs << sparm->nOrientations << endl;
  ofs << sparm->nScales << endl;
  ofs << sparm->nLocalScales << endl;
  ofs << sparm->nScalingCoefficients << endl;

  ofs.precision(9);

  double* smw = sm->w + 1;
  for(int i = 0; i < sm->sizePsi; i++) {
    ofs << smw[i] << endl;
  }

  ofs << endl;
  ofs.close();
}

int         finalize_iteration(double ceps, int cached_constraint,
			       SAMPLE sample, STRUCTMODEL *sm,
			       CONSTSET cset, double *alpha, 
			       STRUCT_LEARN_PARM *sparm)
{
  /* This function is called just before the end of each cutting plane iteration.
     ceps is the amount by which the most violated constraint found in the current iteration was violated.
     cached_constraint is true if the added constraint was constructed from the cache. If the return value is FALSE,
     then the algorithm is allowed to terminate. If it is TRUE, the algorithm will keep iterating even if the desired
     precision sparm->epsilon is already reached. */
  totalLhsXw = totalWPsiGT-totalWPsi;
  double rhs = totalLossMVC/sample.n;
  double lhsXw = totalLhsXw/sample.n;
  double slack = rhs - lhsXw;
  SSVM_PRINT("[finalize_iteration] iteration=%d,ceps=%g,totalLossMVC=%g,n=%d,rhs=totalLossMVC/n=%g,totalLhsXw=%g,totalLhsXw/n=%g,slack=%g,totalWPsi=%g,totalWPsiGT=%g,C=%g\n",
             sparm->iterationId, ceps, totalLossMVC,sample.n,rhs,totalLhsXw,lhsXw,slack,totalWPsi,totalWPsiGT,sparm->C);

  if( (sparm->iterationId>1) &&
      (oldC != sparm->C || (sparm->iterationId%sparm->stepForParameterFiles)==0) ) {
    stringstream sout;
    sout << parameterDir;
    sout << "/iteration_" << sparm->iterationId << ".txt";
    SSVM_PRINT("[finalize_iteration] Writing parameter file %s\n", sout.str().c_str());
    save_parameters(sout.str().c_str(), sparm, sm);
  }

  // reset for next iteration
  totalLoss = 0;
  totalLossMVC = 0;
  totalWPsi = 0;
  totalWPsiGT = 0;

#if VERBOSITY > 3
  // output average and VOC scores
  stringstream soutScores;
  soutScores << "scores/scores_" << sparm->iterationId << ".txt";
  ofstream ofsScores(soutScores.str().c_str());

  ulong* _TPs = new ulong[sparm->nClasses];
  ulong* _FPs = new ulong[sparm->nClasses];
  ulong* _FNs = new ulong[sparm->nClasses];

  for(int c=0; c < sparm->nClasses; c++) {
    _TPs[c] = 0;
    _FPs[c] = 0;
    _FNs[c] = 0;
  }

  for(int sid = 0; sid < sample.n; sid++) {
    for(int c=0; c < sparm->nClasses; c++) {
      _TPs[c] += sample.examples[sid].x.TPs[c];
      _FPs[c] += sample.examples[sid].x.FPs[c];
      _FNs[c] += sample.examples[sid].x.FNs[c];
    }
  }

  for(int c=0; c < sparm->nClasses; c++) {
    ofsScores << _TPs[c] << " " << _FPs[c] << " " << _FNs[c] << endl;
  }

  ofsScores.close();
  delete[] _TPs;
  delete[] _FPs;
  delete[] _FNs;
#endif

  if( (sparm->iterationId>1) &&
      (oldC != sparm->C || (sparm->iterationId%sparm->stepForOutputFiles)==0) ) {

    const bool compress_image = true;

    // switch to BP (or graph-cuts) if sampling is specified
    bool use_sampling = false;
    if(sparm->giType == T_GI_SAMPLING) {
      use_sampling = true;
      sparm->giType = T_GI_LIBDAI;
    }

    //int iterationId = sparm->iterationId-1;
    int iterationId = sparm->iterationId;

    string config_tmp;
    bool exportOverlay = false;
    if(Config::Instance()->getParameter("exportOverlay", config_tmp)) {
      exportOverlay = (atoi(config_tmp.c_str()) == 1);
      SSVM_PRINT("[svm_struct] exportOverlay = %d\n", (int)exportOverlay);
    }

    // output images and compute score for the previous iteration (before updating C)
    stringstream soutWeightFile;
    soutWeightFile << parameterDir;
    soutWeightFile << "/iteration_" << iterationId << ".txt";

    // output directory for images
    stringstream soutOutput;
    if(!useSlice3d) {
      soutOutput.setf(ios::scientific, ios::floatfield);
      soutOutput.precision(0);
      soutOutput << "iteration_" << iterationId << "_";
      soutOutput << "_" << sparm->nGradientLevels << "_" << sparm->nOrientations << "_" << sparm->nScales << "_";
      soutOutput << sparm->startC << "_" << sparm->C;
    }

    EnergyParam param(soutWeightFile.str().c_str());

    //TODO : temporary hack
#if USE_LONG_RANGE_EDGES
    param.nDistances = sparm->nDistances;
#endif

#if VERBOSITY > 4
    if(oldC != sparm->C)
#endif
      {
        // Might not want to do the prediction on training images to save time
        if(predictTrainingImages) {
          string soutOutputDir = inferenceDir_training;
          if(!useSlice3d) {
            soutOutputDir += soutOutput.str() + "/";
            mkdir(soutOutputDir.c_str(), 0777);
          }

          SSVM_PRINT("[SVM_struct] Running inference for training dataset. Output results will be stored in %s\n",
                     soutOutputDir.c_str());

          SSVM_PRINT("[svm_struct] training segmentImage::START %ld\n", clock() - time_0);

          string training_score_filename = scoreDir + "training_score.txt";

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
          for(int sid = 0; sid < sample.n; sid++) {

            int threadId = omp_get_thread_num();
            SSVM_PRINT("[svm_struct] segmentImage Thread %d/%d\n", threadId, omp_get_num_threads());

            stringstream sout;
            sout << soutOutputDir.c_str();
            if(useSlice3d) {
              // add name of the cube
              sout << getNameFromPathWithoutExtension(sample.examples[sid].x.slice->getName());
              sout << "_";
              sout << iterationId;
            }

            if(useSlice3d) {
              segmentImage(sample.examples[sid].x,
                           sout.str().c_str(),
                           sparm->giType,
                           soutWeightFile.str().c_str(),
                           &sparm->labelToClassIdx,
                           training_score_filename,
                           0,0,
                           compress_image,
                           (exportOverlay)?"overlay_training/":0,
                           sparm->metric_type);
            } else {
              for(int c = 0; c < sparm->nClasses; c++) {
                sample.examples[sid].x.TPs[c] = 0;
                sample.examples[sid].x.FPs[c] = 0;
                sample.examples[sid].x.FNs[c] = 0;
                sample.examples[sid].x.count[c] = 0;
              }

              Slice* slice = static_cast<Slice*>(sample.examples[sid].x.slice);
              computeScore(slice,
                           sample.examples[sid].x.feature,
                           sample.examples[sid].x.imgAnnotation,
                           sparm->giType,
                           param,
                           sparm->labelToClassIdx,
                           sparm->classIdxToLabel,
                           true,
                           param.nClasses,
                           sample.examples[sid].x.TPs, sample.examples[sid].x.FPs,
                           sample.examples[sid].x.FNs, sample.examples[sid].x.count,
                           soutOutputDir.c_str(),
                           sample.examples[sid].x.nodeCoeffs,
                           sample.examples[sid].x.edgeCoeffs);
            }
          }

          SSVM_PRINT("[svm_struct] training segmentImage::END %ld\n", clock() - time_0);

          if(!useSlice3d) {
            ulong *TPs = (ulong *)my_malloc(sizeof(ulong)*sparm->nClasses);
            ulong *FPs = (ulong *)my_malloc(sizeof(ulong)*sparm->nClasses);
            ulong *FNs = (ulong *)my_malloc(sizeof(ulong)*sparm->nClasses);
            ulong *count = (ulong *)my_malloc(sizeof(ulong)*sparm->nClasses);
            for(int c = 0; c < sparm->nClasses; c++) {
              TPs[c] = 0;
              FPs[c] = 0;
              FNs[c] = 0;
              count[c] = 0;
            }
            for(int sid = 0; sid < sample.n; sid++) {
              for(int c = 0; c < sparm->nClasses; c++) {
                TPs[c] += sample.examples[sid].x.TPs[c];
                FPs[c] += sample.examples[sid].x.FPs[c];
                FNs[c] += sample.examples[sid].x.FNs[c];
                count[c] += sample.examples[sid].x.count[c];
              }
            }

            outputScore(training_score_filename.c_str(), sparm->classIdxToLabel, TPs, FPs, FNs,
                        count, soutWeightFile.str().c_str());
            free(TPs);
            free(FPs);
            free(FNs);
            free(count);
          }
        }
      }

    if(testDir != "") {
        string soutTestOutputDir = inferenceDir_test;
        if(!useSlice3d) {
          soutTestOutputDir += soutOutput.str() + "/";
          mkdir(soutTestOutputDir.c_str(), 0777);
        }

        SSVM_PRINT("[SVM_struct] Running inference for test dataset containing %ld images\n",
                      nTestExamples); 
        SSVM_PRINT("[SVM_struct] Output results will be stored in %s\n",
                      soutTestOutputDir.c_str());

        clock_t t = clock() - time_0;
        SSVM_PRINT("[svm_struct] test segmentImage::START %ld %f\n", t,
                   t/(float)CLOCKS_PER_SEC);

        string testing_score_filename = scoreDir + "test_score.txt";

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
        for(int sid = 0; sid < nTestExamples; sid++)
          {
            stringstream sout;
            sout << soutTestOutputDir.c_str();
            if(useSlice3d) {
              // add name of the cubeness
              sout << getNameFromPathWithoutExtension(test_examples[sid].x.slice->getName());
              sout << "_";
              sout << iterationId;
            }

            if(useSlice3d) {
              segmentImage(test_examples[sid].x,
                           sout.str().c_str(),
                           sparm->giType,
                           soutWeightFile.str().c_str(),
                           &sparm->labelToClassIdx,
                           testing_score_filename,
                           0,0,
                           compress_image,
                           (exportOverlay)?"overlay_test/":0,
                           sparm->metric_type);
            } else {
              for(int c = 0; c < sparm->nClasses; c++) {
                sample.test_examples[sid].x.TPs[c] = 0;
                sample.test_examples[sid].x.FPs[c] = 0;
                sample.test_examples[sid].x.FNs[c] = 0;
                sample.test_examples[sid].x.count[c] = 0;
              }

              Slice* slice = static_cast<Slice*>(sample.test_examples[sid].x.slice);
              computeScore(slice,
                           sample.test_examples[sid].x.feature,
                           sample.test_examples[sid].x.imgAnnotation,
                           sparm->giType,
                           param,
                           sparm->labelToClassIdx,
                           sparm->classIdxToLabel,
                           true,
                           param.nClasses,
                           sample.test_examples[sid].x.TPs, sample.test_examples[sid].x.FPs,
                           sample.test_examples[sid].x.FNs, sample.test_examples[sid].x.count,
                           soutTestOutputDir.c_str(),
                           test_examples[sid].x.nodeCoeffs,
                           test_examples[sid].x.edgeCoeffs);
            }
          }
        t = clock() - time_0;
        SSVM_PRINT("[svm_struct] test segmentImage::END %ld %f\n", t, t/(float)CLOCKS_PER_SEC);

        if(!useSlice3d) {
            ulong *TPs = (ulong *)my_malloc(sizeof(ulong)*sparm->nClasses);
            ulong *FPs = (ulong *)my_malloc(sizeof(ulong)*sparm->nClasses);
            ulong *FNs = (ulong *)my_malloc(sizeof(ulong)*sparm->nClasses);
            ulong *count = (ulong *)my_malloc(sizeof(ulong)*sparm->nClasses);
          for(int c = 0; c < sparm->nClasses; c++) {
            TPs[c] = 0;
            FPs[c] = 0;
            FNs[c] = 0;
            count[c] = 0;
          }
          for(int sid = 0; sid < nTestExamples; sid++) {
            for(int c = 0; c < sparm->nClasses; c++) {
              TPs[c] += sample.test_examples[sid].x.TPs[c];
              FPs[c] += sample.test_examples[sid].x.FPs[c];
              FNs[c] += sample.test_examples[sid].x.FNs[c];
              count[c] += sample.test_examples[sid].x.count[c];
            }
          }

          outputScore(testing_score_filename.c_str(), sparm->classIdxToLabel, TPs, FPs, FNs,
                      count, soutWeightFile.str().c_str());
            free(TPs);
            free(FPs);
            free(FNs);
            free(count);
        }
      }

      if(validationDir != "" && isDirectory(validationDir.c_str()))
        {
          string soutValidationOutputDir = inferenceDir_validation;
          if(!useSlice3d) {
            soutValidationOutputDir += soutOutput.str() + "/";
            mkdir(soutValidationOutputDir.c_str(), 0777);
          }

          SSVM_PRINT("[SVM_struct] Running inference for validation dataset containing %ld images\n",
                        nValidationExamples); 
          SSVM_PRINT("[SVM_struct] Output results will be stored in %s\n",
                        soutValidationOutputDir.c_str());

	  clock_t t = clock() - time_0;
	  SSVM_PRINT("[svm_struct] validation segmentImage::START %ld %f\n", t,
                     t/(float)CLOCKS_PER_SEC);

          string validation_score_filename = scoreDir + "validation_score.txt";

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
          for(int sid = 0; sid < nValidationExamples; sid++)
            {
              stringstream sout;
              sout << soutValidationOutputDir.c_str();
              if(useSlice3d) {
                sout << getNameFromPathWithoutExtension(validation_examples[sid].x.slice->getName());
                sout << "_";
                sout << iterationId;
              }

              if(useSlice3d) {
                segmentImage(validation_examples[sid].x,
                             sout.str().c_str(),
                             sparm->giType,
                             soutWeightFile.str().c_str(),
                             &sparm->labelToClassIdx,
                             validation_score_filename,
                             0,0,
                             compress_image,
                             (exportOverlay)?"overlay_validation/":0,
                             sparm->metric_type);
              } else {
                for(int c = 0; c < sparm->nClasses; c++) {
                  sample.validation_examples[sid].x.TPs[c] = 0;
                  sample.validation_examples[sid].x.FPs[c] = 0;
                  sample.validation_examples[sid].x.FNs[c] = 0;
                  sample.validation_examples[sid].x.count[c] = 0;
                }

                Slice* slice = static_cast<Slice*>(sample.validation_examples[sid].x.slice);
                computeScore(slice,
                             sample.validation_examples[sid].x.feature,
                             sample.validation_examples[sid].x.imgAnnotation,
                             sparm->giType,
                             param,
                             sparm->labelToClassIdx,
                             sparm->classIdxToLabel,
                             true,
                             param.nClasses,
                             sample.validation_examples[sid].x.TPs, sample.validation_examples[sid].x.FPs,
                             sample.validation_examples[sid].x.FNs, sample.validation_examples[sid].x.count,
                             soutValidationOutputDir.c_str(),
                             validation_examples[sid].x.nodeCoeffs,
                             validation_examples[sid].x.edgeCoeffs);
              }
            }
	  t = clock() - time_0;
	  SSVM_PRINT("[svm_struct] validation segmentImage::END %ld %f\n", t,
                     t/(float)CLOCKS_PER_SEC);

          if(!useSlice3d) {
            ulong *TPs = (ulong *)my_malloc(sizeof(ulong)*sparm->nClasses);
            ulong *FPs = (ulong *)my_malloc(sizeof(ulong)*sparm->nClasses);
            ulong *FNs = (ulong *)my_malloc(sizeof(ulong)*sparm->nClasses);
            ulong *count = (ulong *)my_malloc(sizeof(ulong)*sparm->nClasses);
            for(int c = 0; c < sparm->nClasses; c++) {
              TPs[c] = 0;
              FPs[c] = 0;
              FNs[c] = 0;
              count[c] = 0;
            }
            for(int sid = 0; sid < nValidationExamples; sid++) {
              for(int c = 0; c < sparm->nClasses; c++) {
                TPs[c] += sample.validation_examples[sid].x.TPs[c];
                FPs[c] += sample.validation_examples[sid].x.FPs[c];
                FNs[c] += sample.validation_examples[sid].x.FNs[c];
                count[c] += sample.validation_examples[sid].x.count[c];
              }
            }

            outputScore(validation_score_filename.c_str(), sparm->classIdxToLabel, TPs, FPs, FNs,
                        count, soutWeightFile.str().c_str());
            free(TPs);
            free(FPs);
            free(FNs);
            free(count);
          }
        }

      oldC = sparm->C;
      if(use_sampling) {
        sparm->giType = T_GI_SAMPLING;
      }

    }

  SSVM_PRINT("------------------------------------------------NEW ITERATION\n");
  sparm->iterationId++;

  return(0);
}

void        print_struct_learning_stats(SAMPLE sample, STRUCTMODEL *sm,
					CONSTSET cset, double *alpha, 
					STRUCT_LEARN_PARM *sparm)
{
  /* This function is called after training and allows final touches to
     the model sm. But primarly it allows computing and printing any
     kind of statistic (e.g. training error) you might want. */
}

void        print_struct_testing_stats(SAMPLE sample, STRUCTMODEL *sm,
				       STRUCT_LEARN_PARM *sparm, 
				       STRUCT_TEST_STATS *teststats)
{
  /* This function is called after making all test predictions in
     svm_struct_classify and allows computing and printing any kind of
     evaluation (e.g. precision/recall) you might want. You can use
     the function eval_prediction to accumulate the necessary
     statistics for each prediction. */
}

void        eval_prediction(long exnum, EXAMPLE ex, LABEL ypred, 
			    STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm, 
			    STRUCT_TEST_STATS *teststats)
{
  /* This function allows you to accumlate statistic for how well the
     predicition matches the labeled example. It is called from
     svm_struct_classify. See also the function
     print_struct_testing_stats. */
  if(exnum == 0) { /* this is the first time the function is
		      called. So initialize the teststats */
  }
}

void        write_struct_model(char *file, STRUCTMODEL *sm, 
			       STRUCT_LEARN_PARM *sparm)
{
  /* Writes structural model sm to file file. */
  double* smw = sm->w+1;

#if VERBOSITY > 3
  SSVM_PRINT("[SVM_struct] write_struct_model w=[");
  for(int i = 0; i < sm->sizePsi; i++) {
    SSVM_PRINT(" %g", smw[i]);
  }
  SSVM_PRINT("]\n");
#endif

  SSVM_PRINT("[SVM_struct] Writing struct_model to %s\n", file);
  ofstream ofs(file);
  ofs << sm->sizePsi << endl;
  ofs << sparm->nClasses << endl;
  ofs << sparm->nGradientLevels << endl;
  ofs << sparm->nOrientations << endl;
  ofs << sparm->nScales << endl;
  ofs << sparm->nLocalScales << endl;
  ofs << sparm->nScalingCoefficients << endl;
  for(int i = 0; i < sm->sizePsi; i++) {
    ofs << smw[i] << endl;
  }
  ofs.close();
}

void        write_label(FILE *fp, LABEL y)
{
  /* Writes label y to file handle fp. */
} 

void        free_pattern(SPATTERN x) {
  /* Frees the memory of x. */
  delete x.slice;
  delete x.feature;
  if(x.imgAnnotation != 0) {
    cvReleaseImage(&x.imgAnnotation);
  }
  delete[] x.TPs;
  delete[] x.FPs;
  delete[] x.FNs;
  if(x.count != 0) {
    delete[] x.count;
  }
  if(x.nodeCoeffs) {
    delete x.nodeCoeffs;
  }
}

void        free_label(LABEL y) {
  /* Frees the memory of y. */
  if(!y.cachedNodeLabels) {
    delete[] y.nodeLabels;
  }
}

void        free_struct_model(STRUCTMODEL sm) 
{
  /* Frees the memory of model. */
  /* if(sm.w) free(sm.w); */ /* this is free'd in free_model */
  if(sm.svm_model) {
    free_model(sm.svm_model,1);
  } else {
    if(sm.w) {
      free(sm.w);
    }
  }

  /* add free calls for user defined data here */
  /*
  if(sm.wsum) {
    delete[] sm.wsum;
  }
  */
}

void        free_struct_sample(SAMPLE s)
{
  /* Frees the memory of sample s. */
  for(int i=0;i<s.n;i++) { 
    free_pattern(s.examples[i].x);
    free_label(s.examples[i].y);
  }
  free(s.examples);

  if(s.test_examples) {
    for(int i=0;i<s.nTest;i++) { 
      free_pattern(s.test_examples[i].x);
      free_label(s.test_examples[i].y);
    }
    free(s.test_examples);
  }

  if(s.validation_examples) {
    for(int i=0;i<s.nValidation;i++) { 
      free_pattern(s.validation_examples[i].x);
      free_label(s.validation_examples[i].y);
    }
    free(s.validation_examples);
  }

  for(int il = 0; il < maxBuffers; il++) {
    delete[] tempNodeLabels[il];
  }
  delete[] tempNodeLabels;
}

void        print_struct_help()
{
  /* Prints a help text that is appended to the common help text of
     svm_struct_learn. */
  printf("         --* string  -> custom parameters that can be adapted for struct\n");
  printf("                        learning. The * can be replaced by any character\n");
  printf("                        and there can be multiple options starting with --.\n");
}

void         parse_struct_parameters(STRUCT_LEARN_PARM *sparm)
{
  /* Parses the command line parameters that start with -- */
  int i;

  for(i=0;(i<sparm->custom_argc) && ((sparm->custom_argv[i])[0] == '-');i++) {
    switch ((sparm->custom_argv[i])[2]) 
      { 
      case 'a': i++; /* strcpy(learn_parm->alphafile,argv[i]); */ break;
      case 'e': i++; /* sparm->epsilon=atof(sparm->custom_argv[i]); */ break;
      case 'k': i++; /* sparm->newconstretrain=atol(sparm->custom_argv[i]); */ break;
      default: printf("\nUnrecognized option %s!\n\n",sparm->custom_argv[i]);
	       exit(0);
      }
  }
}

void        print_struct_help_classify()
{
  /* Prints a help text that is appended to the common help text of
     svm_struct_classify. */
  printf("         --* string -> custom parameters that can be adapted for struct\n");
  printf("                       learning. The * can be replaced by any character\n");
  printf("                       and there can be multiple options starting with --.\n");
}

void         parse_struct_parameters_classify(STRUCT_LEARN_PARM *sparm)
{
  /* Parses the command line parameters that start with -- for the
     classification module */
  int i;

  for(i=0;(i<sparm->custom_argc) && ((sparm->custom_argv[i])[0] == '-');i++) {
    switch ((sparm->custom_argv[i])[2]) 
      { 
      /* case 'x': i++; strcpy(xvalue,sparm->custom_argv[i]); break; */
      default: printf("\nUnrecognized option %s!\n\n",sparm->custom_argv[i]);
	       exit(0);
      }
  }
}

void get_best_parameter_vector_id(int& best_idx)
{
  const int score_idx = 9;
  string training_score_filename = scoreDir + "training_score.txt";
  ifstream ifs(training_score_filename.c_str());
  if(ifs.fail()) {
    printf("[svm_struct] Error while loading %s\n", training_score_filename.c_str());
    best_idx = -1;
    return;
    //exit(-1);
  }
  string line;
  getline(ifs, line); // skip first line

  double best_score = 0;
  int idx = 1;
  best_idx = 1;
  while(getline(ifs, line)) {
    vector<string> tokens;
    splitStringUsing(line, tokens, SEPARATOR);
    //assert(tokens.size() >= (score_idx-1));
    if(tokens.size() < (score_idx-1)) {
      printf("[svm_struct] tokens.size()=%ld score_idx=%d\n", tokens.size(), score_idx);
      exit(-1);
    }
    double score = atof(tokens[score_idx].c_str());
     if(score > best_score) {
      best_score = score;
      best_idx = idx;
    }
    ++idx;
  }
  ifs.close();

  string config_tmp;
  Config* config = Config::Instance();
  int stepForOutputFiles = 1;
  if(config->getParameter("stepForOutputFiles", config_tmp)) {
    stepForOutputFiles = atoi(config_tmp.c_str());
    if(stepForOutputFiles == 0) {
      printf("[SVM_struct] Error: stepForOutputFiles = 0\n");
      exit(-1);
    }
  }

  printf("[svm_struct] best_idx %d %d\n", best_idx, stepForOutputFiles);
  best_idx *= stepForOutputFiles;
}

void get_best_parameter_vector(EnergyParam* param)
{
  int best_idx = 0;
  get_best_parameter_vector_id(best_idx);
  if(best_idx != -1) {
    stringstream parameterFile;
    parameterFile << parameterDir << "iteration_";
    parameterFile << best_idx;
    parameterFile << ".txt";
    SSVM_PRINT("[svm_struct] Best parameter vector = %s\n", parameterFile.str().c_str());
    param = new EnergyParam(parameterFile.str().c_str());
  }
}

void finalize()
{
  string outputModel = "model.txt";
  if(Config::Instance()->getParameter("outputModel", outputModel)) {
    SSVM_PRINT("[svm_struct] outputModel = %s\n", outputModel.c_str());
  } else {
    outputModel = "model.txt";
  }

  // pick best model
  int best_idx = 0;
  get_best_parameter_vector_id(best_idx);

  /*
  stringstream cmd;
  cmd << "cp parameter_vector0/iteration_";
  cmd << best_idx;
  cmd << ".txt ";
  cmd << outputModel;
  SSVM_PRINT("[svm_struct] Running command %s\n", cmd.str().c_str());
  int i = system(cmd.str().c_str());
  SSVM_PRINT("[svm_struct] Command executed. Returned code = %d\n", i);
  */

  stringstream src_filename;
  src_filename << "cp parameter_vector0/iteration_";
  src_filename << best_idx;
  src_filename << ".txt ";
  SSVM_PRINT("[svm_struct] Copying %s to %s\n",
             src_filename.str().c_str(), outputModel.c_str());
  copyFile(src_filename.str().c_str(), outputModel.c_str());

}
