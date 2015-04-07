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
// Written and (C) by Aurelien Lucchi                                  //
// Contact <aurelien.lucchi@gmail.com> for comments & bug reports      //
/////////////////////////////////////////////////////////////////////////


#include "Config.h"
#include "inference.h"
#include "energyParam.h"
#include "graphInference.h"
#include "Slice.h"
#include "gi_sampling.h"
#include "gi_max.h"
#include "gi_MF.h"
#include "utils.h"
#include "globalsE.h"

#ifdef _WIN32
#include "direct.h"
#define mkdir(x) _mkdir(x)
#endif

// inference libraries
#if USE_LIBDAI
#include "gi_libDAI.h"
#endif
#if USE_MAXFLOW
#include "gi_maxflow.h"
#endif
#if USE_MULTIOBJ
#include "gi_multiobject.h"
#endif

#include <sstream>
#include <vector>

#ifdef _WIN32
#include <io.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#endif
#include "inference_globals.h" //OPENMP flag

using namespace std;

#define SEPARATOR "\t"

//---------------------------------------------------------------------FUNCTIONS

GraphInference* createGraphInferenceInstance(int algo_type,
                                             Slice_P* slice,
                                             const EnergyParam& param,
                                             Feature* feature)
{
  return createGraphInferenceInstance(algo_type, slice, param, feature, 0, 0);
}

GraphInference* createGraphInferenceInstance(int algo_type,
                                             Slice_P* slice,
                                             const EnergyParam& param,
                                             Feature* feature,
                                             labelType* groundTruthLabels,
                                             double* lossPerLabel)
{
  return createGraphInferenceInstance(algo_type, slice, param, feature,
                                      groundTruthLabels, lossPerLabel, 0);
}

GraphInference* createGraphInferenceInstance(int algo_type,
                                             Slice_P* slice,
                                             const EnergyParam& param,
                                             Feature* feature,
                                             labelType* groundTruthLabels,
                                             double* lossPerLabel,
                                             map<sidType, nodeCoeffType>* _nodeCoeffs)
{
  return createGraphInferenceInstance(algo_type, slice, param, feature,
                                      groundTruthLabels, lossPerLabel,
                                      _nodeCoeffs, 0);
}

GraphInference* createGraphInferenceInstance(int algo_type,
                                             Slice_P* slice,
                                             const EnergyParam& param,
                                             Feature* feature,
                                             labelType* groundTruthLabels,
                                             double* lossPerLabel,
                                             map<sidType, nodeCoeffType>* _nodeCoeffs,
                                             map<sidType, edgeCoeffType>* _edgeCoeffs)
{

#if USE_LONG_RANGE_EDGES
  const int nDistances = DEFAULT_LONG_RANGE_EDGES_DISTANCE;
#endif

  string config_tmp;
  bool useGCForSubModularEnergy = false;
  if(Config::Instance()->getParameter("useGCForSubModularEnergy", config_tmp)) {
    if(atoi(config_tmp.c_str())) {
      useGCForSubModularEnergy = true;
    }
  }
  bool useGC = false;
  double d = 0; // diagonal
  double a = 0; //anti-diagonal
  if(useGCForSubModularEnergy) {
    // For 2 classes, use graph-cuts if pairwise potential is attractive
    if(param.includeLocalEdges && param.nClasses == 2) {
      double* pw = param.weights + param.nUnaryWeights;
      if(param.nGradientLevels == 0) {
        useGC = true;
      } else {
        useGC = true;
#if USE_LONG_RANGE_EDGES
        for(int p = 0; useGC && (p < nDistances); ++p) {
          a = 0;
          d = 0;
#else
          {
#endif
            for(int g = 0; g < param.nGradientLevels; ++g) {
              d += pw[0] + pw[3];
              a += pw[1] + pw[2];
              if(a > d) {
                useGC = false;
                break;
              }
              pw += param.nClasses*param.nClasses;
            }
          }
        }
    }
  }

#ifdef USE_MAXFLOW
  if(useGC && algo_type != T_GI_OPENGM && algo_type != T_GI_MAX && algo_type != T_GI_MULTIOBJ) {
    algo_type = T_GI_MAXFLOW;
  }
#else
  printf("[inference] Warning: Energy function is submodular but maxflow library was not enabled!");
#endif

  if(!param.includeLocalEdges) {
    algo_type = T_GI_MAX;
  }

  GraphInference* gi = 0;
  switch(algo_type)
    {
#if USE_LIBDAI
    case T_GI_LIBDAI:
    case T_GI_LIBDAI_ICM:
    case T_GI_LIBDAI_ICM_QPBO:
      gi = new GI_libDAI(slice,
                         &param,
                         param.weights,
                         groundTruthLabels,
                         lossPerLabel,
                         feature,
                         _nodeCoeffs,
                         _edgeCoeffs);
      break;
#endif

#if USE_MAXFLOW
    case T_GI_MAXFLOW:
      gi = new GI_maxflow(slice,
                          &param,
                          param.weights,
                          groundTruthLabels,
                          lossPerLabel,
                          feature,
                          _nodeCoeffs,
                          _edgeCoeffs);
      break;
#endif

    case T_GI_MF:
      gi = new GI_MF(slice,
                     &param,
                     param.weights,
                     groundTruthLabels,
                     lossPerLabel,
                     feature,
                     _nodeCoeffs);
      break;

    case T_GI_MAX:
      {
        gi = new GI_max(slice,
                        &param,
                        param.weights,
                        groundTruthLabels,
                        lossPerLabel,
                        feature,
                        _nodeCoeffs);
    }
    break;

#if USE_MULTIOBJ
    case T_GI_MULTIOBJ:
      {
        gi = new GI_multiobject(slice,
                                &param,
                                param.weights,
                                groundTruthLabels,
                                lossPerLabel,
                                feature,
                                _nodeCoeffs,
                                _edgeCoeffs);
      }
      break;
#endif

    case T_GI_SAMPLING:      
        gi = new GI_sampling(slice,
                             &param,
                             param.weights,
                             groundTruthLabels,
                             lossPerLabel,
                             feature,
                             _nodeCoeffs,
                             1.0);     
      break;

    default:
      printf("[inference] Unknown algo type %d\n", algo_type);
      exit(-1);
      break;
    }

  return gi; 
}

void postprocessBoundaryLabels(Slice_P* g, labelType* nodeLabels)
{
  // check if there is any boundary label surrounded by foreground
  // and change them to foreground.
  ulong nReplaced = 0;
  const map<sidType, supernode* >& _supernodes = g->getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    int sid = it->first;
    int label = nodeLabels[sid];
    if(label == BOUNDARY) {
      supernode* s = it->second;
      bool onlyForeground = true;
      vector < supernode* >* lNeighbors = &(s->neighbors);
      for(vector < supernode* >::iterator itN = lNeighbors->begin();
          itN != lNeighbors->end(); itN++) {

        // set edges once
        //if(it->first < (*itN)->id) {
        //  continue;
        //}

        supernode* sn = (*itN);
        int labeln = sn->getLabel();
        if(labeln != FOREGROUND) {
          onlyForeground = false;
          break;
        }
      }
      if(onlyForeground) {
        // change label to foreground
        nodeLabels[sid] = FOREGROUND;
        ++nReplaced;
      }
    }
  }
}

labelType* computeLabels(Slice_P* g, Feature* feature, const EnergyParam& param,
                         int algoType, double* energy)
{
  return computeLabels(g, feature, param, algoType, energy, 0, 0, 0, 0);
}

labelType* computeLabels(Slice_P* g, Feature* feature, const EnergyParam& param,
                         int algoType, double* energy,
                         labelType* groundTruthLabels, double* lossPerLabel)
{
  return computeLabels(g, feature, param, algoType, energy, groundTruthLabels,
                       lossPerLabel, 0, 0);
}

labelType* computeLabels(Slice_P* g, Feature* feature, const EnergyParam& param,
                         int algoType, double* energy,
                         labelType* groundTruthLabels, double* lossPerLabel,
                         map<sidType, nodeCoeffType>* _nodeCoeffs,
                         map<sidType, edgeCoeffType>* _edgeCoeffs)
{
  if(algoType == T_GI_SAMPLING) {
    return computeLabels_sampling(g, feature, param, algoType, energy,
                                  groundTruthLabels, lossPerLabel, _nodeCoeffs, _edgeCoeffs);
  } else {
    GraphInference* gi_Inference =
      createGraphInferenceInstance(algoType, g, param, feature, groundTruthLabels,
                                   lossPerLabel, _nodeCoeffs, _edgeCoeffs);
    size_t maxiter = 100;
    ulong nNodes = g->getNbSupernodes();
    labelType* nodeLabels = new labelType[nNodes];
    gi_Inference->run(nodeLabels, 0, maxiter);
    if(param.nClasses == 3) {
      postprocessBoundaryLabels(g, nodeLabels);
    }

    if(energy) {
      *energy = gi_Inference->computeEnergy(nodeLabels);
    }

    delete gi_Inference;
    return nodeLabels;
  }
}

// compute an estimate of the score
// return max among a subset of sampled superpixels and only compute score for class 0
double compute_score(Slice_P* slice, Feature* feature, const EnergyParam& param,
                     labelType* groundTruthLabels, double* lossPerLabel)
{
  const int c = 0; // class 0

  GraphInference gi(slice, (const EnergyParam*)&param, param.weights, feature, 0, 0);

  double max_potential = 0;

  string config_tmp;
  double sampling_rate = 1.0;
  if(Config::Instance()->getParameter("sampling_rate", config_tmp)) {
    sampling_rate = atof(config_tmp.c_str()); 
    printf("[gi_sampling] sampling_rate = %g\n", sampling_rate);
  }
  bool draw_samples = sampling_rate != 1.0;
  ulong nSupernodes = slice->getNbSupernodes();
  int sid = 0;
  for(int i = 0; i < nSupernodes*sampling_rate; ++i) {

    if(draw_samples) {
      // Select a pixel at random
      sid = rand() * ((double)nSupernodes/(double)RAND_MAX);
    } else {
      sid = i;
    }

    supernode* s = slice->getSupernode(sid);
    vector < supernode* >* lNeighbors = &(s->neighbors);

    double potential = gi.computeUnaryPotential(slice, sid, c);
    for(vector < supernode* >::iterator itN = lNeighbors->begin();
        itN != lNeighbors->end(); itN++) {
      const int c2 = rand()*param.nClasses / (double)RAND_MAX;
      double pairwisePotential = gi.computePairwisePotential(slice, s, (*itN),
                                                             c, c2);
      potential += pairwisePotential;
    }

    if(lossPerLabel) {
      if(c != groundTruthLabels[sid]) {
        // add loss of the ground truth label
        potential += lossPerLabel[groundTruthLabels[sid]];
      }
    }

    if(potential > max_potential) {
      max_potential = potential;
    }
  }

  return max_potential;
}

labelType* computeLabels_sampling(Slice_P* g, Feature* feature, const EnergyParam& param,
                                  int algoType, double* energy,
                                  labelType* groundTruthLabels, double* lossPerLabel,
                                  map<sidType, nodeCoeffType>* _nodeCoeffs,
                                  map<sidType, edgeCoeffType>* _edgeCoeffs)
{
  int nParallelChains = 1;
  string config_tmp;
  if(Config::Instance()->getParameter("sampling_nParallelChains", config_tmp)) {
    nParallelChains = atoi(config_tmp.c_str());
  }
  printf("[inference] nParallelChains=%d\n", nParallelChains);

  double temperature_0 = 0.001;
  if(Config::Instance()->getParameter("sampling_initial_temperature", config_tmp)) {
    temperature_0 = atof(config_tmp.c_str()); 
  }

  double score = compute_score(g, feature, param, groundTruthLabels,
                               lossPerLabel);
  printf("[inference] score=%g\n", score);

  // run several chains in parallel
  double *energies = new double[nParallelChains];
  labelType** inferredLabels = new labelType*[nParallelChains];
      
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int iChain = 0; iChain < nParallelChains; ++iChain) {

    GI_sampling* gi_sampling
      = new GI_sampling(g, &param, param.weights, groundTruthLabels, lossPerLabel,
                        feature, _nodeCoeffs, 1.0);

    inferredLabels[iChain] = new labelType[g->getNbSupernodes()];
    gi_sampling->setInitializedLabels(false);

    double temperature = temperature_0/(pow(10.0,iChain));

    energies[iChain] = gi_sampling->runOnce(inferredLabels[iChain],
                                            MVC_MAX_ITER, groundTruthLabels,
                                            0, temperature);

    delete gi_sampling;
  }

  // copy the best chain
  double minEnergy = energies[0];
  int bestChain = 0;
  for(int iChain = 0; iChain < nParallelChains; ++iChain) {
    double temperature = temperature_0/(pow(10.0,iChain));
    printf("[inference] sampling (%d, %g) energy = %g\n",
           iChain, temperature, energies[iChain]);
    if(energies[iChain] < minEnergy) {
      minEnergy = energies[iChain];
      bestChain = iChain;
    }
  }

  if(energy) {
    *energy = minEnergy;
  }

  for(int iChain = 0; iChain < nParallelChains; ++iChain) {
    if(iChain != bestChain) {
      delete[] inferredLabels[iChain];
    }
  }

  double temperature = temperature_0/(pow(10.0,bestChain));
  printf("[inference] best sampling (%d, %g) energy = %g\n",
         bestChain, temperature, energies[bestChain]);

  delete[] energies;
  return inferredLabels[bestChain];
}


labelType* computeCombinedLabels(Slice_P* g, Feature* feature,
                                 labelType* groundTruthLabels, double* lossPerLabel,
                                 int nRFs, vector<double>* alphas, int example_id,
                                 map<sidType, nodeCoeffType>* _nodeCoeffs,
                                 map<sidType, edgeCoeffType>* _edgeCoeffs,
                                 const char* combined_probability_output_file)
{
  const int maxiter = 100;
  int nNodes = g->getNbSupernodes();
  int nClasses = 0;

  labelType** nodeLabels = new labelType*[nRFs];
  for(int i = 0; i < nRFs; ++i) {
    stringstream parameter_vector_dir;
    parameter_vector_dir << "parameter_vector";
    parameter_vector_dir << i;
    string parameter_vector_file_pattern = parameter_vector_dir.str() + "/iteration_";
    int idx = 1; // first file starts with index 1
    string parameter_vector_file_last = findLastFile(parameter_vector_file_pattern, ".txt", &idx);

    printf("[inference] Running inference with %s\n", parameter_vector_file_last.c_str());
    EnergyParam param(parameter_vector_file_last.c_str());
    nClasses = param.nClasses;

    GraphInference* gi_Inference =
      createGraphInferenceInstance(T_GI_LIBDAI, g, param, feature, groundTruthLabels,
                                   lossPerLabel, _nodeCoeffs, _edgeCoeffs);
    nodeLabels[i] = new labelType[nNodes];
    gi_Inference->run(nodeLabels[i], 0, maxiter);
    delete gi_Inference;
  }

  printf("[inference] Combining models\n");
  float* combined_predictions = new float[nNodes];
  float sum_alphas = 0;
  for(int i = 0; i < nRFs; ++i) {
    sum_alphas += alphas[example_id][i];
  }
  for(int n = 0; n < nNodes; ++n) {
    combined_predictions[n] = 0;
    for(int i = 0; i < nRFs; ++i) {
      combined_predictions[n] += alphas[example_id][i]*nodeLabels[i][n];
    }
    combined_predictions[n] /= nRFs;
  }

  if(combined_probability_output_file != 0) {
    g->exportProbabilities(combined_probability_output_file,
                           nClasses, combined_predictions);
  }

  float threshold = sum_alphas/3.0;
  labelType* combined_labels = new labelType[nNodes];
  for(int n = 0; n < nNodes; ++n) {
    combined_labels[n] = (combined_predictions[n] > threshold);
  }

  for(int i = 0; i < nRFs; ++i) {
    delete[] nodeLabels[i];
  }
  delete[] nodeLabels;
  delete[] combined_predictions;
  //delete[] combined_labels;
  return combined_labels;
}

IplImage* getSegmentedImage(Slice* g,
                            Feature* feature,
                            int algoType,
                            const EnergyParam& param,
                            map<labelType, ulong>* labelToClassIdx)
{
  return getSegmentedImage(g, feature, algoType, param,
                           labelToClassIdx, 0, 0);
}

IplImage* getSegmentedImage(Slice* g,
                            Feature* feature,
                            int algoType,
                            const EnergyParam& param,
                            map<labelType, ulong>* labelToClassIdx,
                            map<sidType, nodeCoeffType>* _nodeCoeffs,
                            map<sidType, edgeCoeffType>* _edgeCoeffs)
{
  labelType* nodeLabels = computeLabels(g, feature, param, algoType, 0, 0, 0,
                                        _nodeCoeffs, _edgeCoeffs);

  int nNodes = g->getNbSupernodes();
  IplImage* img = g->getColoredAnnotationImage(param.nClasses,
                                               nodeLabels,
                                               nNodes,
                                               labelToClassIdx);

  delete[] nodeLabels;
  return img;
}

double segmentImage(SPATTERN x,
                    const char* output_file,
                    int algoType,
                    const char* weight_file,
                    map<labelType, ulong>* labelToClassIdx,
                    const string& output_roc_file,
                    labelType* groundTruthLabels,
                    double* lossPerLabel,
                    const bool compress_image,
                    const char* overlay_dir,
                    const int metric_type)
{
  Slice_P* g = x.slice;
  Feature* feature = x.feature;
  map<sidType, nodeCoeffType>* _nodeCoeffs = x.nodeCoeffs;
  map<sidType, edgeCoeffType>* _edgeCoeffs = x.edgeCoeffs;

  EnergyParam param(weight_file);
  bool output_is_dir = isDirectory(output_file);
  string output_dir = getDirectoryFromPath(output_file);
  double energy = 0;
  labelType* nodeLabels = computeLabels(g, feature, param, algoType, 0,
                                        groundTruthLabels, lossPerLabel,
                                        _nodeCoeffs, _edgeCoeffs);

  stringstream soutColoredImage;
  if(output_is_dir) {
    string stemp = getNameFromPathWithoutExtension(x.slice->getName());

    soutColoredImage << output_dir << "/";
    soutColoredImage << getNameFromPathWithoutExtension(x.slice->getName()) << ".png";
  } else {
    soutColoredImage << output_file;
  }

  bool deleteLabelMap = false;
  if(labelToClassIdx == 0) {
    labelToClassIdx = new map<labelType, ulong>;
    deleteLabelMap = true;

    string paramColormap;
    getColormapName(paramColormap);
    getLabelToClassMap(paramColormap.c_str(), *labelToClassIdx);
  }

  // output image
  int nNodes = g->getNbSupernodes();
  g->exportSupernodeLabels(soutColoredImage.str().c_str(),
                           param.nClasses,
                           nodeLabels,
                           nNodes,
                           labelToClassIdx);

  map<labelType, ulong> labelCount;
  g->countSupernodeLabels(nodeLabels, labelCount);
  for(map<labelType, ulong>::iterator it = labelCount.begin();
      it != labelCount.end(); ++it) {
    printf("%d=%ld ", it->first, it->second);
  }
  printf("\n");

  if(overlay_dir != 0) {
    mkdir(overlay_dir, 0777);
    stringstream soutOverlayImage;
    if(output_is_dir) {
      soutOverlayImage << overlay_dir << "/" << getNameFromPathWithoutExtension(x.slice->getName()) << "_overlay";
    } else {
      soutOverlayImage << overlay_dir << "/" << getNameFromPathWithoutExtension(output_file) << "_overlay";
    }
    if(g->getType() != SLICEP_SLICE3D) {
      soutOverlayImage << ".png";
    }
    
    PRINT_MESSAGE("[inference] Exporting overlay to %s\n", soutOverlayImage.str().c_str());
    g->exportOverlay(soutOverlayImage.str().c_str(), nodeLabels);
    
    if(g->getType() == SLICEP_SLICE3D) {
      zipAndDeleteCube(soutOverlayImage.str().c_str());
    }
  }

  eSlicePType sliceType = g->getType();
  if(sliceType == SLICEP_SLICE3D) {
    if(!output_roc_file.empty()) {
      Slice3d* slice3d = static_cast<Slice3d*>(g);
      if(slice3d->isSupernodeLabelsLoaded()) {

        float score;
        float score_bg;
        float score_fg;
        float true_pos;
        float true_neg;
        float false_neg;
        float false_pos;
        ulong total_pos = 0;
        ulong total_neg = 0;

        const bool useColorAnnotations = false;

        switch(metric_type) {
        case METRIC_SUPERNODE_BASED_01:
          compareMultiLabelVolumes(*slice3d,
                                   nodeLabels,
                                   BACKGROUND,
                                   true_neg, true_pos, false_neg, false_pos,
                                   false,
                                   useColorAnnotations,
                                   &total_pos, &total_neg); // normalization
          break;
        case METRIC_NODE_BASED_01:
          compareMultiLabelVolumes_nodeBased(*slice3d,                                             
                                             x.cubeAnnotation,
										     nodeLabels,
                                             BACKGROUND,
                                             true_neg, true_pos, false_neg, false_pos,
                                             false,
                                             useColorAnnotations,
                                             &total_pos, &total_neg); // normalization
          break;
        }

        score_bg = true_pos/(true_pos+false_neg+false_pos);
        score_fg = true_neg/(true_neg+false_neg+false_pos);
        score = (score_fg + score_bg)/2.0;
        float accuracy = (true_pos + true_neg) / (total_pos + total_neg);

        PRINT_MESSAGE("[inference] true_pos=%.2g, false_neg=%.2g, false_pos=%.2g, true_neg=%.2g, score=%f\n",
                      true_pos,false_neg,false_pos,true_neg,score);

        bool roc_file_exists = fileExists(output_roc_file.c_str());
        ofstream ofsRoc(output_roc_file.c_str(), ios::app);
        if(!roc_file_exists) {
          // Add header
          ofsRoc << "Name" << SEPARATOR << "TP" << SEPARATOR << "FN" << SEPARATOR;
          ofsRoc << "FP" << SEPARATOR << "TN" << SEPARATOR << "Total+" << SEPARATOR;
          ofsRoc << "Total-" << SEPARATOR << "Score_F" << SEPARATOR << "Score_B";
          ofsRoc << SEPARATOR << "Score" << SEPARATOR << "TPR" << SEPARATOR << "FNR";
          ofsRoc << SEPARATOR << "FPR" << SEPARATOR << "TNR" << SEPARATOR << "Accuracy" << endl;
        }

        ofsRoc << soutColoredImage.str() << SEPARATOR;
        ofsRoc << true_pos << SEPARATOR << false_neg << SEPARATOR;
        ofsRoc << false_pos << SEPARATOR << true_neg << SEPARATOR;
        ofsRoc << total_pos << SEPARATOR << total_neg << SEPARATOR;
        ofsRoc << score_fg << SEPARATOR << score_bg << SEPARATOR << score << SEPARATOR;

        if(total_pos != 0) {
          true_pos = true_pos*(100.0f/total_pos); // TPR = TP / P
          false_neg = false_neg*(100.0f/total_pos); // FNR = FN / P
        }
        if(total_neg != 0) {
          false_pos = false_pos*(100.0f/total_neg); // FPR = FP / N
          true_neg = true_neg*(100.0f/total_neg); // TNR = TN / N
        }

        ofsRoc << true_pos << SEPARATOR << false_neg << SEPARATOR << false_pos << SEPARATOR << true_neg;
        ofsRoc << SEPARATOR << accuracy << endl;
        ofsRoc.close();
      }

      if(compress_image && g->getType() == SLICEP_SLICE3D) {
        zipAndDeleteCube(soutColoredImage.str().c_str());
      }
    }
  }

  delete[] nodeLabels;
  if(deleteLabelMap)
    delete labelToClassIdx;

  return energy;
}

//-----------------------

void computeScore(Slice* g,
                  Feature* feature,
                  const IplImage* imgAnnotation,
                  int algoType,
                  EnergyParam& param,
                  map<labelType, ulong>& labelToClassIdx,
                  map<ulong, labelType>& classIdxToLabel,
                  bool convertToLabel,       
                  int nClasses,
                  ulong* TPs,
                  ulong* FPs,
                  ulong* FNs,
                  ulong* count,
                  const char* outputDir)
{
  return computeScore(g, feature, imgAnnotation, algoType, param,
                      labelToClassIdx, classIdxToLabel,
                      convertToLabel, nClasses,
                      TPs, FPs, FNs, count, outputDir, 0, 0);
}

void computeScore(Slice* g,
                  Feature* feature,
                  const IplImage* imgAnnotation,
                  int algoType,
                  EnergyParam& param,
                  map<labelType, ulong>& labelToClassIdx,
                  map<ulong, labelType>& classIdxToLabel,
                  bool convertToLabel,       
                  int nClasses,
                  ulong* TPs,
                  ulong* FPs,
                  ulong* FNs,
                  ulong* count,
                  const char* outputDir,
                  map<sidType, nodeCoeffType>* _nodeCoeffs,
                  map<sidType, edgeCoeffType>* _edgeCoeffs)
{
  IplImage* img = getSegmentedImage(g, feature, algoType, param,
                                    &labelToClassIdx, _nodeCoeffs, _edgeCoeffs);

#if 1
  string outputImageName(outputDir);
  outputImageName += "/";
  mkdir(outputImageName.c_str(), 0777);
  outputImageName += getNameFromPathWithoutExtension(g->getName());
  outputImageName += ".png";
  cvSaveImage(outputImageName.c_str(), img);
#endif

  if(imgAnnotation != 0) {
    computeScore(img, imgAnnotation, classIdxToLabel, convertToLabel,
                 nClasses, TPs, FPs, FNs, count);
  }

  cvReleaseImage(&img);
}

void computeScore(const char* image_dir,
                  const char* image_pattern,
                  const char* superpixelDir,
                  const char* maskDir,
                  int nImages,
                  const char* output_file,
                  map<ulong, labelType>& classIdxToLabel,
                  const char* weight_file)
{
  return computeScore(image_dir, image_pattern,
                      superpixelDir, maskDir,
                      nImages, output_file,
                      classIdxToLabel,
                      weight_file, 0, 0);
}

void computeScore(const char* image_dir,
                  const char* image_pattern,
                  const char* superpixelDir,
                  const char* maskDir,
                  int nImages,
                  const char* output_file,
                  map<ulong, labelType>& classIdxToLabel,
                  const char* weight_file,
                  map<sidType, nodeCoeffType>* _nodeCoeffs,
                  map<sidType, edgeCoeffType>* _edgeCoeffs)
{
  int nClasses = classIdxToLabel.size();
  ulong *TPs = new ulong[nClasses];
  ulong *FPs = new ulong[nClasses];
  ulong *FNs = new ulong[nClasses];
  ulong *count = new ulong[nClasses];
  for(int c=0; c < nClasses; c++) {
    TPs[c] = 0;
    FPs[c] = 0;
    FNs[c] = 0;
    count[c] = 0;
  }

  vector<string> lImages;
  if(isDirectory(image_dir)) {
    getFilesInDir(image_dir, lImages,image_pattern, true);
    INFERENCE_PRINT("[inference] Predicting %ld images in directory %s\n", lImages.size(),image_dir);
  } else {
    lImages.push_back(image_dir);
  }

  string paramSuperpixelLabelsDir;
  if(superpixelDir == 0) {
    Config::Instance()->getParameter("superpixel_labels", paramSuperpixelLabelsDir);
    superpixelDir = paramSuperpixelLabelsDir.c_str();
    INFERENCE_PRINT("[inference] superpixelDir=%s\n", superpixelDir);
  }

  if(nImages == -1) {
    nImages = lImages.size();
  }
  INFERENCE_PRINT("[inference] Number of valid images to be processed = %d\n", nImages);

  //bool convertToLabel = (weight_file == 0);
  bool convertToLabel = true;

  string colormapFilename;
  getColormapName(colormapFilename);

  printf("[inference] Colormap=%s\n", colormapFilename.c_str());
  map<labelType, ulong> labelToClassIdx;
  getLabelToClassMap(colormapFilename.c_str(), labelToClassIdx);

  vector<eFeatureType> feature_types;
  if(weight_file != 0) {
    int paramFeatureTypes = 0;
    string config_tmp;
    if(Config::Instance()->getParameter("featureTypes", config_tmp)) {
      paramFeatureTypes = atoi(config_tmp.c_str());
      getFeatureTypes(paramFeatureTypes, feature_types);
    }
  }

  EnergyParam* param = 0;
  if(weight_file != 0) {
    param = new EnergyParam(weight_file);
  }

  /*
#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  */
  for(int iImage = 0; iImage < nImages; iImage++) {
    string image_name = lImages[iImage];

    string groundtruthName;
    bool found_GT = getGroundTruthName(groundtruthName, maskDir, image_name);

    // skip image
    if(!found_GT) { 
      printf("[inference] Could not find %s in %s\n", groundtruthName.c_str(),
             maskDir);
      continue;
    }

    printf("[inference] Found %s in %s\n", groundtruthName.c_str(),
           maskDir);

    IplImage* imgAnnotation = cvLoadImage(groundtruthName.c_str());
    if(!imgAnnotation) {
      printf("[inference] Error: input mask %s was not found\n", groundtruthName.c_str());
      exit(-1);
    }

    IplImage* img = 0;
    if(weight_file != 0) {
      Slice slice(image_name.c_str());
      Feature* feature = Feature::getFeature(&slice, feature_types);

      img = getSegmentedImage(&slice,
                              feature,
                              T_GI_LIBDAI,
                              *param,
                              &labelToClassIdx);
      delete feature;
    } else {
      img = cvLoadImage(image_name.c_str());
    }

    computeScore(img, imgAnnotation, classIdxToLabel, convertToLabel,
                 nClasses, TPs, FPs, FNs, count);
    cvReleaseImage(&img);
    cvReleaseImage(&imgAnnotation);
  }

  outputScore(output_file, classIdxToLabel, TPs, FPs, FNs, count, weight_file);

  if(param) {
    delete param;
  }
  delete[] TPs;
  delete[] FPs;
  delete[] FNs;
  delete[] count;
}

void outputScore_ROC(const char* output_file, map<ulong, labelType>& classIdxToLabel,
                     ulong* TPs, ulong* FPs, ulong* FNs, ulong* count,
                     const char* weight_file)
{
  printf("[inference] Output results for %ld classes in %s\n",
         classIdxToLabel.size(), output_file);
  ofstream ofs(output_file, ios::app);
  if(weight_file != 0) {
    ofs << weight_file << SEPARATOR;
  }
  ofs << "class" << SEPARATOR << "count" << SEPARATOR << "TP" << SEPARATOR << "FP" << SEPARATOR << "FN" << SEPARATOR << "\n";
  for(map<ulong, labelType>::iterator it = classIdxToLabel.begin();
      it != classIdxToLabel.end(); ++it) {
    uint c = it->second;
    ofs << c << SEPARATOR << count[c] << SEPARATOR << TPs[c] << SEPARATOR << FPs[c] << SEPARATOR << FNs[c] << endl;
  }
  ofs.close();
}

// output score in output_file
void outputScore(const char* output_file, map<ulong, labelType>& classIdxToLabel,
                 ulong* TPs, ulong* FPs, ulong* FNs, ulong* count,
                 const char* weight_file)
{
  string config_tmp;
  Config::Instance()->getParameter("computeStats_ROC", config_tmp);
  bool computeStats_ROC = config_tmp.c_str()[0] == '1';
  if(computeStats_ROC) {
    outputScore_ROC(output_file, classIdxToLabel, TPs, FPs, FNs, count, weight_file);
    return;
  }

  string paramMSRC;
  Config::Instance()->getParameter("msrc", paramMSRC);
  bool useMSRC = paramMSRC.c_str()[0] == '1' && (classIdxToLabel.size() >= 21);

  labelType voidLabel = 0;
  labelType moutainLabel = 0;
  labelType horseLabel = 0;
  if(useMSRC) {
    voidLabel = classIdxToLabel[0];
    moutainLabel = classIdxToLabel[4161600];
    horseLabel = classIdxToLabel[8323328];
  }

  printf("[inference] Output results for %ld classes in %s\n",
         classIdxToLabel.size(), output_file);
  bool existingFile = fileExists(output_file);
  ofstream ofs(output_file, ios::app);
  if(!existingFile) {
    // add header to non existing file only
    for(uint c = 0; c < classIdxToLabel.size(); ++c) {
      ofs << c << " ";
    }
    ofs << "global ";
    ofs << "average";
    if(useMSRC) {
      ofs << " averageVOC";
    }
    ofs << endl;
  }

  // add name of parameter file if any
  if(weight_file != 0) {
    ofs << weight_file << " ";
  }

  double avgScore = 0;
  double avgVOCScore = 0;
  double avgScoreForC = 0;
  double globalScore = 0;
  ulong totalCount = 0;
  for(map<ulong, labelType>::iterator it = classIdxToLabel.begin();
      it != classIdxToLabel.end(); ++it) {
    uint c = it->second;
    //SSVM_PRINT("[SVM_struct] class=%d TP=%ld FP=%ld FN=%ld\n",c,x.TPs[c],x.FPs[c],x.FNs[c]);
    if(!useMSRC) {
      if(count[c]>0) {
        avgScoreForC = TPs[c]/(double)count[c];
        avgScore += avgScoreForC;
        globalScore += TPs[c];
        totalCount += count[c];
      }
      float d = TPs[c] + FPs[c] + FNs[c];
      avgVOCScore += TPs[c]/d;
      ofs << c << "=" << (TPs[c]/d)*100 << " & ";      
    } else {
      if(c != voidLabel && c != moutainLabel && c != horseLabel) {
        if(count[c]>0) {
          avgScoreForC = TPs[c]/(double)count[c];
          avgScore += avgScoreForC;
          globalScore += TPs[c];
          totalCount += count[c];
        }
        ofs << c << "=" << avgScoreForC*100 << " & ";
      }   
    }
  }

  int nClasses = classIdxToLabel.size();
  if(useMSRC) {
    avgScore /= (nClasses - 3);
  } else {
    avgScore /= nClasses;
  }
  globalScore /= totalCount;

  ofs << globalScore*100 << " & " << avgScore*100;
  if(useMSRC) {
    ofs << endl;
  } else {
    avgVOCScore /= nClasses;
    ofs << " & " << avgVOCScore*100 << endl;
  }
  ofs.close();
}

void computeScore(const IplImage* img,
                  const IplImage* imgAnnotation,
                  map<ulong, labelType>& classIdxToLabel,
                  bool convertToLabel,       
                  int nClasses,
                  ulong* TPs,
                  ulong* FPs,
                  ulong* FNs,
                  ulong* count)
{
  uchar* ptrM;
  uchar* ptrI;

  //printf("names %s %s %d %d %d %d\n", img_name.c_str(), groundtruthName.c_str(),img->width,img->height,img->nChannels,imgAnnotation->nChannels);

  // Compute TP,FP and FN used to compute VOC score
  for(int u=0; u < img->width; u++) {
    for(int v=0; v < img->height; v++) {
      ptrI = &((uchar*)(img->imageData + img->widthStep*v))[u*img->nChannels];
      ptrM = &((uchar*)(imgAnnotation->imageData + imgAnnotation->widthStep*v))[u*imgAnnotation->nChannels];

      int labelM = 0;
      int labelI = 0;
      if(convertToLabel) {
        ulong classIdx_M = ((ulong)ptrM[2]*65025) + ptrM[1]*255 + ptrM[0];
        ulong classIdx_I = ((ulong)ptrI[2]*65025) + ptrI[1]*255 + ptrI[0];

        /*
        printf("class->label(%d) (%d,%d,%d)=%ld (%d,%d,%d)=%ld\n", nClasses,
               (int)ptrM[0],(int)ptrM[1],(int)ptrM[2],classIdx_M,
               (int)ptrI[0],(int)ptrI[1],(int)ptrI[2],classIdx_I);
        */

        assert(classIdxToLabel.find(classIdx_M) != classIdxToLabel.end());
        assert(classIdxToLabel.find(classIdx_I) != classIdxToLabel.end());
        labelM = classIdxToLabel[classIdx_M];
        labelI = classIdxToLabel[classIdx_I];

        /*
        printf("class->label(%d) (%d,%d,%d)=%ld->%d (%d,%d,%d)=%ld->%d\n", nClasses,
               (int)ptrM[0],(int)ptrM[1],(int)ptrM[2],classIdx_M, labelM,
               (int)ptrI[0],(int)ptrI[1],(int)ptrI[2],classIdx_I, labelI);
        */

      } else {
        labelM = ptrM[0];
        labelI = ptrI[0];

        /*
        printf("class->label(%d) (%d,%d,%d)->%d (%d,%d,%d)->%d\n", nClasses,
               (int)ptrM[0],(int)ptrM[1],(int)ptrM[2], labelM,
               (int)ptrI[0],(int)ptrI[1],(int)ptrI[2], labelI);
        */
      }

      assert(labelI<nClasses);
      assert(labelM<nClasses);

      if(labelM == labelI) {
        // +1 true positive
        TPs[labelM]++;
      } else {
        // +1 false positive for predicted class
        FPs[labelI]++;
        // +1 false negative for groundtruth class
        FNs[labelM]++;
      }
      
      count[labelM]++;
    }
  }
}

void computeVOCLoss(labelType* nodeLabels, labelType* ybar,
                    int nNodes, int nClasses,
                    double& loss, int& nDiff,
                    ulong* TPs, ulong* FPs, ulong* FNs) {
  // Compute TP,FP and FN used to compute VOC score
  loss  = 0;
  nDiff = 0;
  for(int c = 0; c < nClasses; c++) {
    TPs[c] = 0;
    FPs[c] = 0;
    FNs[c] = 0;
  }
  set<uint> lClasses;
  int g; // ground truth
  int p; // prediction
  for(int n = 0; n < nNodes; n++) {
    g = nodeLabels[n];
    p = ybar[n];
    if(p == g) {
      // +1 true positive
      TPs[g]++;
    } else {
      // +1 false positive for class p
      FPs[p]++;
      // +1 false negative for class g
      FNs[g]++;
    }
      
    if(lClasses.count(p) == 0)
      lClasses.insert(p);
    if(lClasses.count(g) == 0)
      lClasses.insert(g);
  }

  double score = 0;
  for(set<uint>::iterator it = lClasses.begin(); it != lClasses.end(); it++) {
    uint c = *it;
    //SSVM_PRINT("[SVM_struct] example=%d class=%d TP=%ld FP=%ld FN=%ld\n",x.id,c,x.TPs[c],x.FPs[c],x.FNs[c]);
    double d = TPs[c] + FPs[c] + FNs[c];
    if(d > 1e-10) {
      score += TPs[c]/d;
    }
  }

  // return loss = 1 - score
  loss = 1.0 - (score/lClasses.size());
}

//------------------------------------------------------------------------------
