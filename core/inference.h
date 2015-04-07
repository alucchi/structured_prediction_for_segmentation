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

#ifndef INFERENCE_H
#define INFERENCE_H

#include <string>
#include <vector>

#include "energyParam.h"
#include "graphInference.h"
#include "svm_struct_api_types.h"

// SliceMe
#include "Feature.h"
#include "Slice_P.h"
#include "Slice.h"

//------------------------------------------------------------------------------

GraphInference* createGraphInferenceInstance(int algo_type,
                                             Slice_P* slice,
                                             const EnergyParam& param,
                                             Feature* _feature);

GraphInference* createGraphInferenceInstance(int algo_type,
                                             Slice_P* slice,
                                             const EnergyParam& param,
                                             Feature* feature,
                                             labelType* groundTruthLabels,
                                             double* lossPerLabel);

GraphInference* createGraphInferenceInstance(int algo_type,
                                             Slice_P* slice,
                                             const EnergyParam& param,
                                             Feature* feature,
                                             labelType* groundTruthLabels,
                                             double* lossPerLabel,
                                             map<sidType, nodeCoeffType>* _nodeCoeffs);

GraphInference* createGraphInferenceInstance(int algo_type,
                                             Slice_P* slice,
                                             const EnergyParam& param,
                                             Feature* feature,
                                             labelType* groundTruthLabels,
                                             double* lossPerLabel,
                                             map<sidType, nodeCoeffType>* _nodeCoeffs,
                                             map<sidType, edgeCoeffType>* _edgeCoeffs);

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
                    const int metric_type);

void computeScore(const char* image_dir,
                  const char* image_pattern,
                  const char* superpixelDir,
                  const char* maskDir,
                  int nImages,
                  const char* output_file,
                  map<ulong, labelType>& classIdxToLabel,
                  const char* weight_file);

void computeScore(const char* image_dir,
                  const char* image_pattern,
                  const char* superpixelDir,
                  const char* maskDir,
                  int nImages,
                  const char* output_file,
                  map<ulong, labelType>& classIdxToLabel,
                  const char* weight_file,
                  map<sidType, nodeCoeffType>* _nodeCoeffs,
                  map<sidType, edgeCoeffType>* _edgeCoeffs);

void computeScore(const IplImage* img,
                  const IplImage* imgAnnotation,
                  map<ulong, labelType>& classIdxToLabel,
                  bool convertToLabel,
                  int nClasses,
                  ulong* TPs,
                  ulong* FPs,
                  ulong* FNs,
                  ulong* count);

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
                  const char* outputDir);

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
                  map<sidType, edgeCoeffType>* _edgeCoeffs);

labelType* computeLabels(Slice_P* g, Feature* feature, const EnergyParam& param,
                         int algoType, double* energy,
                         labelType* groundTruthLabels, double* lossPerLabel,
                         map<sidType, nodeCoeffType>* _nodeCoeffs,
                         map<sidType, edgeCoeffType>* _edgeCoeffs);

labelType* computeLabels(Slice_P* g, Feature* feature, const EnergyParam& param,
                         int algoType, double* energy,
                         labelType* groundTruthLabels, double* lossPerLabel);

labelType* computeLabels(Slice_P* g, Feature* feature, const EnergyParam& param,
                         int algoType, double* energy);

labelType* computeLabels_sampling(Slice_P* g, Feature* feature, const EnergyParam& param,
                                  int algoType, double* energy,
                                  labelType* groundTruthLabels, double* lossPerLabel,
                                  map<sidType, nodeCoeffType>* _nodeCoeffs,
                                  map<sidType, edgeCoeffType>* _edgeCoeffs);

labelType* computeCombinedLabels(Slice_P* g, Feature* feature,
                                 labelType* groundTruthLabels, double* lossPerLabel,
                                 int nRFs, vector<double>* alphas, int example_id,
                                 map<sidType, nodeCoeffType>* _nodeCoeffs,
                                 map<sidType, edgeCoeffType>* _edgeCoeffs,
                                 const char* combined_probability_output_file);

void computeVOCLoss(labelType* nodeLabels, labelType* ybar,
                    int nNodes, int nClasses,
                    double& loss, int& nDiff,
                    ulong* TPs, ulong* FPs, ulong* FNs);

IplImage* getSegmentedImage(Slice* g,
                            Feature* feature,
                            int algoType,
                            const EnergyParam& param,
                            map<labelType, ulong>* labelToClassIdx,
                            map<sidType, nodeCoeffType>* _nodeCoeffs,
                            map<sidType, edgeCoeffType>* _edgeCoeffs);

IplImage* getSegmentedImage(Slice* g,
                            Feature* feature,
                            int algoType,
                            const EnergyParam& param,
                            map<labelType, ulong>* labelToClassIdx);

void outputScore(const char* output_file, map<ulong, labelType>& classIdxToLabel,
                 ulong* TPs, ulong* FPs, ulong* FNs, ulong* count,
                 const char* weight_file);

#endif //INFERENCE_H
