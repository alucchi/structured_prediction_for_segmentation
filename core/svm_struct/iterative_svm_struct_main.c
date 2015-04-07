/***********************************************************************/
/*                                                                     */
/*   svm_struct_main.c                                                 */
/*                                                                     */
/*   Command line interface to the alignment learning module of the    */
/*   Support Vector Machine.                                           */
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


/* the following enables you to use svm-learn out of C++ */
/*
#ifdef __cplusplus
extern "C" {
#endif
#include "../svm_light/svm_common.h"
#include "../svm_light/svm_learn.h"
#ifdef __cplusplus
}
#endif
*/

#include "../svm_light/svm_common.h"
#include "../svm_light/svm_learn.h"

#include "svm_struct_learn.h"
#include "svm_struct_common.h"
#include "../svm_struct_api.h"

#include "../../inference/energyParam.h"
#include "../../inference/graphInference.h"
#include "../../inference/gi_libDAI.h"
#include "../inference/inference.h"

// SliceMe
#include "Config.h"
#include "Feature.h"
#include "F_LoadFromFile.h"
#include "oSVM.h"
#include "Slice3d.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
/* } */

//------------------------------------------------------------------------MACROS

#undef SVM_FEAT_INDEX0
#undef SVM_FEAT_INDEX

// macros used to index w vector
// Caution : Those macros are duplicated in graphInference.h
#define SVM_FEAT_INDEX0 ((sparm->nGradientLevels*sparm->nOrientations*sparm->nClasses*sparm->nClasses) + sparm->nUnaryWeights + (sparm->nGradientLevels==0 && sparm->includeLocalEdges))
#define SVM_FEAT_INDEX(c,f) (SVM_FEAT_INDEX0 + (f*sparm->nUnaryWeights)+c)

#define SEPARATOR "\t"

//---------------------------------------------------------------------FUNCTIONS

#define BUFFER_SIZE 250
char trainfile[BUFFER_SIZE];           /* file with training examples */
char modelfile[BUFFER_SIZE];           /* file for resulting classifier */

void   read_input_parameters(int, char **, char *, char *,long *, long *,
			     STRUCT_LEARN_PARM *, LEARN_PARM *, KERNEL_PARM *,
			     int *);
void   wait_any_key();
void   print_help();

void clearFeatures(vector<osvm_node*>* anchorFeatures);

void loadAnchorPoints(const char* filename, vector<int>* anchorPoints);

void loadAnchorFeaturesFromSVMModel(vector<osvm_node*>* anchorFeatures,
                                    const oSVM& osvm);

void loadAnchorFeatures(vector<osvm_node*>* anchorFeatures,
                        const vector<int>& anchorPoints,
                        const int nNewAnchorPoints,
                        Slice_P* slice,
                        Feature* original_feature,
                        oSVM& original_osvm,
                        bool rescale_fv);

void reloadFeatures(vector<osvm_node*>* anchorFeatures,
                    Slice_P* slice, Feature* input_feature,
                    Feature* output_feature,
                    const osvm_parameter& param,
                    oSVM& original_osvm,
                    oSVM& kernelized_osvm,
                    bool rescale_fv);

void pruneAnchorPoints(vector<int>* anchorPoints,
                       const STRUCTMODEL& structmodel,
                       const STRUCT_LEARN_PARM* sparm,
                       const int nClasses,
                       const int nPointsToRemove,
                       const char* filename);

void computeFeatureWeights(vector< pair<double, int> >* weightedList,
                           vector<int>* anchorPoints,
                           const STRUCTMODEL& structmodel,
                           const STRUCT_LEARN_PARM* sparm,
                           const char* filename);

void selectAnchorPoints(vector<int>* anchorPoints,
                        const vector<int>& anchorPoints_candidates,
                        Slice_P* slice,
                        const double min_square_distance_anchor_points);

string loadFeatureFilename(const char* filename);

void exportList(const vector<int>& anchorPoints, const string& filename);

void exportFeatures(const vector<osvm_node*>& anchorFeatures, const string& filename);

string findLastFile(const string& file_pattern, const string& extension, int idx = 0);

int main (int argc, char* argv[])
{  
  SAMPLE sample;  /* training sample */
  LEARN_PARM learn_parm;
  KERNEL_PARM kernel_parm;
  STRUCT_LEARN_PARM struct_parm;
  STRUCTMODEL structmodel;
  int alg_type;

  svm_struct_learn_api_init(argc,argv);

  read_input_parameters(argc,argv,trainfile,modelfile,&verbosity,
			&struct_verbosity,&struct_parm,&learn_parm,
			&kernel_parm,&alg_type);

  if(struct_verbosity>=1) {
    //verbose = true;
    printf("Reading training examples..."); fflush(stdout);
  }
  /* read the training examples */
  sample=read_struct_examples(trainfile,&struct_parm);
  if(struct_verbosity>=1) {
    printf("done\n"); fflush(stdout);
  }
  
  const char* ssvm_iteration_dir = "ssvm_iteration/";
  mkdir(ssvm_iteration_dir, 0777);
  const char* anchor_points_dir = "anchor_points/";
  mkdir(anchor_points_dir, 0777);
  const char* anchor_features_dir = "anchor_features/";
  mkdir(anchor_features_dir, 0777);

  const char* parameter_vector_dir = "parameter_vector";

  string featureFile;
  if(!Config::Instance()->getParameter("feature_file", featureFile)) {
    printf("[SVM_struct] Error : path for feature file was not set\n");
    exit(-1);
  }
  PRINT_MESSAGE("[SVM_struct] Loading featureFile = %s\n", featureFile.c_str());
  //const char* kernelized_model_file = "/net/cvlabfiler1/home/lucchi/EM2/kernel/features/kernel_hist_rays3d-f7_original";
  //string kernelized_model_file = loadFeatureFilename(featureFile.c_str());

  string config_tmp;
  bool rescale_fv = true;
  if(Config::Instance()->getParameter("rescale_fv", config_tmp)) {
    if(config_tmp[0] == '0') {
      rescale_fv = false;
    }
  }
  PRINT_MESSAGE("[SVM_struct] rescale_fv=%d\n", (int)rescale_fv);

  string original_model_file;
  if(!Config::Instance()->getParameter("original_model_file", original_model_file)) {
    printf("[SVM_struct] Error : path for original_model_file was not set\n");
    exit(-1);
  }

  string kernelized_model_file;
  if(!Config::Instance()->getParameter("kernelized_model_file", kernelized_model_file)) {
    printf("[SVM_struct] Error : path for kernelized_model_file was not set\n");
    exit(-1);
  }

  int max_ssvm_iteration = 5;
  if(Config::Instance()->getParameter("max_ssvm_iteration", config_tmp)) {
    max_ssvm_iteration = atoi(config_tmp.c_str());
  }
  PRINT_MESSAGE("[SVM_struct] max_ssvm_iteration = %d\n", max_ssvm_iteration);

  uint max_anchor_points = 10000;
  if(Config::Instance()->getParameter("max_anchor_points", config_tmp)) {
    max_anchor_points = atoi(config_tmp.c_str());
  }
  PRINT_MESSAGE("[SVM_struct] max_anchor_points = %d\n", max_anchor_points);

  uint max_anchor_points_per_iteration = max_anchor_points/max_ssvm_iteration;
  if(Config::Instance()->getParameter("max_anchor_points_per_iteration", config_tmp)) {
    max_anchor_points_per_iteration = atoi(config_tmp.c_str());
  }
  PRINT_MESSAGE("[SVM_struct] max_anchor_points_per_iteration = %d\n", max_anchor_points_per_iteration);

  double min_square_distance_anchor_points = 500;
  if(Config::Instance()->getParameter("min_square_distance_anchor_points", config_tmp)) {
    min_square_distance_anchor_points = atof(config_tmp.c_str());
  }
  PRINT_MESSAGE("[SVM_struct] min_square_distance_anchor_points = %g\n", min_square_distance_anchor_points);

  // feature type
  vector<eFeatureType> feature_types;
  int paramFeatureTypes = 0;
  if(Config::Instance()->getParameter("originalfeatureTypes", config_tmp)) {
    paramFeatureTypes = atoi(config_tmp.c_str());
    int idxFeature = 0;
    int i = 0;
    while(1) {
      idxFeature = (int)pow(2.0,i);
      if(idxFeature == F_END_FEATURETYPE)
        break;
      
      if(paramFeatureTypes & idxFeature) {
        printf("[Main] Adding feature %d\n", idxFeature);
        feature_types.push_back((eFeatureType)idxFeature);
      }
      i++;
    }
  }  
  if(feature_types.size() == 0) {
    printf("[Main] Error : no features specified in the config file\n");
    exit(-1);
  }

  // Load anchor points and associated features
  vector< vector<int> > anchorPoints;
  vector<osvm_node*> anchorFeatures;
  const bool predict_probability = true;
  PRINT_MESSAGE("[SVM_struct] Loading original model from %s\n",
                original_model_file.c_str());

  oSVM original_osvm(F_UNKNOWN, original_model_file.c_str(), predict_probability);

  PRINT_MESSAGE("[SVM_struct] Loading kernelized model from %s\n",
                kernelized_model_file.c_str());
  // kernelized_osvm is used to rescaled features
  // TODO : a better way to rescale would be to load a .range file computed
  // on the whole dataset
  oSVM kernelized_osvm(F_UNKNOWN, kernelized_model_file.c_str(),
                       predict_probability);

#if 0
  string list_anchor_points;
  if(!Config::Instance()->getParameter("list_anchor_points", list_anchor_points)) {
    printf("[SVM_struct] Error : path for list_anchor_points was not set\n");
    exit(-1);
  }

  vector<int> example_anchorPoints;
  loadAnchorPoints(list_anchor_points.c_str(), &example_anchorPoints);
  anchorPoints.push_back(example_anchorPoints);

  PRINT_MESSAGE("[SVM_struct] Loading anchor features from %s\n",
                original_model_file.c_str());
  loadAnchorFeaturesFromSVMModel(&anchorFeatures, original_osvm);
  assert(example_anchorPoints.size() == anchorFeatures.size());

  struct_parm.ssvm_iteration = 0;
#else
  // Load features from file
  string anchor_points_dir_last;
  int idx_dir = 0;
  do {
    anchor_points_dir_last = anchor_points_dir;
    anchor_points_dir_last += StringPrintf("anchor_points_", idx_dir);
    if(isDirectory(anchor_points_dir_last)) {
      ++idx_dir;
    } else {
      break;
    }
  } while(1);
  --idx_dir;
  anchor_points_dir_last = anchor_points_dir;
  anchor_points_dir_last += StringPrintf("anchor_points_", idx_dir);
  anchor_points_dir_last += "/";
  struct_parm.ssvm_iteration = idx_dir + 1;

  vector<string> files;
  getFilesInDir(anchor_points_dir_last.c_str(), files, "txt", true);
  PRINT_MESSAGE("[SVM_struct] Loading anchor points from %ld files in %s\n",
                files.size(), anchor_points_dir_last.c_str());

  if(files.size() == 0) {
    PRINT_MESSAGE("[SVM_struct] Loading anchor features from %s\n",
                  original_model_file.c_str());

    string list_anchor_points;
    if(!Config::Instance()->getParameter("list_anchor_points", list_anchor_points)) {
      printf("[SVM_struct] Error : path for list_anchor_points was not set\n");
      exit(-1);
    }

    vector<int> example_anchorPoints;
    loadAnchorPoints(list_anchor_points.c_str(), &example_anchorPoints);
    anchorPoints.push_back(example_anchorPoints);

    loadAnchorFeaturesFromSVMModel(&anchorFeatures, original_osvm);
    assert(example_anchorPoints.size() == anchorFeatures.size());

    struct_parm.ssvm_iteration = 0;

  } else {
    for(uint example_id = 0; example_id < files.size(); ++example_id) {
      string list_anchor_points = anchor_points_dir_last + StringPrintf("/example_", example_id);
      list_anchor_points += ".txt";
      PRINT_MESSAGE("[SVM_struct] Loading anchor points from %s\n",
                    list_anchor_points.c_str());

      Slice_P* slice = sample.examples[example_id].x.slice;
      Feature* input_feature = Feature::getFeature(*slice, feature_types);

      // load anchor points
      vector<int> example_anchorPoints;
      loadAnchorPoints(list_anchor_points.c_str(), &example_anchorPoints);
      PRINT_MESSAGE("[SVM_struct] Loaded %ld anchor points from %s\n",
                    example_anchorPoints.size(), list_anchor_points.c_str());
      anchorPoints.push_back(example_anchorPoints);

      // load anchor features
      loadAnchorFeatures(&anchorFeatures, example_anchorPoints,
                         example_anchorPoints.size(), slice, input_feature,
                         original_osvm, rescale_fv);
    }

    // reload example features
    for(int i = 0; i < sample.n; ++i) {
      Slice_P* slice = sample.examples[i].x.slice;
      Feature* output_feature = sample.examples[i].x.feature;
      Feature* input_feature = Feature::getFeature(*slice, feature_types);

      reloadFeatures(&anchorFeatures, slice, input_feature, output_feature,
                     original_osvm.model->param, original_osvm, kernelized_osvm,
                     rescale_fv);
    }

    for(int i = 0; i < sample.nTest; ++i) {
      Slice_P* slice = sample.test_examples[i].x.slice;
      Feature* output_feature = sample.test_examples[i].x.feature;
      Feature* input_feature = Feature::getFeature(*slice, feature_types);

      reloadFeatures(&anchorFeatures, slice, input_feature, output_feature,
                     original_osvm.model->param, original_osvm, kernelized_osvm,
                     rescale_fv);
    }

    for(int i = 0; i < sample.nValidation; ++i) {
      Slice_P* slice = sample.validation_examples[i].x.slice;
      Feature* output_feature = sample.validation_examples[i].x.feature;
      Feature* input_feature = Feature::getFeature(*slice, feature_types);

      reloadFeatures(&anchorFeatures, slice, input_feature, output_feature,
                     original_osvm.model->param, original_osvm, kernelized_osvm,
                     rescale_fv);
    }

    // sanity check
    string anchorFeatures_filename = anchor_features_dir;
    anchorFeatures_filename += StringPrintf("anchor_features_", struct_parm.ssvm_iteration-1);
    mkdir(anchorFeatures_filename.c_str(), 0777);
    anchorFeatures_filename += StringPrintf("/reloaded_example", 0);
    printf("[Main] Exporting anchor features to %s\n", anchorFeatures_filename.c_str());
    exportFeatures(anchorFeatures, anchorFeatures_filename);
  }

  PRINT_MESSAGE("[SVM_struct] Starting at iteration %d\n", struct_parm.ssvm_iteration);
#endif
  PRINT_MESSAGE("[SVM_struct] Loaded %ld anchor features\n", anchorFeatures.size());

  // check if parameter vector files exist
  string parameter_vector_dir_last = findLastFile(parameter_vector_dir, "");
  string parameter_vector_file_pattern = parameter_vector_dir_last + "/iteration_";
  string parameter_vector_file_last = findLastFile(parameter_vector_file_pattern, ".txt", 1);
  printf("Checking parameter vector file %s\n", parameter_vector_file_last.c_str());
  bool skip_learning = false;
  if(fileExists(parameter_vector_file_last)) {
    skip_learning = true;
    EnergyParam e(parameter_vector_file_last.c_str());
    e.print();
    structmodel.svm_model = 0;
    structmodel.sizePsi = e.sizePsi;
    structmodel.w = new double[structmodel.sizePsi+1];
    structmodel.w[0] = 0;
    for(long l = 0; l < structmodel.sizePsi; ++l) {
      structmodel.w[l+1] = e.weights[l];
    }
  }

  ofstream ofs_stats("stats.txt");
  ofs_stats << "Iteration" << SEPARATOR << "Misclassified" << SEPARATOR
            << "Candidates" << SEPARATOR << "Added" << SEPARATOR
            << "Removed" << SEPARATOR << "AnchorPoints" << endl;
  ofs_stats << "0" << SEPARATOR << "-1" << SEPARATOR << "-1" << SEPARATOR
            << "-1" << SEPARATOR << anchorFeatures.size() << endl;

#define ALTERNATING_SCHEME 1

#if ALTERNATING_SCHEME
  bool alternating_do_pruning = true;
  bool reachedMaxNumberOfAnchorPoints = false;
#endif

  uint nNewAnchorPoints = 0;
  const int algo_type = T_GI_LIBDAI;
  const int max_iteration_BP = 100;
  PRINT_MESSAGE("[Main] Learning...iteration %d\n", struct_parm.ssvm_iteration);
  do {

    PRINT_MESSAGE("-----------------------------------------ITERATIVE_SSVM-%d-START\n",
                  struct_parm.ssvm_iteration);

    // add separation line between iterations
    ofstream ofs("min.txt", ios::app);
    ofs << "-----------------------------------------ITERATIVE_SSVM "
        << struct_parm.ssvm_iteration << endl;
    ofs.close();

    if(skip_learning) {
      skip_learning = false;
    } else {
      /* Do the learning and return structmodel. */
      if(alg_type == 0)
        svm_learn_struct(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,NSLACK_ALG);
      else if(alg_type == 1)
        svm_learn_struct(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,NSLACK_SHRINK_ALG);
      else if(alg_type == 2)
        svm_learn_struct_joint(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,ONESLACK_PRIMAL_ALG);
      else if(alg_type == 3)
        svm_learn_struct_joint(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,ONESLACK_DUAL_ALG);
      else if(alg_type == 4)
        svm_learn_struct_joint(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,ONESLACK_DUAL_CACHE_ALG);
      else if(alg_type == 9)
        svm_learn_struct_joint_custom(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel);
      else
        exit(1);
    }

    /* Warning: The model contains references to the original data 'docs'.
       If you want to free the original data, and only keep the model, you 
       have to make a deep copy of 'model'. */
    char _modelfile[BUFFER_SIZE];
    sprintf(_modelfile, "%s_%d", modelfile, struct_parm.ssvm_iteration);
    printf("[Main] Writing learned model to %s\n", _modelfile);
    write_struct_model(_modelfile, &structmodel, &struct_parm);

    //TODO : Loop over all the training examples : sample.n
    const int example_id = 0;

    //-----------------------------------------
    /* Run inference on training cube and add misclassified points to the training set */
    EnergyParam param(_modelfile);

    Slice_P* slice = sample.examples[example_id].x.slice;
    Feature* feature = sample.examples[example_id].x.feature;

    // Allocate memory
    int nNodes = slice->getNbSupernodes();
    labelType* inferredLabels = new labelType[nNodes];

    GraphInference* gi_Inference =
      createGraphInferenceInstance(algo_type, slice, param, feature);
    gi_Inference->run(inferredLabels,
                      0,
                      max_iteration_BP);
    delete gi_Inference;

    stringstream sout;
    sout << ssvm_iteration_dir;
    sout << slice->getName();
    sout << "_";
    sout << struct_parm.ssvm_iteration;
    slice->exportSupernodeLabels(sout.str().c_str(),
                                 param.nClasses,
                                 inferredLabels, nNodes, 0);

    //-----------------------------------------
    // Update anchor points and associated features
    if((struct_parm.ssvm_iteration + 1) < max_ssvm_iteration) {
      labelType* groundTruthLabels = sample.examples[example_id].y.nodeLabels;
      vector<int> misclassifiedNodes;
      for(int n = 0; n < nNodes; ++n) {
        if(groundTruthLabels[n] != inferredLabels[n]) {
          misclassifiedNodes.push_back(n);
        }
      }
      printf("[Main] %ld examples are misclassified at iteration %d\n",
             misclassifiedNodes.size(), struct_parm.ssvm_iteration);

      vector<int> anchorPoints_candidates;
      selectAnchorPoints(&anchorPoints_candidates, misclassifiedNodes,
                         slice, min_square_distance_anchor_points);

      printf("[Main] %ld examples are potential candidates\n",
             anchorPoints_candidates.size());

      Feature* input_feature = Feature::getFeature(*slice, feature_types);

      nNewAnchorPoints = anchorPoints_candidates.size();
      if(nNewAnchorPoints > max_anchor_points_per_iteration) {
        // add as many points as we can for current iteration
        nNewAnchorPoints = max_anchor_points_per_iteration;
      }
      if(anchorFeatures.size() > max_anchor_points) {
        //nNewAnchorPoints = 0;
#if ALTERNATING_SCHEME
        alternating_do_pruning = true;
#endif
      }
      printf("[Main] nNewAnchorPoints = %d/%ld\n", nNewAnchorPoints, anchorPoints_candidates.size());

      string featureWeights_filename = StringPrintf("featureWeights_", struct_parm.ssvm_iteration);
      int nPointsToRemove = 0;

      if(nNewAnchorPoints + anchorFeatures.size() > max_anchor_points) {

#if ALTERNATING_SCHEME
        reachedMaxNumberOfAnchorPoints = true;
#endif

        // remove anchor points with small weights
        if(anchorFeatures.size() > max_anchor_points) {
          nPointsToRemove = (anchorFeatures.size()+nNewAnchorPoints) - max_anchor_points;        
        } else {
          nPointsToRemove = min(max_anchor_points_per_iteration,
                                (uint)((anchorFeatures.size()+nNewAnchorPoints) - max_anchor_points));
        }

#if ALTERNATING_SCHEME
        if(reachedMaxNumberOfAnchorPoints && !alternating_do_pruning) {
          printf("[Main] alternating_do_pruning=false. Set number of points to remove to 0\n");
          nPointsToRemove = 0;
        }
#endif
        printf("[Main] Pruning %d anchor points\n", nPointsToRemove);
        if(nPointsToRemove > 0) {
          pruneAnchorPoints(&anchorPoints[example_id], structmodel, &struct_parm,
                            struct_parm.nClasses, nPointsToRemove,
                            featureWeights_filename.c_str());

          printf("[Main] %ld anchor points were kept after pruning\n",
                 anchorPoints[example_id].size());

          // reload anchor features
          clearFeatures(&anchorFeatures);
          loadAnchorFeatures(&anchorFeatures, anchorPoints[example_id],
                             anchorPoints[example_id].size(), slice,
                             input_feature, original_osvm, rescale_fv);
        }

      } else {
        printf("[Main] Computing weights %ld %ld\n",
               anchorPoints[example_id].size(), anchorFeatures.size());
        vector< pair<double, int> > weightedList;
        computeFeatureWeights(&weightedList, &anchorPoints[example_id], structmodel, &struct_parm,
                              featureWeights_filename.c_str());
      }

      printf("[Main] nNewAnchorPoints = %d/%ld anchorFeatures.size()=%ld max_anchor_points=%d\n",
             nNewAnchorPoints, anchorPoints_candidates.size(),
             anchorFeatures.size(), max_anchor_points);
      // add as many points as we can
      if(anchorFeatures.size() < max_anchor_points) {
          nNewAnchorPoints = max_anchor_points - anchorFeatures.size();
      } else {
        if(nNewAnchorPoints + anchorFeatures.size() > max_anchor_points) {
          nNewAnchorPoints = anchorFeatures.size() - max_anchor_points;
        }
      }
      if(nNewAnchorPoints > anchorPoints_candidates.size()) {
        nNewAnchorPoints = anchorPoints_candidates.size();
      }

#if ALTERNATING_SCHEME
      if(reachedMaxNumberOfAnchorPoints && alternating_do_pruning) {
        nNewAnchorPoints = 0;
      }
#endif

      printf("[Main] nNewAnchorPoints = %d/%ld\n", nNewAnchorPoints, anchorPoints_candidates.size());
      for(uint k = 0; k < nNewAnchorPoints; ++k) {
        anchorPoints[example_id].push_back(anchorPoints_candidates[k]);
      }

      string anchorPoints_filename = anchor_points_dir;
      anchorPoints_filename += StringPrintf("anchor_points_", struct_parm.ssvm_iteration);
      mkdir(anchorPoints_filename.c_str(), 0777);
      anchorPoints_filename += StringPrintf("/example_", example_id);
      anchorPoints_filename += ".txt";
      printf("[Main] Exporting anchor points to %s\n", anchorPoints_filename.c_str());
      exportList(anchorPoints[example_id], anchorPoints_filename);

      printf("[Main] Loading %d new anchor features\n", nNewAnchorPoints);
      loadAnchorFeatures(&anchorFeatures, anchorPoints_candidates,
                         nNewAnchorPoints, slice, input_feature, original_osvm,
                         rescale_fv);

      string anchorFeatures_filename = anchor_features_dir;
      anchorFeatures_filename += StringPrintf("anchor_features_", struct_parm.ssvm_iteration);
      mkdir(anchorFeatures_filename.c_str(), 0777);
      anchorFeatures_filename += StringPrintf("/example", example_id);
      printf("[Main] Exporting anchor features to %s\n", anchorFeatures_filename.c_str());
      exportFeatures(anchorFeatures, anchorFeatures_filename);

      // ---------------------------------
      // Reload features
      printf("[Main] Reloading features for training set\n");
      reloadFeatures(&anchorFeatures, slice, input_feature, feature,
                     original_osvm.model->param, original_osvm, kernelized_osvm,
                     rescale_fv);

      printf("[Main] Reloading features for test and validation sets\n");
      for(int i = 0; i < sample.nTest; ++i) {
        Slice_P* test_slice = sample.test_examples[i].x.slice;
        Feature* test_feature = sample.test_examples[i].x.feature;
        Feature* test_input_feature = Feature::getFeature(*test_slice, feature_types);

        reloadFeatures(&anchorFeatures, test_slice, test_input_feature, test_feature,
                       original_osvm.model->param, original_osvm, kernelized_osvm,
                       rescale_fv);
      }

      for(int i = 0; i < sample.nValidation; ++i) {
        Slice_P* validation_slice = sample.validation_examples[i].x.slice;
        Feature* validation_feature = sample.validation_examples[i].x.feature;
        Feature* validation_input_feature = Feature::getFeature(*validation_slice, feature_types);

        reloadFeatures(&anchorFeatures, validation_slice, validation_input_feature, validation_feature,
                       original_osvm.model->param, original_osvm, kernelized_osvm,
                       rescale_fv);
      }

      // ---------------------------------
      // Log results
      ofs_stats << struct_parm.ssvm_iteration + 1 << SEPARATOR;
      ofs_stats << misclassifiedNodes.size() << SEPARATOR;
      ofs_stats << anchorPoints_candidates.size() << SEPARATOR;
      ofs_stats << nNewAnchorPoints << SEPARATOR;
      ofs_stats << nPointsToRemove << SEPARATOR;
      ofs_stats << anchorPoints[example_id].size() << endl;
    }

#if ALTERNATING_SCHEME
    alternating_do_pruning = !alternating_do_pruning;
#endif

    delete[] inferredLabels;
    ++struct_parm.ssvm_iteration;

#if ALTERNATING_SCHEME
  } while ((struct_parm.ssvm_iteration < max_ssvm_iteration));
#else
  } while ((struct_parm.ssvm_iteration < max_ssvm_iteration) && (nNewAnchorPoints>0));
#endif

  ofs_stats.close();

  printf("[Main] Cleaning\n");
  clearFeatures(&anchorFeatures);
  free_struct_sample(sample);
  free_struct_model(structmodel);

  svm_struct_learn_api_exit();

  return 0;
}

string loadFeatureFilename(const char* filename)
{
  ifstream ifs(filename);
  if(ifs.fail()) {
    printf("[Main] Error while loading %s\n", filename);
    exit(-1);
  }
  
  // Load path containing files
  string featurePath;
  getline(ifs, featurePath);

  string trainingDir;
  Config::Instance()->getParameter("trainingDir", trainingDir);
  string cubeName = getLastDirectoryFromPath(trainingDir) + "/";

  // Load name of the files containing the features
  string fn;
  getline(ifs,fn);
  ifs.close();

  return featurePath + cubeName + fn;
}

void loadAnchorPoints(const char* filename, vector<int>* anchorPoints)
{
  ifstream ifs(filename);
  if(ifs.fail()) {
    printf("[Main] Error while loading %s\n", filename);
    exit(-1);
  }
  string line;
  while(getline(ifs, line)) {
    //printf("%s\n", line.c_str());
    anchorPoints->push_back(atoi(line.c_str()));
  }
  ifs.close();
}

void loadAnchorFeaturesFromSVMModel(vector<osvm_node*>* anchorFeatures,
                                    const oSVM& osvm)
{
  osvm_model* model = osvm.model;
  PRINT_MESSAGE("[Main] Loaded kernelized model containing %d vectors\n",
                model->l);

  // Copy sparse vectors
  anchorFeatures->resize(model->l);
  for(int n = 0; n < model->l; ++n) {
    int feature_size = 0;
    for(; model->SV[n][feature_size].index != -1; ++feature_size);
    oSVM::initSVMNode((*anchorFeatures)[n], feature_size);
    for(int i = 0; model->SV[n][i].index != -1; ++i) {
      (*anchorFeatures)[n][i].index = model->SV[n][i].index;
      (*anchorFeatures)[n][i].value = model->SV[n][i].value;
    }
  }

  //printf("CHECKING ANCHOR FEATURES\n");
  //oSVM::print((*anchorFeatures)[0],"SV");
}

/**
 * @param nNewAnchorPoints is the maximum number of anchor points that can be accepted
 */
void loadAnchorFeatures(vector<osvm_node*>* anchorFeatures,
                        const vector<int>& anchorPoints,
                        const int nNewAnchorPoints,
                        Slice_P* slice,
                        Feature* original_feature,
                        oSVM& original_osvm,
                        bool rescale_fv)
{
  // Copy sparse vectors
  const int feature_size = original_feature->getSizeFeatureVector();
  PRINT_MESSAGE("[Main] Considering %ld candidate anchor features of size %d\n",
                anchorPoints.size(), feature_size);
  int nAcceptedPoints = 0;
  for(vector<int>::const_iterator it = anchorPoints.begin();
      it != anchorPoints.end(); ++it) {
    osvm_node* x = 0;
    oSVM::initSVMNode(x, feature_size);
    original_feature->getFeatureVector(x, *slice, *it);

    if(rescale_fv) {
      original_osvm.rescale(x);
    }

    anchorFeatures->push_back(x);
    ++nAcceptedPoints;
    if(nAcceptedPoints >= nNewAnchorPoints) {
      break;
    }
  }

  PRINT_MESSAGE("[Main] Added %d anchor features. Total number of anchor features = %ld\n",
                nAcceptedPoints, anchorFeatures->size());
}


void clearFeatures(vector<osvm_node*>* anchorFeatures)
{
  for(vector<osvm_node*>::iterator it = anchorFeatures->begin();
      it != anchorFeatures->end(); ++it) {
    delete[] *it;
  }
  anchorFeatures->clear();
}

void reloadFeatures(vector<osvm_node*>* anchorFeatures,
                    Slice_P* slice, Feature* input_feature,
                    Feature* output_feature,
                    const osvm_parameter& param,
                    oSVM& original_osvm,
                    oSVM& kernelized_osvm,
                    bool rescale_fv)
{
  // Initialize feature vectors
  osvm_node* x = 0;
  const int input_feature_size = input_feature->getSizeFeatureVector(); 
  oSVM::initSVMNode(x, input_feature_size);
  osvm_node* x_k = 0;
  const int output_feature_size = anchorFeatures->size(); 
  oSVM::initSVMNode(x_k, output_feature_size);

  eFeatureType featureType = output_feature->getFeatureType();
  if(featureType != F_LOADFROMFILE) {
    printf("[Main] featureType != F_LOADFROMFILE\n");
    exit(-1);
  }
  F_LoadFromFile* fLoadFromFile = static_cast<F_LoadFromFile*>(output_feature);
  fLoadFromFile->clearFeatures();
  map<ulong, fileFeatureType*>* output_feature_map = fLoadFromFile->getMutableFeatures();
  fLoadFromFile->setFeatureSize(output_feature_size);

  sidType sid;
  const map<sidType, supernode* >& _supernodes = slice->getSupernodes();
  PRINT_MESSAGE("[Main] Creating new features for %ld supernodes\n",
                _supernodes.size());
  PRINT_MESSAGE("[Main] Mapping %d->%d\n", input_feature_size, output_feature_size);
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    sid = it->first;
    // Create feature vector
    input_feature->getFeatureVector(x, *slice, sid);

    // DEBUG
    if(sid == 0) {
      //printf("svm_param %g %g\n", svm_param.C, svm_param.gamma);
      printf("svm_param %g\n", param.gamma);
      oSVM::print((*anchorFeatures)[anchorFeatures->size()-1],
                  StringPrintf("SV_", anchorFeatures->size()-1).c_str());
      oSVM::print(x,"originalFeatureVector");
    }
    // DEBUG

    // TODO : should compute rescale parameters from the whole datasets
    original_osvm.rescale(x);

    // Kernalize vector
    // Dimension of the output vector is the number of anchor points
    for(int j=0; j < output_feature_size; j++) {
      x_k[j].value =
        Kernel::k_function(x, (*anchorFeatures)[j], param);
    }

    // DEBUG
    if(sid == 0) {
      oSVM::print(x,"originalFeatureVector_rescaled");
      oSVM::printNonZeros(x_k,"FeatureVector");
      //oSVM::print(x_k,"FeatureVector");
    }

    if(rescale_fv) {
      kernelized_osvm.rescale(x_k);
    }

    if(sid == 0) {
      oSVM::printNonZeros(x_k,"rescaledFeatureVector");
      //oSVM::print(x_k,"rescaledFeatureVector");
    }

    // copy output vector
    (*output_feature_map)[sid] = new fileFeatureType[output_feature_size];
    for(int i = 0; i < output_feature_size; ++i) {
      (*output_feature_map)[sid][i] = 0;
    }
    for(int i = 0; i < output_feature_size; ++i) {
      (*output_feature_map)[sid][i] = x_k[i].value;
    }
    /*
    for(int i = 0; x_k[i].index != -1; ++i) {
      (*output_feature_map)[sid][x_k[i].index-1] = x_k[i].value;
    }
    */

    /*
    // sparse vector : would need to change F_LoadFromFile to store osvm_node
    int feature_size = 0;
    for(; x_k[feature_size].index != -1; ++feature_size);
    oSVM::initSVMNode(output_feature[sid], feature_size);
    for(int i = 0; x_k[i].index != -1; ++i) {
      output_feature[sid].index = x_k[i].index;
      output_feature[sid].value = x_k[i].value;
    }
    */
  }

  delete[] x;
  delete[] x_k;
}

void exportList(const vector<int>& anchorPoints, const string& filename)
{
  ofstream ofs(filename.c_str());
  for(vector<int>::const_iterator it = anchorPoints.begin();
      it != anchorPoints.end(); ++it) {
    ofs << *it << endl;
  }
  ofs.close();
}

void exportFeatures(const vector<osvm_node*>& anchorFeatures, const string& filename)
{
  ofstream ofs(filename.c_str());
  for(vector<osvm_node*>::const_iterator it = anchorFeatures.begin();
      it != anchorFeatures.end(); ++it) {
    osvm_node* n = *it;
    for(int i = 0; n[i].index != -1; ++i) {
      ofs << n[i].index << ":" << n[i].value << " ";
    }
    ofs << endl;
  }
  ofs.close();
}

void selectAnchorPoints(vector<int>* anchorPoints_out,
                        const vector<int>& anchorPoints_in,
                        Slice_P* slice,
                        const double min_square_distance_anchor_points)
{
  PRINT_MESSAGE("[Main] Considering %ld candidate anchor points\n",
                anchorPoints_in.size());
  double total_avg_dist = 0;
  int nAcceptedPoints = 0;
  for(vector<int>::const_iterator it = anchorPoints_in.begin();
      it != anchorPoints_in.end(); ++it) {

    bool addFeature = true;
    double avg_dist = 0;
    int np = 0;

    // Compute euclidean distance to all anchor points
    supernode* sc = slice->getSupernode(*it); // candidate
    node center_c;
    sc->getCenter(center_c);
    node center_p;
    float dist = 0;
    for(vector<int>::iterator itPoint = anchorPoints_out->begin();
        itPoint != anchorPoints_out->end(); ++itPoint) {
      supernode* sp = slice->getSupernode(*itPoint);
      sp->getCenter(center_p);
      dist = node_square_distance(center_c, center_p);
      if(dist < min_square_distance_anchor_points) {
        addFeature = false;
        break;
      }
      avg_dist += dist;
      ++np;
    }
    if(np !=0) {
      avg_dist /= np;
    }

    total_avg_dist += avg_dist;

    if(addFeature) {
      anchorPoints_out->push_back(*it);
      ++nAcceptedPoints;
    }
  }
  if(anchorPoints_in.size() != 0) {
    total_avg_dist /= anchorPoints_in.size();
  }

  PRINT_MESSAGE("[Main] Average distance = %g. %d anchor points accepted\n",
                total_avg_dist, nAcceptedPoints);
}

void pruneAnchorPoints(vector<int>* anchorPoints,
                       const STRUCTMODEL& structmodel,
                       const STRUCT_LEARN_PARM* sparm,
                       const int nClasses,
                       const int nPointsToRemove,
                       const char* filename)
{
  vector< pair<double, int> > weightedList;
  computeFeatureWeights(&weightedList, anchorPoints, structmodel, sparm,
                        filename);
  // copy new list
  anchorPoints->clear();
  for(uint i = nPointsToRemove; i < weightedList.size(); ++i) {
    anchorPoints->push_back(weightedList[i].second);
  }
}

void computeFeatureWeights(vector< pair<double, int> >* weightedList,
                           vector<int>* anchorPoints,
                           const STRUCTMODEL& structmodel,
                           const STRUCT_LEARN_PARM* sparm,
                           const char* filename)
{
  if(structmodel.sizePsi < (int)anchorPoints->size()) {
    printf("[Main] Error : structmodel.sizePsi < anchorPoints->size() %ld %ld\n",
           structmodel.sizePsi, anchorPoints->size());
    exit(-1);
  }
  // Loop over all anchor points and compute their weight over classes.
  int idx = 0;
  double t = 0;
  double w_norm = 0;
  for(vector<int>::iterator it = anchorPoints->begin();
      it != anchorPoints->end(); ++it) {
    double weight = 0;
    for(int c =0; c < sparm->nUnaryWeights; ++c) {
      //printf("%d %d %d\n", idx,c,SVM_FEAT_INDEX(c, idx));
      t = structmodel.w[SVM_FEAT_INDEX(c, idx)];
      weight += t*t;
    }
    w_norm += weight;
    weightedList->push_back(make_pair(weight, idx));
    ++idx;
  }
  w_norm = sqrt(w_norm);
  printf("[Main] pruneAnchorPoints : w_norm = %g\n", w_norm);

  std::sort(weightedList->begin(), weightedList->end());
  printf("[Main] pruneAnchorPoints : minimum weight = %g, maximum weight = %g\n",
         (*weightedList)[0].first, (*weightedList)[weightedList->size()-1].first);

  if(filename) {
    ofstream ofs(filename);
    ofs << structmodel.sizePsi << " " << anchorPoints->size() << endl << endl;
    for(uint i = 0; i < weightedList->size(); ++i) {
      ofs << (*weightedList)[i].second << " " << (*weightedList)[i].first << endl;
    }  
    ofs.close();
  }
}

string findLastFile(const string& file_pattern, const string& extension, int idx)
{
  string last_file;
  do {
    last_file = file_pattern;
    last_file += VarToString(idx);
    last_file += extension;
    if(fileExists(last_file)) {
    //if(isDirectory(last_file)) {
      ++idx;
    } else {
      break;
    }
  } while(1);
  --idx;
  last_file = file_pattern;
  last_file += VarToString(idx);
  last_file += extension;
  return last_file; 
}

/*---------------------------------------------------------------------------*/

void read_input_parameters(int argc,char *argv[],char *trainfile,
			   char *modelfile,
			   long *verbosity,long *struct_verbosity, 
			   STRUCT_LEARN_PARM *struct_parm,
			   LEARN_PARM *learn_parm, KERNEL_PARM *kernel_parm,
			   int *alg_type)
{
  long i;
  char type[100];
  
  /* set default */
  (*alg_type)=DEFAULT_ALG_TYPE;
  struct_parm->C=-0.01;
  struct_parm->maxC=1.0;
  struct_parm->giType=0;
  struct_parm->slack_norm=1;
  struct_parm->epsilon=DEFAULT_EPS;
  struct_parm->custom_argc=0;
  struct_parm->loss_function=DEFAULT_LOSS_FCT;
  struct_parm->loss_type=DEFAULT_RESCALING;
  struct_parm->newconstretrain=100;
  struct_parm->ccache_size=5;
  struct_parm->batch_size=100;

  strcpy (modelfile, "svm_struct_model");
  strcpy (learn_parm->predfile, "trans_predictions");
  strcpy (learn_parm->alphafile, "");
  (*verbosity)=0;/*verbosity for svm_light*/
  (*struct_verbosity)=1; /*verbosity for struct learning portion*/
  learn_parm->biased_hyperplane=1;
  learn_parm->remove_inconsistent=0;
  learn_parm->skip_final_opt_check=0;
  learn_parm->svm_maxqpsize=10;
  learn_parm->svm_newvarsinqp=0;
  learn_parm->svm_iter_to_shrink=-9999;
  learn_parm->maxiter=100000;
  learn_parm->kernel_cache_size=40;
  learn_parm->svm_c=99999999;  /* overridden by struct_parm->C */
  learn_parm->eps=0.001;       /* overridden by struct_parm->epsilon */
  learn_parm->transduction_posratio=-1.0;
  learn_parm->svm_costratio=1.0;
  learn_parm->svm_costratio_unlab=1.0;
  learn_parm->svm_unlabbound=1E-5;
  learn_parm->epsilon_crit=0.001;
  learn_parm->epsilon_a=1E-10;  /* changed from 1e-15 */
  learn_parm->compute_loo=0;
  learn_parm->rho=1.0;
  learn_parm->xa_depth=0;
  kernel_parm->kernel_type=0;
  kernel_parm->poly_degree=3;
  kernel_parm->rbf_gamma=1.0;
  kernel_parm->coef_lin=1;
  kernel_parm->coef_const=1;
  strcpy(kernel_parm->custom,"empty");
  strcpy(type,"c");

  for(i=1;(i<argc) && ((argv[i])[0] == '-');i++) {
    switch ((argv[i])[1]) 
      { 
      case '?': print_help(); exit(0);
      case 'a': i++; strcpy(learn_parm->alphafile,argv[i]); break;
      case 'c': i++; struct_parm->C=atof(argv[i]); break;
      case 'x': i++; struct_parm->maxC=atof(argv[i]); break; //al
      case 'i': i++; struct_parm->giType=atoi(argv[i]); break; //al
      case 'p': i++; struct_parm->slack_norm=atol(argv[i]); break;
      case 'e': i++; struct_parm->epsilon=atof(argv[i]); break;
      case 'k': i++; struct_parm->newconstretrain=atol(argv[i]); break;
      case 'h': i++; learn_parm->svm_iter_to_shrink=atol(argv[i]); break;
      case '#': i++; learn_parm->maxiter=atol(argv[i]); break;
      case 'm': i++; learn_parm->kernel_cache_size=atol(argv[i]); break;
      case 'w': i++; (*alg_type)=atol(argv[i]); break;
      case 'o': i++; struct_parm->loss_type=atol(argv[i]); break;
      case 'n': i++; learn_parm->svm_newvarsinqp=atol(argv[i]); break;
      case 'q': i++; learn_parm->svm_maxqpsize=atol(argv[i]); break;
      case 'l': i++; struct_parm->loss_function=atol(argv[i]); break;
      case 'f': i++; struct_parm->ccache_size=atol(argv[i]); break;
      case 'b': i++; struct_parm->batch_size=atof(argv[i]); break;
      case 't': i++; kernel_parm->kernel_type=atol(argv[i]); break;
      case 'd': i++; kernel_parm->poly_degree=atol(argv[i]); break;
      case 'g': i++; kernel_parm->rbf_gamma=atof(argv[i]); break;
      case 's': i++; kernel_parm->coef_lin=atof(argv[i]); break;
      case 'r': i++; kernel_parm->coef_const=atof(argv[i]); break;
      case 'u': i++; strcpy(kernel_parm->custom,argv[i]); break;
      case '-': strcpy(struct_parm->custom_argv[struct_parm->custom_argc++],argv[i]);i++; strcpy(struct_parm->custom_argv[struct_parm->custom_argc++],argv[i]);break; 
      case 'v': i++; (*struct_verbosity)=atol(argv[i]); break;
      case 'y': i++; (*verbosity)=atol(argv[i]); break;
      default: printf("\nUnrecognized option %s!\n\n",argv[i]);
	       print_help();
	       exit(0);
      }
  }
  if(i>=argc) {
    printf("\nNot enough input parameters!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  strcpy (trainfile, argv[i]);
  if((i+1)<argc) {
    strcpy (modelfile, argv[i+1]);
  }
  if(learn_parm->svm_iter_to_shrink == -9999) {
    learn_parm->svm_iter_to_shrink=100;
  }

  if((learn_parm->skip_final_opt_check) 
     && (kernel_parm->kernel_type == SVM_LINEAR)) {
    printf("\nIt does not make sense to skip the final optimality check for linear kernels.\n\n");
    learn_parm->skip_final_opt_check=0;
  }    
  if((learn_parm->skip_final_opt_check) 
     && (learn_parm->remove_inconsistent)) {
    printf("\nIt is necessary to do the final optimality check when removing inconsistent \nexamples.\n");
    wait_any_key();
    print_help();
    exit(0);
  }    
  if((learn_parm->svm_maxqpsize<2)) {
    printf("\nMaximum size of QP-subproblems not in valid range: %ld [2..]\n",learn_parm->svm_maxqpsize); 
    wait_any_key();
    print_help();
    exit(0);
  }
  if((learn_parm->svm_maxqpsize<learn_parm->svm_newvarsinqp)) {
    printf("\nMaximum size of QP-subproblems [%ld] must be larger than the number of\n",learn_parm->svm_maxqpsize); 
    printf("new variables [%ld] entering the working set in each iteration.\n",learn_parm->svm_newvarsinqp); 
    wait_any_key();
    print_help();
    exit(0);
  }
  if(learn_parm->svm_iter_to_shrink<1) {
    printf("\nMaximum number of iterations for shrinking not in valid range: %ld [1,..]\n",learn_parm->svm_iter_to_shrink);
    wait_any_key();
    print_help();
    exit(0);
  }
  if(struct_parm->C<0) {
    printf("\nYou have to specify a value for the parameter '-c' (C>0)!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if(((*alg_type) < 0) || (((*alg_type) > 5) && ((*alg_type) != 9))) {
    printf("\nAlgorithm type must be either '0', '1', '2', '3', '4', or '9'!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if(learn_parm->transduction_posratio>1) {
    printf("\nThe fraction of unlabeled examples to classify as positives must\n");
    printf("be less than 1.0 !!!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if(learn_parm->svm_costratio<=0) {
    printf("\nThe COSTRATIO parameter must be greater than zero!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if(struct_parm->epsilon<=0) {
    printf("\nThe epsilon parameter must be greater than zero!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if((struct_parm->ccache_size<=0) && ((*alg_type) == 4)) {
    printf("\nThe cache size must be at least 1!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if(((struct_parm->batch_size<=0) || (struct_parm->batch_size>100))  
     && ((*alg_type) == 4)) {
    printf("\nThe batch size must be in the interval ]0,100]!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if((struct_parm->slack_norm<1) || (struct_parm->slack_norm>2)) {
    printf("\nThe norm of the slacks must be either 1 (L1-norm) or 2 (L2-norm)!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if((struct_parm->loss_type != SLACK_RESCALING) 
     && (struct_parm->loss_type != MARGIN_RESCALING)) {
    printf("\nThe loss type must be either 1 (slack rescaling) or 2 (margin rescaling)!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if(learn_parm->rho<0) {
    printf("\nThe parameter rho for xi/alpha-estimates and leave-one-out pruning must\n");
    printf("be greater than zero (typically 1.0 or 2.0, see T. Joachims, Estimating the\n");
    printf("Generalization Performance of an SVM Efficiently, ICML, 2000.)!\n\n");
    wait_any_key();
    print_help();
    exit(0);
  }
  if((learn_parm->xa_depth<0) || (learn_parm->xa_depth>100)) {
    printf("\nThe parameter depth for ext. xi/alpha-estimates must be in [0..100] (zero\n");
    printf("for switching to the conventional xa/estimates described in T. Joachims,\n");
    printf("Estimating the Generalization Performance of an SVM Efficiently, ICML, 2000.)\n");
    wait_any_key();
    print_help();
    exit(0);
  }

  parse_struct_parameters(struct_parm);
}

void wait_any_key()
{
  printf("\n(more)\n");
  (void)getc(stdin);
}

void print_help()
{
  printf("\nSVM-struct learning module: %s, %s, %s\n",INST_NAME,INST_VERSION,INST_VERSION_DATE);
  printf("   includes SVM-struct %s for learning complex outputs, %s\n",STRUCT_VERSION,STRUCT_VERSION_DATE);
  printf("   includes SVM-light %s quadratic optimizer, %s\n",VERSION,VERSION_DATE);
  copyright_notice();
  printf("   usage: svm_struct_learn [options] example_file model_file\n\n");
  printf("Arguments:\n");
  printf("         example_file-> file with training data\n");
  printf("         model_file  -> file to store learned decision rule in\n");

  printf("General Options:\n");
  printf("         -?          -> this help\n");
  printf("         -v [0..3]   -> verbosity level (default 1)\n");
  printf("         -y [0..3]   -> verbosity level for svm_light (default 0)\n");
  printf("Learning Options:\n");
  printf("         -c float    -> C: trade-off between training error\n");
  printf("                        and margin (default 0.01)\n");
  printf("         -p [1,2]    -> L-norm to use for slack variables. Use 1 for L1-norm,\n");
  printf("                        use 2 for squared slacks. (default 1)\n");
  printf("         -o [1,2]    -> Rescaling method to use for loss.\n");
  printf("                        1: slack rescaling\n");
  printf("                        2: margin rescaling\n");
  printf("                        (default %d)\n",DEFAULT_RESCALING);
  printf("         -l [0..]    -> Loss function to use.\n");
  printf("                        0: zero/one loss\n");
  printf("                        ?: see below in application specific options\n");
  printf("                        (default %d)\n",DEFAULT_LOSS_FCT);
  printf("Optimization Options (see [2][5]):\n");
  printf("         -w [0,..,9] -> choice of structural learning algorithm (default %d):\n",(int)DEFAULT_ALG_TYPE);
  printf("                        0: n-slack algorithm described in [2]\n");
  printf("                        1: n-slack algorithm with shrinking heuristic\n");
  printf("                        2: 1-slack algorithm (primal) described in [5]\n");
  printf("                        3: 1-slack algorithm (dual) described in [5]\n");
  printf("                        4: 1-slack algorithm (dual) with constraint cache [5]\n");
  printf("                        9: custom algorithm in svm_struct_learn_custom.c\n");
  printf("         -e float    -> epsilon: allow that tolerance for termination\n");
  printf("                        criterion (default %f)\n",DEFAULT_EPS);
  printf("         -k [1..]    -> number of new constraints to accumulate before\n"); 
  printf("                        recomputing the QP solution (default 100) (-w 0 and 1 only)\n");
  printf("         -f [5..]    -> number of constraints to cache for each example\n");
  printf("                        (default 5) (used with -w 4)\n");
  printf("         -b [1..100] -> percentage of training set for which to refresh cache\n");
  printf("                        when no epsilon violated constraint can be constructed\n");
  printf("                        from current cache (default 100%%) (used with -w 4)\n");
  printf("SVM-light Options for Solving QP Subproblems (see [3]):\n");
  printf("         -n [2..q]   -> number of new variables entering the working set\n");
  printf("                        in each svm-light iteration (default n = q). \n");
  printf("                        Set n < q to prevent zig-zagging.\n");
  printf("         -m [5..]    -> size of svm-light cache for kernel evaluations in MB\n");
  printf("                        (default 40) (used only for -w 1 with kernels)\n");
  printf("         -h [5..]    -> number of svm-light iterations a variable needs to be\n"); 
  printf("                        optimal before considered for shrinking (default 100)\n");
  printf("         -# int      -> terminate svm-light QP subproblem optimization, if no\n");
  printf("                        progress after this number of iterations.\n");
  printf("                        (default 100000)\n");
  printf("Kernel Options:\n");
  printf("         -t int      -> type of kernel function:\n");
  printf("                        0: linear (default)\n");
  printf("                        1: polynomial (s a*b+c)^d\n");
  printf("                        2: radial basis function exp(-gamma ||a-b||^2)\n");
  printf("                        3: sigmoid tanh(s a*b + c)\n");
  printf("                        4: user defined kernel from kernel.h\n");
  printf("         -d int      -> parameter d in polynomial kernel\n");
  printf("         -g float    -> parameter gamma in rbf kernel\n");
  printf("         -s float    -> parameter s in sigmoid/poly kernel\n");
  printf("         -r float    -> parameter c in sigmoid/poly kernel\n");
  printf("         -u string   -> parameter of user defined kernel\n");
  printf("Output Options:\n");
  printf("         -a string   -> write all alphas to this file after learning\n");
  printf("                        (in the same order as in the training set)\n");
  printf("Application-Specific Options:\n");
  print_struct_help();
  wait_any_key();

  printf("\nMore details in:\n");
  printf("[1] T. Joachims, Learning to Align Sequences: A Maximum Margin Aproach.\n");
  printf("    Technical Report, September, 2003.\n");
  printf("[2] I. Tsochantaridis, T. Joachims, T. Hofmann, and Y. Altun, Large Margin\n");
  printf("    Methods for Structured and Interdependent Output Variables, Journal\n");
  printf("    of Machine Learning Research (JMLR), Vol. 6(Sep):1453-1484, 2005.\n");
  printf("[3] T. Joachims, Making Large-Scale SVM Learning Practical. Advances in\n");
  printf("    Kernel Methods - Support Vector Learning, B. Schölkopf and C. Burges and\n");
  printf("    A. Smola (ed.), MIT Press, 1999.\n");
  printf("[4] T. Joachims, Learning to Classify Text Using Support Vector\n");
  printf("    Machines: Methods, Theory, and Algorithms. Dissertation, Kluwer,\n");
  printf("    2002.\n");
  printf("[5] T. Joachims, T. Finley, Chun-Nam Yu, Cutting-Plane Training of Structural\n");
  printf("    SVMs, Machine Learning Journal, to appear.\n");
}



