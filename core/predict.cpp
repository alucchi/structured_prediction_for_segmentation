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

#include <cv.h>
#include <highgui.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
//#include <argp.h>
#include <vector>
#include <sys/stat.h>
#include <sys/types.h>

//argp replacement
#ifdef _WIN32
#include "getopt.h"
#else
#include "unistd.h"
#include <getopt.h>
#endif

// SliceMe
#include "Config.h"
#include "Slice.h"
#include "colormap.h"
#include "utils.h"
#include "globalsE.h"
#include "globals.h"
#include "inference.h"
#include "Feature.h"
#include "F_Combo.h"

#include "gi_libDAI.h"

#include "energyParam.h"
#include "svm_struct_api_types.h"
#include "svm_struct_api.h"

#ifdef _WIN32
#include "direct.h"
#define mkdir(x,y) _mkdir(x)
#endif

using namespace std;

/* Variables for the arguments.*/
const char *argp_program_version =
  "predict 0.1";
const char *argp_program_bug_address =
  "<aurelien.lucchi@epfl.ch>";
/* Program documentation. */
static char doc[] =
  "Inference";

/* A description of the arguments we accept. */
static char args_doc[] = "args";

/* Program options */
static struct option long_options[] = {
  {"all", no_argument, 0, 'a'}, //"export marginals and also run inference using unary potentials only (useful for debugging)"},
  {"config_file", required_argument, 0, 'c'}, //"config_file"},
  {"image_dir", required_argument, 0, 'i'}, //"input directory"},
  {"algo_type", required_argument, 0, 'g'}, //"algo_type"},
  {"image_pattern", required_argument, 0, 'k'}, //"image_pattern"},
  {"superpixel_labels", required_argument, 0, 'l'}, //"path of the file containing the labels for the superpixels (labels have be ordered by rows)"},
  {"mask_dir", required_argument, 0, 'm'}, //"mask directory"},
  {"nImages", required_argument, 0, 'n'}, //"number of images to process"},
  {"output_dir", required_argument, 0, 'o'}, //"output filename"},
  {"superpixelStepSize", required_argument, 0, 's'}, //"superpixel step size"},
  {"dataset_type", required_argument, 0, 't'}, //"type (0=training, 1=test)"},
  {"verbose", no_argument, 0, 'v'}, //"verbose"},
  {"weight_file", required_argument, 0, 'w'}, //"weight_file"},
  {"overlay", required_argument, 0, 'y'}, //"overlay directory"},
  {"help", no_argument, 0, 'h'}, //"usage"},
  { 0, 0, 0, 0}
};


struct arguments
{
  bool export_all;
  char* image_dir;
  char* superpixel_labels;
  char* output_dir;
  char* weight_file;
  char* image_pattern;
  char* mask_dir;
  int nImages;
  int superpixelStepSize;
  int algo_type;
  char* config_file;
  char* overlay_dir;
  int dataset_type;
};

arguments args;


//------------------------------------------------------------------------------
void print_usage(){
 printf(
  "usage: \n \
  predict.exe -c config.txt -w model.txt \n \
  -a all: export marginals and also run inference using unary potentials only (useful for debugging) \n \
  -c config_file \n \
  -i image_dir input directory \n \
  -g algo_type \n \
  -k image_pattern \n \
  -l superpixel_labels : path of the file containing the labels for the superpixels (labels have be ordered by rows) \n \
  -m mask_dir : mask directory \n \
  -n nImages : number of images to process \n \
  -o output_dir : output filename \n \
  -s superpixelStepSize : superpixel step size \n \
  -t dataset_type : type (0=training, 1=test) \n \
  -v : verbose \n \
  -w weight_file : model obtained from training \n \
  -y : overlay directory\n");
}

/* Parse a single option. */
static int parse_opt (int key, char *arg, struct arguments *argments)
{
  switch (key)
    {
    case 'a':
      //TODO change argument from required_argument to no_argument with flag
      argments->export_all = true;
      break;
    case 'c':
      argments->config_file = arg;
      break;
    case 'g':
      if(arg!=0)
        argments->algo_type = atoi(arg);
      break;
    case 'i':
      argments->image_dir = arg;
      break;
    case 'k':
      if(arg!=0)
        argments->image_pattern = arg;
      break;
    case 'l':
      argments->superpixel_labels = arg;
      break;
    case 'm':
      argments->mask_dir = arg;
      break;
    case 'n':
      if(arg!=0)
        argments->nImages = atoi(arg);
      break;
    case 'o':
      argments->output_dir = arg;
      break;
    case 's':
      argments->superpixelStepSize = atoi(arg);
      break;
    case 't':
      if(arg!=0)
        argments->dataset_type = atoi(arg);
      break;
    case 'v':
      //TODO change argument from no_argument to no_argument with flag
      verbose = true;
      break;
    case 'w':
      argments->weight_file = arg;
      break;
    case 'y':
      argments->overlay_dir = arg;
      break;
    case 'h':
      print_usage();
      return 1;
      break;
    default:
      std::cout << "some option was wrong or missing." << std::endl;
      print_usage();
      return -1;
    }
  return 0;
}

//------------------------------------------------------------------------------

int main(int argc,char* argv[])
{
  args.image_dir = 0;
  args.mask_dir = (char*)"";
  args.nImages = -1;
  args.superpixel_labels = 0;
  args.output_dir = (char*)"./inference/";
  args.algo_type = T_GI_MULTIOBJ;
  args.weight_file = 0;
  verbose = false;
  args.image_pattern = (char*)"png";
  args.superpixelStepSize = SUPERPIXEL_DEFAULT_STEP_SIZE;
  args.config_file = 0;
  args.overlay_dir = 0;
  args.export_all = false;
  args.dataset_type = 0;
  const bool compress_image = false;

  int option_index = 0;
  int key;
  int parsing_output;

  if(argc < 2){
     fprintf(stderr, "Insufficient number of arguments. Missing configuration and model file.\n Example: predict -c config.txt -w model.txt\n usage with -h");
     exit(EXIT_FAILURE);
  }

  while((key = getopt_long(argc, argv, "ac:g:i:k:l:m:n:o:s:t:vw:y:h", long_options, &option_index)) != -1){
      parsing_output = parse_opt(key, optarg, &args);
      if(parsing_output == -1){
          fprintf(stderr, "Wrong argument. Parsing failed.");
          exit(EXIT_FAILURE);
      }else if(parsing_output == 1){
          exit(EXIT_SUCCESS);
      }
  } 

  if(args.overlay_dir == 0) {
    args.overlay_dir = args.output_dir;
  }

  string config_tmp;
  Config* config = new Config(args.config_file);
  Config::setInstance(config);

  set_default_parameters(config);

  mkdir(args.output_dir, 0777);

  string imageDir;
  string maskDir;
  if(args.image_dir != 0) {
    imageDir = args.image_dir;
  } else {
    if(args.dataset_type == 0) {
      Config::Instance()->getParameter("trainingDir", imageDir);
      Config::Instance()->getParameter("maskTrainingDir", maskDir);
    } else {
      Config::Instance()->getParameter("testDir", imageDir);
      Config::Instance()->getParameter("maskTestDir", maskDir);
    }
  }

  vector<eFeatureType> feature_types;
  int paramFeatureTypes = 0;
  if(config->getParameter("featureTypes", config_tmp)) {
    paramFeatureTypes = atoi(config_tmp.c_str());
    getFeatureTypes(paramFeatureTypes, feature_types);
  }
  
  EnergyParam param(args.weight_file);
  
  if(config->getParameter("giType", config_tmp)) {
    args.algo_type = atoi(config_tmp.c_str());
    printf("[SVM_struct] giType = %d\n", args.algo_type);
    if(args.algo_type == T_GI_MULTIOBJ && param.nClasses == 2) {
      args.algo_type = T_GI_MAXFLOW;
    }
    if(args.algo_type == T_GI_MULTIOBJ && param.nClasses > 3) {
      args.algo_type = T_GI_LIBDAI;
    }    
  }

  // TODO: Find a better way to change the labels!
  if(param.nClasses == 3) {
    printf("[svm_struct] Set class labels\n");
    BACKGROUND = 0;
    BOUNDARY = 1;
    FOREGROUND = 2;
  }

  Slice_P* slice = 0;
  Feature* feature = 0;
  int featureSize = 0;
  loadDataAndFeatures(imageDir, maskDir, config, slice, feature, &featureSize);

  // rescale features
  bool rescale_features = true;
  if(Config::Instance()->getParameter("rescale_features", config_tmp)) {
    rescale_features = config_tmp.c_str()[0] == '1';
  }
  if(rescale_features) {
    const char* scale_filename = "scale.txt";
    printf("[Main] Rescaling features\n");
    slice->rescalePrecomputedFeatures(scale_filename);
  }

  if( (args.weight_file == 0) || !fileExists(args.weight_file)) {
    printf("[Main] No parameter file provided. Loading vector of 1's\n");
    // load vector of 1's for unary term
    param.nClasses = 2;
    param.sizePsi = featureSize;
    param.nUnaryWeights = featureSize;
    param.nScalingCoefficients = featureSize;
    param.nScales = 1;
    param.weights = new double[featureSize];
    for(int f = 0; f < featureSize; ++f) {
      param.weights[f] = 1;
    }
  }

  string colormapFilename;
  getColormapName(colormapFilename);
  printf("[Main] Colormap=%s\n", colormapFilename.c_str());
  map<labelType, ulong> labelToClassIdx;
  getLabelToClassMap(colormapFilename.c_str(), labelToClassIdx);

  labelType* groundTruthLabels = 0;
  double* lossPerLabel = 0;
  const int nExamples = 1;
  const string score_filename = "";
  for(int sid = 0; sid < nExamples; sid++) {

    SPATTERN p;
    p.id = 0;
    p.slice = slice;
    p.feature = feature;

    segmentImage(p,
                 args.output_dir,
                 args.algo_type,
                 args.weight_file,
                 &labelToClassIdx,
                 score_filename,
                 groundTruthLabels, lossPerLabel,
                 compress_image,
                 args.overlay_dir,
                 METRIC_SUPERNODE_BASED_01);

    if(args.export_all) {

      if(args.algo_type != T_GI_MAX) {
        printf("[Main] Running inference with unary potentials only\n");
        // run inference with unary potentials only (algo_type = T_GI_MAX)
        stringstream sout_output;
        sout_output << args.output_dir;
        sout_output << "/unary/";
        mkdir(sout_output.str().c_str(), 0777);

        stringstream sout_overlay;
        sout_overlay << args.output_dir;
        sout_overlay << "/unary/";
        mkdir(sout_overlay.str().c_str(), 0777);

        segmentImage(p,
                     args.output_dir,
                     T_GI_MAX,
                     args.weight_file,
                     &labelToClassIdx,
                     score_filename,
                     groundTruthLabels, lossPerLabel,
                     compress_image,
                     args.overlay_dir,
                     METRIC_SUPERNODE_BASED_01);
      }

      // export marginals
      printf("[Main] Exporting marginals\n");
      GI_libDAI* gi_Inference = new GI_libDAI(slice,
                                              &param,
                                              param.weights,
                                              groundTruthLabels,
                                              lossPerLabel,
                                              feature,
                                              0, 0);

      int nNodes = slice->getNbSupernodes();
      float* marginals = new float[nNodes];
      GI_libDAI* libDAI = gi_Inference;

      for(int l = 0; l < param.nClasses; ++l) {
        libDAI->getMarginals(marginals, l);
        
        stringstream sout;
        sout << args.output_dir;
        sout << "/marginals_" << l << "/";
        mkdir(sout.str().c_str(), 0777);
        sout << getNameFromPathWithoutExtension(slice->getName());
        sout << ".png";
        printf("Exporting %s\n", sout.str().c_str());
        slice->exportProbabilities(sout.str().c_str(), param.nClasses, marginals);
      }
      
    }
  }

  printf("[Main] Cleaning\n");
  printf("[Main] Done\n");
  return 0;
}
