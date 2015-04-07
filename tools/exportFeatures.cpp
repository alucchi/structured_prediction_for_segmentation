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

//-------------------------------------------------------------------- INCLUDES

#include <argp.h>
#include <cv.h>
#include <highgui.h>
#include <set>
#include <sstream>
#include <sys/stat.h>

// SliceMe
#include "Slice3d.h"
#include "utils.h"
#include "globals.h"
#include "Config.h"
#include "F_Precomputed.h"

using namespace std;

//--------------------------------------------------------------------- GLOBALS

struct arguments
{
  char* output_filename;
  char* config_file;
};

struct arguments a_args;

//----------------------------------------------------------------------- PARSER

/* Program documentation. */
static char doc[] =
  "Export features";

/* A description of the arguments we accept. */
static char args_doc[] = "";

/* The options we understand. */
static struct argp_option options[] = {
  {"config_file",'c',  "config_file",0, "config_file"},
  {"output_filename",'o',  "output_filename", 0, "output filename"},
  {"verbose",'v',  "verbose",0, "verbose"},
  { 0 }
};


/* Parse a single option. */
static error_t
parse_opt (int key, char *arg, struct argp_state *state)
{
  /* Get the input argument from argp_parse, which we
     know is a pointer to our arguments structure. */
  struct arguments *argments = (arguments*)state->input;

  switch (key)
    {
    case 'c':
      argments->config_file = arg;
      break;
    case 'o':
      argments->output_filename = arg;
      break;
    case 'v':
      if(*arg == '1')
        verbose = true;
      else
        verbose = false;
      break;

    case ARGP_KEY_ARG:
      switch(state->arg_num)
        {
          //case 0:
          //argments->image = arg;
          //break;
        default:
          // Too many arguments
          printf("Too many arguments %s\n", arg);
          argp_usage (state);
          break;
        }
      break;

    case ARGP_KEY_END:
      /* Not enough arguments. */
      //if (state->arg_num < 2)
      //  argp_usage (state);
      break;

    default:
      return ARGP_ERR_UNKNOWN;
    }
  return 0;
}

/* Our argp parser. */
static struct argp argp = { options, parse_opt, args_doc, doc };

//-------------------------------------------------------------------- FUNCTIONS

void exportFeatures(string imageDir,
                    string maskDir,
                    Config* config,
                    char* _outputFilename)
{
  Slice3d* slice3d = new Slice3d(imageDir.c_str());
  if(!slice3d) {
    printf("[SVM_struct] Error while loading %s\n", imageDir.c_str());
    exit(-1);
  }

  string config_tmp;

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

  string outputFilename;
  if(_outputFilename == 0) {
    stringstream sout_feature_filename;
    sout_feature_filename << slice3d->inputDir << "/features_";
    sout_feature_filename << slice3d->getSupernodeStep() << "_" << slice3d->getCubeness();
    sout_feature_filename << "_" << paramFeatureTypes;
    sout_feature_filename << "_" << DEFAULT_FEATURE_DISTANCE;
    outputFilename = sout_feature_filename.str();
  } else {
    outputFilename = _outputFilename;
  }

  printf("[SVM_struct] Checking if %s exists\n", outputFilename.c_str());
  Feature* feature = 0;
  bool featuresLoaded = false;
  if(fileExists(outputFilename)) {
    printf("[SVM_struct] Loading features from %s\n", outputFilename.c_str());
    int featureSize = -1;
    if(slice3d->loadFeatures(outputFilename.c_str(), &featureSize)) {
      featuresLoaded = true;
      feature = new F_Precomputed(slice3d->getPrecomputedFeatures(), featureSize/DEFAULT_FEATURE_DISTANCE);
      printf("[SVM_struct] Features Loaded succesfully\n");
    } else {
      printf("[SVM_struct] Features not loaded succesfully\n");
    }
  } else {
    printf("[SVM_struct] File %s does not exist\n", outputFilename.c_str());
  }

  if(!featuresLoaded) {
    feature = Feature::getFeature(slice3d, feature_types);

    int featureSize = feature->getSizeFeatureVector();
    printf("[SVM_struct] Feature size = %d\n", featureSize);
    slice3d->precomputeFeatures(feature);

    // Dump features
    if(feature) {
      feature->save(*slice3d, outputFilename.c_str());
    }

  }

}

//------------------------------------------------------------------------- MAIN

int main(int argc,char*argv[])
{
  a_args.output_filename = 0;

  printf("[Main] Parsing arguments\n");
  argp_parse (&argp, argc, argv, 0, 0, &a_args);

  Config* config = new Config(a_args.config_file);
  Config::setInstance(config);
  set_default_parameters(config);

  string imageDir;
  Config::Instance()->getParameter("trainingDir", imageDir);

  string maskDir;
  Config::Instance()->getParameter("maskTrainingDir", maskDir);

  exportFeatures(imageDir, maskDir, config, a_args.output_filename);

  printf("[Main] Cleaning\n");
  delete config;
  return 0;
}
