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

/*
 * Create superpixel features from features extracted from Overfeat  
 */

//-------------------------------------------------------------------- INCLUDES

#include <argp.h>
#include <cv.h>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

// SliceMe
#include "Slice3d.h"
#include "utils.h"
#include "globals.h"
#include "Config.h"

using namespace std;

//--------------------------------------------------------------------- GLOBALS

struct arguments
{
  char* output_filename;
  char* config_file;
  char* feature_typename;
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
  {"feature_typename",'f',  "feature_typename", 0, "feature typename"},
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
    case 'f':
      argments->feature_typename = arg;
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

string floatptr_to_string(float* ptr, int d) {
  stringstream sout;
  sout << "1:" << ptr[0];
  for(int i = 1; i < d; ++i) {
    sout << " " << (i+1) << ":" << ptr[i];
  }
  return sout.str();
}

void exportFeatures(string imageDir,
                    string maskDir,
                    Config* config,
                    string featureDir,
                    string featureTypename,
                    string outputDir)
{
  string imagePattern = "bmp";
  vector<string> lImages;
  getFilesInDir(imageDir.c_str(), lImages, imagePattern.c_str(), false);
  if(lImages.size() == 0) {
    printf("[Main] No files with extension %s in %s\n", imagePattern.c_str(), imageDir.c_str());
    exit(-1);
  }

  // load colormap
  map<ulong, labelType> classIdxToLabel;
  string colormapFilename;
  config->getParameter("colormapFilename", colormapFilename);
  printf("[SVM_struct] colormapFilename=%s\n", colormapFilename.c_str());
  getClassToLabelMap(colormapFilename.c_str(), classIdxToLabel);

  // create output directory if it doesn't exist
  if(!isDirectory(outputDir.c_str())) {
    mkdir(outputDir.c_str(), 0777);
  }

  for(uint iImage = 0; iImage < lImages.size(); iImage++) {

    // Create graph based on image data
    string imagePath = imageDir + lImages[iImage];
    PRINT_MESSAGE("[Main] Loading slice %s\n",imagePath.c_str());
    Slice* slice = new Slice(imagePath.c_str(), "");
    if(!slice) {
      printf("[Main] Error while loading slice %s\n", imagePath.c_str());
      exit(-1);
    }

    // load groundtruth image
    string baseName = getNameFromPathWithoutExtension(lImages[iImage]);
    string groundtruthName = maskDir + lImages[iImage];
    if(!fileExists(groundtruthName.c_str())) {
      // check if ground truth file for MSRC exists ?
      groundtruthName =  maskDir + baseName + "_GT.bmp";
    }
    printf("[Main] Loading %s\n", groundtruthName.c_str());
    slice->generateSupernodeLabelsFromMultiClassMaskImage(groundtruthName.c_str(),
                                                          classIdxToLabel);

    // export superpixels
    printf("img %d %d %d %d\n", slice->img_width, slice->img_height, sizeof(sidType), sizeof(sidType)*slice->img_width*slice->img_height);
    string superpixelFilename = baseName + ".dat";
    slice->exportSuperpixels(superpixelFilename.c_str());

    // load deep learning features
    string featureFilename = featureDir + baseName + "_" + featureTypename;
    cout << "Loading " << featureFilename << endl;
    ifstream ifs(featureFilename.c_str());
    string line;
    getline(ifs,line);
    vector<string> tokens0;
    splitString(line, tokens0);
    int d = atoi(tokens0[0].c_str()); // first is the number of features
    int h = atoi(tokens0[1].c_str()); // the second is the number of rows (h)
    int w = atoi(tokens0[2].c_str()); // and the last is the number of columns (w)
    int wh = w*h;
    cout << d << " " << w << " " << h << endl;

    // allocate memory
    float** features = new float*[wh];
    for(int i = 0; i < wh; ++i) {
      features[i] = new float[d];
    }

    // copy from file -> memory
    getline(ifs,line);
    vector<string> tokens;
    splitString(line, tokens);
    assert((wh*d) == tokens.size());
    int idx = 0;
    // From Overfeat:
    // "The feature is the first dimension
    // (so that to obtain the next feature, you must add w*h to your index), followed by the
    // row (to obtain the next row, add w to your index)."
    for(int j = 0; j < d; ++j) {
      for(int i = 0; i < wh; ++i) {
        features[i][j] = atof(tokens[idx].c_str());
        ++idx;
      }
    }

    ifs.close();

    string outputFilename = outputDir + baseName + ".txt";
    ofstream ofs(outputFilename.c_str());
    PRINT_MESSAGE("[Main] Exporting %s\n", outputFilename.c_str());
    // Extract features for each superpixel
    const map<sidType, supernode* >& _supernodes = slice->getSupernodes();
    for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
        it != _supernodes.end(); it++) {
      node center;
      supernode* s = it->second;
      s->getCenter(center);
      int x = center.x * w/((float)slice->img_width);
      int y = center.y * h/((float)slice->img_height);
      int feature_idx = y*w + x;
      assert(feature_idx < wh);
      ofs << (int)s->getLabel();
      ofs << " ";
      ofs << floatptr_to_string(features[feature_idx], d) << endl;
    }
    ofs.close();

    delete slice;
    for(int i = 0; i < wh; ++i) {    
      delete[] features[i];
    }
    delete[] features;

  }
}

//------------------------------------------------------------------------- MAIN

int main(int argc,char*argv[])
{
  a_args.output_filename = 0;
  a_args.feature_typename = "feat_4";

  printf("[Main] Parsing arguments\n");
  argp_parse (&argp, argc, argv, 0, 0, &a_args);

  Config* config = new Config(a_args.config_file);
  Config::setInstance(config);
  set_default_parameters(config);

  string imageDir;
  Config::Instance()->getParameter("trainingDir", imageDir);

  string maskDir;
  Config::Instance()->getParameter("maskTrainingDir", maskDir);

  string featureDir = "/home/alucchi/tmp/features/";

  string outputDir(a_args.output_filename);
  string featureTypename(a_args.feature_typename);
  cout << "featureFilename " << featureTypename << endl;
  exportFeatures(imageDir, maskDir, config, featureDir, featureTypename, outputDir);

  printf("[Main] Cleaning\n");
  delete config;
  return 0;
}
