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
// Contact aurelien.lucchi (at) gmail.com                              // 
// for comments & bug reports                                          //
/////////////////////////////////////////////////////////////////////////

#include <mex.h>
#include <stdio.h>

#ifdef _WIN32
#include <windows.h>
#endif

// SliceMe
#include <Config.h>
#include <Feature.h>
#include <Slice.h>

#include <energyParam.h>
#include <graphInference.h>
#include <inference.h>
#include <svm_struct_api.h>

#include <globals.h>


void mexFunction(int nlhs,       mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
  char *pParamFilename;
  uchar *pResult;
  const mwSize *dim_array;
    
  /* Check for proper number of input and output arguments */    
  if (nrhs < 2) {
    mexErrMsgTxt("Usage: mex_predict config_file parameter_file dataset_id(0=training, 1=testing).\n");
  }

  char* config_file;
  int n0 = mxGetN(prhs[0]);
  config_file = (char *)mxCalloc(n0+1,sizeof(char));
  mxGetString(prhs[0],config_file,n0+1);
  mexPrintf("Loading config file %s\n", config_file);
  Config* config = new Config(config_file);
  Config::setInstance(config);
  set_default_parameters(config);
  mxFree(config_file);

  int n1 = mxGetN(prhs[1]);
  pParamFilename = (char *)mxCalloc(n1+1,sizeof(char));
  mxGetString(prhs[1],pParamFilename,n1+1);

  mexPrintf("Loading parameter file %s\n", pParamFilename);
  EnergyParam param(pParamFilename);

  string imageDir;
  string maskDir;
  if(nrhs > 2) {
    char* input_file;
    int n = mxGetN(prhs[2]);
    input_file = (char *)mxCalloc(n+1,sizeof(char));
    mxGetString(prhs[2],input_file,n+1);

    imageDir = input_file;
    mxFree(input_file);
  } else {
    if(1) {
      Config::Instance()->getParameter("trainingDir", imageDir);
      Config::Instance()->getParameter("maskTrainingDir", maskDir);
    } else {
      Config::Instance()->getParameter("testDir", imageDir);
      Config::Instance()->getParameter("maskTestDir", maskDir);
    }

    string config_tmp;
    Config::Instance()->getParameter("slice3d", config_tmp);
    bool useSlice3d = config_tmp.c_str()[0] == '1';
    if(!useSlice3d) {
      vector<string> files;
      getFilesInDir(imageDir.c_str(),files,"png");
      if(files.size() > 0) {
        imageDir += files[0];
      }
    }
  }

  verbose = true;

  mexPrintf("Opening %s\n", imageDir.c_str());

  Slice_P* slice = 0;
  Feature* feature = 0;
  int featureSize = 0;
  loadDataAndFeatures(imageDir, maskDir, config, slice, feature, &featureSize, 0);

  STRUCT_LEARN_PARM sparm;
  energyParamToSparm(param, &sparm);

  string config_tmp;
  int giType = T_GI_LIBDAI;
  if(config->getParameter("giType", config_tmp)) {
    giType = atoi(config_tmp.c_str());
  }
  bool useGC = isSubmodular(&sparm, param.weights);
  if(useGC) {
    giType = T_GI_MAXFLOW;
  }

  labelType* labels = computeLabels(slice, feature, param,
                                    giType, 0);

  // Create output matrix
  const mwSize dims[] = {slice->getWidth(), slice->getHeight(), slice->getDepth()};
  plhs[0] = mxCreateNumericArray(3,dims,mxUINT8_CLASS,mxREAL);
  pResult = (uchar*)mxGetData(plhs[0]);

  node n;
  ulong idx = 0;
  labelType label;
  supernode* s;
  unsigned long width = slice->getWidth();
  unsigned long sliceSize = slice->getWidth()*slice->getHeight();
  const map<sidType, supernode* >& _supernodes = slice->getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    s = it->second;
    label = labels[it->first];

    nodeIterator ni = s->getIterator();
    ni.goToBegin();
    while(!ni.isAtEnd()) {
      ni.get(n);
      ni.next();
      idx = ((ulong)n.z*sliceSize) + n.y*width + n.x;
      pResult[idx] = (label/(float)param.nClasses)*255;
    }
  }

  delete[] labels;
  delete slice;
  delete feature;
}
