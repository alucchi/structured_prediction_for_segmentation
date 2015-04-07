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

#include "../svm_light/svm_common.h"
#include "../svm_light/svm_learn.h"

# include "svm_struct_learn.h"
# include "svm_struct_common.h"
# include "../svm_struct_api.h"

#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

// SliceMe
#include "globals.h"
#include "Config.h"

#include "../../inference/energyParam.h"
#include "../../inference/graphInference.h"
#include "../../inference/inference.h"

//------------------------------------------------------------------------------

#define BUFFER_SIZE 250

char trainfile[200];           /* file with training examples */
char modelfile[200];           /* file for resulting classifier */

const char* loss_dir = "loss_update/";

enum eUpdateType
{
  UPDATE_NODE_ONLY = 0,
  UPDATE_NODE_EDGE,
  UPDATE_01
};

#define USE_COMBINED_MODELS 1

//------------------------------------------------------------------------------

void   read_input_parameters(int, char **, char *, char *,long *, long *,
			     STRUCT_LEARN_PARM *, LEARN_PARM *, KERNEL_PARM *,
			     int *);
void   wait_any_key();
void   print_help();

//------------------------------------------------------------------------------

void _exportNodeCoeffImage(Slice_P* slice, map<sidType, nodeCoeffType>& nodeCoeffs, const char* outputFilename)
{
  if(slice->getType() == SLICEP_SLICE) {
    IplImage* img = cvCreateImage(cvSize(slice->getWidth(), slice->getHeight()),IPL_DEPTH_32F,1);
    node n;
    float* pfData = 0;
    const map<sidType, supernode* >& _supernodes = slice->getSupernodes();
    for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
        it != _supernodes.end(); it++) {
      supernode* s = it->second;
      nodeIterator ni = s->getIterator();
      ni.goToBegin();
      while(!ni.isAtEnd()) {
        ni.get(n);
        ni.next();
        pfData = (((float*)(img->imageData + n.y*img->widthStep)+n.x*img->nChannels));
        *pfData = nodeCoeffs[it->first];
      }
    }
    //saveFloatImage(outputFilename, img);
    IplImage* u_img = 0;
    float2ucharImage(img, u_img);

    IplImage* c_img = cvLoadImage(slice->getName().c_str(), CV_LOAD_IMAGE_COLOR);
    uchar* puData;
    uchar* pcData;
    for(int y=0;y<c_img->height;y++) {
      for(int x=0;x<c_img->width;x++) {
        pcData = &((uchar*)(c_img->imageData + c_img->widthStep*y))[x*c_img->nChannels];
        puData = &((uchar*)(u_img->imageData + u_img->widthStep*y))[x*u_img->nChannels];
        *pcData = *puData;
      }
    }

    cvSaveImage(outputFilename, c_img);

    string tmp = getDirectoryFromPath(outputFilename) + getNameFromPathWithoutExtension(outputFilename);
    tmp += "_tmp.png";
    printf("tmp %s\n", tmp.c_str());
    cvSaveImage(tmp.c_str(), u_img);

    cvReleaseImage(&img);
    cvReleaseImage(&u_img);
    cvReleaseImage(&c_img);
  } else {
    sizeSliceType width = slice->getWidth();
    sizeSliceType sliceSize = width*slice->getHeight();
    sizeSliceType cubeSize = sliceSize*slice->getDepth();
    labelType* labelCube = new labelType[cubeSize];
    memset(labelCube,0,cubeSize*sizeof(labelType));

    node n;
    sizeSliceType idx;
    const map<sidType, supernode* >& _supernodes = slice->getSupernodes();
    for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
        it != _supernodes.end(); it++) {
      supernode* s = it->second;
      nodeIterator ni = s->getIterator();
      ni.goToBegin();
      while(!ni.isAtEnd()) {
        ni.get(n);
        ni.next();
        idx = n.z*sliceSize + n.y*width + n.x;
        labelCube[idx] = nodeCoeffs[it->first];
      }
    }

    exportCube(labelCube, outputFilename, slice->getDepth(), slice->getHeight(), slice->getWidth());
    delete[] labelCube;
  }
}

void exportEdgeCoeffImage(Slice_P* slice, map<sidType, edgeCoeffType>& edgeCoeffs, const char* outputFilename, double maxEdgeCoeff)
{
  if(slice->getType() == SLICEP_SLICE) {
    Slice* _slice = static_cast<Slice*>(slice);
    IplImage* _img = cvCloneImage(_slice->img);

    int sid;
    int nsid;
    node center;
    node center2;
    CvPoint pt,pt2;
    uchar r = 0;
    ulong edgeId = 0;
    map<sidType, edgeCoeffType>::iterator iEdgeWeights;
    const map<int, supernode* >& _supernodes = slice->getSupernodes();
    for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
        it != _supernodes.end(); it++) {
      supernode* s = it->second;
      sid = it->first;
      s->getCenter(center);
      for(vector < supernode* >::iterator itN = s->neighbors.begin();
          itN != s->neighbors.end(); itN++) {
        nsid = (*itN)->id;

        if(sid < nsid) {
          continue;
        }

        (*itN)->getCenter(center2);

        iEdgeWeights = edgeCoeffs.find(edgeId);
        if(iEdgeWeights == edgeCoeffs.end()) {
          printf("[Main] Error: edge key %ld for nodes (%d, %d) not found\n", edgeId, sid, nsid);
          return;
        }
        if(edgeCoeffs[edgeId] > maxEdgeCoeff) {
          r = 255;
        } else {
          r = edgeCoeffs[edgeId]*(255.0/maxEdgeCoeff);
        }

        // Draw a line between the 2 centers
        pt.x = center.x;
        pt.y = center.y;
        pt2.x = center2.x;
        pt2.y = center2.y;
        cvLine(_img, pt, pt2, CV_RGB(r,r,r),1);

        ++edgeId;
      }
    }

    printf("[Main] Exporting edge coefficient image%s\n", outputFilename);
    cvSaveImage(outputFilename, _img);
    cvReleaseImage(&_img);
  }
}

void updateCoeffs(SAMPLE& sample, EnergyParam& param,
                  const STRUCT_LEARN_PARM& struct_parm, eUpdateType updateType,
                  vector<double>* alphas, vector<double>* alphas_edge, int nRFs)
{
  long nExamples = sample.n;
  int ssvm_iteration = struct_parm.ssvm_iteration;
  int algo_type = struct_parm.giType;
  const int max_iteration_BP = 100;
  labelType* groundTruthLabels_inference = 0;
  double* lossPerLabel = 0;

  for(int i = 0; i < nExamples; i++) { /*** example loop ***/

    Slice_P* slice = sample.examples[i].x.slice;
    Feature* feature = sample.examples[i].x.feature;
    map<sidType, nodeCoeffType>* nodeCoeffs = sample.examples[i].x.nodeCoeffs;
    map<sidType, edgeCoeffType>* edgeCoeffs = sample.examples[i].x.edgeCoeffs;
    labelType* groundTruthLabels = sample.examples[i].y.nodeLabels;

    // Allocate memory
    int nNodes = slice->getNbSupernodes();
    int nEdges = slice->getNbEdges();
    labelType* inferredLabels = new labelType[nNodes];
    labelType* inferredLabels_unaryOnly = 0;

    if(updateType == UPDATE_NODE_EDGE && param.includeLocalEdges) {
      param.includeLocalEdges = false;

#if USE_COMBINED_MODELS
      inferredLabels_unaryOnly =
        computeCombinedLabels(slice, feature, groundTruthLabels_inference,
                              lossPerLabel, nRFs, alphas, i,
                              nodeCoeffs, edgeCoeffs, 0);
#else
      inferredLabels_unaryOnly = new labelType[nNodes];

      // Run inference
      GraphInference* gi_Inference =
        createGraphInferenceInstance(algo_type, slice, param, feature, 0, 0,
                                     nodeCoeffs);
      gi_Inference->run(inferredLabels_unaryOnly, 0, max_iteration_BP);
      delete gi_Inference;
#endif

      stringstream sout;
      sout << loss_dir;
      sout << getNameFromPathWithoutExtension(slice->getName());
      sout << "_";
      sout << setw(5) << setfill('0') << ssvm_iteration;
      sout << "_noEdges.tif";
      printf("[Main] Exporting %s\n", sout.str().c_str());
      slice->exportOverlay(sout.str().c_str(), inferredLabels_unaryOnly);

      stringstream soutBW;
      soutBW << loss_dir;
      soutBW << getNameFromPathWithoutExtension(slice->getName());
      soutBW << "_";
      soutBW << setw(5) << setfill('0') << ssvm_iteration;
      soutBW << "_BW_noEdges.tif";
      printf("[SVM_struct_custom] Exporting %s\n", soutBW.str().c_str());        
      slice->exportSupernodeLabels(soutBW.str().c_str(),
                                   param.nClasses,
                                   inferredLabels_unaryOnly, nNodes,
                                   &(struct_parm.labelToClassIdx));

      param.includeLocalEdges = true;
    }

#if USE_COMBINED_MODELS
    stringstream soutPb;
    soutPb << loss_dir;
    soutPb << "pb_";
    soutPb << getNameFromPathWithoutExtension(slice->getName());
    soutPb << "_";
    soutPb << setw(5) << setfill('0') << ssvm_iteration;
    soutPb << ".tif";
    printf("[Main] Exporting %s\n", soutPb.str().c_str());

    double* lossPerLabel = 0;
    inferredLabels =
      computeCombinedLabels(slice, feature, groundTruthLabels_inference,
                            lossPerLabel, nRFs, alphas, i,
                            nodeCoeffs, edgeCoeffs, soutPb.str().c_str());

#else
    // Run inference
    GraphInference* gi_Inference =
      createGraphInferenceInstance(algo_type, slice, param, feature, 0, 0,
                                   nodeCoeffs, edgeCoeffs);
    gi_Inference->run(inferredLabels, 0, max_iteration_BP);
    delete gi_Inference;
#endif

    // compute epsilon
    double epsilon = 0;
    double Z_beforeUpdate = 0;
    for(int n = 0; n < nNodes; ++n) {
      if(groundTruthLabels[n] != inferredLabels[n]) {
        epsilon += (*nodeCoeffs)[n];
      }
      Z_beforeUpdate += (*nodeCoeffs)[n];
    }
    epsilon /= Z_beforeUpdate; //nNodes;

    double alpha = 0.5*(log((1-epsilon)/epsilon));
    alphas[i].push_back(alpha);
    double alpha_pos = 1;
    double alpha_neg = 0;
    double Z = 0; //normalizing constant

    switch(updateType) {
    case UPDATE_01:
      {
        // Increase loss for misclassified nodes
        alpha_pos = 1;
        alpha_neg = 0;
        for(int n = 0; n < nNodes; ++n) {
          if(groundTruthLabels[n] != inferredLabels[n]) {
            (*nodeCoeffs)[n] = alpha_pos;
          } else {
            (*nodeCoeffs)[n] = alpha_neg;
          }
          Z += (*nodeCoeffs)[n];
        }

        double node_norm = nNodes/Z;
        for(int n = 0; n < nNodes; ++n) {
          (*nodeCoeffs)[n] *= node_norm;
        }
      }
      break;
    case UPDATE_NODE_ONLY:
      {
        // Increase loss for misclassified nodes
        alpha_pos = exp(alpha);
        alpha_neg = exp(-alpha);
        for(int n = 0; n < nNodes; ++n) {
          if(groundTruthLabels[n] != inferredLabels[n]) {
            (*nodeCoeffs)[n] = (*nodeCoeffs)[n]*(alpha_pos);
          } else {
            (*nodeCoeffs)[n] = (*nodeCoeffs)[n]*(alpha_neg);
          }
          Z += (*nodeCoeffs)[n];
        }

        double node_norm = nNodes/Z;
        for(int n = 0; n < nNodes; ++n) {
          (*nodeCoeffs)[n] *= node_norm;
        }
      }
      break;
    case UPDATE_NODE_EDGE:
      {
      // Increase loss for misclassified nodes
      alpha_pos = exp(alpha);
      alpha_neg = exp(-alpha);
      for(int n = 0; n < nNodes; ++n) {
        if(groundTruthLabels[n] != inferredLabels_unaryOnly[n]) {
          if(groundTruthLabels[n] != inferredLabels[n]) {
            (*nodeCoeffs)[n] = (*nodeCoeffs)[n]*(alpha_pos);
          } else {
            (*nodeCoeffs)[n] = (*nodeCoeffs)[n]*(alpha_neg);
          }
        } else {
          (*nodeCoeffs)[n] = (*nodeCoeffs)[n]*(alpha_neg);
        }
        Z += (*nodeCoeffs)[n];
      }

      // normalize node coefficients
      double node_norm = nNodes/Z;
      for(int n = 0; n < nNodes; ++n) {
        (*nodeCoeffs)[n] *= node_norm;
      }

      // compute epsilon for edges
      ulong edgeId = 0;
      sidType sid;
      sidType sidn;
      const map<int, supernode* >& _supernodes = slice->getSupernodes();
      double epsilon_edge = 0;
      double Ze_beforeUpdate = 0;
      for(map<sidType, supernode* >::const_iterator itNode = _supernodes.begin();
          itNode != _supernodes.end(); itNode++) {
        sid = itNode->first;
        for(vector<supernode*>::iterator itNode2 = itNode->second->neighbors.begin();
            itNode2 != itNode->second->neighbors.end(); itNode2++) {
          // set edges once
          sidn = (*itNode2)->id;
          if(sid < sidn) {
            continue;
          }

          if( (groundTruthLabels[sid] == inferredLabels_unaryOnly[sid]) &&
              (groundTruthLabels[sidn] == inferredLabels_unaryOnly[sidn]) &&
              ( (groundTruthLabels[sid] != inferredLabels[sid]) ||
                (groundTruthLabels[sidn] != inferredLabels[sidn]))) {
            // both unary predictions are correct but full model
            // (unary + pairwise) made a mistake on one of the nodes.
            epsilon_edge += (*edgeCoeffs)[edgeId];
          }

          Ze_beforeUpdate += (*edgeCoeffs)[edgeId];
          ++edgeId;
        }
      }
      epsilon_edge /= Ze_beforeUpdate; //nEdges;

      double alpha_edge = 0.5*(log((1-epsilon_edge)/epsilon_edge));
      alphas_edge[i].push_back(alpha_edge);
      double alpha_edge_pos = 1;
      double alpha_edge_neg = 0;

      double Ze = 0; //normalizing constant for edge coefficients
      edgeId = 0;
      for(map<sidType, supernode* >::const_iterator itNode = _supernodes.begin();
          itNode != _supernodes.end(); itNode++) {
        sid = itNode->first;
        for(vector<supernode*>::iterator itNode2 = itNode->second->neighbors.begin();
            itNode2 != itNode->second->neighbors.end(); itNode2++) {
          // set edges once
          sidn = (*itNode2)->id;
          if(sid < sidn) {
            continue;
          }

          if( (groundTruthLabels[sid] == inferredLabels_unaryOnly[sid]) &&
              (groundTruthLabels[sidn] == inferredLabels_unaryOnly[sidn]) &&
              ( (groundTruthLabels[sid] != inferredLabels[sid]) ||
                (groundTruthLabels[sidn] != inferredLabels[sidn]))) {
            // both unary predictions are correct but full model
            // (unary + pairwise) made a mistake on one of the nodes.
            (*edgeCoeffs)[edgeId] = (*edgeCoeffs)[edgeId]*(alpha_edge_pos);
          } else {
            (*edgeCoeffs)[edgeId] = (*edgeCoeffs)[edgeId]*(alpha_edge_neg);
          }

          Ze += (*edgeCoeffs)[edgeId];
          ++edgeId;
        }
      }

      // normalize edge coefficients
      double edge_norm = nEdges/Ze;
      double maxEdgeCoeff = 0;
      double minEdgeCoeff = 0;
      for(int e = 0; e < nEdges; ++e) {
        (*edgeCoeffs)[e] *= edge_norm;
        if(maxEdgeCoeff < (*edgeCoeffs)[e]) {
          maxEdgeCoeff = (*edgeCoeffs)[e];
        }
        if(minEdgeCoeff > (*edgeCoeffs)[e]) {
          minEdgeCoeff = (*edgeCoeffs)[e];
        }
      }
      printf("[Main] Ze=%g, maxEdgeCoeff=%g, minEdgeCoeff=%g\n", Ze, maxEdgeCoeff, minEdgeCoeff);

      stringstream soutEC;
      soutEC << loss_dir;
      soutEC << "edge_weights_";
      soutEC << getNameFromPathWithoutExtension(slice->getName());
      soutEC << "_";
      soutEC << setw(5) << setfill('0') << ssvm_iteration;
      soutEC << ".tif";
      exportEdgeCoeffImage(slice, *edgeCoeffs, soutEC.str().c_str(), maxEdgeCoeff);

      ofstream ofse("weights_edge.txt", ios::app);
      ofse << epsilon_edge << " " << alpha_edge << " " << alpha_edge_pos << " " << alpha_edge_neg << " " << Ze_beforeUpdate << " " << Ze << endl;
      ofse.close();
    }
    break;
    }

    ofstream ofs("weights.txt", ios::app);
    ofs << epsilon << " " << alpha << " " << alpha_pos << " " << alpha_neg << " " << Z_beforeUpdate << " " << Z << endl;
    ofs.close();

    stringstream sout;
    sout << loss_dir;
    sout << getNameFromPathWithoutExtension(slice->getName());
    sout << "_";
    sout << setw(5) << setfill('0') << ssvm_iteration;
    sout << ".tif";
    printf("[Main] Exporting %s\n", sout.str().c_str());
    slice->exportOverlay(sout.str().c_str(), inferredLabels);

    stringstream soutBW;
    soutBW << loss_dir;
    soutBW << getNameFromPathWithoutExtension(slice->getName());
    soutBW << "_";
    soutBW << setw(5) << setfill('0') << ssvm_iteration;
    soutBW << "_BW.tif";
    printf("[Main] Exporting %s\n", soutBW.str().c_str());        
    slice->exportSupernodeLabels(soutBW.str().c_str(),
                                 param.nClasses,
                                 inferredLabels, nNodes,
                                 &(struct_parm.labelToClassIdx));

    stringstream soutNC;
    soutNC << loss_dir;
    soutNC << "weights_";
    soutNC << getNameFromPathWithoutExtension(slice->getName());
    soutNC << "_";
    soutNC << setw(5) << setfill('0') << ssvm_iteration;
    soutNC << ".tif";
    _exportNodeCoeffImage(slice, *nodeCoeffs, soutNC.str().c_str());

    delete[] inferredLabels;
    if(inferredLabels_unaryOnly != 0) {
      delete[] inferredLabels_unaryOnly;
    }
  }
}

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

  string config_tmp;
  bool update_loss_function = false;
  if(Config::Instance()->getParameter("update_loss_function", config_tmp)) {
    update_loss_function = config_tmp.c_str()[0] == '1';
  }
  printf("[Main] update_loss_function=%d\n", (int)update_loss_function);
  if(!update_loss_function) {
    printf("update_loss_function should be true\n");
    exit(-1);
  }  

  eUpdateType updateType = UPDATE_NODE_EDGE;
  if(Config::Instance()->getParameter("update_type", config_tmp)) {
    updateType = (eUpdateType)atoi(config_tmp.c_str());
  }
  printf("[Main] update_type=%d\n", (int)updateType);

  mkdir(loss_dir, 0777);

  // check if parameter vector files exist
  const char* parameter_vector_dir = "parameter_vector";
  int idx = 0;
  string parameter_vector_dir_last = findLastFile(parameter_vector_dir, "", &idx);
  string parameter_vector_file_pattern = parameter_vector_dir_last + "/iteration_";
  int idx_1 = 1;
  string parameter_vector_file_last = findLastFile(parameter_vector_file_pattern, ".txt", &idx_1);
  printf("[Main] Checking parameter vector file %s\n", parameter_vector_file_last.c_str());

  // vectors used to store RF weights
  vector<double>* alphas = new vector<double>[sample.n];
  vector<double>* alphas_edge = 0;
  if(updateType == UPDATE_NODE_EDGE) {
    alphas_edge = new vector<double>[sample.n];
  }
  int alphas_idx = 0;

  if(fileExists("alphas.txt") && fileExists(parameter_vector_file_last)) {

    // Loading alpha coefficients
    ifstream ifs("alphas.txt");
    string line;
    int lineIdx = 0;
    while(lineIdx < sample.n && getline(ifs, line)) {
      vector<string> tokens;
      splitString(line, tokens);
      for(vector<string>::iterator it = tokens.begin();
          it != tokens.end(); ++it) {
        alphas[lineIdx].push_back(atoi(it->c_str()));
      }
      ++lineIdx;
    }
    ifs.close();
    if(lineIdx > 0) {
      alphas_idx = alphas[0].size();
    }

    // Loading parameters
    printf("[Main] Found parameter vector file %s\n", parameter_vector_file_last.c_str());
    struct_parm.ssvm_iteration = idx + 1;
    update_output_dir(struct_parm.ssvm_iteration);

    //EnergyParam param(parameter_vector_file_last.c_str());
    //updateCoeffs(sample, param, struct_parm, updateType, alphas, alphas_idx);
    //alphas_idx = 1;

  } else {
    struct_parm.ssvm_iteration = 0;

    // insert alpha coefficients for first iteration
    for(int i = 0; i < sample.n; ++i) {
      alphas[i].push_back(1);
    }

    ofstream ofs("alphas.txt", ios::app);
    int i = 0;
    for(; i < sample.n - 1; ++i) {
      ofs << alphas[i][alphas_idx] << " ";
    }
    if(i < sample.n) {
      ofs << alphas[i][alphas_idx];
    }
    ofs << endl;
    ofs.close();

    // edges
    for(int i = 0; i < sample.n; ++i) {
      alphas_edge[i].push_back(1);
    }

    ofstream ofse("alphas_edge.txt", ios::app);
    i = 0;
    for(; i < sample.n - 1; ++i) {
      ofse << alphas_edge[i][alphas_idx] << " ";
    }
    if(i < sample.n) {
      ofse << alphas_edge[i][alphas_idx];
    }
    ofse << endl;
    ofse.close();

    ++alphas_idx;
  }

  const int nMaxIterations = 5;
  bool bIterate = true;
  do {

    printf("-----------------------------------------SSVM-ITERATION-%d-START\n",
           struct_parm.ssvm_iteration);

    struct_parm.iterationId = 1;

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

    char _modelfile[BUFFER_SIZE];
    //sprintf(_modelfile, "%s_%d", modelfile, struct_parm.ssvm_iteration);
    sprintf(_modelfile, "%s_%d", modelfile, struct_parm.ssvm_iteration);
    printf("[Main] Writing learned model to %s\n", _modelfile);
    write_struct_model(_modelfile, &structmodel, &struct_parm);

    // Run inference on training data and increase loss for misclassified points
    printf("[Main] Loading learned model to %s\n", _modelfile);
    EnergyParam param(_modelfile);

    updateCoeffs(sample, param, struct_parm, updateType, alphas, alphas_edge,
                 struct_parm.ssvm_iteration + 1);

    ofstream ofs("alphas.txt", ios::app);
    int i = 0;
    for(; i < sample.n - 1; ++i) {
      ofs << alphas[i][alphas_idx] << " ";
    }
    if(i < sample.n) {
      ofs << alphas[i][alphas_idx];
    }
    ofs << endl;
    ofs.close();

    ofstream ofse("alphas_edge.txt", ios::app);
    i = 0;
    for(; i < sample.n - 1; ++i) {
      ofse << alphas_edge[i][alphas_idx] << " ";
    }
    if(i < sample.n) {
      ofse << alphas_edge[i][alphas_idx];
    }
    ofse << endl;
    ofse.close();

    ++alphas_idx;

    printf("-----------------------------------------SSVM-ITERATION-%d-END\n",
           struct_parm.ssvm_iteration);

    ++struct_parm.ssvm_iteration;

    bIterate = (nMaxIterations == -1 || struct_parm.ssvm_iteration < nMaxIterations);

  } while(bIterate);

  // Output final segmentation for all examples
  long nExamples = sample.n;
  int nRFs = struct_parm.ssvm_iteration;
  double* lossPerLabel = 0;
  labelType* groundTruthLabels = 0;

  for(int i = 0; i < nExamples; i++) { /*** example loop ***/

    Slice_P* slice = sample.examples[i].x.slice;
    Feature* feature = sample.examples[i].x.feature;
    //map<sidType, nodeCoeffType>* nodeCoeffs = sample.examples[i].x.nodeCoeffs;
    //map<sidType, edgeCoeffType>* edgeCoeffs = sample.examples[i].x.edgeCoeffs;
    map<sidType, nodeCoeffType>* nodeCoeffs = 0;
    map<sidType, edgeCoeffType>* edgeCoeffs = 0;
    int nNodes = slice->getNbSupernodes();

    stringstream soutPb;
    soutPb << loss_dir;
    soutPb << "pb_";
    soutPb << getNameFromPathWithoutExtension(slice->getName());
    soutPb << "_";
    soutPb << "combined";
    //soutPb << setw(5) << setfill('0') << ssvm_iteration;
    soutPb << ".tif";
    printf("[Main] Exporting %s\n", soutPb.str().c_str());

    labelType* inferredLabels =
      computeCombinedLabels(slice, feature, groundTruthLabels,
                            lossPerLabel, nRFs, alphas, i,
                            nodeCoeffs, edgeCoeffs, soutPb.str().c_str());

    stringstream sout;
    sout << loss_dir;
    sout << getNameFromPathWithoutExtension(slice->getName());
    sout << "_";
    sout << "combined";
    //sout << setw(5) << setfill('0') << ssvm_iteration;
    sout << ".tif";
    printf("[Main] Exporting %s\n", sout.str().c_str());
    slice->exportOverlay(sout.str().c_str(), inferredLabels);

    stringstream soutBW;
    soutBW << loss_dir;
    soutBW << getNameFromPathWithoutExtension(slice->getName());
    soutBW << "_";
    soutBW << "combined";
    //soutBW << setw(5) << setfill('0') << ssvm_iteration;
    soutBW << "_BW.tif";
    printf("[Main] Exporting %s\n", soutBW.str().c_str());        
    slice->exportSupernodeLabels(soutBW.str().c_str(),
                                 struct_parm.nClasses,
                                 inferredLabels, nNodes,
                                 &(struct_parm.labelToClassIdx));

    delete[] inferredLabels;
  }

  free_struct_sample(sample);
  free_struct_model(structmodel);

  svm_struct_learn_api_exit();

  return 0;
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



