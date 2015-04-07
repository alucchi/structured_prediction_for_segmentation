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

#include "gi_libDAI.h"

// SliceMe
#include "Config.h"
#include "utils.h"

#include "inference_globals.h"

//------------------------------------------------------------------------------

#define LIBDAI_24 1

#define MAX_POTENTIAL 5.0
//#define MAX_POTENTIAL 1.0

using namespace dai;

//------------------------------------------------------------------------------

GI_libDAI::GI_libDAI(Slice_P* _slice, 
                     const EnergyParam* _param,
                     double* _smw,
                     labelType* _groundTruthLabels,
                     double* _lossPerLabel,
                     Feature* _feature,
                     map<sidType, nodeCoeffType>* _nodeCoeffs,
                     map<sidType, edgeCoeffType>* _edgeCoeffs)
{
  GraphInference::init();
  slice = _slice;
  param = _param;
  smw = _smw;
  lossPerLabel = _lossPerLabel;
  groundTruthLabels = _groundTruthLabels;
  feature = _feature;
  nodeCoeffs = _nodeCoeffs;
  edgeCoeffs = _edgeCoeffs;
  maxPotential = 0;
  unaryPotentials = 0;
  nUnaryPotentials = 0;
  nEdgePotentials = 0;
  allocateGraph();
  buildSSVMGraph();
}

void GI_libDAI::allocateGraph()
{
  // Create variables.
  const ulong nSupernodes = slice->getNbSupernodes();
  ulong nFactors = nSupernodes;

  if(param->includeLocalEdges) {
    nFactors += slice->getNbEdges();
  }

  vars.clear();
  vars.reserve(nSupernodes);
  uint varId = 0;
  for(; varId < nSupernodes; varId++) {
    vars.push_back( Var(varId, param->nClasses) );
  }

  factors.clear();
  factors.reserve(nFactors);
}

/**
 * Run BP on a given factor graph
 * @param nodeLabelsBP node inferred by BP
 */
double GI_libDAI::run(labelType* inferredLabels,
                      int id,
                      size_t maxiter,
                      labelType* nodeLabelsGroundTruth,
                      bool computeEnergyAtEachIteration,
                      double* _loss)
{
  double energy = 0;
  double loss = 0;

  string paramMSRC;
  Config::Instance()->getParameter("msrc", paramMSRC);
  bool useMSRC = paramMSRC.c_str()[0] == '1';
  bool replaceVoidMSRC = false;
  if(useMSRC) {
    Config::Instance()->getParameter("msrc_replace_void", paramMSRC);
    replaceVoidMSRC = paramMSRC.c_str()[0] == '1';
  }

  // Set some constants
  Real   tol = 1e-7;
  //Real   tol = 1e-9;
  size_t verb = 0;
  //size_t verb = 3;

  // Store the constants in a PropertySet object
  PropertySet opts;
#if LIBDAI_24
  opts.Set("maxiter",maxiter);  // Maximum number of iterations
  opts.Set("tol",tol);          // Tolerance for convergence
  opts.Set("verbose",verb);     // Verbosity (amount of output generated)
  opts.Set("logdomain",false);
  opts.Set("inference",string("MAXPROD"));
  opts.Set("updates",string("SEQFIX"));
#else
  opts.set("maxiter",maxiter);  // Maximum number of iterations
  opts.set("tol",tol);          // Tolerance for convergence
  opts.set("verbose",verb);     // Verbosity (amount of output generated)
  opts.set("logdomain",false);
  opts.set("inference",string("MAXPROD"));
  opts.set("updates",string("SEQFIX"));
#endif

  FactorGraph fg( factors.begin(), factors.end(), vars.begin(), vars.end(), factors.size(), vars.size() );

  // Construct a BP (belief propagation) object from the FactorGraph fg
  // using the parameters specified by opts and two additional properties,
  // specifying the type of updates the BP algorithm should perform and
  // whether they should be done in the real or in the logdomain
  //BP bp(fg, opts("updates",string("SEQFIX"))("logdomain",true));
  //BP bp(fg, opts("updates",string("SEQFIX"))("logdomain",false));
  //BP bp(fg, opts("updates",string("SEQFIX"))("logdomain",false)("inference",string("MAXPROD")));
  BP bp(fg, opts);
  // Initialize belief propagation algorithm
  bp.init();

  vector<std::size_t> labels;

  // Run belief propagation algorithm
  if(computeEnergyAtEachIteration) {

#if OUTPUT_ENERGY
      // one file per example
      stringstream sEnergyMaxFile;
      sEnergyMaxFile << "x_" << id;
      sEnergyMaxFile << "_energyBPMax.txt";

      ofstream ofsEnergyMax(sEnergyMaxFile.str().c_str(),ios::app);
#endif

      int i = maxiter;
      //for(int i = 1; i <= (int)maxiter; i+=10) // loop to see how energy evolves
        {
          //INFERENCE_PRINT("[gi_libDAI] BP Iteration %d\n", i);
#if LIBDAI_24
          opts.Set("maxiter",(size_t)i);  // Maximum number of iterations
#else
          opts.set("maxiter",(size_t)i);  // Maximum number of iterations
#endif
          bp.setProperties(opts);
          //Real maxDiff = bp.run();
          bp.run();

          labels = bp.findMaximum();

          if(replaceVoidMSRC) {
            if(lossPerLabel == 0) {
              copyLabels_MSRC(labels, inferredLabels, bp, fg);                
            } else {
              copy(labels.begin(),labels.end(),inferredLabels);
            }
          } else {
            copy(labels.begin(),labels.end(),inferredLabels);
          }

          energy = GraphInference::computeEnergy(inferredLabels);

#if OUTPUT_ENERGY
          ofsEnergyMax << energy << " " << loss << endl;
#endif

          //if( maxDiff < tol )
          //  break;
        }

#if OUTPUT_ENERGY
        ofsEnergyMax.close();
#endif

  } else {
    bp.run();
    labels = bp.findMaximum();

    if(replaceVoidMSRC) {
      if(lossPerLabel == 0) {
	copyLabels_MSRC(labels,
			inferredLabels,
			bp,
			fg);
      } else {
	copy(labels.begin(),labels.end(),inferredLabels);
      }
    } else {
      copy(labels.begin(),labels.end(),inferredLabels);
    }

    energy = GraphInference::computeEnergy(inferredLabels);
  }

  if(_loss) {
    *_loss = loss;
  }

  return energy;
}

// use bp.findmaximum instead
void GI_libDAI::getLabels(BP& bp,
                          FactorGraph& fg,
                          vector<std::size_t>& labels)
{
  float maxProb;
  int label;
  ofstream ofs("beliefs2");
  ofs << "------------------------------- UNARY -------------------------------\n";
  for(int sid = 0; sid < slice->getNbSupernodes(); sid++) {
    Factor f = bp.belief(fg.var(sid));
    maxProb = f[0];
    label = 0;
    ofs << sid;
    for(int i=1; i < (int)param->nClasses; i++) {
      ofs << " " << i << ":" << f[i];
      if(f[i] > maxProb) {
        maxProb = f[i];
        label = i;
      }
    }
    ofs << endl;
    labels[sid] = label;
  }
  ofs.close();
}

void GI_libDAI::getMarginals(float* marginals,
                             const int label)
{
  // Set some constants
  Real   tol = 1e-7;
  //Real   tol = 1e-9;
  size_t verb = 0;
  //size_t verb = 3;
  size_t maxiter = 100;

  // Store the constants in a PropertySet object
  PropertySet opts;
#if LIBDAI_24
  opts.Set("maxiter",maxiter);  // Maximum number of iterations
  opts.Set("tol",tol);          // Tolerance for convergence
  opts.Set("verbose",verb);     // Verbosity (amount of output generated)
  opts.Set("logdomain",false);
  opts.Set("inference",string("MAXPROD"));
  opts.Set("updates",string("SEQFIX"));
#else
  opts.set("maxiter",maxiter);  // Maximum number of iterations
  opts.set("tol",tol);          // Tolerance for convergence
  opts.set("verbose",verb);     // Verbosity (amount of output generated)
  opts.set("logdomain",false);
  opts.set("inference",string("MAXPROD"));
  opts.set("updates",string("SEQFIX"));
#endif

  FactorGraph fg( factors.begin(), factors.end(), vars.begin(), vars.end(), factors.size(), vars.size() );
  BP bp(fg, opts);
  bp.init();
  bp.run();

  for(int sid = 0; sid < slice->getNbSupernodes(); sid++) {
    Factor f = bp.belief(fg.var(sid));
    marginals[sid] = f[label];
  }
}

void GI_libDAI::computeNodePotentials()
{
  Real *buf = new Real[param->nClasses];
  factors.clear();

  // allocate memory to store features
  int fvSize = feature->getSizeFeatureVector();
  int max_index = fvSize + 1;
  osvm_node* n = new osvm_node[max_index];
  int i = 0;
  for(i = 0;i < max_index-1; i++)
    n[i].index = i+1;
  n[i].index = -1;

  bool useLossFunction = lossPerLabel!=0;

  string config_tmp;
  int loss_function = 0;
  if(Config::Instance()->getParameter("loss_function", config_tmp)) {
    loss_function = atoi(config_tmp.c_str());
  }

  int sid = 0;
  const map<int, supernode* >& _supernodes = slice->getSupernodes();
  nUnaryPotentials = _supernodes.size();
  unaryPotentials = new Real*[nUnaryPotentials];
  for (uint k = 0; k < nUnaryPotentials; ++k) {
    unaryPotentials[k] = new Real[param->nClasses];
  }

  maxPotential = 0;
  for(map<int, supernode* >::const_iterator its = _supernodes.begin();
      its != _supernodes.end(); its++) {
    sid = its->first;

    if(param->nClasses != 2) {
      for(int i = 0; i < (int)param->nClasses; i++) {
        double unaryPotential = computeUnaryPotential(slice, sid, i);

        if(nodeCoeffs) {
          unaryPotential *= (*nodeCoeffs)[sid];
        }
        buf[i] = unaryPotential;

        if (fabs(buf[i]) > maxPotential) {
          maxPotential = fabs(buf[i]);
        }
      }
    } else {
      // Only 2 classes.
      buf[T_FOREGROUND] = 0;
      const int i = T_BACKGROUND;
      double unaryPotential = computeUnaryPotential(slice, sid, i);

      if(nodeCoeffs) {
        unaryPotential *= (*nodeCoeffs)[sid];
      }

      buf[i] = unaryPotential;

      if (fabs(buf[i]) > maxPotential) {
        maxPotential = fabs(buf[i]);
      }
    }

    if(useLossFunction) {
      // loss function
      for(int s = 0; s < (int)param->nClasses; s++) {
        if(s != groundTruthLabels[sid]) {
          // add loss of the ground truth label

          int _loss = 0;
          if(loss_function == LOSS_NODE_BASED) {
            _loss = lossPerLabel[sid];
          } else {
            _loss = lossPerLabel[groundTruthLabels[sid]];
          }

          if(nodeCoeffs) {
            _loss *= (*nodeCoeffs)[sid];
          }
          buf[s] += _loss;

          if (fabs(buf[s]) > maxPotential) {
            maxPotential = fabs(buf[s]);
          }
        }
        //INFERENCE_PRINT("3.%d:%e ", s, buf[s]);
      }

      // loss function : -1 if correct label (same as adding +1 to all incorrect labels)
      //buf[groundTruthLabels[sid]] = buf[groundTruthLabels[sid]]-1;
    }

    for(int c = 0; c < (int)param->nClasses; c++) {
      unaryPotentials[sid][c] = buf[c];
    }
    //INFERENCE_PRINT("\n");
  }
  delete[] n;
  delete[] buf;
}

void GI_libDAI::addUnaryNodes()
{
  Real scale = 1;
  if (maxPotential != 0) {
    scale = MAX_POTENTIAL/maxPotential;
  }
  INFERENCE_PRINT("[gi_libDAI] maxPotential=%g, scale=%g\n", maxPotential, scale);
  Real *buf = new Real[param->nClasses];
  for (uint sid = 0; sid < nUnaryPotentials; ++sid) {
    for(int c = 0; c < (int)param->nClasses; c++) {
      buf[c] = std::exp(unaryPotentials[sid][c]*scale);
    }
    factors.push_back(Factor(vars[sid], buf));
  }

  for (uint k = 0; k < nUnaryPotentials; ++k) {
    delete[] unaryPotentials[k];
  }
  delete[] unaryPotentials;
  delete[] buf;
}

void GI_libDAI::precomputePotentials()
{
  computeNodePotentials();
  if(param->includeLocalEdges) {
    computeEdgePotentials();
  }
}

void GI_libDAI::addNodes()
{
  addUnaryNodes();
  if(param->includeLocalEdges) {
    addLocalEdges();
  }
}

void GI_libDAI::buildSSVMGraph()
{
  precomputePotentials();
  addNodes();
}

void GI_libDAI::computeEdgePotentials()
{
  // Pairwise factors
  if(param->nGradientLevels > 0) {
    int nPairwiseStates = param->nClasses*param->nClasses;
    int gradientIdx;
    int orientationIdx;
    int idx;
    int oidx;
    supernode *s;
    const map<int, supernode* >& _supernodes = slice->getSupernodes();
    nEdgePotentials = slice->getNbEdges();
    INFERENCE_PRINT("[gi_libDAI] Allocating memory for %d edges and %d states\n",
                    nEdgePotentials, nPairwiseStates);
    for (uint k = 0; k < nEdgePotentials; ++k) {
      edgePotentials[k] = new Real[nPairwiseStates];
    }

    uint edgeId = 0;
    for(map<int, supernode* >::const_iterator its = _supernodes.begin();
        its != _supernodes.end(); its++) {
      s = its->second;
      vector < supernode* >* lNeighbors = &(s->neighbors);
      for(vector < supernode* >::iterator itN = lNeighbors->begin();
          itN != lNeighbors->end(); itN++) {

        // set edges once
        if(its->first < (*itN)->id) {
          continue;
        }

        gradientIdx = slice->getGradientIdx(its->first,(*itN)->id);
        orientationIdx = slice->getOrientationIdx(its->first,(*itN)->id);
        oidx = orientationIdx*param->nClasses*param->nClasses;
#if USE_LONG_RANGE_EDGES
        int distanceIdx = slice->getDistanceIdx(s->id, (*itN)->id);
        oidx += distanceIdx*param->nGradientLevels*param->nClasses*param->nClasses*param->nOrientations;
#endif

        map<uint, Real*>::iterator lookup = edgePotentials.find(edgeId);
        if(lookup == edgePotentials.end()) {
          edgePotentials[edgeId] = new Real[nPairwiseStates];
        }

        for(int p = 0; p < nPairwiseStates; p++ ) {
          // w[0..nClasses-1] contains the unary weights
          double w_sum = 0;
          for(int i = 0; i <= gradientIdx; i++) {
            idx = (i*param->nClasses*param->nClasses*param->nOrientations) + oidx + p;
            w_sum += smw[idx+param->nUnaryWeights]; // param->nUnaryWeights is the offset due to unary terms
          }

          if(edgeCoeffs) {
            w_sum *= (*edgeCoeffs)[edgeId];
          }

          edgePotentials[edgeId][p] = w_sum;

          if (fabs(w_sum) > maxPotential) {
            maxPotential = fabs(w_sum);
          }
        }
        ++edgeId;
      }
    }
  } else {
    // potts model
    double edgeWeight = fabs(smw[param->nUnaryWeights]);
    /*
    if(edgeCoeffs) {
    //TODO : Iterate over all edges and take max ?
    }
    */
    if(edgeWeight > maxPotential) {
      maxPotential = edgeWeight;
    }
  }
}

void GI_libDAI::addLocalEdges()
{
  Real scale = 1;
  if (maxPotential != 0) {
    scale = MAX_POTENTIAL/maxPotential;
  }
  INFERENCE_PRINT("[gi_libDAI] maxPotential=%g, scale=%g\n", maxPotential, scale);

  if(param->nGradientLevels > 0) {
    supernode *s;
    uint edgeId = 0;
    int nPairwiseStates = param->nClasses*param->nClasses;
    const map<int, supernode* >& _supernodes = slice->getSupernodes();
    for(map<int, supernode* >::const_iterator its = _supernodes.begin();
        its != _supernodes.end(); its++) {
      s = its->second;

      vector < supernode* >* lNeighbors = &(s->neighbors);
      for(vector < supernode* >::iterator itN = lNeighbors->begin();
          itN != lNeighbors->end(); itN++) {

        // set edges once
        if(its->first < (*itN)->id) {
          continue;
        }

        Factor fac( VarSet( vars[its->first], vars[(*itN)->id] ), 1.0 );
        for(int p = 0; p < nPairwiseStates; p++ ) {
#if LIBDAI_24
          fac[p] = std::exp(edgePotentials[edgeId][p]*scale);
#else
          fac.set(p, std::exp(edgePotentials[edgeId][p]*scale));
#endif
        }
        factors.push_back(fac);
        ++edgeId;
      }
    }

    for (uint k = 0; k < nEdgePotentials; ++k) {
      delete[] edgePotentials[k];
    }
    edgePotentials.clear();
  } else {
    // potts model
    double gamma = smw[param->nUnaryWeights];
    if(maxPotential != 0) {
      gamma *= MAX_POTENTIAL/maxPotential;
    }
    printf("gi_libDAI potts w[%d]=%g %g\n", param->nUnaryWeights, smw[param->nUnaryWeights], gamma);
    ulong edgeId = 0;
    const map<int, supernode* >& _supernodes = slice->getSupernodes();
    for(map<int, supernode* >::const_iterator its = _supernodes.begin();
        its != _supernodes.end(); its++) {
      vector < supernode* >* lNeighbors = &(its->second->neighbors);
      for(vector < supernode* >::iterator itN = lNeighbors->begin();
          itN != lNeighbors->end(); itN++)  {
        assert(its->first != (*itN)->id);
        // set edges once
        if(its->first < (*itN)->id) {
          continue;
        }

        //Factor fac = createFactorPotts( vars[its->first], vars[(*itN)->id],gamma);
        //Factor fac = createFactorIsing( vars[its->first], vars[(*itN)->id],gamma);
        Factor fac( VarSet(vars[its->first], vars[(*itN)->id]), 1.0 );
        fac[0] = std::exp(gamma);
        if(edgeCoeffs) {
          fac[0] = fac[0]*(*edgeCoeffs)[edgeId];
        }
        fac[3] = fac[0];
        factors.push_back(fac);
        ++edgeId;
      }
    }
  }
}

/**
 * The output image can be used to compare with the ground-truth data
 * use BP::findMaximum()
 */
void GI_libDAI::createColoredAnnotationImageMax(vector<std::size_t>& labels,
                                                const char* outputFilename,
                                                map<labelType, ulong>* labelToClassIdx)
{
  IplImage* img = cvCreateImage(cvSize(slice->getWidth(), slice->getHeight()),
                                IPL_DEPTH_8U,3);

  supernode* s;
  uchar *pData;
  node n;
  uchar r,g,b;
  ulong classIdx = 0;
  int label;
  int *counter = new int[param->nClasses];
  for(int i = 0; i < (int)param->nClasses; i++)
    counter[i]=0;

  const map<int, supernode* >& _supernodes = slice->getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++)
    {
      label = labels[it->first];
      counter[label]++;
      //INFERENCE_PRINT("%d:%d ", it->first,label);

      if(labelToClassIdx == 0)
        {
          // black and white images
          if(label == FOREGROUND)
            {
              r = g = b = 255;
            }
          else
            {
              if(label == BOUNDARY)
                {
                  r = g = b = 128;
                }
              else
                {
                  r = g = b = 0;
                }
            }
        }
      else
        {
          classIdx = (*labelToClassIdx)[label];
          classIdxToRGB(classIdx,r,g,b);
        }

      s = it->second;
      nodeIterator ni = s->getIterator();
      ni.goToBegin();
      while(!ni.isAtEnd())
        {
          ni.get(n);
          ni.next();
          pData = (((uchar*)(img->imageData + n.y*img->widthStep)+n.x*img->nChannels));
          pData[0] = b;
          pData[1] = g;
          pData[2] = r;
        }
    }

  INFERENCE_PRINT("[gi_libDAI] ");
  for(int i = 0; i < (int)param->nClasses; i++)
    INFERENCE_PRINT("counter[%d]=%d ",i,counter[i]);
  INFERENCE_PRINT("\n");

  INFERENCE_PRINT("[gi_libDAI] Saving image to %s\n", outputFilename);
  cvSaveImage(outputFilename,img);
  cvReleaseImage(&img);

  delete[] counter;
}


void GI_libDAI::buildPottsModel(double gamma)
{
  // Unary factors
  Real *buf = new Real[param->nClasses];

  int sid = 0;
  const map<int, supernode* >& _supernodes = slice->getSupernodes();
  for(map<int, supernode* >::const_iterator its = _supernodes.begin();
      its != _supernodes.end(); its++)
    {
      for(int i = 0; i < param->nClasses; i++) {
        buf[i] = min(1.0,slice->getProb(its->first,i)+DELTA_PB);
      }

      factors.push_back( Factor(vars[its->first], &buf[0]));
      sid++;
    }

  // pairwise factors
  for(map<int, supernode* >::const_iterator its = _supernodes.begin();
      its != _supernodes.end(); its++) {
    vector < supernode* >* lNeighbors = &(its->second->neighbors);
    for(vector < supernode* >::iterator itN = lNeighbors->begin();
        itN != lNeighbors->end(); itN++) {
      // set edges once
      if(its->first < (*itN)->id)
        continue;

      Factor fac = createFactorPotts( vars[its->first], vars[(*itN)->id],gamma);
      factors.push_back(fac);
    }
  }
  delete[] buf;
}

void GI_libDAI::copyLabels_MSRC(vector<std::size_t>& labels,
                                labelType* nodeLabels,
                                BP& bp,
                                FactorGraph& fg)
{
  labelType voidLabel = classIdxToLabel[0];
  labelType moutainLabel = classIdxToLabel[4161600];
  labelType horseLabel = classIdxToLabel[8323328];

  INFERENCE_PRINT("[gi_libDAI] copyLabels_MSRC void=%d, moutain=%d, horse=%d\n",
                  (int)voidLabel, (int)moutainLabel, (int)horseLabel);

  int label;
  int n = slice->getNbSupernodes();
  for(int sid = 0; sid < n; sid++) {
    label = labels[sid];

    //if(lossPerLabel == 0 && label > 20) // void label
    if(label == voidLabel || label == moutainLabel || label == horseLabel) {
      //Factor f = bp.beliefV(sid);
      Factor f = bp.belief(fg.var(sid));
      double maxProb = -1;
      for(int i = 0; i < (int)f.states(); i++) {
        if(i != voidLabel && i != moutainLabel && i != horseLabel)
          if(maxProb < f[i]) {
            maxProb = f[i];
            label = i;
          }
      }
      //INFERENCE_PRINT("label=%d\n", label);
    }
    nodeLabels[sid] = label;
  }
}

