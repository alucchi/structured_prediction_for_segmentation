/***********************************************************************/
/*                                                                     */
/*   svm_struct_api_types.h                                            */
/*                                                                     */
/*   Definition of API for attaching implementing SVM learning of      */
/*   structures (e.g. parsing, multi-label classification, HMM)        */ 
/*                                                                     */
/*   Author: Thorsten Joachims                                         */
/*   Date: 13.10.03                                                    */
/*                                                                     */
/*   Copyright (c) 2003  Thorsten Joachims - All rights reserved       */
/*                                                                     */
/*   This software is available for non-commercial use only. It must   */
/*   not be modified and distributed without prior permission of the   */
/*   author. The author is not responsible for implications from the   */
/*   use of this software.                                             */
/*                                                                     */
/***********************************************************************/

#ifndef svm_struct_api_types
#define svm_struct_api_types

#include "svm_light/svm_common.h"
#include "svm_light/svm_learn.h"

// SliceMe
#include "Feature.h"
#include "Slice_P.h"

// nodeCoeff, edgeCoeff types
#include "inference_globals.h"

#include <map>

# define INST_NAME          "Generic and empty API"
# define INST_VERSION       "V0.00"
# define INST_VERSION_DATE  "??.??.??"

/* default precision for solving the optimization problem */
# define DEFAULT_EPS         0.1 
/* default loss rescaling method: 1=slack_rescaling, 2=margin_rescaling */
# define DEFAULT_RESCALING   2
/* default loss function: */
# define DEFAULT_LOSS_FCT    0
/* default optimization algorithm to use: */
# define DEFAULT_ALG_TYPE    9
/* store Psi(x,y) (for ALG_TYPE 1) instead of recomputing it every time: */
# define USE_FYCACHE         1
/* decide whether to evaluate sum before storing vectors in constraint
   cache: 
   0 = NO, 
   1 = YES (best, if sparse vectors and long vector lists), 
   2 = YES (best, if short vector lists),
   3 = YES (best, if dense vectors and long vector lists) */
# define COMPACT_CACHED_VECTORS 1
/* minimum absolute value below which values in sparse vectors are
   rounded to zero. Values are stored in the FVAL type defined in svm_common.h 
   RECOMMENDATION: assuming you use FVAL=float, use 
     10E-15 if COMPACT_CACHED_VECTORS is 1 
     10E-10 if COMPACT_CACHED_VECTORS is 2 or 3 
*/
# define COMPACT_ROUNDING_THRESH 10E-15

#define MAX_SIZE_EXT 10

typedef struct pattern {
  /* this defines the x-part of a training example, e.g. the structure
     for storing a natural language sentence in NLP parsing */
  //int add_your_variables_here;
  int id;
  Slice_P* slice;
  Feature* feature;
  IplImage* imgAnnotation; // ground-truth image (null for 3d cubes)
  uchar* cubeAnnotation; // ground-truth cube (null for 2d images)
  int nEdges;
  //string imageFilename;

  // variables used to store TP,FP and FN for each example.
  // Those results are written to a file after each iteration (see finalize_iteration)
  ulong* TPs;
  ulong* FPs;
  ulong* FNs;
  ulong* count;

  // structure used to score coefficients for each node.
  map<sidType, nodeCoeffType>* nodeCoeffs;

  // structure used to score coefficients for each edge.
  map<sidType, edgeCoeffType>* edgeCoeffs;

  pattern() {
    id = 0; slice = 0; feature = 0; imgAnnotation = 0;
    cubeAnnotation = 0; nEdges = 0;
    //imageFilename = "";
    TPs = 0; FPs = 0; FNs = 0; count = 0;
    nodeCoeffs = 0; edgeCoeffs = 0; }

} SPATTERN;

typedef struct label {
  /* this defines the y-part (the label) of a training example,
     e.g. the parse tree of the corresponding sentence. */
  //int add_your_variables_here;
  int nNodes; // total number of nodes
  labelType* nodeLabels;
  bool cachedNodeLabels;

  // structure used to score coefficients for each node.
  map<sidType, nodeCoeffType>* nodeCoeffs;

  label()
  {
    cachedNodeLabels = false;
  }

} LABEL;

typedef struct structmodel {
  double *w;          /* pointer to the learned weights */
  MODEL  *svm_model;  /* the learned SVM model */
  long   sizePsi;     /* maximum number of weights in w */
  double walpha;
  /* other information that is needed for the stuctural model can be
     added here, e.g. the grammar rules for NLP parsing */
  double *wsum; // weighted sum (used for MIRA)
} STRUCTMODEL;



typedef struct struct_learn_parm {
  double epsilon;              /* precision for which to solve
				  quadratic program */
  double newconstretrain;      /* number of new constraints to
				  accumulate before recomputing the QP
				  solution (used in w=1 algorithm) */
  int    ccache_size;          /* maximum number of constraints to
				  cache for each example (used in w=4
				  algorithm) */
  double batch_size;           /* size of the mini batches in percent
				  of training set size (used in w=4
				  algorithm) */
  double C;                    /* trade-off between margin and loss */
  char   custom_argv[50][300]; /* storage for the --* command line options */
  int    custom_argc;          /* number of --* command line options */
  int    slack_norm;           /* norm to use in objective function
                                  for slack variables; 1 -> L1-norm, 
				  2 -> L2-norm */
  int    loss_type;            /* selected loss type from -r
				  command line option. Select between
				  slack rescaling (1) and margin
				  rescaling (2) */
  int    loss_function;        /* select between different loss
				  functions via -l command line
				  option */
  /* further parameters that are passed to init_struct_model() */
  int alg_type;
  double* lossPerLabel;
  double lossScale;
  double maxC; // maximum value of C
  double startC;    // starting value of C
  int nClasses;
  int nDistances;
  int nUnaryWeights;
  int nGradientLevels;
  int nOrientations;
  int nLocalScales;
  int nScales;
  int nScalingCoefficients;
  bool includeLocalEdges;
  int nMaxIterations;
  int stepForOutputFiles;
  int stepForParameterFiles;
  int featureSize;
  map<labelType, ulong> labelToClassIdx;
  map<ulong, labelType> classIdxToLabel;
  char fileExtension[MAX_SIZE_EXT]; // 9 characters should be enough
  int giType;
  int ssvm_iteration;
  int iterationId;
  double sampling_temperature_0;
  double sampling_rate;
  int maxSqEdgeDistance;
  int metric_type;
} STRUCT_LEARN_PARM;

typedef struct struct_test_stats {
  /* you can add variables for keeping statistics when evaluating the
     test predictions in svm_struct_classify. This can be used in the
     function eval_prediction and print_struct_testing_stats. */
} STRUCT_TEST_STATS;

#endif
