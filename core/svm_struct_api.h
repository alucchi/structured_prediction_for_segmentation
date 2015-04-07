/***********************************************************************/
/*                                                                     */
/*   svm_struct_api.h                                                  */
/*                                                                     */
/*   Definition of API for attaching implementing SVM learning of      */
/*   structures (e.g. parsing, multi-label classification, HMM)        */ 
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

#include "svm_struct_api_types.h"
#include "svm_struct/svm_struct_common.h"

#ifndef svm_struct_api
#define svm_struct_api

#include "Config.h"
#include "energyParam.h"

void        svm_struct_learn_api_init(int argc, char* argv[]);
void        svm_struct_learn_api_exit();
void        svm_struct_classify_api_init(int argc, char* argv[]);
void        svm_struct_classify_api_exit();
SAMPLE      read_struct_examples(char *file, STRUCT_LEARN_PARM *sparm);
void        init_struct_model(SAMPLE sample, STRUCTMODEL *sm, 
			      STRUCT_LEARN_PARM *sparm, LEARN_PARM *lparm, 
			      KERNEL_PARM *kparm);
CONSTSET    init_struct_constraints(SAMPLE sample, STRUCTMODEL *sm, 
				    STRUCT_LEARN_PARM *sparm);
LABEL       find_most_violated_constraint_slackrescaling(SPATTERN x, LABEL y, 
                                                         const STRUCTMODEL *sm, 
                                                         const STRUCT_LEARN_PARM *sparm);
LABEL       find_most_violated_constraint_marginrescaling(SPATTERN x, LABEL y, 
                                                          const STRUCTMODEL *sm, 
                                                          const STRUCT_LEARN_PARM *sparm);
LABEL       classify_struct_example(SPATTERN x, STRUCTMODEL *sm, 
				    STRUCT_LEARN_PARM *sparm);
int         empty_label(LABEL y);
SVECTOR     *psi(SPATTERN x, LABEL y, const STRUCTMODEL *sm, 
                 const STRUCT_LEARN_PARM *sparm);
double      loss(LABEL y, LABEL ybar, const STRUCT_LEARN_PARM *sparm);
int         finalize_iteration(double ceps, int cached_constraint,
			       SAMPLE sample, STRUCTMODEL *sm,
			       CONSTSET cset, double *alpha, 
			       STRUCT_LEARN_PARM *sparm);
void        print_struct_learning_stats(SAMPLE sample, STRUCTMODEL *sm,
					CONSTSET cset, double *alpha, 
					STRUCT_LEARN_PARM *sparm);
void        print_struct_testing_stats(SAMPLE sample, STRUCTMODEL *sm,
				       STRUCT_LEARN_PARM *sparm,
				       STRUCT_TEST_STATS *teststats);
void        eval_prediction(long exnum, EXAMPLE ex, LABEL prediction, 
			    STRUCTMODEL *sm, STRUCT_LEARN_PARM *sparm,
			    STRUCT_TEST_STATS *teststats);
void        write_struct_model(char *file,STRUCTMODEL *sm, 
			       STRUCT_LEARN_PARM *sparm);
void        write_label(FILE *fp, LABEL y);
void        free_pattern(SPATTERN x);
void        free_label(LABEL y);
void        free_struct_model(STRUCTMODEL sm);
void        free_struct_sample(SAMPLE s);
void        print_struct_help();
void        parse_struct_parameters(STRUCT_LEARN_PARM *sparm);
void        print_struct_help_classify();
void        parse_struct_parameters_classify(STRUCT_LEARN_PARM *sparm);
void        svm_learn_struct_joint_custom(SAMPLE sample, 
				   STRUCT_LEARN_PARM *sparm,
				   LEARN_PARM *lparm, KERNEL_PARM *kparm, 
				   STRUCTMODEL *sm);

// Al

void computeLoss(LABEL& y, labelType* ybar, const STRUCT_LEARN_PARM *sparm,
                 double& loss, int& nDiff);

void computeLoss(labelType* y, labelType* ybar, int nNodes, const STRUCT_LEARN_PARM *sparm,
                 double& loss, int& nDiff);

void computeVOCLoss(LABEL& y, labelType* ybar,
                    int nClasses,
                    double& loss, int& nDiff, SPATTERN x);

void computeVOCLoss(LABEL& y, labelType* ybar,
                    int nClasses,
                    double& loss, int& nDiff);

void computeVOCLoss(LABEL& y, labelType* ybar,
                    int nClasses,
                    double& loss, int& nDiff,
                    ulong* TPs, ulong* FPs, ulong* FNs);

SWORD* computePsi(SWORD* words, SPATTERN x, LABEL y, const STRUCTMODEL *sm,
                 const STRUCT_LEARN_PARM *sparm,
                 double* _score = 0);

double computeScore(const STRUCTMODEL *sm, SWORD* fy);

void finalize();

void get_best_parameter_vector_id(int& best_idx);

void get_best_parameter_vector(EnergyParam* param);

void initLossFunction(EXAMPLE  *examples, const long nExamples,
                      double*& lossPerLabel,
                      int nClasses,
                      int loss_function,
                      STRUCT_LEARN_PARM *sparm = 0);

void initLossFunction_classBased(EXAMPLE  *examples, const long nExamples,
                                 double*& lossPerLabel,
                                 int nClasses,
                                 int loss_function,
                                 STRUCT_LEARN_PARM *sparm = 0);

void initLossFunction_nodeBased(EXAMPLE  *examples, const long nExamples,
                                double*& lossPerLabel,
                                int nClasses,
                                int loss_function,
                                STRUCT_LEARN_PARM *sparm = 0);

bool isSubmodular(const STRUCT_LEARN_PARM *sparm, double* smw);

void load_2d_dataset(string imageDir,
                     int nImages,
                     string maskDir,
                     STRUCT_LEARN_PARM *sparm,
                     EXAMPLE*& examples,
                     long* nExamples,
                     uint* maxNbNodes,
		     int* featureSize,
                     Config* config);

void load_3d_dataset(string imageDir,
		     string maskDir,
		     STRUCT_LEARN_PARM *sparm,
                     EXAMPLE*& examples,
                     long* nExamples,
                     uint* maxNbNodes,
		     int* featureSize,
                     Config* config);

void load_learn_parm(STRUCT_LEARN_PARM *sparm, Config* config);

void runInference(SPATTERN x, LABEL y, 
                  const STRUCTMODEL *sm, 
                  const STRUCT_LEARN_PARM *sparm,
                  LABEL& ybar, const int threadId, bool labelFound, int cacheId);

void update_output_dir(const int iteration);

void save_parameters(const char* output_filename, STRUCT_LEARN_PARM *sparm,
                     STRUCTMODEL *sm);

void sparmToEnergyParam(const STRUCT_LEARN_PARM& sparm, EnergyParam* param);

void energyParamToSparm(const EnergyParam& param, STRUCT_LEARN_PARM* sparm);

#endif
