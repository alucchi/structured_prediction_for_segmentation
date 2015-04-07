////////////////////////////////////////////////////////////////////////
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

#include "svm_light/svm_common.h"
#include "svm_light/svm_learn.h"

#include "svm_struct/svm_struct_learn.h"
#include "svm_struct/svm_struct_common.h"
#include "svm_struct_api.h"
#include "graphInference.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>

// SliceMe
#include "globals.h"

void   init_parameters(long *, long *,
                       STRUCT_LEARN_PARM *, LEARN_PARM *, KERNEL_PARM *);

void train(char* config);

int main (int argc, char* argv[])
{
  if(argc < 2){
    fprintf(stderr, "Missing configuration file argument\n");
    fprintf(stderr, "Usage: train configuration_file\n");
    return EXIT_FAILURE;
  }
  train(argv[1]);
  return EXIT_SUCCESS;
}

void train(char* config)
{
  SAMPLE sample;  /* training sample */
  LEARN_PARM learn_parm;
  KERNEL_PARM kernel_parm;
  STRUCT_LEARN_PARM struct_parm;
  STRUCTMODEL structmodel;

  init_parameters(&verbosity,
                  &struct_verbosity,&struct_parm,&learn_parm,
                  &kernel_parm);

  if(struct_verbosity>=1) {
    //verbose = true;
    printf("[main] Reading training examples...\n"); fflush(stdout);
  }

  /* read the training examples */
  sample = read_struct_examples(config, &struct_parm);
  if(struct_verbosity>=1) {
    printf("done\n"); fflush(stdout);
  }

  clock_t time_0 = clock();
  
  /* Do the learning and return structmodel. */
  if(struct_parm.alg_type == 0)
    svm_learn_struct(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,NSLACK_ALG);
  else if(struct_parm.alg_type == 1)
    svm_learn_struct(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,NSLACK_SHRINK_ALG);
  else if(struct_parm.alg_type == 2)
    svm_learn_struct_joint(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,ONESLACK_PRIMAL_ALG);
  else if(struct_parm.alg_type == 3)
    svm_learn_struct_joint(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,ONESLACK_DUAL_ALG);
  else if(struct_parm.alg_type == 4)
    svm_learn_struct_joint(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel,ONESLACK_DUAL_CACHE_ALG);
  else if(struct_parm.alg_type == 9)
    svm_learn_struct_joint_custom(sample,&struct_parm,&learn_parm,&kernel_parm,&structmodel);
  else
    exit(EXIT_FAILURE);

  clock_t t = clock() - time_0;
  printf("[main] Running time = %ld clocks = %f s\n", t, t/(float)CLOCKS_PER_SEC); fflush(stdout);


  /* Warning: The model contains references to the original data 'docs'.
     If you want to free the original data, and only keep the model, you 
     have to make a deep copy of 'model'. */
  if(struct_verbosity>=1) {
    printf("[main] Writing learned model...");fflush(stdout);
  }

  char modelfile[200];           /* file for resulting classifier */
  strcpy (modelfile, "svm_struct_model");
  write_struct_model(modelfile,&structmodel,&struct_parm);
  if(struct_verbosity>=1) {
    printf("done\n");fflush(stdout);
  }

  finalize();

  free_struct_sample(sample);
  free_struct_model(structmodel);

  svm_struct_learn_api_exit();
}

/*---------------------------------------------------------------------------*/

void init_parameters(long *verbosity,long *struct_verbosity, 
                     STRUCT_LEARN_PARM *struct_parm,
                     LEARN_PARM *learn_parm, KERNEL_PARM *kernel_parm)
{
  /* set default */
  struct_parm->C=0.01;
  struct_parm->maxC=1.0;
  struct_parm->giType=T_GI_MULTIOBJ;
  struct_parm->slack_norm=1;
  struct_parm->epsilon=DEFAULT_EPS;
  struct_parm->custom_argc=0;
  struct_parm->loss_function=DEFAULT_LOSS_FCT;
  struct_parm->loss_type=DEFAULT_RESCALING;
  struct_parm->newconstretrain=100;
  struct_parm->ccache_size=5;
  struct_parm->batch_size=100;

  strcpy (learn_parm->predfile, "trans_predictions");
  strcpy (learn_parm->alphafile, "");
  (*verbosity)=0;/*verbosity for svm_light*/
  (*struct_verbosity)=1; /*verbosity for struct learning portion*/
  learn_parm->biased_hyperplane=1;
  learn_parm->remove_inconsistent=0;
  learn_parm->skip_final_opt_check=0;
  learn_parm->svm_maxqpsize=10;
  learn_parm->svm_newvarsinqp=0;
  learn_parm->svm_iter_to_shrink=100; //-9999;
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

  parse_struct_parameters(struct_parm);
}
