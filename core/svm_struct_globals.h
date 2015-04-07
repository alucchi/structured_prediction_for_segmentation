#include "globalsE.h" // VOC...

// From 1 to 5
// 1=output ground truth and predicted labels
// 2=output inference images
// 3=output groundtruth images + images for most violated constraints at each iteration
// 4=output psiGT and psiMVC files
#define VERBOSITY 3

#define SSVM_VERBOSE 1
#define SSVM_PRINT(format, ...) if(SSVM_VERBOSE) printf (format, ## __VA_ARGS__)

// comment following line to disable openmp
#define USE_OPENMP 1
#define NTHREADS 8

//---------------------------------------------------------------------FUNCTIONS

