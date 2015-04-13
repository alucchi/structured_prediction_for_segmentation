#ifndef __GETOPT_H__
#define __GETOPT_H__

#ifdef __cplusplus
extern "C" {
#endif

extern int opterr;		/* if error message should be printed */
extern int optind;		/* index into parent argv vector */
extern int optopt;		/* character checked for validity */
extern int optreset;		/* reset getopt */
extern char *optarg;		/* argument associated with option */

struct option
{
  const char *name;
  int has_arg;
  int *flag;
  int val;
};

#define no_argument       0
#define required_argument 1
#define optional_argument 2

#ifdef __cplusplus
int getopt(int, char* const*, const char*);
#else
int getopt(int, char**, char*);
#endif
int getopt_long(int, char**, char*, struct option*, int*);

#ifdef __cplusplus
}
#endif

#endif /* __GETOPT_H__ */
