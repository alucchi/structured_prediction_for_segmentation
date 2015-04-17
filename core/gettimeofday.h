#ifndef GETTIMEOFDAY_H_
#define GETTIMEOFDAY_H_
#define NOTIME 

#if 0
#define _WINSOCKAPI_    // stops windows.h including winsock.h
#include <winsock2.h>
#include <time.h>

struct timezone 
{
  int  tz_minuteswest; /* minutes W of Greenwich */
  int  tz_dsttime;     /* type of dst correction */
};

#ifdef __cplusplus
extern "C" {
#endif
int gettimeofday(struct timeval *tv, struct timezone *tz);
#ifdef __cplusplus
}
#endif
#endif
#endif //GETTIMEOFDAY_H_
