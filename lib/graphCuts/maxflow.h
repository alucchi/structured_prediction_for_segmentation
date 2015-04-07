/* maxflow.cpp */


#ifndef MAXFLOW_CPP
#define MAXFLOW_CPP

#include <stdio.h>
#include "kgraph.h"
//#include "instances.inc"


namespace maxflow
{

/*
	special constants for node->parent
*/
#define TERMINAL ( (arc *) 1 )		/* to terminal */
#define ORPHAN   ( (arc *) 2 )		/* orphan */


#define INFINITE_D ((int)(((unsigned)-1)/2))		/* infinite distance to the terminal */

/***********************************************************************/

#include "maxflow.tpp"

}
#endif
