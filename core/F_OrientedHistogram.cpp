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

#include "Config.h"
#include "F_OrientedHistogram.h"
#include "Histogram.h"
#include "utils.h"

F_OrientedHistogram::F_OrientedHistogram(int _nOrientation)
{
  nOrientation = _nOrientation;
}

int F_OrientedHistogram::getSizeFeatureVector()
{
  return nOrientation*nOrientation;
}

bool F_OrientedHistogram::getFeatureVector(osvm_node *x, Slice* slice, int sid, int nsid)
{
  //getSizeFeatureVector() should return nOrientation
  return false; // not tested yet, see supervoxel method
}

bool F_OrientedHistogram::getFeatureVector(osvm_node *x, Slice3d* slice3d, int sid, int nsid)
{
  node cs;
  node cn;
  node ct;
  float v1[3];
  float v2[3];
  float v2p[3];

  // B=(v1,vn,vn2) basis used to project vectors
  float B[9];
  int angleXY;
  int angleXZ;
  int idx = 0;

  supernode* s = slice3d->getSupernode(sid);
  supernode* n = slice3d->getSupernode(nsid);
  s->getCenter(cs);
  n->getCenter(cn);

  // compute vector sn
  v1[0] = cs.x - cn.x;
  v1[1] = cs.y - cn.y;
  v1[2] = cs.z - cn.z;

  // Compute B basis
  B[0] = v1[0];
  B[1] = v1[1];
  B[2] = v1[2];
  B[3] = -v1[1];
  B[4] = v1[0];
  B[5] = 0;
  crossProduct(&B[0], &B[3], &B[6]);

  Histogram hist(nOrientation*nOrientation);
  float angleToIdx = (float)nOrientation/360.0f;

  vector < supernode* >* lNeighbors = &(s->neighbors);
  supernode* t;
  for(vector < supernode* >::iterator itN = lNeighbors->begin();
      itN != lNeighbors->end(); itN++)
    {
      //if(t->id == n->id)
      //  continue;

      t = *itN;
      t->getCenter(ct);

      // compute vector st
      v2[0] = cs.x - ct.x;
      v2[1] = cs.y - ct.y;
      v2[2] = cs.z - ct.z;

      // Project v2 in B
      matMulVec_3(B, v2, v2p);
      // compute angle in x-y plane
      angleXY = atan2(v2p[1], v2p[0])*180.0/PI;
      angleXY += 180.0;
      // compute angle in x-z plane
      angleXZ = atan2(v2p[2], v2p[0])*180.0/PI;
      angleXZ += 180.0;

      if(angleXY > 360.0 || angleXZ > 360.0)
        {
          printf("[F_OrientedHistogram] Error : angleXY=%d > 360.0 || angleXZ=%d > 360.0\n", angleXY, angleXZ);
          exit(-1);
        }

      idx = (int)(angleXY*angleToIdx)*nOrientation
        +(int)(angleXZ*angleToIdx);

      if(idx >= hist.nBins)
        idx = hist.nBins - 1;

      if(hist.histData[idx] < 0.1)
        hist.histData[idx] = slice3d->getAvgIntensity(t->id);
      else
        hist.histData[idx] = (hist.histData[idx]+slice3d->getAvgIntensity(t->id))/2;

      if(t->id == n->id)
        {
          printf("[F_OrientedHistogram] v2p=(%f,%f,%f)\n", v2p[0],v2p[1],v2p[2]);
          printf("[F_OrientedHistogram] angleXY=%d angleXZ=%d angleToIdx=%f idx=%d bin[idx]=%f\n",
                 angleXY,angleXZ,angleToIdx,idx, hist.histData[idx]);
        }
    }

  for(int i = 0;i < hist.nBins; i++)
    x[i].value = hist.histData[i];

  return true;
}
