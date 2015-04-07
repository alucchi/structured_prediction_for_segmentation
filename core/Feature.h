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

#ifndef FEATURE_H
#define FEATURE_H

#include "Slice_P.h"
#include "Slice.h"
#include "Slice3d.h"
#include "Slice_P.h"

#define FEATURE_FIELD_SEPARATOR ':'

//-----------------------------------------------------------------------------

class Feature
{
 public:	

  Feature();

  ~Feature();

  static void deleteFeature(Slice_P* slice_p, Feature* _feature);
  static void deleteFeature(Feature* _feature);

  inline int getSizeFeatureVector() {
    int fvSize = getSizeFeatureVectorForOneSupernode();
    if(includeNeighbors) {
      fvSize *= DEFAULT_FEATURE_DISTANCE;
    }
    return fvSize;
  }

  //-------------------------------------

  virtual int getSizeFeatureVectorForOneSupernode() { return 0; }

  virtual bool getFeatureVectorForOneSupernode(osvm_node *x,
                                               Slice_P* slice_p,
                                               const int supernodeId);

  virtual bool getFeatureVectorForOneSupernode(osvm_node *x,
                                               Slice3d* slice,
                                               const int supernodeId)
  { return false; }

  virtual bool getFeatureVectorForOneSupernode(osvm_node *x,
                                               Slice* slice,
                                               const int supernodeId)
  { return false; }

  virtual bool getFeatureVectorForOneSupernode(osvm_node *n,
                                               const int x,
                                               const int y)
  { return false; }

  virtual bool getFeatureVectorForOneSupernode(osvm_node *n,
                                               Slice_P* slice_p,
                                               const int x,
                                               const int y,
                                               const int z)
  { return false; }

  //-------------------------------------

  inline bool getFeatureVector(osvm_node *x,
                                Slice_P* slice_p,
                               const int supernodeId) {
    bool ok_flag = false;
    if(includeNeighbors) {
      ok_flag = getFeatureVectorGivenDistance(x, slice_p, supernodeId);
    } else {
      switch(slice_p->getType())
        {
        case SLICEP_SLICE:
          {
            Slice* slice = static_cast<Slice*>(slice_p);
            ok_flag = getFeatureVectorForOneSupernode(x, slice, supernodeId);
          }
          break;
        case SLICEP_SLICE3D:
          {
            Slice3d* slice3d = static_cast<Slice3d*>(slice_p);
            ok_flag = getFeatureVectorForOneSupernode(x, slice3d, supernodeId);
          }
          break;
        default:
          break;
        }
    }

    return ok_flag;
  }

  inline bool getFeatureVector(osvm_node *n,
                               Slice_P* slice_p,
                               const int x,
                               const int y)
  {
    sidType sid = slice_p->getSid(x, y, 0);
    return getFeatureVector(n, slice_p, sid);
  }

  inline bool getFeatureVector(osvm_node *n,
                               Slice_P* slice_p,
                               const int x,
                               const int y,
                               const int z)
  {
    sidType sid = slice_p->getSid(x, y, z);
    return getFeatureVector(n, slice_p, sid);
  }

  //-------------------------------------

  bool getFeatureVectorGivenDistance(osvm_node *x, Slice_P* slice, int supernodeId);

  static Feature* getFeature(Slice_P* slice_p,
                             std::vector<eFeatureType>& feature_types);

  static Feature* getFeature(Slice* slice,
                             eFeatureType feature_type);

  static Feature* getFeature(Slice* slice,
                             std::vector<eFeatureType>& feature_types);

  static Feature* getFeature(Slice3d* slice3d,
                             eFeatureType feature_type);

  static Feature* getFeature(Slice3d* slice3d,
                             std::vector<eFeatureType>& feature_types);

  virtual eFeatureType getFeatureType() { return F_UNKNOWN; }

  static void initSVMNode(osvm_node*& x, int d);

  static void precomputeFeatures(Slice_P* slice, Feature* feature, float**& output);

  static void precomputeFeatures(Slice_P* slice, Feature* feature, map<ulong, float*> output);

  static void rescaleCache(Slice_P* slice);

  virtual void rescale(Slice_P* slice) { ; }

  void save(Slice_P& slice, const char* filename);

  void saveCube(Slice_P& slice, const char* filename, int feature_dimension);

  void saveSparse(Slice_P& slice, const char* filename);

  void updateFeatureStats(Slice_P* slice);

 protected:

  osvm_node* mean;
  osvm_node* variance;

  // TODO(lucchi) : Delete features !
  static map<ulong, map<ulong, Feature*> > feature_cache;

 private:
  bool includeNeighbors;

};

#endif // FEATURE
