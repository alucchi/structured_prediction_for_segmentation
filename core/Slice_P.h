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

#ifndef SLICE_P_H
#define SLICE_P_H

// standard libraries
#include <map>
#include <cmath>
#include <math.h>
#include <string>
#include <vector>

// SliceMe
#include "globalsE.h"
#include "Supernode.h"
#include "oSVM_types.h"

using namespace std;

//------------------------------------------------------------------------------

enum eSlicePType
  {
    SLICEP_SLICE = 0,
    SLICEP_SLICE3D,
    SLICEP_SLICEP
  };

class Feature;

//------------------------------------------------------------------------------

class Slice_P
{
 public:

  Slice_P();

  ~Slice_P();

  void addLongRangeEdges();

  void addLongRangeEdges_supernodeBased_generic(int nDistances);

  void addLongRangeEdges_supernodeBased(int nDistances);

  int angleToIdx(int angle) {
    int idx = 0;
    if(angle > 45 && angle < 135) {
      idx = 1; // on top
    } else {
      if(angle > 225 && angle < 315) {
          idx = 2; // bottom
      }
    }
    return idx;
  }

  int computeDistanceIdx(supernode* s, supernode* sn, int _nDistances);

  int computeGradientIdx(int sid1, int sid2, int nGradientLevels);

  int computeOrientationIdx(supernode* s, supernode* sn, int _nOrientations);

  virtual labelType computeSupernodeLabel(supernode* s, const uchar* _labels);

  void countSupernodeLabels(labelType* labels, map<labelType, ulong>& labelCount);
  void countSupernodeLabels(map<labelType, ulong>& labelCount);

  /**
   * Export image with mask as an overlay
   * @param filename is the name of the output file
   */
  virtual void exportOverlay(const char* filename) = 0;

  virtual void exportOverlay(const char* filename, labelType* labels) = 0;

  virtual void exportProbabilities(const char* filename, int nClasses,
                                   float* pbs) = 0;

  virtual void exportSupernodeLabels(const char* filename, int nClasses,
				     labelType* labels,
				     int nLabels,
				     const map<labelType, ulong>* labelToClassIdx) = 0;

  static ulong generateId();

  /**
   * Get centers of all the supernodes
   * Caller is responsible for destroying the dynamically allocated list
   * returned by this function
   */
  vector<node>* getCenters();

  /**
   * Return a pointer to the raw image data (pixels or voxels).
   */
  virtual uchar* getRawData() = 0;

  virtual int getIntensity(int x, int y, int z = 0) = 0;

  /**
   * Get average intensity of all the supernodes in the slice
   */
  float getAvgIntensity();

  /**
   * Get average intensity of a given supernode
   * Assume gray image was loaded
   */
  virtual float getAvgIntensity(sidType supernodeId) = 0;

  virtual float getAvgIntensity(int supernodeId, int& r, int &g, int &b) = 0;

  int getCubeness() { return cubeness; }

  ulong getId();

  void getMinMaxIntensity(int& minIntensity, int& maxIntensity);

  void getMinMaxGradient(int& minGradient, int& maxGradient);

  virtual string getName() = 0;

  virtual int getNbChannels() = 0;

  virtual ulong getNbEdges() = 0;
  virtual ulong getNbNodes() = 0;
  virtual ulong getNbSupernodes() = 0;

  virtual ulong getNbUndirectedEdges();

  int getNbLabels() { return nLabels; }

  int getSupernodeStep() { return supernode_step; }

  // TODO : this function should be generic.
  virtual probType getProb(int sid, int label, int scale = 0) = 0;

  /**
   * Returns cost = -log(p)
   */
  inline double getCost(int sid, int label, int scale = 0) {
    return -std::log(min(1.0,getProb(sid,label,scale)+DELTA_PB));
  }

  virtual sidType getSid(int x, int y, int z) = 0;

  inline ulong getDirectedEdgeId(sidType sid1, sidType sid2) {
    return sid1*getNbSupernodes() + sid2;
  }

  inline ulong getEdgeId(sidType sid1, sidType sid2) {
    sidType min_sid = min(sid1, sid2);
    sidType max_sid = max(sid1, sid2);
    return min_sid*getNbSupernodes() + max_sid;
  }

  /**
   * Returns cost = +log(p)
   */
  inline double getScore(int sid, int label, int scale = 0) {
    return std::log(min(1.0,getProb(sid,label,scale)+DELTA_PB));
  }

  inline int getDistanceIdx(int sid1, int sid2) {
    ulong edgeId = getEdgeId(sid1, sid2);
    return distanceIdxs[edgeId];
  }

  inline osvm_node* getFeature(sidType sid) {
    return features[sid];
  }

  inline int getFeatureSize() { return feature_size; }

#if USE_SPARSE_VECTORS
  inline int getFeatureSize(int id) { return feature_sizes[id]; }
#endif

  /**
   * This function assumes that the gradient was already computed with
   * precomputeGradientIndices
   */
  inline int getGradientIdx(int sid1, int sid2) {
    ulong edgeId = getEdgeId(sid1, sid2);
    return gradientIdxs[edgeId];
  }

  inline int getOrientationIdx(int sid1, int sid2) {
    if(orientationIdxs.size() == 0) {
      return 0;
    } else {
      ulong edgeId = getDirectedEdgeId(sid1, sid2);
      return orientationIdxs[edgeId];
    }
  }

  virtual sizeSliceType getWidth() = 0;
  virtual sizeSliceType getHeight() = 0;
  virtual sizeSliceType getDepth() { return 1; }
  virtual ulong getSize() = 0;

  virtual eSlicePType getType() {return SLICEP_SLICEP; } 

  virtual map<sidType, supernode* >* getMutableSupernodes() = 0;
  virtual const map<sidType, supernode* >& getSupernodes() = 0;

  virtual supernode* getSupernode(sidType sid) {
    map<sidType, supernode* >* _supernodes = getMutableSupernodes();
    return (*_supernodes)[sid];
  }

  labelType getSupernodeLabel(sidType sid);

  virtual void generateSupernodeLabels(const char* fn_annotation,
                                       bool includeBoundaryLabels,
                                       bool useColorImages) = 0;

  void getSupernodeLabelsFromBuffer(map<sidType, uchar>& supernodeLabels,
                                    const uchar* labels,
                                    int _nLabels,
                                    bool includeBoundaryLabels,
                                    bool useColorImages);

  virtual void generateSupernodeLabelsFromBuffer(const uchar* labels,
                                                 int _nLabels,
                                                 bool includeBoundaryLabels,
                                                 bool useColorImages);

  map<sidType, osvm_node*>* getPrecomputedFeatures() { return &features; }

  inline bool isFeatureComputed(sidType sid) {
    return (features.find(sid) != features.end());
  }

  bool loadFeatures(const char* filename, int* featureSize);

  // First, compute mean and variance of all precomputed features.
  // All the features get the mean subtracted and get divided by the variance.
  void rescalePrecomputedFeatures(const char* scale_filename = 0);

  void precomputeDistanceIndices(int _nDistances);

  void precomputeFeatures(Feature* feature);

  void precomputeGradientIndices(int _nGradientLevels);

  void precomputeOrientationIndices(int _nOrientations);

  void printDistanceIndicesCount(int nDistances);

  void printEdgeStats();

  void readSupernodeLabels_Multiscale(vector<string>& prediction_filenames,
                                      int nLabels);

  /**
   * Read supernode labels from a file following the LIBSVM file format,
   * i.e "label probability for each class".
   */
  int readSupernodeLabels(const char* prediction_filename,
                          vector<int>& labels,
                          int prob_offset = 0);

  /**
   * Read supernode labels from a file with numbers written with exponential format
   */
  int readSupernodeLabels_e(const char* prediction_filename,
                                   vector<int>& labels,
                                   const int prob_offset = 0);

  int readSupernodeLabels_float(const char* prediction_filename,
                                vector<int>& labels,
                                const int prob_offset = 0);

  int readSupernodeProbs(const char* prediction_filename,
                         vector<int>& labels,
                         int prob_offset);

  void setId(int _id) { id = _id; }

  void setNbLabels(const int _nLabels) {nLabels = _nLabels;}

 protected:

  int supernode_step;

  int cubeness;

  int id;

  /**
   * Number of labels
   */
  int nLabels;

  ulong nbEdges;

  ulong max_distance;

  int feature_size;

  float minPercentToAssignLabel;
  bool includeOtherLabel;

  // precomputed quantities for nodes
  map<sidType, osvm_node*> features;

#if USE_SPARSE_VECTORS
  map<sidType, int> feature_sizes;
#endif

  // precomputed quantities for edges
  map<ulong, int> gradientIdxs;
  map<ulong, int> orientationIdxs;
  map<ulong, int> distanceIdxs;

 public:
  string inputDir;

};

//------------------------------------------------------------------------------

#endif // SLICE_P_H
