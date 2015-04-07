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

#ifndef F_LoadFromFile_H
#define F_LoadFromFile_H

#include <vector>

// SliceMe
#include "Slice.h"
#include "Slice3d.h"
#include "Feature.h"

//-------------------------------------------------------------------------TYPES

// TODO : should really use char for 3d volumes !!!
//typedef char fileFeatureType;
typedef float fileFeatureType;

#define USE_SPARSE_STRUCTURE 0

#if USE_SPARSE_STRUCTURE
typedef map<ulong, fileFeatureType*> featureType;
#else
typedef fileFeatureType** featureType;
#endif

//-------------------------------------------------------------------------CLASS

class F_LoadFromFile : public Feature
{
 public:	

  F_LoadFromFile();

  ~F_LoadFromFile();

  void clearFeatures();

  string getAbsoluteFeaturePath(const string& featureFilename,
                                const string& inputDir);

  int getSizeFeatureVectorForOneSupernode();

  bool getFeatureVector(osvm_node *n,
                        const int x,
                        const int y);
  /**
   * Extract a feature vector for a given supernode in a 2d slice
   */
  bool getFeatureVectorForOneSupernode(osvm_node *x,
                                       Slice* slice,
                                       const int supernodeId);

  /**
   * Extract a feature vector for a given supernode in a 3d volume
   */
  bool getFeatureVectorForOneSupernode(osvm_node *x,
                                       Slice3d* slice3d,
                                       const int supernodeId);

  bool getFeatureVector(osvm_node *x,
                        Slice3d* slice3d,
                        const int gx,
                        const int gy,
                        const int gz);

  string getFeaturePath() { return featurePath; }

  eFeatureType getFeatureType() { return F_LOADFROMFILE; }

  /**
   * Load all the features in a given cube.
   */
  void init(Slice& slice,
            const char* filename);
  void init(Slice3d& slice3d,
            const char* filename);
  void init(Slice_P& slice,
            const char* filename);

  void init(Slice3d& slice3d, const char* filename,
            map<sidType, sidType>& sid_mapping);

  /**
   * Load all the features in a given sub-cube.
   */
  void init(Slice3d& slice3d, const char* filename,
            const node& start, const node& end);

  /**
   * Load features for given list of nodes.
   */
  void init(Slice3d& slice3d,
            const char* filename,
            std::vector<sidType>& lNodes);

  void loadFeatureFilenames(const char* filename,
                            vector<string>* lFeatureFilenames);

  /**
   * Load features from a text file that contain features
   * for each supervoxel.
   */
  void loadTextFeatures(Slice_P& slice,
                        const vector<string>& lFeatureFilenames);

  /**
   * Load features from a cube, e.g a TIF file contain features for
   * each voxel.
   */
  void loadVoxelBasedFeaturesFromTIF(Slice3d& slice,
                                     const vector<string>& lFeatureFilenames);

  void loadSupervoxelBasedFeaturesFromTIF(Slice3d& slice3d,
                                          const vector<string>& lFeatureFilenames,
                                          vector<sidType>& lNodes);

  void loadSupervoxelBasedFeaturesFromBinary(Slice3d& slice3d,
                                             const vector<string>& lFeatureFilenames,
                                             vector<sidType>& lNodes);

  void loadSupervoxelBasedFeaturesFromBinary(Slice3d& slice3d,
                                             const vector<string>& lFeatureFilenames,
                                             vector<sidType>& lNodes,
                                             map<sidType, sidType>& sid_mapping);

  void loadSupervoxelBasedFeaturesFromSetOfBinaries(Slice3d& slice3d,
                                                    const vector<string>& lFeatureFilenames,
                                                    vector<sidType>& lNodes);

  featureType* getMutableFeatures() { return &features; }

  void rescale();

  void rescale(Slice_P* slice);

  void setFeatures(featureType _features, int _nFeatures, int _featureSize) {
    features = _features;
    nFeatures = _nFeatures;
    featureSize = _featureSize;
  }

  void setFeatureSize(const int _featureSize) { featureSize = _featureSize; }

private:
  featureType features;
  int nFeatures; // number of feature vectors
  int featureSize; // size of each feature vector
  string featurePath;
  bool initialized;
};

#endif // F_LoadFromFile_H
