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

#include <fstream>

// SliceMe
#include "Config.h"
#include "F_LoadFromFile.h"
#include "oSVM.h"
#include "utils.h"

//------------------------------------------------------------------------------

#define OUTPUT_FEATURES_TO_TXT_FILE 0

#define UPSIDE_DOWN_FEATURES 0

//------------------------------------------------------------------------------

F_LoadFromFile::F_LoadFromFile() {
  featureSize = 0;
  nFeatures = 0;
  featurePath = "";
  initialized = false;

#if UPSIDE_DOWN_FEATURES
  printf("[F_LoadFromFile] UPSIDE_DOWN_FEATURES is set to 1\n");
#endif
}

void F_LoadFromFile::loadFeatureFilenames(const char* filename,
                                          vector<string>* lFeatureFilenames)
{
  ifstream ifs(filename);
  if(ifs.fail()) {
    printf("[F_LoadFromFile] Error while loading %s\n", filename);
    exit(-1);
  }
  
  // Load path containing files
  getline(ifs, featurePath);

  // Load name of the files containing the features
  string fn;
  while(!ifs.eof()) {
    getline(ifs,fn);
    if(!fn.empty()) {
      lFeatureFilenames->push_back(fn);
    }
  }
  ifs.close();
}

string F_LoadFromFile::getAbsoluteFeaturePath(const string& featureFilename,
                                              const string& inputDir)
{
  string fullpath = featureFilename;
  string inputDir_rel = inputDir;
  if(inputDir_rel[inputDir_rel.length()-1] == '/')
    inputDir_rel = inputDir_rel.substr(0,inputDir_rel.length()-1);
  inputDir_rel = inputDir_rel.substr(inputDir_rel.rfind('/'));
  fullpath += inputDir_rel;
  fullpath += "/";
  return fullpath;
}

void F_LoadFromFile::init(Slice3d& slice3d, const char* filename,
                          vector<sidType>& lNodes)
{
  if(initialized) {
    printf("[F_LoadFromFile] Warning: Features were already loaded.\n");
    return;
  }

  vector<string> lFeatureFilenames;
  loadFeatureFilenames(filename, &lFeatureFilenames);

  if(lFeatureFilenames.size() == 0) {
    // no feature filename specified. Check if there is a file with a name
    // matching filename
    string featName = getNameFromPathWithoutExtension(slice3d.getName());
    string fullpath = featurePath + "/";
    fullpath += featName;
    printf("[F_LoadFromFile] Checking binary file %s\n", fullpath.c_str());
    if(fileExists(fullpath)) {
      lFeatureFilenames.push_back(featName);
      loadTextFeatures(slice3d, lFeatureFilenames);
    } else {
      string fullpath_with_ext = fullpath + ".bin";
      featName += ".bin";
      if(fileExists(fullpath_with_ext)) {
        printf("[F_LoadFromFile] Loading features from binary file %s\n", fullpath_with_ext.c_str());
        vector<string> lFiles;
        lFiles.push_back(featName);
        loadSupervoxelBasedFeaturesFromBinary(slice3d, lFiles, lNodes);
      } else {
        printf("[F_LoadFromFile] No features to be loaded in %s\n", filename);
      }
    }
  } else {
    string ext = getExtension(lFeatureFilenames[0]);
    if(ext == "tif") {
      PRINT_MESSAGE("[F_LoadFromFile] Loading features from TIF file %s\n", lFeatureFilenames[0].c_str());
      loadSupervoxelBasedFeaturesFromTIF(slice3d, lFeatureFilenames, lNodes);
    } else {
      if(ext == "bin") {
        printf("[F_LoadFromFile] Loading features from binary file %s\n", lFeatureFilenames[0].c_str());
        loadSupervoxelBasedFeaturesFromBinary(slice3d, lFeatureFilenames, lNodes);
      } else {
        loadTextFeatures(slice3d, lFeatureFilenames);
      }
    }
  }

  initialized = true;
}

void F_LoadFromFile::loadSupervoxelBasedFeaturesFromTIF(Slice3d& slice3d,
                                                        const vector<string>& lFeatureFilenames,
                                                        vector<sidType>& lNodes)
{
  featureSize = lFeatureFilenames.size();
  nFeatures = lNodes.size();
  PRINT_MESSAGE("[F_LoadFromFile] Allocating memory for %ld nodes and %d features\n",
                lNodes.size(), featureSize);
#if !USE_SPARSE_STRUCTURE
  features = new fileFeatureType*[nFeatures];
#endif
  for(vector<sidType>::iterator itNode = lNodes.begin();
      itNode != lNodes.end(); itNode++) {
    features[*itNode] = new fileFeatureType[featureSize];
  }

  // store features in memory
  fileFeatureType value;
  ulong idx;
  int fileId = 0;
  supernode* s;
  node center;
  for(vector<string>::const_iterator itFile = lFeatureFilenames.begin();
      itFile != lFeatureFilenames.end(); itFile++)
    {

#if OUTPUT_FEATURES_TO_TXT_FILE
      stringstream sout;
      sout << "features_" << fileId << ".txt";
      ofstream ofsFeat(sout.str().c_str());
#endif

      string fullpath = getAbsoluteFeaturePath(featurePath, slice3d.inputDir);
      fullpath += *itFile;
      PRINT_MESSAGE("[F_LoadFromFile] Loading %s\n", fullpath.c_str());
      Slice3d* inputCube = new Slice3d(fullpath.c_str());
      for(vector<sidType>::iterator itNode = lNodes.begin();
          itNode != lNodes.end(); itNode++)
        {
          s = slice3d.getSupernode(*itNode);
          s->getCenter(center);

#if UPSIDE_DOWN_FEATURES
          // Bug fix : Feture cubes from Yunpeng are upside-down
          center.y = slice3d.height - center.y;
#endif

          idx = slice3d.getIndex(center.x,center.y,center.z);
          // read data
          value = inputCube->raw_data[idx];
          features[*itNode][fileId] = (fileFeatureType)value;

#if OUTPUT_FEATURES_TO_TXT_FILE
          ofsFeat << (double) value << endl;
#endif          
        }
      delete inputCube;

#if OUTPUT_FEATURES_TO_TXT_FILE
      ofsFeat.close();
#endif

      fileId++;
    }
  PRINT_MESSAGE("[F_LoadFromFile] All feature files are now loaded in memory\n");
}

void F_LoadFromFile::loadSupervoxelBasedFeaturesFromSetOfBinaries(Slice3d& slice3d,
                                                                  const vector<string>& lFeatureFilenames,
                                                                  vector<sidType>& lNodes)
{
  featureSize = lFeatureFilenames.size();
  nFeatures = lNodes.size();
  PRINT_MESSAGE("[F_LoadFromFile] Allocating memory for %ld nodes and %d features\n",
                lNodes.size(), featureSize);
#if !USE_SPARSE_STRUCTURE
  features = new fileFeatureType*[nFeatures];
#endif
  for(vector<sidType>::iterator itNode = lNodes.begin();
      itNode != lNodes.end(); itNode++) {
    features[*itNode] = new fileFeatureType[featureSize];
  }

  // store features in memory
  fileFeatureType value;
  ulong idx;
  int fileId = 0;
  supernode* s;
  node center;
  for(vector<string>::const_iterator itFile = lFeatureFilenames.begin();
      itFile != lFeatureFilenames.end(); itFile++)
    {
      ifstream ifsF(itFile->c_str(), ios::binary);
      if(ifsF.fail()) {
        printf("[F_LoadFromFile] Failed to load %s\n", itFile->c_str());
        exit(-1);
      }
      for(vector<sidType>::iterator itNode = lNodes.begin();
          itNode != lNodes.end(); itNode++) {
        s = slice3d.getSupernode(*itNode);
        s->getCenter(center);
        idx = slice3d.getIndex(center.x,center.y,center.z);
        ifsF.seekg (idx, ios::beg);
        // read data
        ifsF.read((char*)&value,sizeof(fileFeatureType));
        features[*itNode][fileId] = value;
        //printf("f %ld %d %d %d %d\n",idx,(int)value,*itNode,fileId,(int)features[*itNode][fileId]);
      }
      ifsF.close();
      fileId++;
    }
  PRINT_MESSAGE("[F_LoadFromFile] All feature files are now loaded in memory\n");
}


void F_LoadFromFile::loadSupervoxelBasedFeaturesFromBinary(Slice3d& slice3d,
                                                           const vector<string>& lFeatureFilenames,
                                                           vector<sidType>& lNodes)
{
  uint featureSizePerFile = 0;
  string config_tmp;
  if(Config::Instance()->getParameter("featureSizePerFile", config_tmp)) {
    featureSizePerFile = atoi(config_tmp.c_str());
  }

  featureSize = featureSizePerFile;
  nFeatures = lNodes.size();
  PRINT_MESSAGE("[F_LoadFromFile] Allocating memory for %ld nodes and %d features\n",
                lNodes.size(), featureSize);
#if !USE_SPARSE_STRUCTURE
  features = new fileFeatureType*[nFeatures];
#endif
  for(vector<sidType>::iterator itNode = lNodes.begin();
      itNode != lNodes.end(); itNode++) {
    features[*itNode] = new fileFeatureType[featureSize];
  }

  for(vector<string>::const_iterator itFile = lFeatureFilenames.begin();
      itFile != lFeatureFilenames.end(); itFile++)
    {
      printf("[F_LoadFromFile] Loading %s\n", itFile->c_str());

      /*
       * File format :
       * <unsigned int numRows>
       * <unsigned int numCols>
       * <numRows * numCols float values (32-bit float data)>
       * The data is in column-major format, and each row belongs to a given supervoxel.
       */

      string fulllpath = featurePath + *itFile;
      ifstream ifsF(fulllpath.c_str(), ios::binary);
      if(ifsF.fail()) {
        printf("[F_LoadFromFile] Failed to load %s\n", itFile->c_str());
        exit(-1);
      }
      uint nRows;
      uint nCols;
      ifsF.read((char*)&nRows, sizeof(uint));
      ifsF.read((char*)&nCols, sizeof(uint));
      ulong n = nRows*nCols;
      printf("[F_LoadFromFile] nRows = %d, nCols = %d, nSupernodes = %ld, featureSize = %d. Need %f Mb\n",
             nRows, nCols, lNodes.size(), featureSize, n/(1024.0*1024.0));
      assert(nRows >= lNodes.size());
      assert(nCols >= (uint)featureSize);
      bool useSparseFeatures = lNodes.size() < nRows;
      float* _features = new float[n];
      ifsF.read((char*)_features, sizeof(float)*n);

#if 1
      if(useSparseFeatures) {
        for(int i = 0; i < featureSize; ++i) {
          ulong idx = i*nRows;
          for(vector<sidType>::iterator itNode = lNodes.begin();
              itNode != lNodes.end(); itNode++) {
            ulong fidx = idx + *itNode;
            if(fidx > n) {
              printf("f %d/%d, %d/%d, %ld %ld %ld\n", i, featureSize, *itNode, nRows, idx, fidx, n); 
              exit(-1);
            }
            assert(fidx < n);
            features[*itNode][i] = _features[fidx];
          }
        }
      } else {
        ulong fidx = 0;
        for(int j = 0; j < featureSize; ++j) {
          for(uint i = 0; i < nRows; ++i) {
            features[i][j] = _features[fidx];
            ++fidx;
          }
        }
      }
#else
      if(useSparseFeatures) {
        for(vector<sidType>::iterator itNode = lNodes.begin();
            itNode != lNodes.end(); itNode++) {
          ulong idx = (*itNode)*nCols;
          for(int i = 0; i < featureSize; ++i) {
            ulong fidx = idx + i;
            if(fidx > n) {
              printf("f %d/%d, %d/%d, %ld %ld %ld\n", i, featureSize, *itNode, nRows, idx, fidx, n); 
              exit(-1);
            }
            assert(fidx < n);
            features[*itNode][i] = _features[fidx];
          }
        }
      } else {
        ulong fidx = 0;
        for(uint i = 0; i < nRows; ++i) {
          for(int j = 0; j < featureSize; ++j) {
            features[i][j] = _features[fidx];
            ++fidx;
          }
        }
      }
#endif

      ifsF.close();
      delete[] _features;
    }

  PRINT_MESSAGE("[F_LoadFromFile] Rescaling features...\n");
  // rescale features
  string range_filename = slice3d.getName() + ".range";
  ofstream ofs_range(range_filename.c_str());
  for(int i = 0; i < featureSize; ++i) {
    fileFeatureType min_value = features[lNodes[0]][i];
    fileFeatureType max_value = features[lNodes[0]][i];
    for(int featIdx = 0; featIdx < nFeatures; ++featIdx) {
      fileFeatureType* _feat = features[featIdx];
      if(_feat[i] < min_value) {
        min_value = _feat[i];
      }
      if(_feat[i] > max_value) {
        max_value = _feat[i];
      }
    }
    PRINT_MESSAGE("[F_LoadFromFile] Feature %d : (min,max)=(%g,%g)\n", i,
                  min_value, max_value);

    ofs_range << min_value << " " << max_value << endl;

    fileFeatureType df = max_value - min_value;
    for(int featIdx = 0; featIdx < nFeatures; ++featIdx) {
      fileFeatureType* _feat = features[featIdx];
      _feat[i] = (_feat[i]-min_value)/df;
    }
  }
  ofs_range.close();

  PRINT_MESSAGE("[F_LoadFromFile] All feature files are now loaded in memory\n");
}

void F_LoadFromFile::init(Slice3d& slice3d, const char* filename,
                          map<sidType, sidType>& sid_mapping)
{
  if(initialized) {
    printf("[F_LoadFromFile] Warning: Features were already loaded.\n");
    return;
  }

  vector<string> lFeatureFilenames;
  loadFeatureFilenames(filename, &lFeatureFilenames);

  if(lFeatureFilenames.size() == 0) {
    // no feature filename specified. Check if there is a file with a name
    // matching filename
    string featName = getNameFromPathWithoutExtension(slice3d.getName());
    featName += ".bin";
    string fullpath = featurePath + "/";
    fullpath += featName;
    printf("[F_LoadFromFile] Checking binary file %s\n", fullpath.c_str());
    if(fileExists(fullpath)) {
      printf("[F_LoadFromFile] Loading features from binary file %s\n", fullpath.c_str());
      vector<string> lFiles;
      lFiles.push_back(featName);
      vector<sidType> lSamples;
      for(map<sidType, supernode* >::iterator it = slice3d.mSupervoxels->begin();
          it != slice3d.mSupervoxels->end(); it++) {
        lSamples.push_back(it->first);
      }
      loadSupervoxelBasedFeaturesFromBinary(slice3d, lFiles, lSamples, sid_mapping);
    } else {
      printf("[F_LoadFromFile] No features to be loaded in %s\n", filename);
    }
  } else {
    string ext = getExtension(lFeatureFilenames[0]);
    if(ext == "tif") {
      //PRINT_MESSAGE("[F_LoadFromFile] Loading feature cube %s\n", lFeatureFilenames[0].c_str());
      //loadVoxelBasedFeaturesFromTIF(slice3d, lFeatureFilenames);
      printf("[F_LoadFromFile] Not implemented yet\n");
      exit(-1);
    } else {
      if(ext == "bin") {
        PRINT_MESSAGE("[F_LoadFromFile] Loading features from binary file %s\n", lFeatureFilenames[0].c_str());
        vector<sidType> lSamples;
        for(map<sidType, supernode* >::iterator it = slice3d.mSupervoxels->begin();
            it != slice3d.mSupervoxels->end(); it++) {
          lSamples.push_back(it->first);
        }
        loadSupervoxelBasedFeaturesFromBinary(slice3d, lFeatureFilenames, lSamples, sid_mapping);
      } else {
        PRINT_MESSAGE("[F_LoadFromFile] Loading features from text file %s\n", lFeatureFilenames[0].c_str());
        loadTextFeatures(slice3d, lFeatureFilenames);
      }
    }
  }

  initialized = true;
}

void F_LoadFromFile::loadSupervoxelBasedFeaturesFromBinary(Slice3d& slice3d,
                                                           const vector<string>& lFeatureFilenames,
                                                           vector<sidType>& lNodes,
                                                           map<sidType, sidType>& sid_mapping)
{
  uint featureSizePerFile = 0;
  string config_tmp;
  if(Config::Instance()->getParameter("featureSizePerFile", config_tmp)) {
    featureSizePerFile = atoi(config_tmp.c_str());
  }

  featureSize = featureSizePerFile;
  nFeatures = lNodes.size();
  PRINT_MESSAGE("[F_LoadFromFile] Allocating memory for %ld nodes and %d features\n",
                lNodes.size(), featureSize);
#if !USE_SPARSE_STRUCTURE
  features = new fileFeatureType*[nFeatures];
#endif
  for(vector<sidType>::iterator itNode = lNodes.begin();
      itNode != lNodes.end(); itNode++) {
    features[*itNode] = new fileFeatureType[featureSize];
  }

  for(vector<string>::const_iterator itFile = lFeatureFilenames.begin();
      itFile != lFeatureFilenames.end(); itFile++)
    {
      printf("[F_LoadFromFile] Loading %s\n", itFile->c_str());

      /*
       * File format :
       * <unsigned int numRows>
       * <unsigned int numCols>
       * <numRows * numCols float values (32-bit float data)>
       * The data is in column-major format, and each row belongs to a given supervoxel.
       */

      string fulllpath = featurePath + *itFile;
      ifstream ifsF(fulllpath.c_str(), ios::binary);
      if(ifsF.fail()) {
        printf("[F_LoadFromFile] Failed to load %s\n", itFile->c_str());
        exit(-1);
      }
      uint nRows;
      uint nCols;
      ifsF.read((char*)&nRows, sizeof(uint));
      ifsF.read((char*)&nCols, sizeof(uint));
      ulong n = nRows*nCols;
      printf("[F_LoadFromFile] nRows = %d, nCols = %d, nSupernodes = %ld, featureSize = %d. Need %f Mb\n",
             nRows, nCols, lNodes.size(), featureSize, n/(1024.0*1024.0));
      assert(nRows >= lNodes.size());
      assert(nCols >= (uint)featureSize);
      bool useSparseFeatures = lNodes.size() < nRows;
      assert(useSparseFeatures);
      float* _features = new float[n];
      ifsF.read((char*)_features, sizeof(float)*n);

#if 1
      for(int i = 0; i < featureSize; ++i) {
        ulong idx = i*nRows;
        for(vector<sidType>::iterator itNode = lNodes.begin();
            itNode != lNodes.end(); itNode++) {
          ulong fidx = idx + sid_mapping[*itNode];
          if(fidx > n) {
            printf("f %d/%d, %d/%d, %ld %ld %ld\n", i, featureSize, *itNode, nRows, idx, fidx, n); 
            exit(-1);
          }
          assert(fidx < n);
          features[*itNode][i] = _features[fidx];
        }
      }
#else
      for(vector<sidType>::iterator itNode = lNodes.begin();
          itNode != lNodes.end(); itNode++) {
        ulong idx = sid_mapping[*itNode]*nCols;
        for(int i = 0; i < featureSize; ++i) {
          ulong fidx = idx + i;
          if(fidx > n) {
            printf("f %d/%d, %d/%d, %ld %ld %ld\n", i, featureSize, *itNode, nRows, idx, fidx, n); 
            exit(-1);
          }
          assert(fidx < n);
          features[*itNode][i] = _features[fidx];
        }
      }
#endif

      ifsF.close();
      delete[] _features;
    }

  PRINT_MESSAGE("[F_LoadFromFile] Rescaling features...\n");
  // rescale features
  string range_filename = slice3d.getName() + ".range";
  ofstream ofs_range(range_filename.c_str());
  for(int i = 0; i < featureSize; ++i) {
    fileFeatureType min_value = features[lNodes[0]][i];
    fileFeatureType max_value = features[lNodes[0]][i];
    for(int featIdx = 0; featIdx < nFeatures; ++featIdx) {
      fileFeatureType* _feat = features[featIdx];
      if(_feat[i] < min_value) {
        min_value = _feat[i];
      }
      if(_feat[i] > max_value) {
        max_value = _feat[i];
      }
    }
    PRINT_MESSAGE("[F_LoadFromFile] Feature %d : (min,max)=(%g,%g)\n", i,
                  min_value, max_value);

    ofs_range << min_value << " " << max_value << endl;

    fileFeatureType df = max_value - min_value;
    for(int featIdx = 0; featIdx < nFeatures; ++featIdx) {
      fileFeatureType* _feat = features[featIdx];
      _feat[i] = (_feat[i]-min_value)/df;
    }
  }
  ofs_range.close();

  PRINT_MESSAGE("[F_LoadFromFile] All feature files are now loaded in memory\n");
}

void F_LoadFromFile::init(Slice_P& slice_p, const char* filename)
{
  switch(slice_p.getType()) {
  case SLICEP_SLICE:
    {
      Slice* slice = static_cast<Slice*>(&slice_p);
      init(*slice, filename);
    }
    break;
  case SLICEP_SLICE3D:
    {
      Slice3d* slice3d = static_cast<Slice3d*>(&slice_p);
      init(*slice3d, filename);
    }
    break;
  default:
    break;
  }
}

void F_LoadFromFile::init(Slice3d& slice3d, const char* filename)
{
  vector<string> lFeatureFilenames;
  loadFeatureFilenames(filename, &lFeatureFilenames);

  if(lFeatureFilenames.size() == 0) {
    // no feature filename specified. Check if there is a file with a name
    // matching filename
    string featName = getNameFromPathWithoutExtension(slice3d.getName());
    featName += ".bin";
    string fullpath = featurePath + "/";
    fullpath += featName;
    printf("[F_LoadFromFile] Checking binary file %s\n", fullpath.c_str());
    if(fileExists(fullpath)) {
      printf("[F_LoadFromFile] Loading features from binary file %s\n", fullpath.c_str());
      vector<string> lFiles;
      lFiles.push_back(featName);
      vector<sidType> lSamples;
      for(map<sidType, supernode* >::iterator it = slice3d.mSupervoxels->begin();
          it != slice3d.mSupervoxels->end(); it++) {
        lSamples.push_back(it->first);
      }
      loadSupervoxelBasedFeaturesFromBinary(slice3d, lFiles, lSamples);
    } else {
      printf("[F_LoadFromFile] No features to be loaded in %s\n", filename);
    }
  } else {
    string ext = getExtension(lFeatureFilenames[0]);
    if(ext == "tif" || ext == "nrrd") {
      PRINT_MESSAGE("[F_LoadFromFile] Loading feature cube %s\n", lFeatureFilenames[0].c_str());
      //loadVoxelBasedFeaturesFromTIF(slice3d, lFeatureFilenames);
      vector<sidType> lSamples;
      for(map<sidType, supernode* >::iterator it = slice3d.mSupervoxels->begin();
          it != slice3d.mSupervoxels->end(); it++) {
        lSamples.push_back(it->first);
      }
      loadSupervoxelBasedFeaturesFromTIF(slice3d, lFeatureFilenames, lSamples);
    } else {
      if(ext == "bin") {
        PRINT_MESSAGE("[F_LoadFromFile] Loading features from binary file %s\n", lFeatureFilenames[0].c_str());
        vector<sidType> lSamples;
        for(map<sidType, supernode* >::iterator it = slice3d.mSupervoxels->begin();
            it != slice3d.mSupervoxels->end(); it++) {
          lSamples.push_back(it->first);
        }
        loadSupervoxelBasedFeaturesFromBinary(slice3d, lFeatureFilenames, lSamples);
      } else {
        PRINT_MESSAGE("[F_LoadFromFile] Loading features from text file %s\n", lFeatureFilenames[0].c_str());
        loadTextFeatures(slice3d, lFeatureFilenames);
      }
    }
  }
}

void F_LoadFromFile::loadVoxelBasedFeaturesFromTIF(Slice3d& slice3d,
                                                   const vector<string>& lFeatureFilenames)
{
  featureSize = lFeatureFilenames.size();
  ulong cubeSize = slice3d.getSize();
  nFeatures = cubeSize;
#if !USE_SPARSE_STRUCTURE
  features = new fileFeatureType*[nFeatures];
#endif
  for(ulong i = 0; i < cubeSize; ++i) {
    features[i] = new fileFeatureType[featureSize];
  }

  // store features in memory
  fileFeatureType value;
  ulong idx;
  int fileId = 0;
  node center;
  for(vector<string>::const_iterator itFile = lFeatureFilenames.begin();
      itFile != lFeatureFilenames.end(); itFile++)
    {
      //printf("[F_LoadFromFile] Loading %s\n", itFile->c_str());

#if OUTPUT_FEATURES_TO_TXT_FILE
      stringstream sout;
      sout << "features_" << fileId << ".txt";
      ofstream ofsFeat(sout.str().c_str());
#endif

      string fullpath = getAbsoluteFeaturePath(featurePath, slice3d.inputDir);
      fullpath += *itFile;
      PRINT_MESSAGE("[F_LoadFromFile] Loading %s\n", fullpath.c_str());
      Slice3d* inputCube = new Slice3d(fullpath.c_str());
      PRINT_MESSAGE("[F_LoadFromFile] Loading %s. Dimension = (%d,%d,%d)\n",
                    fullpath.c_str(), inputCube->getWidth(), inputCube->getHeight(), inputCube->getDepth());
      ulong featIdx = 0;
      for(int z = 0; z < inputCube->depth; ++z) {
        for(int y = 0; y < inputCube->height; ++y) {
          for(int x = 0; x < inputCube->width; ++x) {

#if UPSIDE_DOWN_FEATURES
            // Bug fix : Feture cubes from Yunpeng are upside-down
            idx = slice3d.getIndex(x, slice3d.height-y, z);
#else
            idx = slice3d.getIndex(x, y, z);
#endif

            // read data
            value = inputCube->raw_data[idx];
            features[featIdx][fileId] = (fileFeatureType)value;
            ++featIdx;
#if OUTPUT_FEATURES_TO_TXT_FILE
            ofsFeat << (double) value << endl;
#endif
          //printf("f %ld %d %d %d %d\n",idx,(int)value,*itNode,fileId,(int)features[*itNode][fileId]);
          }
        }
      }
      delete inputCube;

#if OUTPUT_FEATURES_TO_TXT_FILE
      ofsFeat.close();
#endif

      fileId++;
    }
}

void F_LoadFromFile::init(Slice3d& slice3d, const char* filename,
                          const node& start, const node& end)
{
  vector<string> lFeatureFilenames;
  loadFeatureFilenames(filename, &lFeatureFilenames);

  featureSize = lFeatureFilenames.size();
  ulong cubeSize = slice3d.getSize();
  nFeatures = cubeSize;
#if !USE_SPARSE_STRUCTURE
  features = new fileFeatureType*[nFeatures];
#endif
  for(ulong i = 0; i < cubeSize; ++i) {
    features[i] = new fileFeatureType[featureSize];
  }

  // store features in memory
  fileFeatureType value;
  ulong idx;
  int fileId = 0;
  node center;
  for(vector<string>::iterator itFile = lFeatureFilenames.begin();
      itFile != lFeatureFilenames.end(); itFile++)
    {
      //printf("[F_LoadFromFile] Loading %s\n", itFile->c_str());

#if OUTPUT_FEATURES_TO_TXT_FILE
      stringstream sout;
      sout << "features_" << fileId << ".txt";
      ofstream ofsFeat(sout.str().c_str());
#endif

      string fullpath = getAbsoluteFeaturePath(featurePath, slice3d.inputDir);
      fullpath += *itFile;
      PRINT_MESSAGE("[F_LoadFromFile] Loading %s\n", fullpath.c_str());
      Slice3d* inputCube = new Slice3d(fullpath.c_str());
      ulong featIdx = 0;
      for(int z = start.z; z < end.z; ++z) {
        for(int y = start.y;  y < end.y; ++y) {
          for(int x = start.x; x < end.x; ++x) {
#if UPSIDE_DOWN_FEATURES
            // Bug fix : Feture cubes from Yunpeng are upside-down
            idx = slice3d.getIndex(x, slice3d.height-y, z);
#else
            idx = slice3d.getIndex(x, y, z);
#endif
            // read data
            value = inputCube->raw_data[idx];
            features[featIdx][fileId] = (fileFeatureType)value;
            ++featIdx;
#if OUTPUT_FEATURES_TO_TXT_FILE
            ofsFeat << (double) value << endl;
#endif
          //printf("f %ld %d %d %d %d\n",idx,(int)value,*itNode,fileId,(int)features[*itNode][fileId]);
          }
        }
      }
      delete inputCube;

#if OUTPUT_FEATURES_TO_TXT_FILE
      ofsFeat.close();
#endif

      fileId++;
    }
}

F_LoadFromFile::~F_LoadFromFile()
{
  clearFeatures();
}

int F_LoadFromFile::getSizeFeatureVectorForOneSupernode()
{
  return featureSize;
}

bool F_LoadFromFile::getFeatureVector(osvm_node *n,
                                      const int x,
                                      const int y)
{
  return false;
}

bool F_LoadFromFile::getFeatureVectorForOneSupernode(osvm_node *x, Slice* slice,
                                                     int supernodeId)
{
  for(int i = 0; i < featureSize; i++) {
    x[i].value = (double)(features[supernodeId][i]);
  }
  return true;
}

bool F_LoadFromFile::getFeatureVectorForOneSupernode(osvm_node *x, Slice3d* slice3d,
                                                     int supernodeId)
{
  for(int i = 0; i < featureSize; i++) {
    x[i].value = (double)(features[supernodeId][i]);
  }
  return true;
}

bool F_LoadFromFile::getFeatureVector(osvm_node *x,
                                      Slice3d* slice3d,
                                      const int gx,
                                      const int gy,
                                      const int gz)
{
  int idx = slice3d->getIndex(gx,gy,gz);
  for(int i = 0; i < featureSize; i++) {
    x[i].value = features[idx][i];
  }
  return true;
}


void F_LoadFromFile::init(Slice& slice, const char* filename)
{
  vector<string> lFeatureFilenames;
 
  // TODO(al) : quick fix for the ECCV12 deadline...
  bool isDir = isDirectory(filename);
  if(isDir) {
    featurePath = string(filename);
  } else {
    loadFeatureFilenames(filename, &lFeatureFilenames);
  }

  if(lFeatureFilenames.size() == 0) {
    string fullpath = getNameFromPathWithoutExtension(slice.getName());
    fullpath += ".txt";
    lFeatureFilenames.push_back(fullpath);
  }

  loadTextFeatures(slice, lFeatureFilenames);
}

// Load text file
void F_LoadFromFile::loadTextFeatures(Slice_P& slice,
                                      const vector<string>& lFeatureFilenames)
{
  int label_offset = 1;
  // Open first file to get feature size
  string fullFeaturePath;
  string fullpath;
  switch(slice.getType())
    {
    case SLICEP_SLICE:
      {
        fullFeaturePath = featurePath + "/";
        fullpath = fullFeaturePath;
        fullpath += getNameFromPathWithoutExtension(slice.getName());
        fullpath += ".txt";
        break;
      }
    case SLICEP_SLICE3D:
      {
        fullFeaturePath = featurePath + "/";
        fullFeaturePath += slice.getName();
        fullFeaturePath += "/";
        fullpath = fullFeaturePath + lFeatureFilenames[0].c_str();
        break;
      }
    default:
      printf("[F_LoadFromFile] Unknown type of feature %d\n", (int)slice.getType());
      exit(-1);
      break;
  }

  printf("[F_LoadFromFile] featurePath = %s, slice.getName = %s, fullpath = %s\n", featurePath.c_str(), slice.getName().c_str(), fullpath.c_str());
  ifstream ifsF(fullpath.c_str());
  string tmp_line;
  if(ifsF.fail()) {
    printf("[F_LoadFromFile] Failed to load first text file %s\n", fullpath.c_str());
    exit(-1);
  }

  uint featureSizePerFile = 0;
  string config_tmp;
  if(Config::Instance()->getParameter("featureSizePerFile", config_tmp)) {
    featureSizePerFile = atoi(config_tmp.c_str());
  }

  vector<string> tmp_tokens;
  if(featureSizePerFile == 0) {
    // Estimate size by sampling lines in the first file
    // Read the first 100 lines. Feature vector size is assumed to be the maximum
    // vector size found among those 100 lines.
    const int nTries = 10;
    for(int t = 0; t < nTries; ++t) {
      if(!getline(ifsF, tmp_line)) {
        break;
      }
      tmp_tokens.clear();
      splitString(tmp_line, tmp_tokens);
      for(uint i = 0; i < tmp_tokens.size()-1; ++i) {
        string field = tmp_tokens[i+label_offset];
        size_t index_ch = field.find(FEATURE_FIELD_SEPARATOR);
        assert(index_ch != string::npos);
        string str_field_index = field.substr(0,index_ch);
        int field_index = atoi(str_field_index.c_str());
        if(featureSizePerFile < (uint)field_index) {
          featureSizePerFile = field_index;
        }
      }
    }
  }

  featureSize = (featureSizePerFile*lFeatureFilenames.size());
  ifsF.close();

  const map<sidType, supernode* >& _supernodes = slice.getSupernodes();
  nFeatures = _supernodes.size();
  PRINT_MESSAGE("[F_LoadFromFile] Allocating memory to store %ld vectors of size %d\n",
		_supernodes.size(), featureSize);
#if !USE_SPARSE_STRUCTURE
  features = new fileFeatureType*[nFeatures];
#endif
  for(uint i = 0; i < _supernodes.size(); ++i) {
    features[i] = new fileFeatureType[featureSizePerFile];
  }

  bool node_based_features = false;
  if(Config::Instance()->getParameter("node_based_features", config_tmp)) {
    if(config_tmp[0] == '1') {
      node_based_features = true;
    }
  }
  PRINT_MESSAGE("[F_LoadFromFile] node_based_features = %d\n", (int)node_based_features);

  ulong nNodes = slice.getNbNodes();
  fileFeatureType** node_features = 0;
  if(node_based_features) {
    node_features = new fileFeatureType*[nNodes];
    for(uint i = 0; i < nNodes; ++i) {
      node_features[i] = new fileFeatureType[featureSizePerFile];
    }

    label_offset = 0; // no labels
  }

  // store features in memory
  int fileIdx = 0;
  node center;
  for(vector<string>::const_iterator itFile = lFeatureFilenames.begin();
      itFile != lFeatureFilenames.end(); itFile++) {
    string fullpath = fullFeaturePath + *itFile;
    PRINT_MESSAGE("[F_LoadFromFile] Loading %s\n", fullpath.c_str());
    ifstream ifsF(fullpath.c_str());
    if(ifsF.fail()) {
      printf("[F_LoadFromFile] Failed to load text file %s\n", fullpath.c_str());
      exit(-1);
    }
    
    uint featIdx = fileIdx*featureSizePerFile;

    if(node_based_features) {

      string line;
      vector<string> tokens;
      ulong nodeId = 0;
      while(getline(ifsF, line)) {
        tokens.clear();
        splitString(line, tokens);
        // check that size of the feature vector is equal to the number of tokens-1
        // -1 is due to the ground truth label stored in the text file
        if(featureSizePerFile < (tokens.size()-1)) {
          printf("[F_LoadFromFile] featureSizePerFile < tokens.size()-1 : %d != %ld\n",
                 featureSizePerFile, tokens.size()-1);
          exit(-1);
        }
        for(uint i = 0; i < featureSizePerFile; ++i) {
          node_features[nodeId][i] = 0;
        }
        for(uint i = 0; i < tokens.size()-1; ++i) {
          string field = tokens[i+label_offset];
          size_t index_ch = field.find(FEATURE_FIELD_SEPARATOR);
          int field_index = i + FEATURE_FIRST_INDEX;
          float value = 0;
          if(index_ch != string::npos) {
            assert(0);
            string str_field_index = field.substr(0,index_ch);
            field_index = atoi(str_field_index.c_str());
            assert(field_index <= (int)featureSizePerFile);
            assert(field_index >= FEATURE_FIRST_INDEX);
            string str_field_value = field.substr(index_ch+1);
            value = atof(str_field_value.c_str());
          } else {
            value = atof(field.c_str());
          }
          node_features[nodeId][field_index - FEATURE_FIRST_INDEX] = value;
        }
        ++nodeId;
      }

      if(nodeId < nNodes) {
         printf("[F_LoadFromFile] Error: Feature file contains %ld lines but %ld are expected\n", nodeId, nNodes);
         exit(-1);
      }

      // aggregate values
      fileFeatureType *feature_values = new fileFeatureType[featureSizePerFile];
      const map<sidType, supernode* >& _supernodes = slice.getSupernodes();
      for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
          it != _supernodes.end(); it++)
        {

          for(uint i = 0; i < featureSizePerFile; ++i) {
            feature_values[i] = 0;
          }

          supernode* s = it->second;
          node n;
          nodeIterator ni = s->getIterator();
          ni.goToBegin();
          while(!ni.isAtEnd()) {
            ni.get(n);

            ulong nodeId = (n.z * slice.getDepth()) + (n.y*slice.getWidth()) + n.x;

            for(uint i = 0; i < featureSizePerFile; ++i) {
              feature_values[i] += node_features[nodeId][i];
            }

            ni.next();
          }

          for(uint i = 0; i < featureSizePerFile; ++i) {
            features[it->first][featIdx+i] = feature_values[i]/s->size();
          }
      }
      delete [] feature_values;
    } else {
      string line;
      vector<string> tokens;
      for(uint sid = 0; sid < _supernodes.size(); ++sid) {
        if(!getline(ifsF, line)) {
          printf("[F_LoadFromFile] Error : Feature file contains %d lines but %ld are expected\n", sid, _supernodes.size());
          exit(-1);
        }
        tokens.clear();
        splitString(line, tokens);
        // check that size of the feature vector is equal to the number of tokens-1
        // -1 is due to the ground truth label stored in the text file
        if(featureSizePerFile < (tokens.size()-1)) {
          printf("[F_LoadFromFile] featureSizePerFile < tokens.size()-1 : %d != %ld\n",
                 featureSizePerFile, tokens.size()-1);
          exit(-1);
        }
        //uint featIdx = fileIdx*featureSizePerFile;
        for(uint i = 0; i < featureSizePerFile; ++i) {
          features[sid][featIdx+i] = 0;
        }
        for(uint i = 0; i < tokens.size()-1; ++i) {
          string field = tokens[i+label_offset];
          size_t index_ch = field.find(FEATURE_FIELD_SEPARATOR);
          assert(index_ch != string::npos);
          string str_field_index = field.substr(0,index_ch);
          int field_index = atoi(str_field_index.c_str());
          assert(field_index <= (int)featureSizePerFile);
          assert(field_index >= FEATURE_FIRST_INDEX);
          string field_value = field.substr(index_ch+1);
          features[sid][featIdx+field_index - FEATURE_FIRST_INDEX] = atof(field_value.c_str());
        }
      }
    }

    ifsF.close();
    ++fileIdx;
  }

  if(node_based_features) {
    for(uint i = 0; i < nNodes; ++i) {
      delete[] node_features[i];
    }
    delete[] node_features;
  }

}

void F_LoadFromFile::clearFeatures()
{
  for(int featIdx = 0; featIdx < nFeatures; ++featIdx) {
    delete[] features[featIdx];
  }
#if !USE_SPARSE_STRUCTURE
  delete[] features;
#endif

  initialized = false;
}

//todo: delete this
void F_LoadFromFile::rescale(Slice_P* slice)
{
  printf("[F_LoadFromFile] rescaling features\n");
  for(int i = 0; i < nFeatures; ++i) {
    for(int j = 0; j < featureSize; ++j) {
      features[i][j] *= 1e3;
    }
  }
}

void F_LoadFromFile::rescale()
{
  int max_index = featureSize + 1;

  osvm_node* mean;
  oSVM::initSVMNode(mean, max_index);
  osvm_node* variance;
  oSVM::initSVMNode(variance, max_index);

  for(int i = 0; i < featureSize; ++i) {
    mean[i].value = 0;
  }
  for(int i = 0; i < featureSize; ++i) {
    variance[i].value = 0;
  }

  for(int n = 1; n <= nFeatures; ++n) {

    fileFeatureType* x = features[n-1];

    // use running average to avoid overflow
    for(int i = 0; i < featureSize; ++i) {
      mean[i].value = ((n-1.0)/n*mean[i].value) + x[i]/n;
    }

    // New variance V(n) = (n-1*V(n-1) + (x(n)-M(n))^2)/(n)
    for(int i = 0; i < featureSize; ++i) {
      variance[i].value = (((n-1.0)*variance[i].value) + (x[i]-mean[i].value)*(x[i]-mean[i].value))/n;
    }

  }

  printf("Mean:");
  for(int i = 0; i < featureSize; ++i) {
    printf("%g ", mean[i].value);
  }
  printf("\n");

  printf("Variance:");
  for(int i = 0; i < featureSize; ++i) {
    printf("%g ", variance[i].value);
  }
  printf("\n");

  // prevent division by 0
  for(int i = 0; i < featureSize; ++i) {
    if(variance[i].value == 0) {
      variance[i].value = 1.0;
    }
  }

  const int sid_to_print = 100;

  for(int n = 0; n < nFeatures; ++n) {

    fileFeatureType* x = features[n];

    if(n == sid_to_print) {
      printf("x:");
      for(int i = 0; i < featureSize; ++i) {
        printf("%g ", x[i]);
      }
      printf("\n");
    }

    for(int i = 0; i < featureSize; ++i) {
      x[i] -= mean[i].value;
      x[i] /= variance[i].value;
    }

    if(n == sid_to_print) {
      printf("x:");
      for(int i = 0; i < featureSize; ++i) {
        printf("%g ", x[i]);
      }
      printf("\n");
    }

  }

  delete[] mean;
  delete[] variance;
}
