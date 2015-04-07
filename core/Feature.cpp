#include "Config.h"
#include "F_Bias.h"
#include "F_Combo.h"
#include "F_Dft.h"
#include "F_Gaussian.h"
#include "F_Glcm.h"
#include "F_ColorHistogram.h"
#include "F_Histogram.h"
#include "F_LoadFromFile.h"
#include "F_Position.h"
#include "oSVM.h"

#include <deque>
#include <map>
#include <set>

#ifdef USE_ITK
#include "F_Filter.h"
#include "F_GradientStats.h"
#endif

#ifdef USE_SIFT
#include "F_Sift.h"
#endif

map<ulong, map<ulong, Feature*> > Feature::feature_cache;

using namespace std;

struct search_item {
  supernode* s;
  int distance;
};

//------------------------------------------------------------------------------

Feature::Feature()
{
  mean = 0;
  variance = 0;
  includeNeighbors = true;
  string config_tmp;
  if(Config::Instance()->getParameter("include_neighbors", config_tmp)) {
    if(config_tmp.c_str()[0] == '1') {
      includeNeighbors = true;
    } else {
      includeNeighbors = false;
    }
  }
}

Feature::~Feature()
{
  if(mean != 0) {
    delete[] mean;
  }
  if(variance != 0) {
    delete[] variance;
  }
}

Feature* Feature::getFeature(Slice_P* slice_p,
                             vector<eFeatureType>& feature_types)
{
  switch(slice_p->getType())
    {
    case SLICEP_SLICE:
      {
        Slice* slice = static_cast<Slice*>(slice_p);
        return getFeature(slice, feature_types);
      }
      break;
      case SLICEP_SLICE3D:
        {
          Slice3d* slice3d = static_cast<Slice3d*>(slice_p);
          return getFeature(slice3d, feature_types);
        }
        break;
    default:
      break;
    }
  return 0;
}

Feature* Feature::getFeature(Slice* slice,
                             std::vector<eFeatureType>& feature_types)
{
  ulong featureId = getFeatureTypeId(feature_types);
  map<ulong, map<ulong, Feature*> >::iterator lookup = feature_cache.find(slice->getId());
  if(lookup != feature_cache.end()) {
    map<ulong, Feature*>::iterator lookup_f = lookup->second.find(featureId);
    if(lookup_f != lookup->second.end()) {
      PRINT_MESSAGE("[Feature] Re-use feature from cache : (%p,%ld)\n", slice, featureId);
      return lookup_f->second;
    }
  }

  Feature* feat = 0;
  if(feature_types.size() == 1) {
    // Function will insert feature in cache if necessary
    feat = Feature::getFeature(slice,feature_types[0]);
  } else {
    PRINT_MESSAGE("[Feature] Inserting combo feature in cache : (%p,%ld)\n", slice, featureId);
    feat = new F_Combo(feature_types, slice);
    feature_cache[slice->getId()][featureId] = feat;
  }
  return feat;
}

Feature* Feature::getFeature(Slice* slice,
                             eFeatureType feature_type)
{
  ulong featureId = (ulong)feature_type;
  map<ulong, map<ulong, Feature*> >::iterator lookup = feature_cache.find(slice->getId());
  if(lookup != feature_cache.end()) {
    map<ulong, Feature*>::iterator lookup_f = lookup->second.find(featureId);
    if(lookup_f != lookup->second.end()) {
      PRINT_MESSAGE("[Feature] Re-use feature from cache : (%p,%ld)\n", slice, featureId);
      return lookup_f->second;
    }
  }

  Feature* _feature = 0; // caller is responsible for deleting returned feature

  switch(feature_type)
    {
      // DO NOT FORGET to add a break when adding new features!

    case F_BIAS:
      _feature = new F_Bias;
      break;

    case F_GAUSSIAN:
      _feature = new F_Gaussian;
      break;

    case F_GLCM:
      _feature = new F_Glcm;
      break;

    case F_POSITION:
      _feature = new F_Position;
      break;      

    case F_LOADFROMFILE:
      {
        F_LoadFromFile* feat_File = new F_LoadFromFile;
        string featureFile;
        if(!Config::Instance()->getParameter("feature_file", featureFile)) {
          printf("[Feature] Error : path for feature file was not set\n");
          exit(-1);
        }

        if(isDirectory(featureFile)) {
          string dirName = getLastDirectoryFromPath(slice->getName());
          featureFile += dirName;
          PRINT_MESSAGE("[Feature] Loading features from directory %s\n", featureFile.c_str());
          feat_File->init(*slice, featureFile.c_str());
        } else {
          if(fileExists(featureFile)) {
            PRINT_MESSAGE("[Feature] Loading features from file %s\n", featureFile.c_str());
            // Load features for each supernode. Another option is to load features
            // for a given list of points :
            // feat->init(*slice,featureFile.c_str(),lSamples);
            feat_File->init(*slice, featureFile.c_str());
          } else {
            PRINT_MESSAGE("[Feature] WARNING : Feature file %s does not exist\n", featureFile.c_str());
          }
        }
        _feature = feat_File;
        break;
      }

#ifdef USE_ITK
    case F_FILTER:
      {
        _feature = new F_Filter(*((Slice_P*)slice));
        break;
      }
#endif

#ifdef USE_SIFT
    case F_SIFT:
      {
        int octaves = 3;
        string paramOctaves;
        if(Config::Instance()->getParameter("sift_octaves", paramOctaves)) {
          octaves = atoi(paramOctaves.c_str());
        }
        PRINT_MESSAGE("[Feature] sift_octaves=%d\n",octaves);

        int levels = 2;
        string paramLevels;
        if(Config::Instance()->getParameter("sift_levels", paramLevels)) {
          levels = atoi(paramLevels.c_str());
          PRINT_MESSAGE("[Feature] sift_levels=%d\n",levels);
        }

        int omin = 0;

        _feature = new F_Sift(slice->getName().c_str(),octaves,levels,omin);
        break;
      }
#endif

    case F_DFT:
      {
        _feature = new F_Dft(slice->getName().c_str(), 800, 600);
        break;
      }

    case F_HISTOGRAM:
      {
      string config_tmp;
      eHistogramType hist_type = INCLUDE_NEIGHBORS_IN_SEPARATE_BINS;
      if(Config::Instance()->getParameter("hist_type", config_tmp)) {
	hist_type = (eHistogramType)atoi(config_tmp.c_str());
      }
      printf("[Feature] Histogram: %d bins, hist_type = %d\n", NB_BINS, (int)hist_type);

      _feature = new F_Histogram(NB_BINS,
                                 MAX_PIXEL_VALUE,
                                 hist_type,
                                 true,
                                 slice->img);

      if(((F_Histogram*)_feature)->useColorImage) {
        slice->generateColorImage();
      }
      }
      break;

    case F_COLOR_HISTOGRAM:
      {
        int histogramType = INCLUDE_NEIGHBORS_PLUS_LOCATION;
        string config_tmp;
        if(Config::Instance()->getParameter("histogram_type", config_tmp)) {
          histogramType = atoi(config_tmp.c_str());
        }

        _feature = new F_ColorHistogram(NB_BINS,
                                        MAX_PIXEL_VALUE,
                                        (eHistogramType)histogramType,
                                        //INCLUDE_NEIGHBORS_IN_SEPARATE_BINS,
                                        //INCLUDE_NEIGHBORS_PLUS_LOCATION,
                                        true,
                                        slice->img);

      slice->generateColorImage();
      }
      break;

    default:
      printf("[Feature] Error in getFeature(Slice): unknown feature type %d\n", feature_type);
      exit(-1);
      break;
    }

  PRINT_MESSAGE("[Feature] Inserting feature of size %d in cache : (%ld,%p,%ld)\n",
                _feature->getSizeFeatureVector(), slice->getId(), slice, featureId);
  feature_cache[slice->getId()][featureId] = _feature;

  return _feature;
}

Feature* Feature::getFeature(Slice3d* slice3d,
                             std::vector<eFeatureType>& feature_types)
{
  ulong featureId = getFeatureTypeId(feature_types);
  map<ulong, map<ulong, Feature*> >::iterator lookup = feature_cache.find(slice3d->getId());
  if(lookup != feature_cache.end()) {
    map<ulong, Feature*>::iterator lookup_f = lookup->second.find(featureId);
    if(lookup_f != lookup->second.end()) {
      PRINT_MESSAGE("[Feature] Re-use feature from cache : (%p,%ld)\n", slice3d, featureId);
      return lookup_f->second;
    }
  }

  Feature* feat = 0;
  if(feature_types.size() == 1) {
    // Function will insert feature in cache if necessary
    feat = Feature::getFeature(slice3d,feature_types[0]);
  } else {
    feat = new F_Combo(feature_types, slice3d);
    PRINT_MESSAGE("[Feature] Inserting feature of size %d in cache : (%ld,%p,%ld)\n",
                  feat->getSizeFeatureVector(), slice3d->getId(), slice3d, featureId);
    feature_cache[slice3d->getId()][featureId] = feat;
  }
  return feat;
}

Feature* Feature::getFeature(Slice3d* slice3d,
                             eFeatureType feature_type)
{
  ulong featureId = (ulong)feature_type;
  map<ulong, map<ulong, Feature*> >::iterator lookup = feature_cache.find(slice3d->getId());
  if(lookup != feature_cache.end()) {
    map<ulong, Feature*>::iterator lookup_f = lookup->second.find(featureId);
    if(lookup_f != lookup->second.end()) {
      PRINT_MESSAGE("[Feature] Re-use feature from cache : (%p,%ld)\n", slice3d, featureId);
      return lookup_f->second;
    }
  }

  Feature* feat = 0; // caller is responsible for deleting returned feature
  switch(feature_type)
    {
    case F_BIAS:
      feat = new F_Bias;
      break;

    case F_HISTOGRAM:
      {
        string config_tmp;
        eHistogramType hist_type = INCLUDE_NEIGHBORS_IN_SEPARATE_BINS;
        if(Config::Instance()->getParameter("hist_type", config_tmp)) {
          hist_type = (eHistogramType)atoi(config_tmp.c_str());
        }
        printf("[Feature] Histogram: %d bins, hist_type = %d\n", NB_BINS, (int)hist_type);

        feat = new F_Histogram(NB_BINS, MAX_PIXEL_VALUE, hist_type, false);
      }
      break;

    case F_GAUSSIAN:
      feat = new F_Gaussian;
      break;

    case F_LOADFROMFILE:
      {
      F_LoadFromFile* feat_File = new F_LoadFromFile;

      string featureFile;
      if(!Config::Instance()->getParameter("feature_file", featureFile)) {
        printf("[Main] Error : path for feature file was not set\n");
        exit(-1);
      }
      PRINT_MESSAGE("[Feature] Loading featureFile = %s\n", featureFile.c_str());
      // Load features for each supernode. Another option is to load features
      // for a given list of points :
      // feat->init(*slice3d,featureFile.c_str(),lSamples);
      feat_File->init(*slice3d,featureFile.c_str());
      feat = feat_File;
      }
      break;

    case F_GLCM:
      feat = new F_Glcm;
      break;

#ifdef USE_ITK

    case F_FILTER:
      {
        feat = new F_Filter(*((Slice_P*)slice3d));
        break;
      }

    case F_GRADIENTSTATS:
      {
        feat = new F_GradientStats(*slice3d, GRADIENT_MAG);
        break;
      }

#endif

    default:
      printf("[Feature] Error in getFeature(Slice3d): unknown feature type %d\n", feature_type);
      exit(-1);
      break;
    }

  PRINT_MESSAGE("[Feature] Inserting feature of size %d in cache: (%ld,%p,%ld)\n",
                feat->getSizeFeatureVector(), slice3d->getId(), slice3d, featureId);
  feature_cache[slice3d->getId()][featureId] = feat;
  return feat;
}

bool Feature::getFeatureVectorForOneSupernode(osvm_node *x,
                                              Slice_P* slice_p,
                                              const int supernodeId)
{
  bool ok_flag = false;
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
        Slice3d* slice3d = (Slice3d*)slice_p;
        ok_flag = getFeatureVectorForOneSupernode(x, slice3d, supernodeId);
      }
      break;
    default:
      break;
    }
  return ok_flag;
}

void Feature::initSVMNode(osvm_node*& x, int d)
{
  x = new osvm_node[d+1];
  int i = 0;
  for(i = 0; i < d; i++)
    x[i].index = i+1;
  x[i].index = -1;
}

void Feature::saveCube(Slice_P& slice, const char* filename, int feature_dimension)
{
  int fvSize = getSizeFeatureVector();
  PRINT_MESSAGE("[Feature] Saving feature dimension %d to %s\n", feature_dimension, filename);
  ulong idx;
  node n;
  osvm_node* fn;
  initSVMNode(fn, fvSize);

  ulong sizeSlice = slice.getHeight() * slice.getWidth();
  ulong size = slice.getSize();
  ulong width = slice.getWidth();
  uchar* data = new uchar[size];

  const map<sidType, supernode* >& _supernodes = slice.getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    getFeatureVector(fn, &slice, it->first);

    nodeIterator ni = it->second->getIterator();
    ni.goToBegin();
    while(!ni.isAtEnd()) {
      ni.get(n);
      ni.next();
      idx = n.z*sizeSlice + n.y*width + n.x;
      data[idx] = fn[feature_dimension].value;
    }
  }

#if USE_ITK
  exportTIFCube(data, filename, slice.getDepth(), slice.getHeight(), slice.getWidth());
#endif

  delete[] data;
}

void Feature::save(Slice_P& slice, const char* filename)
{
  int fvSize = getSizeFeatureVector();
  osvm_node* n;
  initSVMNode(n, fvSize);

  ofstream ofs(filename);
  const map<sidType, supernode* >& _supernodes = slice.getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    ofs << (int) it->second->getLabel();
    getFeatureVector(n, &slice, it->first);
    for(int i = 0; i < fvSize; ++i) {
      ofs << " " << i+1 << FEATURE_FIELD_SEPARATOR << n[i].value;
    }
    ofs << " " << endl;
  }
  ofs.close();
  delete[] n;
}

void Feature::saveSparse(Slice_P& slice, const char* filename)
{
  int fvSize = getSizeFeatureVector();
  osvm_node* n;
  initSVMNode(n, fvSize);

  ofstream ofs(filename);
  const map<sidType, supernode* >& _supernodes = slice.getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    ofs << (int) it->second->getLabel();
    getFeatureVector(n, &slice, it->first);
    for(int i = 0; i < fvSize; ++i) {
      if(n[i].value != 0) {
        ofs << " " << i+1 << FEATURE_FIELD_SEPARATOR << n[i].value;
      }
    }
    ofs << " " << endl;
  }
  ofs.close();
  delete[] n;
}

void Feature::deleteFeature(Slice_P* slice_p, Feature* _feature)
{
  if(slice_p->getType() == SLICEP_SLICE) {
    delete _feature;
  } else {
    // do not delete features for slice3d
  }
}

void Feature::deleteFeature(Feature* _feature)
{
  // feature_cache is only for 3d slices and only contain a few features
  // so iterating over the cache should be fast.
  bool found_feature = false;
  for(map<ulong, map<ulong, Feature*> >::iterator it = feature_cache.begin();
      it != feature_cache.begin(); ++it) {
    for(map<ulong, Feature*>::iterator itFeat = it->second.begin();
        itFeat != it->second.end(); ++itFeat) {
      if(itFeat->second == _feature) {
        found_feature = true; 
        break;
      }
    }
    if(found_feature){
      break;
    }
  }
  if(!found_feature){
    delete _feature;
  }
}

void Feature::rescaleCache(Slice_P* slice)
{
  printf("[Feature] rescaling features in the cache: %ld\n", feature_cache.size());
  // feature_cache is only for 3d slices and only contain a few features
  // so iterating over the cache should be fast.
  for(map<ulong, map<ulong, Feature*> >::iterator it = feature_cache.begin();
      it != feature_cache.end(); ++it) {
    printf("[Feature] rescaling feature %ld in the cache\n", it->first);
    for(map<ulong, Feature*>::iterator itFeat = it->second.begin();
        itFeat != it->second.end(); ++itFeat) {
      printf("[Feature] rescaling feature %ld in the cache\n", itFeat->first);
      Feature* feat = itFeat->second;
      feat->rescale(slice);
    }
  }
}

bool Feature::getFeatureVectorGivenDistance(osvm_node *x, Slice_P* slice_p, int supernodeId)
{
  supernode* s = slice_p->getSupernode(supernodeId);

  int nDistances = DEFAULT_FEATURE_DISTANCE;
  int *countPerDistance = new int[nDistances];
  for(int c = 0; c < nDistances; ++c) {
    countPerDistance[c] = 0;
  }

  // add supernode s to a queue of supernodes to be visited
  deque<search_item> pending_list;
  map<sidType, int> distances;
  set<sidType> visited;
  distances[s->id] = 0;
  countPerDistance[0] = 1;
  search_item si;
  si.s = s;
  si.distance = 0;
  pending_list.push_front(si);

  // iterate over all supernodes in the queue
  while(!pending_list.empty()) {

    // retrieve supernode and associated distance from s
    search_item sqi = pending_list.front();
    supernode* sq = sqi.s;
    int dist = sqi.distance;

    // remove sq from the queue and mark as visited
    pending_list.pop_front();

    // add neighbors to the queue
    if(dist < (nDistances-1)) {
      // iterate over neighbors of sq and add them to the queue if there are not already in it
      for(vector<supernode*>::iterator itNq = sq->neighbors.begin();
          itNq != sq->neighbors.end(); ++itNq) {
        supernode* snq = *itNq;

        set<sidType>::iterator lookup_visited = visited.find(snq->id);
        if(lookup_visited == visited.end()) {
          search_item snqi;
          snqi.s = snq;
          snqi.distance = dist + 1;
          pending_list.push_back(snqi);
          visited.insert(sq->id);
          distances[snq->id] = dist + 1;
          //assert(distances[snq->id] < nDistances);
          ++countPerDistance[distances[snq->id]];
        }

      }
    }
  }

  int sizeFV_full = getSizeFeatureVector();
  for(int i = 0; i < sizeFV_full; ++i) {
    x[i].index = i+1; // probably don't need this
    x[i].value = 0;
  }

  // Go over all the supernodes stored in the list and get corresponding features
  int sizeFV = getSizeFeatureVectorForOneSupernode();
  osvm_node* xt = new osvm_node[sizeFV];
  for(map<sidType, int>::iterator it = distances.begin(); it != distances.end(); ++it) {
    sidType sidn = it->first;
    supernode* sn = slice_p->getSupernode(sidn);
    int distance = it->second;

    if(slice_p->isFeatureComputed(sidn)) {
      // Feature already exists. Copy relevant part of the vector.
      osvm_node* xn = slice_p->getFeature(sidn);
      for(int i = 0; i < sizeFV; ++i) {
        x[distance*sizeFV+i].value += xn[i].value;
      }
    } else {
      getFeatureVectorForOneSupernode(xt, slice_p, sn->id);
      for(int i = 0; i < sizeFV; ++i) {
        x[distance*sizeFV+i].value += xt[i].value;
      }
    }
  }
  delete[] xt;

  // normalize
  for(int d = 0; d < nDistances; ++d) {
    if(countPerDistance[d] != 0) {
      for(int i = 0;i < sizeFV; i++) {
        x[sizeFV*d + i].value /= countPerDistance[d];
      }
    }
  }

  /*
  // todo: also add l1 and l2 normalization
  if(normalize_l1_norm) {
    for(int d = 0; d < nDistances; ++d) {
      double count = 0;
      for(int i = 0;i < sizeFV; i++) {
        count += x[sizeFV*d + i].value;
      }
      if(count != 0) {
        for(int i = 0;i < sizeFV; i++) {
          x[sizeFV*d + i].value /= count;
        }
      }
    }
  }
  */

  delete [] countPerDistance;
  return true;
}

void Feature::updateFeatureStats(Slice_P* slice)
{
  PRINT_MESSAGE("[Feature] Update feature mean and variance\n");

  int fvSize = getSizeFeatureVector();
  int max_index = fvSize + 1;

  osvm_node* x;
  oSVM::initSVMNode(x, max_index);
  if(mean != 0) {
    delete[] mean;
  }
  oSVM::initSVMNode(mean, max_index);
  if(variance != 0) {
    delete[] variance;
  }
  oSVM::initSVMNode(variance, max_index);

  for(int i = 0;x[i].index != -1; i++) {
    mean[i].value = 0;
  }
  for(int i = 0;x[i].index != -1; i++) {
    variance[i].value = 0;
  }

  const map<sidType, supernode* >& _supernodes = slice->getSupernodes();
  double n = 1;
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {

    this->getFeatureVector(x, slice, it->first);

    // use running average to avoid overflow
    for(int i = 0;x[i].index != -1; i++) {
      mean[i].value = ((n-1.0)/n*mean[i].value) + x[i].value/n;
    }

    // New variance V(n) = (n-1*V(n-1) + (x(n)-M(n))^2)/(n)
    for(int i = 0;x[i].index != -1; i++) {
      variance[i].value = (((n-1.0)*variance[i].value) + (x[i].value-mean[i].value)*(x[i].value-mean[i].value))/n;
    }

    ++n;
  }

  delete[] x;
}

void Feature::precomputeFeatures(Slice_P* slice, Feature* feature, float**& output)
{
  int fvSize = feature->getSizeFeatureVector();
  int max_index = fvSize + 1;

  // allocate memory
  ulong nSupernodes = slice->getNbSupernodes();
  output = new float*[nSupernodes];
  for(ulong s = 0; s < nSupernodes; ++s) {
    output[s] = new float[fvSize];
  }

  osvm_node* n = new osvm_node[max_index];
  int i = 0;
  for(i = 0;i < max_index-1; i++)
    n[i].index = i+1;
  n[i].index = -1;

  const map<sidType, supernode* >& _supernodes = slice->getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {

    feature->getFeatureVector(n, slice, it->first);
    for(i = 0;i < max_index-1; i++) {
      output[it->first][i] = n[i].value;
    }
  }

  delete[] n;
}

void Feature::precomputeFeatures(Slice_P* slice, Feature* feature, map<ulong, float*> output)
{
  int fvSize = feature->getSizeFeatureVector();
  int max_index = fvSize + 1;

  osvm_node* n = new osvm_node[max_index];
  int i = 0;
  for(i = 0;i < max_index-1; i++)
    n[i].index = i+1;
  n[i].index = -1;

  const map<sidType, supernode* >& _supernodes = slice->getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {

    feature->getFeatureVector(n, slice, it->first);

    output[it->first] = new float[fvSize];
    for(i = 0;i < max_index-1; i++) {
      output[it->first][i] = n[i].value;
    }
  }
  delete[] n;
}
