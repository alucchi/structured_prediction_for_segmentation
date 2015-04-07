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

#include "Slice_P.h"
#include "Feature.h"
#include "Config.h"
#include "utils.h"
#include "globalsE.h"
#include "oSVM.h"

#include <fstream>
#include <deque>
#include <stdlib.h>

//------------------------------------------------------------------------------

ulong Slice_P::generateId()
{
  static ulong id = 0;
  return id++;
}

//------------------------------------------------------------------------------

Slice_P::Slice_P()
{
  max_distance = -1;
  id = Slice_P::generateId();
}

Slice_P::~Slice_P()
{
  for(map<sidType, osvm_node*>::iterator it = features.begin();
      it != features.end(); ++it) {
    delete[] it->second;
  }
}

ulong Slice_P::getId()
{
  if(id != -1) {
    return id;
  } else {
    //return static_cast<ulong>(this);
    assert(0);
  }
}

labelType Slice_P::getSupernodeLabel(sidType sid)
{
  supernode* s = getSupernode(sid);
  if(s) {
    return s->getLabel();
  } else {
    return -1;
  }
}

float Slice_P::getAvgIntensity()
{
  //if(avgIntensity == -1)
  float avgIntensity;
    {
      avgIntensity = 0;
      int n = 1;
      const map<sidType, supernode* >& _supernodes = getSupernodes();
      for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
          it != _supernodes.end(); it++) {
        // New mean M(n+1) = (n*M(n) + M(n+1))/(n+1)
        avgIntensity = ((avgIntensity*n) + getAvgIntensity(it->first))/(n+1.0);
        n++;
      }
    }
  return avgIntensity;
}

void Slice_P::getMinMaxIntensity(int& minIntensity, int& maxIntensity)
{
  int avgIntensity = 0;
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  minIntensity = _supernodes.begin()->first;
  maxIntensity = minIntensity;
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    avgIntensity = getAvgIntensity(it->first);
    if(minIntensity > avgIntensity) {
      minIntensity = avgIntensity;
    }
    if(maxIntensity < avgIntensity) {
      maxIntensity = avgIntensity;
    }
  }
}

void Slice_P::getMinMaxGradient(int& minGradient, int& maxGradient)
{
  int avgIntensity1 = 0;
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  minGradient = 0;
  maxGradient = 0;
  bool initialized = false;

  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    avgIntensity1 = getAvgIntensity(it->first);

    for(vector < supernode* >::iterator itN = it->second->neighbors.begin();
        itN != it->second->neighbors.end(); itN++) {
      float gradient = abs(avgIntensity1 - getAvgIntensity((*itN)->id));

      if(minGradient > gradient || (!initialized)) {
        minGradient = gradient;
      }
      if(maxGradient < gradient || (!initialized)) {
        maxGradient = gradient;
      }
      initialized = true;
    }
  }
}

void Slice_P::precomputeDistanceIndices(int _nDistances)
{
  /*
  // compute max distance
    if(max_distance == -1) {
    ulong w = getWidth();
    ulong h = getHeight();
    max_distance = w*w + h*h;
    printf("max_distance=%ld\n", max_distance);
    }
  */

  if(distanceIdxs.size() == 0) {
    const map<sidType, supernode* >& _supernodes = getSupernodes();
    for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
        it != _supernodes.end(); it++) {
      for(vector < supernode* >::iterator itN = it->second->neighbors.begin();
          itN != it->second->neighbors.end(); itN++) {

        if(it->first < (*itN)->id) {
          continue;
        }

        computeDistanceIdx(it->second, *itN, _nDistances);
      }
    }
  } else {
    printf("[Slice_P]::precomputeDistanceIndices Distance indices were already precomputed\n");
  }
}

int Slice_P::computeDistanceIdx(supernode* s, supernode* sn, int _nDistances)
{
  ulong edgeId = getEdgeId(s->id, sn->id);
  int distanceIdx = 0;

  if(distanceIdxs.find(edgeId) != distanceIdxs.end()) {
    distanceIdx = distanceIdxs[edgeId];
  } else {
    node cs;
    node cn;
    s->getCenter(cs);
    sn->getCenter(cn);
    double sq_distance = node_square_distance(cs, cn);
    double dist_ratio = sq_distance / MAX_SQ_DISTANCE_LONG_RANGE_EDGES;
    distanceIdx = (int)(dist_ratio*_nDistances+0.5);
    if(distanceIdx >= _nDistances) {
      distanceIdx = _nDistances-1;
    }
    distanceIdxs[edgeId] = distanceIdx;
    //printf("[Slice_P]::precomputeDistanceIndices (%d,%d) %d -> (%g,%d)\n",
    //       s->id, sn->id, _nDistances, sq_distance, distanceIdx);
  }

  return distanceIdx;
}

void Slice_P::precomputeGradientIndices(int _nGradientLevels)
{
  if(gradientIdxs.size() == 0) {
    const map<sidType, supernode* >& _supernodes = getSupernodes();

    /*
    // update MAX_INTENSITY_GRADIENT
    if(getNbChannels() == 3) {
      MAX_INTENSITY_GRADIENT = 0;

      for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
          it != _supernodes.end(); it++) {
        for(vector < supernode* >::iterator itN = it->second->neighbors.begin();
            itN != it->second->neighbors.end(); itN++) {

          if(it->first < (*itN)->id) {
            continue;
          }

          int r1,g1,b1;
          int r2,g2,b2;
          getAvgIntensity(it->first,r1,g1,b1);
          getAvgIntensity((*itN)->id,r2,g2,b2);

          float gradient = abs(r1-r2) + abs(g1-g2) + abs(b1-b2);
          if(gradient > MAX_INTENSITY_GRADIENT) {
            MAX_INTENSITY_GRADIENT = gradient;
          }
        }
      }
      printf("[Slice_P] MAX_INTENSITY_GRADIENT = %d\n", MAX_INTENSITY_GRADIENT);
    }
    */

    for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
        it != _supernodes.end(); it++) {
      for(vector < supernode* >::iterator itN = it->second->neighbors.begin();
          itN != it->second->neighbors.end(); itN++) {

        if(it->first < (*itN)->id) {
          continue;
        }

        computeGradientIdx(it->first, (*itN)->id, _nGradientLevels);
      }
    }
  } else {
    printf("[Slice_P]::precomputeGradientIndices Gradient indices were already precomputed\n");
  }
}

int Slice_P::computeGradientIdx(int sid1, int sid2, int nGradientLevels)
{
  ulong edgeId = getEdgeId(sid1, sid2);
  int gradientIdx = 0;

  if(gradientIdxs.find(edgeId) != gradientIdxs.end()) {
    gradientIdx = gradientIdxs[edgeId];
  } else {
    if(getNbChannels() == 3) {
      int r1,g1,b1;
      int r2,g2,b2;
      getAvgIntensity(sid1,r1,g1,b1);
      getAvgIntensity(sid2,r2,g2,b2);

      float gradient = abs(r1-r2) + abs(g1-g2) + abs(b1-b2);
      //gradient /= (3*MAX_INTENSITY);
      gradient /= MAX_INTENSITY_GRADIENT;
      gradientIdx = (int)(gradient*nGradientLevels+0.5);

      /*
      // 3d gradient : too many parameters...
      float dr = abs(r1-r2)/(float)MAX_INTENSITY;
      float dg = abs(g1-g2)/(float)MAX_INTENSITY;
      float db = abs(b1-b2)/(float)MAX_INTENSITY;
      gradientIdx = dr*nGradientLevels*nGradientLevels + dg*nGradientLevels + db;
      */

      if(gradientIdx >= nGradientLevels) {
        gradientIdx = nGradientLevels-1;
      }
    } else {
      float gradient = abs(getAvgIntensity(sid1) - getAvgIntensity(sid2));
      gradient /= MAX_INTENSITY_GRADIENT;
      gradientIdx = (int)(gradient*nGradientLevels+0.5);
      if(gradientIdx >= nGradientLevels) {
        gradientIdx = nGradientLevels-1;
      }
    }
    gradientIdxs[edgeId] = gradientIdx;
  }

  return gradientIdx;
}

void Slice_P::precomputeOrientationIndices(int _nOrientations)
{
  if(orientationIdxs.size() == 0 && _nOrientations > 1) {
    const map<sidType, supernode* >& _supernodes = getSupernodes();
    for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
        it != _supernodes.end(); it++) {
      for(vector < supernode* >::iterator itN = it->second->neighbors.begin();
          itN != it->second->neighbors.end(); itN++) {

        // process all edges
        //if(it->first < (*itN)->id) {
        //  continue;
        //}

        computeOrientationIdx(it->second, *itN, _nOrientations);
      }
    }
  } else {
    printf("[Slice_P]::precomputeOrientationIndices Orientation indices were already precomputed OR no orientation specified\n");
  }
}

int Slice_P::computeOrientationIdx(supernode* s, supernode* sn, int _nOrientations)
{
  ulong edgeId = getDirectedEdgeId(s->id, sn->id);
  int orientationIdx = 0;

  if(orientationIdxs.find(edgeId) != orientationIdxs.end()) {
    orientationIdx = orientationIdxs[edgeId];
  } else {
    node cs;
    node cn;
    float v[3];
    int angleXY;

    s->getCenter(cs);
    sn->getCenter(cn);
    // compute vector sn
    if(s->id < sn->id) {
      v[0] = cn.x - cs.x;
      v[1] = cn.y - cs.y;
      v[2] = cn.z - cs.z;
    } else {
      v[0] = cs.x - cn.x;
      v[1] = cs.y - cn.y;
      v[2] = cs.z - cn.z;
    }

    angleXY = atan2(v[1], v[0])*180.0/PI;
    angleXY += 180.0;
    orientationIdx = angleToIdx(angleXY);
    orientationIdxs[edgeId] = orientationIdx;
  }
  return orientationIdx;
}

void Slice_P::precomputeFeatures(Feature* feature)
{
  if(features.size() == 0) {

    int fvSize = feature->getSizeFeatureVector();
    int max_index = fvSize + 1;
    feature_size = fvSize;

#if USE_SPARSE_VECTORS
    int upper_index = fvSize;
    int nScales = 1;
    string config_tmp;
    if(Config::Instance()->getParameter("nScales", config_tmp)) {
      nScales = atoi(config_tmp.c_str());
    }
    printf("[Slice_P] Using sparse vectors with nScales=%d\n", nScales);
    upper_index -= nScales*21;
#endif

    const map<sidType, supernode* >& _supernodes = getSupernodes();
    printf("[Slice_P] precomputing features for %ld nodes\n", _supernodes.size());
    for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
        it != _supernodes.end(); it++) {
      //printf("-");
      osvm_node* n = new osvm_node[max_index];
      int i = 0;
      for(i = 0;i < max_index-1; i++)
        n[i].index = i+1;
      n[i].index = -1;

      feature->getFeatureVector(n, this, it->first);

#if USE_SPARSE_VECTORS
      osvm_node* n_sparse = 0;
      feature_sizes[it->first] = oSVM::create_sparse_vector(n, n_sparse, upper_index);
      delete[] n;
      features[it->first] = n_sparse;
#else
      features[it->first] = n;
#endif
    }
  } else {
    printf("[Slice_P]::precomputeFeatures : Features were already precomputed\n");
  }
  printf("-\n");
}

void Slice_P::rescalePrecomputedFeatures(const char* scale_filename)
{
  if(features.size() == 0) {
    printf("[Slice_P]::rescalePrecomputedFeatures: Features were not precomputed\n");
    return;
  }

  // get feature dimension
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  osvm_node* x1 = features[_supernodes.begin()->first];
  int fvSize = 0;
  for(int i = 0;x1[i].index != -1; i++) {
    ++fvSize;
  }
  int max_index = fvSize + 1;

  printf("[Slice_P] Rescaling features of dimension %d\n", fvSize);

  osvm_node* mean = 0; // E_x
  osvm_node* variance = 0;
  oSVM::initSVMNode(mean, fvSize);
  oSVM::initSVMNode(variance, fvSize);

  if(scale_filename != 0 && (fileExists(scale_filename))) {

    printf("[Slice_P] Loading scale of features from file %s\n", scale_filename);

    ifstream ifs(scale_filename);
    string line;
    vector<string> tokens;

    // mean
    getline(ifs, line);
    splitString(line, tokens);
    printf("tokens.size %ld\n", tokens.size());
    for(int i = 0; i < fvSize; ++i) {
      mean[i].value = atof(tokens[i].c_str());
    }

    // variance
    getline(ifs, line);
    tokens.clear();
    splitString(line, tokens);
    printf("tokens.size %ld\n", tokens.size());
    for(int i = 0; i < fvSize; ++i) {
      variance[i].value = atof(tokens[i].c_str());
    }

    ifs.close();
  } else {

    printf("[Slice_P] Compute mean and variance\n");

    osvm_node* E_x2 = 0;
    oSVM::initSVMNode(E_x2, fvSize);

    for(int i = 0; i < fvSize; ++i) {
      mean[i].value = 0;
    }

    for(int i = 0; i < fvSize; ++i) {
      E_x2[i].value = 0;
    }

    for(int i = 0; i < fvSize; ++i) {
      variance[i].value = 0;
    }

    // compute mean
    //double n = 1;
    for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
        it != _supernodes.end(); it++) {

      osvm_node* x = features[it->first];

      // use running average to avoid overflow
      for(int i = 0;x[i].index != -1; i++) {
        //mean[i].value = ((n-1.0)/n*mean[i].value) + (x[i].value/n);
        mean[i].value += x[i].value;
      }

      for(int i = 0;x[i].index != -1; i++) {
        //E_x2[i].value = ((n-1.0)/n*E_x2[i].value) + ((x[i].value*x[i].value)/n);
        E_x2[i].value += x[i].value*x[i].value;
      }

      //++n;
    }

    double n = (double) _supernodes.size();
    for(int i = 0; i < fvSize; ++i) {
      mean[i].value /= n;
      E_x2[i].value /= n;
    }

    // compute variance
    for(int i = 0; i < fvSize; ++i) {
      variance[i].value = E_x2[i].value - (mean[i].value*mean[i].value);
    }

    printf("[Slice_P] Output scale to %s\n", scale_filename);
    ofstream ofs(scale_filename);
    for(int i = 0; i < fvSize - 1; ++i) {
      ofs << mean[i].value << " ";
    }
    ofs << mean[fvSize - 1].value << endl;
    for(int i = 0; i < fvSize - 1; ++i) {
      ofs << variance[i].value << " ";
    }
    ofs << variance[fvSize - 1].value << endl;

    ofs.close();

  }

  printf("Mean:");
  for(int i = 0; mean[i].index != -1; i++) {
    printf("%d:%g ", mean[i].index, mean[i].value);
  }
  printf("\n");

  printf("Variance:");
  for(int i = 0; variance[i].index != -1; i++) {
    printf("%d:%g ", variance[i].index, variance[i].value);
  }
  printf("\n");

  // prevent division by 0
  for(int i = 0; variance[i].index != -1; i++) {
    if(variance[i].value == 0) {
      variance[i].value = 1.0;
    }
  }

  const int sid_to_print = 100;

  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {

    osvm_node* x = features[it->first];

    if(it->first == sid_to_print) {
      printf("x_100 (before rescaling):");
      for(int i = 0; x[i].index != -1; i++) {
        printf("%g ", x[i].value);
      }
      printf("\n");
    }

    for(int i = 0;x[i].index != -1; i++) {
      x[i].value -= mean[i].value;
      x[i].value /= sqrt(variance[i].value);
    }

    if(it->first == sid_to_print) {
      printf("x_100 (after rescaling):");
      for(int i = 0; x[i].index != -1; i++) {
        printf("%g ", x[i].value);
      }
      printf("\n");
    }

  }

  delete[] mean;
  delete[] variance;
}

bool Slice_P::loadFeatures(const char* filename, int* featureSize)
{
  ifstream ifsF(filename);
  if(ifsF.fail()) {
    printf("[Slice_P] Failed to load text file %s\n", filename);
    return false;
  }

  // check that the file contains the correct number of lines
  vector<string> tokens;
  ulong nLines = 0;
  string line;
  if(getline(ifsF, line)) {
    ++nLines;
    splitString(line, tokens);
    *featureSize = tokens.size() - 1;
    printf("[Slice_P] featureSize = %d\n", *featureSize);
    while(getline(ifsF, line)) {
      ++nLines;
    }
  }
  if(nLines != getNbSupernodes()) {
    printf("[Slice_P] Failed to load text file %s. Incorrect number of lines %ld. Expected %ld lines\n",
           filename, nLines, getNbSupernodes());
    return false;
  }
  ifsF.clear();
  ifsF.seekg(0, ios::beg);

#if USE_SPARSE_VECTORS
  int upper_index = *featureSize;
  int nScales = 1;
  string config_tmp;
  if(Config::Instance()->getParameter("nScales", config_tmp)) {
    nScales = atoi(config_tmp.c_str());
  }
  printf("[Slice_P] nScales=%d\n", nScales);
  upper_index -= nScales*21;
#endif

  const int label_offset = 1;
  int max_index = *featureSize + 1;
  ulong nodeId = 0;
  while(getline(ifsF, line)) {
    tokens.clear();
    splitString(line, tokens);

    osvm_node* n = new osvm_node[max_index];
    int i = 0;
    for(i = 0;i < max_index-1; i++)
      n[i].index = i+1;
    n[i].index = -1;

    for(uint i = 0; i < tokens.size()-1; ++i) {
      string field = tokens[i+label_offset];
      size_t index_ch = field.find(FEATURE_FIELD_SEPARATOR);
      int field_index = i + FEATURE_FIRST_INDEX;
      float value = 0;
      if(index_ch != string::npos) {
        string str_field_index = field.substr(0,index_ch);
        field_index = atoi(str_field_index.c_str());
        assert(field_index <= (int)*featureSize);
        assert(field_index >= FEATURE_FIRST_INDEX);
        string str_field_value = field.substr(index_ch+1);
        value = atof(str_field_value.c_str());
      } else {
        value = atof(field.c_str());
      }
      n[field_index - FEATURE_FIRST_INDEX].value = value;
    }

#if USE_SPARSE_VECTORS
    osvm_node* n_sparse = 0;
    feature_sizes[nodeId] = oSVM::create_sparse_vector(n, n_sparse, upper_index);
    delete[] n;
    features[nodeId] = n_sparse;
#else
    features[nodeId] = n;
#endif

    ++nodeId;

  }

  printf("[Slice_P] Loaded %ld features\n", nodeId);

  ifsF.close();
  return true;
}

vector<node>* Slice_P::getCenters()
{
  vector < node >* lCenters = new vector < node >;
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  lCenters->reserve(_supernodes.size());
  node center;
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    it->second->getCenter(center);
    lCenters->push_back(center);
  }
  return lCenters;
}

// this function is more generic and can add neighbors at any given distance
void Slice_P::addLongRangeEdges_supernodeBased(int nDistances)
{
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  ulong nSupernodes = _supernodes.size();

  //#ifdef WITH_OPENMP
  //#pragma omp parallel for
  //#endif
  for(ulong sid = 0; sid < nSupernodes; ++sid) {
    map<sidType, supernode* >::const_iterator it = _supernodes.find(sid);
    supernode* s = it->second;

    // add supernode s to a queue of supernodes to be visited
    //deque<search_item> pending_list;
    deque<supernode*> pending_list;
    map<sidType, int> distances;
    set<sidType> visited;
    distances[s->id] = -1;
    pending_list.push_front(s);

    for(vector<supernode*>::iterator itN = s->neighbors.begin();
        itN != s->neighbors.end(); ++itN) {
        // store distance to avoid online computation
      ulong edgeId = getEdgeId(s->id, (*itN)->id);
      //assert(distanceIdxs.find(edgeId) == distanceIdxs.end());
      distanceIdxs[edgeId] = 0;      
      distances[(*itN)->id] = 0;
    }

    // iterate over all supernodes in the queue
    while(!pending_list.empty()) {

      // retrieve supernode and associated distance from s
      supernode* sq = pending_list.front();
      int dist = distances[sq->id];

      // remove sq from the queue and mark as visited
      pending_list.pop_front();
      visited.insert(sq->id);

      // insert sq in the list of neighbors
      if(s->id != sq->id) {
        // check if retrieved supernode is already in the list of neighbors
        vector<supernode*>::iterator itN = find(s->neighbors.begin(), s->neighbors.end(), sq); // slow
        if(itN == s->neighbors.end()) {
          // did not find neighbors, add it
          s->neighbors.push_back(sq);

          // store distance to avoid online computation
          ulong edgeId = getEdgeId(s->id, sq->id);
          //assert(distanceIdxs.find(edgeId) == distanceIdxs.end());
          distanceIdxs[edgeId] = dist;
        }
        /*
        else {
          if(distanceIdxs.find(edgeId) == distanceIdxs.end()) {
            printf("[Slice_P] Error %d %d %ld %d\n", s->id, sq->id, edgeId, dist);
            exit(-1);
          }
        }
        */
      }

      // add neighbors to the queue
      if(dist < (nDistances-1)) { //-1 as distance indices start at 0
        // iterate over neighbors of sq and add them to the queue if there are not already in it
        for(vector<supernode*>::iterator itNq = sq->neighbors.begin();
            itNq != sq->neighbors.end(); ++itNq) {

          // check that distance is 0 (i.e. immediate neighbor)
          ulong edgeId = getEdgeId(sq->id, (*itNq)->id);
          if(distanceIdxs.count(edgeId) > 0 && distanceIdxs[edgeId] != 0) {
            continue;
          }

          supernode* snq = *itNq;
          deque<supernode*>::iterator lookup = find(pending_list.begin(), pending_list.end(), snq);
          if(lookup == pending_list.end()) {
            set<sidType>::iterator lookup_visited = visited.find(snq->id);
            if(lookup_visited == visited.end()) {
              pending_list.push_back(snq);
              distances[snq->id] = dist + 1;
            }
          }
        }
      }
    }
  }

  // make sure it's symmetric
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    supernode* s = it->second;
    for(vector<supernode*>::iterator itN = s->neighbors.begin();
        itN != s->neighbors.end(); ++itN) {
      supernode* sn = *itN;
      // check if retrieved supernode is already in the list of neighbors
      vector<supernode*>::iterator lookup = find(sn->neighbors.begin(), sn->neighbors.end(), s); // slow
      if(lookup == sn->neighbors.end()) {
        sn->neighbors.push_back(s);
      }
    }
  }

  // count edges
  nbEdges = 0;
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    nbEdges += it->second->neighbors.size();
  }
}

#if 0
void Slice_P::addLongRangeEdges_supernodeBased(int nDistances)
{
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  ulong nSupernodes = _supernodes.size();

  //#ifdef WITH_OPENMP
  //#pragma omp parallel for
  //#endif
  for(int sid = 0; sid < nSupernodes; ++sid) {
    map<sidType, supernode* >::const_iterator it = _supernodes.find(sid);
    supernode* s = it->second;
    printf("[Slice_P] %d/%d\n", s->id, _supernodes.size());

    // add first layer
    for(vector<supernode*>::iterator itN = s->neighbors.begin();
        itN != s->neighbors.end(); ++itN) {
      supernode* sn = *itN;
        // store distance to avoid online computation
      ulong edgeId = getEdgeId(s->id, sn->id);
      //assert(distanceIdxs.find(edgeId) == distanceIdxs.end());
      distanceIdxs[edgeId] = 0;      
    }

    // add second layer
    set<supernode*> layer2;
    for(vector<supernode*>::iterator itN = s->neighbors.begin();
        itN != s->neighbors.end(); ++itN) {
      supernode* sn = *itN;
      for(vector<supernode*>::iterator itNN = sn->neighbors.begin();
          itNN != sn->neighbors.end(); ++itNN) {
        supernode* snn = *itNN;
        // insert sq in the list of neighbors
        if(s->id != snn->id) {
          // check if retrieved supernode is already in the list of neighbors
          vector<supernode*>::iterator lookup = find(s->neighbors.begin(), s->neighbors.end(), snn); // slow
          if(lookup == s->neighbors.end()) {
            // did not find neighbors, add it
            layer2.insert(snn);
          }
        }
      }
    }
    printf("layer2 %ld\n", layer2.size());
    for(set<supernode*>::iterator itlayer2 = layer2.begin(); itlayer2 != layer2.end(); ++itlayer2) {
      supernode* snn = *itlayer2;
      s->neighbors.push_back(snn);
      // store distance to avoid online computation
      ulong edgeId = getEdgeId(s->id, snn->id);
      //assert(distanceIdxs.find(edgeId) == distanceIdxs.end());
      distanceIdxs[edgeId] = 1;
    }

  }    

  // count edges
  nbEdges = 0;
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    nbEdges += it->second->neighbors.size();
  }
}
#endif

void Slice_P::addLongRangeEdges()
{
  vector<node>* centers = getCenters();
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  ulong nSupernodes = _supernodes.size();

#ifdef WITH_OPENMP
  #pragma omp parallel for
#endif
  for(long sid = 0; sid < nSupernodes; ++sid) {
    map<sidType, supernode* >::const_iterator it = _supernodes.find(sid);
    supernode* s = it->second;
    node cs = (*centers)[s->id];

    supernode* sn = 0;
    node cn;
    for(map<sidType, supernode* >::const_iterator it2 = _supernodes.begin();
        it2 != _supernodes.end(); it2++) {
      if(s->id != it2->first) {
        sn = it2->second;
        cn = (*centers)[sn->id];
        double sq_distance = node_square_distance(cs, cn);
        if(sq_distance < MAX_SQ_DISTANCE_LONG_RANGE_EDGES) {
          s->neighbors.push_back(sn);
        }
      }
    }
  }
  delete centers;

  // count edges
  nbEdges = 0;
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    nbEdges += it->second->neighbors.size();
  }
}

void Slice_P::printEdgeStats()
{
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  ulong mean_edges = 0;
  ulong max_edges = 0;
  ulong min_edges = _supernodes.size()*_supernodes.size();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    mean_edges += it->second->neighbors.size();
    if(max_edges < it->second->neighbors.size()) {
      max_edges = it->second->neighbors.size();
    }
    if(min_edges > it->second->neighbors.size()) {
      min_edges = it->second->neighbors.size();
    }
  }
  mean_edges /= _supernodes.size();
  printf("[Slice_P] Edges: mean = %ld, min = %ld, max = %ld\n", mean_edges, min_edges, max_edges);
}

void Slice_P::printDistanceIndicesCount(int nDistances)
{
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  ulong nSupernodes = _supernodes.size();
  map<int, ulong> counts;
  for(int i = 0; i < nDistances; ++i) {
    counts[i] = 0;
  }

  for(ulong sid = 0; sid < nSupernodes; ++sid) {
    map<sidType, supernode* >::const_iterator it = _supernodes.find(sid);
    supernode* s = it->second;

    supernode* sn = 0;
    for(vector < supernode* >::iterator it2 = s->neighbors.begin();
        it2 != s->neighbors.end(); it2++) { 
      sn = *it2;

      // TODO : edge Id or < ?
      int idx = getDistanceIdx(s->id, sn->id);
      ++counts[idx];
    }
  }

  for(map<int, ulong>::iterator it = counts.begin();
      it != counts.end(); ++it) {
    printf("%d -> %ld\n", it->first, it->second);
  }
}

ulong Slice_P::getNbUndirectedEdges()
{
  ulong nUndirectedEdges = 0;
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    for(vector < supernode* >::iterator itN = it->second->neighbors.begin();
        itN != it->second->neighbors.end(); itN++) {

      if(it->first < (*itN)->id) {
        continue;
      }

      ++nUndirectedEdges;
    }
  }
  return nUndirectedEdges;
}

/*
ulong Slice_P::getNbUndirectedEdges2()
{
  ulong nUndirectedEdges = 0;
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    for(vector < supernode* >::iterator itN = it->second->neighbors.begin();
        itN != it->second->neighbors.end(); itN++) {

      if(it->first > (*itN)->id) {
        continue;
      }

      ++nUndirectedEdges;
    }
  }
  return nUndirectedEdges;
}

ulong Slice_P::getNbUndirectedEdges3()
{
  ulong nUndirectedEdges = 0;
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    for(vector < supernode* >::iterator itN = it->second->neighbors.begin();
        itN != it->second->neighbors.end(); itN++) {
      ++nUndirectedEdges;
    }
  }
  return nUndirectedEdges;
}
*/

void Slice_P::generateSupernodeLabelsFromBuffer(const uchar* labels,
                                                int _nLabels,
                                                bool includeBoundaryLabels,
                                                bool useColorImages)
{
  nLabels = _nLabels;
  supernode* s;
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    s = it->second;
    s->setLabel(computeSupernodeLabel(s, labels));
  }
}

void Slice_P::getSupernodeLabelsFromBuffer(map<sidType, uchar>& supernodeLabels,
                                           const uchar* labels,
                                           int _nLabels,
                                           bool includeBoundaryLabels,
                                           bool useColorImages)
{
  nLabels = _nLabels;
  supernode* s;
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    s = it->second;
    supernodeLabels[it->first] = computeSupernodeLabel(s, labels);
  }
}


inline labelType Slice_P::computeSupernodeLabel(supernode* s, const uchar* _labels)
{
  labelType label = 0;
  node n;
  // add +1 to deal with both 0 and 1 based indices.
  int *counts = new int[nLabels + 1];
  for(int l = 0; l < nLabels + 1; ++l) {
    counts[l] = 0;
  }

  ulong width = getWidth();
  ulong sizeSlice = getWidth()*getHeight();
  nodeIterator ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd()) {
    ni.get(n);
    ni.next();
    ++counts[_labels[n.z*sizeSlice + n.y*width + n.x]];
  }

  int maxCount = 0;
  for(int l = 0; l < nLabels + 1; ++l) {
    if(maxCount < counts[l]) {
      maxCount = counts[l];
      label = l;
    }
  }

  delete[] counts;
  return label;
}


void Slice_P::countSupernodeLabels(map<labelType, ulong>& labelCount)
{
  labelType label;
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    label = it->second->getLabel();
    if(labelCount.count(label)==0)
      labelCount[label] = 1;
    else
      labelCount[label]++;
  }
}

void Slice_P::countSupernodeLabels(labelType* labels, map<labelType, ulong>& labelCount)
{
  labelType label;
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    label = labels[it->first];
    if(labelCount.count(label)==0)
      labelCount[label] = 1;
    else
      labelCount[label]++;
  }
}

//-------------------------------------------------------------LOADING FUNCTIONS

void Slice_P::readSupernodeLabels_Multiscale(vector<string>& prediction_filenames,
                                           int nLabels)
{
  int nScales = prediction_filenames.size();
  // allocate memory to store probabilities for the given number of scales
  int probSize = nLabels*nScales;
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    supernode_data* d = new supernode_data(probSize, 0);
    it->second->data = d;
  }

  int prob_offset = 0;
  for(int i=0; i < nScales; i++) {
    vector<int> labels; // not used
    readSupernodeLabels(prediction_filenames[i].c_str(),
                        labels, prob_offset);
    prob_offset += nLabels;
  }

  // Sanity check to make sure we don't have negative probabilities
  double pb;
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    for(int i=0; i < nScales; i++) {
      for(int l=0; l < nLabels; l++) {
        pb = this->getProb(it->first,l,i);
        if(pb < 0) {
          printf("[Slice_P] readSupernodeLabels_Multiscale. Neagative probability detected at scale=%d/%d and label=%d/%d\n", i,nScales,l,nLabels);
          exit(-1);
        }
      }
    }
  }
}

int Slice_P::readSupernodeLabels(const char* prediction_filename,
                                 vector<int>& labels,
                                 int prob_offset)
{
  int ret = 0;
  ifstream predict(prediction_filename);
  if(predict.fail()) {
    printf("[Slice_P] Error while loading prediction file %s\n",prediction_filename);
    return -1;
  }
  char line[SVM_MAX_LINE_LENGTH];
  predict.getline(line, SVM_MAX_LINE_LENGTH); // first line contains labels

  if(strncmp(line,"labels",6) == 0) {
    // load second line
    predict.getline(line, SVM_MAX_LINE_LENGTH);
  }
  predict.close();

  // parse line to detect type of file
  int len = strlen(prediction_filename);
  //printf("prediction_filename+len-6 %s\n", prediction_filename+len-6);
  //if(len > 6 && strcmp(prediction_filename+len-6,".local")==0)
  if(len > 6 && strstr(prediction_filename,".local")!=0)
    {
      //printf("[Slice] readSupernodeProbs\n");
      ret = readSupernodeProbs(prediction_filename,labels, prob_offset);
    }
  else
    {
      if(strstr(line,"e+"))
        {
          //printf("[Slice] readSupernodeLabels_e\n");
          ret = readSupernodeLabels_e(prediction_filename,labels, prob_offset);
        }
      else
        {
          //printf("[Slice] readSupernodeLabels_float\n");
          ret = readSupernodeLabels_float(prediction_filename,labels, prob_offset);
        }
    }

  return ret;
}


// Import file format : label cost
// where label is 1-based index and cost=-log(p)
int Slice_P::readSupernodeLabels_e(const char* prediction_filename,
                                 vector<int>& labels,
                                 int prob_offset)
{
  int idx;
  char* ps;
  char* pe;
  int l;
  int nLabels;

  ifstream predict(prediction_filename);
  if(predict.fail()) {
    printf("[Slice] Error while loading prediction file %s\n",prediction_filename);
    return -1;
  }
  char line[SVM_MAX_LINE_LENGTH];
  predict.getline(line, SVM_MAX_LINE_LENGTH); // first line contains labels

  string paramVOC;
  Config::Instance()->getParameter("voc", paramVOC);
  bool useVOC = paramVOC.c_str()[0] == '1';

  if(strncmp(line,"labels",6) == 0)
    {
      // read labels
      ps = strchr(line,' ');
      if(ps != 0)
        {
          ps++;
          while((pe=strchr(ps,' ')))
            {
              *pe = 0;
              l = atoi(ps);
              labels.push_back(l);
              ps = pe+1;
            }
          l = atoi(ps);
          labels.push_back(l);
        }

      nLabels = labels.size();
    }
  else
    {
      // count number of elements
      nLabels = 0;
      char* ps = line;
      char* pe = line + strlen(line);
      //printf("line %s %c %c\n",line,*ps,*pe);
      while(ps != pe)
        {
          //printf("p %x %c %x\n",ps,*ps,pe);
          while(ps != pe && *ps == ' ') ps++;
          while(ps != pe && *ps != ' ') ps++;

          if(ps != pe)
            {
              if(!useVOC)
                labels.push_back(nLabels);
              nLabels++;
            }
        }

      if(useVOC)
        {
          // HACK VOC
          for(int i = 1 ; i < 21; i++)
            labels.push_back(i);
          labels.push_back(0);
        }

      /*
      while((pLine=strchr(pLine+1,' ')))
        {
          labels.push_back(nLabels);
          nLabels++;
        }
      */
      predict.seekg(0, ios::beg);
      PRINT_MESSAGE("[Slice] Warning : no labels provided. %d labels detected\n",nLabels);
    }

  this->setNbLabels(nLabels);

  PRINT_MESSAGE("[Slice] Index/Labels :");
  idx = 0;
  for(vector<int>::iterator itLabels = labels.begin();
      itLabels != labels.end(); itLabels++, idx++)
    PRINT_MESSAGE(" %d/%d", idx,labels[idx]);
  PRINT_MESSAGE("\n");

  // Read file content
  PRINT_MESSAGE("[Slice] Loading probabilities from prediction file %s\n", prediction_filename);

  const map<sidType, supernode* >& _supernodes = getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++)
    {
      predict.getline(line, SVM_MAX_LINE_LENGTH);
      //printf("line %d:%s\n",idx,line);
      // skip blank characters
      ps = line;
      while(*ps == ' ') ps++;
      //ps=strchr(line,' ');
      if(ps != 0 && *ps != '\n' && *ps != 0)
        {
          // prediction file contains both labels + probabilities
          //*ps = 0;

          pe=strchr(ps,' ');
          *pe = 0;

          // load data
          supernode_data* d = it->second->data;
          if(!d)
            {
              d = new supernode_data(nLabels, 0);
              it->second->data = d;
            }

          //d->label = atoi(line);

          if(useVOC)
            d->label = (int)atof(ps);
          else
            d->label = (int)atof(ps)-1; // subtract 1

          //while(*ps != ' ') ps++;
          ps = pe+1;

          // skip blank characters
          while(*ps == ' ') ps++;

          idx = 0;
          while((pe=strchr(ps,' ')))
            {
              *pe = 0;
              d->prob_estimates[prob_offset+labels[idx]] = exp(-atof(ps)); // file contains -log(p)
              //printf("pb %f %g\n",atof(ps),d->prob_estimates[labels[idx]]);
              ps = pe+1;
              // skip blank characters
              while(*ps == ' ') ps++;
              idx++;
            }
          d->prob_estimates[prob_offset+labels[idx]] = exp(-atof(ps));
          //printf("pb %f\n",d->prob_estimates[labels[idx]]);

          // normalize probabilities between 0 and 1
          double totalP = 0;
          for(int i=0; i < nLabels; i++)
            totalP += d->prob_estimates[prob_offset+i];
          for(int i=0; i < nLabels; i++)
            d->prob_estimates[prob_offset+i] /= totalP;
          
          /*
          totalP = 0;
          for(int i=0; i < nLabels; i++)
            {
              printf("%d:%g ",i,d->prob_estimates[i]);
              totalP += d->prob_estimates[i];
            }
          printf("\ntotalP=%g\n",totalP);
          */

          //if(it->second->data)
          //  delete it->second->data;
          //it->second->data = d;
        }
      else
        {
          // no probabilities in the prediction file
          supernode_data* d = new supernode_data(nLabels, 0);
          d->label = atoi(line);
          //printf("label %d\n",d->label);
          for(int i=0;i<nLabels;i++)
            d->prob_estimates[i] = 0;

          if(it->second->data)
            delete it->second->data;
          it->second->data = d;
        }
    }

  predict.close();
  return 0;
}

int Slice_P::readSupernodeLabels_float(const char* prediction_filename,
                                         vector<int>& labels,
                                         int prob_offset)
{
  int idx;
  char* ps;
  char* pe;
  int l;
  int nLabels;

  ifstream predict(prediction_filename);
  if(predict.fail())
    {
      printf("[Slice] Error while loading prediction file %s\n",prediction_filename);
      return -1;
    }
  char line[SVM_MAX_LINE_LENGTH];
  predict.getline(line, SVM_MAX_LINE_LENGTH); // first line contains labels

  if(strncmp(line,"labels",6) == 0)
    {
      // read labels
      ps = strchr(line,' ');
      if(ps != 0)
        {
          ps++;
          while((pe=strchr(ps,' ')))
            {
              *pe = 0;
              l = atoi(ps);
              labels.push_back(l);
              ps = pe+1;
            }
          l = atoi(ps);
          labels.push_back(l);
        }

      nLabels = labels.size();
    }
  else
    {
      // count number of elements
      nLabels = 0;
      char* ps = line;
      char* pe = line + strlen(line);
      //printf("line %s %c %c\n",line,*ps,*pe);
      while(ps != pe)
        {
          //printf("p %x %c %x\n",ps,*ps,pe);
          while(ps != pe && *ps == ' ') ps++;
          while(ps != pe && *ps != ' ') ps++;

          if(ps != pe)
            {
              labels.push_back(nLabels);
              nLabels++;
            }
        }
      predict.seekg(0, ios::beg);
      printf("[Slice] Warning : no labels provided. %d labels detected\n",nLabels);
    }

  this->setNbLabels(nLabels);

  PRINT_MESSAGE("[Slice] Index/Labels :");
  idx = 0;
  for(vector<int>::iterator itLabels = labels.begin();
      itLabels != labels.end(); itLabels++, idx++)
    PRINT_MESSAGE(" %d/%d", idx,labels[idx]);
  PRINT_MESSAGE("\n");

  // Read file content
  PRINT_MESSAGE("[Slice] Loading probabilities from prediction file %s\n", prediction_filename);
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++)
    {
      predict.getline(line, SVM_MAX_LINE_LENGTH);
      //printf("line %s\n",line);
      // skip blank characters
      ps = line;
      while(*ps == ' ') ps++;
      if(ps != 0)
        {
          // prediction file contains both labels + probabilities

          pe=strchr(ps,' ');
          *pe = 0;

          // load data
          supernode_data* d = it->second->data;
          if(!d)
            {
              d = new supernode_data(nLabels, 0);
              it->second->data = d;
            }

          d->label = (int)atof(ps);
          //printf("p label %d\n",d->label);

          //while(*ps != ' ') ps++;
          ps = pe+1;

          // skip blank characters
          while(*ps == ' ') ps++;

          idx = 0;
          while((pe=strchr(ps,' ')))
            {
              *pe = 0;
              //printf("pb %f\n",atof(ps));
              d->prob_estimates[prob_offset+labels[idx]] = atof(ps);
              ps = pe+1;
              // skip blank characters
              while(*ps == ' ') ps++;
              idx++;
            }
          //printf("pb %f\n",atof(ps));
          d->prob_estimates[prob_offset+labels[idx]] = atof(ps);
          //if(it->second->data)
          //  delete it->second->data;
          //it->second->data = d;
        }
      else
        {
          // no probabilities in the prediction file
          supernode_data* d = new supernode_data(nLabels, 0);
          d->label = atoi(line);
          //printf("label %d\n",d->label);
          for(int i=0;i<nLabels;i++)
            d->prob_estimates[i] = 0;

          if(it->second->data)
            delete it->second->data;
          it->second->data = d;
        }
    }

  predict.close();
  return 0;
}


int Slice_P::readSupernodeProbs(const char* prediction_filename,
                              vector<int>& labels,
                              int prob_offset)
{
  int idx;
  char* ps;
  char* pe;
  int l;
  int nLabels;

  ifstream predict(prediction_filename);
  if(predict.fail()) {
    printf("[Slice] Error while loading prediction file %s\n",prediction_filename);
    return -1;
  }
  char line[SVM_MAX_LINE_LENGTH];
  predict.getline(line, SVM_MAX_LINE_LENGTH); // first line contains labels

  string paramVOC;
  Config::Instance()->getParameter("voc", paramVOC);
  bool useVOC = paramVOC.c_str()[0] == '1';

  if(strncmp(line,"labels",6) == 0)
    {
      // read labels
      ps = strchr(line,' ');
      if(ps != 0)
        {
          ps++;
          while((pe=strchr(ps,' ')))
            {
              *pe = 0;
              l = atoi(ps);
              labels.push_back(l);
              ps = pe+1;
            }
          l = atoi(ps);
          labels.push_back(l);
        }

      nLabels = labels.size();
    }
  else
    {
      // count number of elements
      nLabels = 0;
      char* ps = line;
      char* pe = line + strlen(line);
      //printf("line %s %c %c\n",line,*ps,*pe);
      while(ps != pe) {
        //printf("p %x %c %x\n",ps,*ps,pe);
        while(ps != pe && *ps == ' ') ps++;
        while(ps != pe && *ps != ' ') ps++;

        if(ps != pe) {
          if(!useVOC) {
            labels.push_back(nLabels);
          }
          nLabels++;
        }
      }

      if(useVOC) {
        // HACK VOC
        for(int i = 1 ; i < 21; i++)
          labels.push_back(i);
        labels.push_back(0);
      } else {
        labels.push_back(nLabels);
      }

      nLabels++;

      predict.seekg(0, ios::beg);
      printf("[Slice] Warning : no labels provided. %d labels detected\n",nLabels);
    }

  setNbLabels(nLabels);

  PRINT_MESSAGE("[Slice] Index/Labels :");
  idx = 0;
  for(vector<int>::iterator itLabels = labels.begin();
      itLabels != labels.end(); itLabels++, idx++)
    PRINT_MESSAGE(" %d/%d", idx,labels[idx]);
  PRINT_MESSAGE("\n");

  // Read file content
  PRINT_MESSAGE("[Slice] Loading probabilities from prediction file %s\n", prediction_filename);

  float maxProb;
  float pb;
  int label = 0;
  
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++)
    {
      supernode* s = it->second;
      predict.getline(line, SVM_MAX_LINE_LENGTH);
      //printf("line %d/%d:%s\n",sid,this->mSupernodes.size(),line);

      supernode_data* d = it->second->data;
      if(!d) {
        d = new supernode_data(nLabels, 0);
        it->second->data = d;
      }

      //printf("p label %d\n",d->label);

      ps = line;
      // skip blank characters
      while(*ps == ' ') ps++;

      idx = 0;
      maxProb = 0;
      while((pe=strchr(ps,' ')))
        {
          *pe = 0;
          //pb = atof(ps);
          pb = exp(-atof(ps));
          if(pb < 0)
            {
              printf("[SliceData] sid=%d idx=%d pb=%g\n",s->id,idx,pb);
              exit(-1);
            }
          //printf("%d:%f ",idx, atof(ps));
          d->prob_estimates[prob_offset+labels[idx]] = pb;
          ps = pe+1;
          while(*ps == ' ') ps++;

          if(maxProb < pb)
            {
              maxProb = pb;
              label = labels[idx];
            }

          idx++;
        }
      //printf("%d:%f, label=%d\n",idx, atof(ps),label);
      //d->prob_estimates[labels[idx]] = atof(ps);
      pb = exp(-atof(ps));
      if(pb < 0)
        {
          printf("[SliceData] sid=%d idx=%d pb=%g\n",s->id,idx,pb);
          exit(-1);
        }
      d->prob_estimates[prob_offset+labels[idx]] = pb;
      if(maxProb < pb)
        {
          //maxProb = pb;
          label = labels[idx];
        }

      d->label = label;

      //if(s->data)
      //  delete s->data;
      //s->data = d;
    }

  //PRINT_MESSAGE("[SliceData] Reading prediction file : DONE\n");

  predict.close();
  return 0;
}

