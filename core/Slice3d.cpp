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


// standard libraries
#include <sstream>
#include <time.h>

// Third-party libraries
#include "LKM.h"

// SliceMe
#include "Slice3d.h"
#include "globalsE.h"
#include "utils.h"

#define USE_RUN_LENGTH_ENCODING

Slice3d::Slice3d(unsigned char* a_raw_data,
                 int awidth, int aheight,
                 int adepth,
                 sizeSliceType vstep,
                 int anChannels,
                 bool _loadNeighbors)
{
  init();
  delete_raw_data = false; // caller is responsible for freeing memory
  width = awidth;
  height = aheight;
  depth = adepth;
  raw_data = a_raw_data;
  nChannels = anChannels;
  loadNeighbors = _loadNeighbors;

  if(vstep > depth) {
    supernode_step = depth;
  } else {
    supernode_step = vstep;
  }

  sliceSize = width*height;
}

Slice3d::Slice3d(const char* input_dir,
                 int vstep,
                 int nImgs,
                 bool _loadNeighbors)
{
  width = UNITIALIZED_SIZE;
  height = UNITIALIZED_SIZE;
  depth = UNITIALIZED_SIZE;
  nChannels = 1;
  loadNeighbors = _loadNeighbors;

  inputDir = string(input_dir);
  init();
  loadFromDir(input_dir, nImgs);

  sliceSize = width*height;

  if(vstep > depth) {
    supernode_step = depth;
  } else {
    supernode_step = vstep;
  }
}

Slice3d::Slice3d(const char* input_dir,
                 int awidth,
                 int aheight,
                 int adepth,
                 int vstep,
                 bool _loadNeighbors)
{
  init();
  width = awidth;
  height = aheight;
  depth = adepth;
  nChannels = 1;
  loadNeighbors = _loadNeighbors;

  inputDir = string(input_dir);
  loadFromDir(input_dir, depth);

  sliceSize = width*height;

  if(vstep > depth) {
    supernode_step = depth;
  } else {
    supernode_step = vstep;
  }
}

Slice3d::Slice3d(const char* input_dir,
                 const node& _start,
                 const node& _end,
                 int vstep,
                 bool _loadNeighbors)
{
  init();
  width = _end.x - _start.x;
  height = _end.y - _start.y;
  depth = _end.z - _start.z;
  start_x = _start.x;
  start_y = _start.y;
  start_z = _start.z;
  nChannels = 1;
  loadNeighbors = _loadNeighbors;

  inputDir = string(input_dir);
  loadFromDir(input_dir, depth);

  sliceSize = width*height;

  if(vstep > depth) {
    supernode_step = depth;
  } else {
    supernode_step = vstep;
  }
}

Slice3d::Slice3d() {
  init();
}

void Slice3d::init()
{
  mSupervoxels = 0;
  nbEdges = 0;
  supernodeLabelsLoaded = false;
  maxDegree = -1;
  minPercentToAssignLabel = MIN_PERCENT_TO_ASSIGN_LABEL;
  nLabels = 0;
  raw_data = 0;
  includeOtherLabel = true;
  delete_raw_data = true;
  cubeness = SUPERVOXEL_DEFAULT_CUBENESS;

  start_x = 0;
  start_y = 0;
  start_z = 0;

#ifdef USE_REVERSE_INDEXING
  klabels = 0;
#endif
}

Slice3d::~Slice3d()
{
  if(mSupervoxels) {
    for(map< sidType, supernode* >::iterator it = mSupervoxels->begin();
        it != mSupervoxels->end();it++) {
      delete it->second;
    }
    delete mSupervoxels;
  }
  
  if(delete_raw_data && raw_data) {
    delete[] raw_data;
  }
}

uchar Slice3d::at(int x, int y, int z)
{
  return raw_data[z*sliceSize+y*width+x];
}

int Slice3d::raw2Double(double**& ptr_data)
{
  // supervoxel library needs a cube made of ints so we have to convert the cube
  // ask for enough memory for the texels and make sure we got it before proceeding
  ptr_data = new double* [depth];
  if (ptr_data == 0) {
    printf("[Slice3d] Error while allocating memory for 3d volume\n");
    return -1;
  }

  int i = 0;
  for(int z = 0;z < depth; z++) {
    ptr_data[z] = new double[sliceSize];
    for(int xy = 0; xy < sliceSize; xy++) {
      ptr_data[z][xy] = (double)raw_data[i];
      i++;
    }
  }
  return 0;
}

int Slice3d::raw2RGB(unsigned int**& ptr_data)
{
  // supervoxel library needs a cube made of ints so we have to convert the cube
  // ask for enough memory for the texels and make sure we got it before proceeding
  ptr_data = new unsigned int* [depth];
  if (ptr_data == 0) {
    printf("[Slice3d] Error while allocating memory for 3d volume\n");
    return -1;
  }

  int i = 0;
  char c;
  for(int z=0;z<depth;z++) {
    ptr_data[z] = new unsigned int[sliceSize];
    for(int xy = 0;xy<sliceSize;xy++) {
      c = raw_data[i];
      ptr_data[z][xy] = c | (c << 8) | (c << 16);
      i++;
    }
  }
  return 0;
}

void Slice3d::loadSupervoxels(const char* imageDir)
{
  cubeness = SUPERVOXEL_DEFAULT_CUBENESS;
  loadSupervoxels(imageDir, DEFAULT_VOXEL_STEP, SUPERVOXEL_DEFAULT_CUBENESS);
}

void Slice3d::loadSupervoxels(const char* imageDir, const int voxel_step, const float _cubeness)
{
  cubeness = _cubeness;
  supernode_step = voxel_step;
  stringstream soutSupervoxels_nrrd;
  soutSupervoxels_nrrd << imageDir << "supervoxels_" << voxel_step << "_" << cubeness << ".nrrd";
  if(fileExists(soutSupervoxels_nrrd.str().c_str())) {
    uint* outputData = 0;
    int width;
    int height;
    int depth;
    PRINT_MESSAGE("[Slice3d] Loading supervoxels from nrrd file %s\n", soutSupervoxels_nrrd.str().c_str());
#ifdef USE_ITK
    importNRRDCube_uint(soutSupervoxels_nrrd.str().c_str(), outputData, width, height, depth);
#else
    PRINT_MESSAGE("[Slice3d] Set USE_ITK to true to import NRRD cubes\n");
    exit(-1);
#endif
    importSupervoxelsFromBuffer((const uint*)outputData, width, height, depth);
    delete[] outputData;
  } else {
    stringstream soutSupervoxels;
    soutSupervoxels << imageDir << "supervoxels_" << voxel_step << "_" << cubeness;
    if(fileExists(soutSupervoxels.str().c_str())) {
      PRINT_MESSAGE("[Slice3d] Loading supervoxels from %s\n", soutSupervoxels.str().c_str());
      importSupervoxelsFromBinaryFile(soutSupervoxels.str().c_str());
    } else {
      PRINT_MESSAGE("[Slice3d] Generating supervoxels (cubeness=%d)\n", cubeness);
      generateSupervoxels(cubeness);
      if(strlen(imageDir) > 0) {
        PRINT_MESSAGE("[Slice3d] Exporting supervoxels to %s\n", soutSupervoxels.str().c_str());
        exportSupervoxelsToBinaryFile(soutSupervoxels.str().c_str());
      }
    }
  }
}

void Slice3d::generateSupervoxels(const double _cubeness)
{
  cubeness = _cubeness;
  int slice_size = width*height;

  // voxel step should not be greater than the number of slices
  if(supernode_step > depth) {
    supernode_step = depth;
  }

  // klabels is a 2d array indexed by z coordinates. Slices are ordered by yx.
#ifndef USE_REVERSE_INDEXING
  sidType** klabels;
#endif

  PRINT_MESSAGE("[Slice3d] Generating supervoxels. vol_size=(%d, %d, %d). voxel_step=%d. cubeness=%d, %fMb needed\n",
                width,height,depth,supernode_step,cubeness,slice_size*depth/(1024.0*1024.0));

  if(cubeness == -1)
    {
      // use cubical supervoxels
      PRINT_MESSAGE("[Slice3d] Uniform sampling...\n");
      klabels = new sidType*[depth];
      for(int z=0;z<depth;z++) {
        klabels[z] = new sidType[slice_size];
      }

      int sid = 0;
      for(int z=0;z<depth; z += supernode_step)
        {
          for(int x=0;x<width; x+= supernode_step)
            for(int y=0;y<height; y += supernode_step)
              {
                // do not use voxel+step
                for(int sz=z;sz<min(depth,(sizeSliceType)z+supernode_step);sz++)
                  {
                    for(int sx=x;sx<min(width,(sizeSliceType)x+supernode_step);sx++)
                      for(int sy=y;sy<min(height,(sizeSliceType)y+supernode_step);sy++)
                        klabels[sz][sy*width+sx] = sid;
                  }
                sid++;
              }
        }
      nLabels = sid;
      PRINT_MESSAGE("[Slice3d] Uniform sampling done. %d labels created\n",nLabels);
    }
  else  
    {
      double** ptr_data;
      raw2Double(ptr_data);
      LKM* lkm = new LKM(false); // do not free memory
      lkm->DoSupervoxelSegmentationForGrayVolume(ptr_data,
                                                 (int)width,(int)height,(int)depth,
                                                 klabels,
                                                 nLabels,
                                                 (int)supernode_step,
                                                 cubeness);

      PRINT_MESSAGE("[Slice3d] Supervoxelization done\n");
      for(int z=0;z<depth;z++) {
        delete[] ptr_data[z];
      }
      delete[] ptr_data;

      delete lkm;
    }

  createIndexingStructures(klabels);
}


void Slice3d::createIndexingStructures(sidType** _klabels, bool force)
{
  if(mSupervoxels !=0) {
    if(force) {
      for(map< sidType, supernode* >::iterator it = mSupervoxels->begin();
          it != mSupervoxels->end();it++) {
        delete it->second;
      }
      delete mSupervoxels;
    } else {
      printf("[Slice3d] Error in createIndexingStructures : structures already existing\n");
      return;
    }
  }

  ulong slice_size = width*height;
  // Creating indexation structure
  PRINT_MESSAGE("[Slice3d] Creating indexing structure. %fMb needed\n",
                (sizeof(supernode)*slice_size*depth/(supernode_step*supernode_step)
                 + sizeof(node) * slice_size*depth)/(1024.0*1024.0));

  mSupervoxels = new map< sidType, supernode* >;
  supernode* s;
  map<sidType,supernode*>::iterator itVoxel;

  PRINT_MESSAGE("[Slice3d] Cube size = (%d,%d,%d)=%ld voxels\n", width, height, depth,slice_size*depth);

  sidType sid;

#ifdef USE_RUN_LENGTH_ENCODING
  sidType previousSid;
  lineContainer* line;

  for(int d=0;d<depth;d++) { 
    for(int y=0;y<height;y++) {

      // create new line
      line = new lineContainer;
      line->coord.x = 0;
      line->coord.y = y;
      line->coord.z = d;
      sid = _klabels[d][y*width]; // x=0

      // add line to supernode
      itVoxel = mSupervoxels->find(sid);
      if(itVoxel == mSupervoxels->end()) {
        // Create new supernode and add it to the list
        s = new supernode;
        s->id = sid;
        (*mSupervoxels)[sid] = s;
      } else {
        // Supernode already exists
        s = itVoxel->second;
      }

      s->addLine(line);

      previousSid = sid;

      for(int x = 0;x < width;x++) {
        sid = _klabels[d][y*width+x];
        if(sid == previousSid) {
          line->length++;
        } else {
          // create new line
          line = new lineContainer;
          line->coord.x = x;
          line->coord.y = y;
          line->coord.z = d;
          line->length = 1;
          previousSid = sid;

          // add line to supernode
          itVoxel = mSupervoxels->find(sid);
          if(itVoxel == mSupervoxels->end()) {
            // Create new supernode and add it to the list
            s = new supernode;
            s->id = sid;
            (*mSupervoxels)[sid] = s;
          } else {
            // Supernode already exists
            s = itVoxel->second;
          }

          s->addLine(line);
        }
      }
    }
  }
#else
  for(int d=0;d<depth;d++)
    for(int y=0;y<height;y++)
      for(int x=0;x<width;x++) {
        sid = _klabels[d][y*width+x];
        itVoxel = mSupervoxels->find(sid);
        if(itVoxel == mSupervoxels->end()) {
          // Create new supernode and add it to the list
          s = new supernode;
          s->id = sid;
          (*mSupervoxels)[sid] = s;
        } else {
          // Supernode already exists
          s = itVoxel->second;
        }

        // Add node to the supernode
        node* p = new node;
        p->z = d;
        p->y = y;
        p->x = x;
        s->addNode(p);
      }
#endif

  PRINT_MESSAGE("[Slice3d] %d supervoxels created\n", (int)mSupervoxels->size());

  /*
#if USE_LONG_RANGE_EDGES
  PRINT_MESSAGE("[Slice3d] Adding long range edges...\n");
  addLongRangeEdges();
  return;
#endif
  */

  if(loadNeighbors) {
    PRINT_MESSAGE("[Slice3d] Indexing neighbors...\n");

    stringstream sout_neighbors;
    sout_neighbors << inputDir << "neighbors_" << supernode_step << "_" << cubeness;
    if(fileExists(sout_neighbors.str().c_str())) {
      PRINT_MESSAGE("[Slice3d] Loading neighbors from %s\n", sout_neighbors.str().c_str());
      ifstream ifs(sout_neighbors.str().c_str());
      string line;
      while(getline(ifs,line)) {
        vector<string> tokens;
        splitString(line, tokens);
        int sid = atoi(tokens[0].c_str());
        supernode* s = (*mSupervoxels)[sid];
        vector<supernode*>* ptrNeighbors = &(s->neighbors);
        for(int i = 1; i < tokens.size(); ++i) {
          int nsid = atoi(tokens[i].c_str());
          supernode* sn = (*mSupervoxels)[nsid];
          s->neighbors.push_back(sn);
        }
      }
      ifs.close();
    } else {
      const int nh_size = 1; // neighborhood size
      sidType nsid;
      vector<supernode*>* ptrNeighbors;
      supernode* sn;
      bool existingSupernode;
      nbEdges = 0;
      for(int z = nh_size; z < depth - nh_size; z++) {
        for(int x = nh_size; x < width - nh_size; x++) {
          for(int y = nh_size; y < height - nh_size; y++) {
            sid = _klabels[z][y*width+x];
            for(int nx = x-nh_size; nx <= x+nh_size; nx++) {
              for(int ny = y-nh_size; ny <= y+nh_size; ny++) {
                for(int nz = z-nh_size; nz <= z+nh_size; nz++) {
                  nsid = _klabels[nz][ny*width+nx];
                  if(sid > nsid) {
                    s = (*mSupervoxels)[sid];
                    sn = (*mSupervoxels)[nsid];
                    if(sn == 0) {
                      printf("[Slice3d] Error : supernode %d is null (coordinate=(%d,%d,%d))\n",nsid,nx,ny,nz);
                      exit(-1);
                    }
                    ptrNeighbors = &(s->neighbors);
                    existingSupernode = false;
                    // Searching in a hash map would be faster but it would also consume
                    // more memory ...
                    for(vector<supernode*>::iterator itN = ptrNeighbors->begin();
                        itN != ptrNeighbors->end(); itN++) {
                      if((*itN)->id == nsid) {
                        existingSupernode = true;
                        break;
                      }
                    }
                    if(!existingSupernode) {
                      s->neighbors.push_back(sn);
                      sn->neighbors.push_back(s);
                      nbEdges++;
                    }
                  }
                }
              }
            }
          }
        }
      }

      PRINT_MESSAGE("Exporting neighbors to %s\n", sout_neighbors.str().c_str());
      ofstream ofs(sout_neighbors.str().c_str());
      for(map<sidType, supernode* >::iterator it = mSupervoxels->begin();
          it != mSupervoxels->end(); it++) {
        s = it->second;
        stringstream sout;
        sout << s->id;
        for(vector < supernode* >::iterator itN = s->neighbors.begin();
            itN != s->neighbors.end();itN++) {
          supernode* ns = *itN;
          sout << " " << ns->id;
        }
        ofs << sout.str() << endl;
      }
      ofs.close();

    }

    maxDegree = -1;
    int d;
    for(map<sidType, supernode* >::iterator it = mSupervoxels->begin();
        it != mSupervoxels->end(); it++) {
      d = it->second->neighbors.size();
      if(maxDegree < d) {
        maxDegree = d;
      }
    }    
    PRINT_MESSAGE("[Slice3d] %ld undirected edges created. Maximum degree = %d\n", nbEdges, maxDegree);
  }
}

uchar* Slice3d::createNodeLabelVolume()
{
  ulong sliceSize = width*height;
  ulong volSize = sliceSize*depth;
  uchar* labelVolume = new uchar[volSize];
  memcpy(labelVolume,raw_data,sizeof(char)*volSize);

  // Loop through supernodes and their neighbors
  supernode* s;
  node n;
  int superpixelLabel;

  for(map<sidType, supernode* >::iterator it = mSupervoxels->begin();
      it != mSupervoxels->end(); it++) {
    s = it->second;
    superpixelLabel = s->getLabel();

    nodeIterator ni = s->getIterator();
    ni.goToBegin();
    while(!ni.isAtEnd()) {
      ni.get(n);
      ni.next();
      labelVolume[n.z*sliceSize
                  +n.y*width
                  +n.x] = superpixelLabel*(255.0f/NUMBER_TYPE);
    }  
  }

  return labelVolume;
}

int Slice3d::getIntensity(int x, int y, int z)
{
  ulong sliceSize = width*height;
  return raw_data[z*sliceSize + y*width + x];
}

// TODO : Using the first channel is only OK for grey images
float Slice3d::getAvgIntensity(sidType supernodeId)
{
  map<sidType, supernode* >::iterator iLabel = mSupervoxels->find(supernodeId);
  if(iLabel == mSupervoxels->end()) {
    printf("[Slice3d] Error for sid %d : iLabel==g->mSupernodes.end()\n",supernodeId);
    return -1;
  }

  float intensity = 0;
  ulong idx = 0;
  supernode* s = iLabel->second;
  node n;
  nodeIterator ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd())
    {
      ni.get(n);
      ni.next();
      idx = ((ulong)n.z)*width*height + (n.y*width) + n.x;
      intensity += (float) raw_data[idx];
    }
  intensity /= s->size();
  return intensity;
}

void Slice3d::importSupervoxelsFromBinaryFile(const char* filename)
{
#ifndef USE_REVERSE_INDEXING
  sidType** klabels = 0;
#endif

  if(klabels != 0) {
    printf("[Slice3d] Error in importSupervoxels : supervoxels have already been generated\n");
    return;
  }
  PRINT_MESSAGE("[Slice3d] Importing supervoxel labels from binary file %s. depth=%d, height=%d, width=%d, supernode_step=%d\n",
                filename,depth,height,width,supernode_step);

  klabels = new sidType*[depth];
  ulong sliceSize = width*height;

  ifstream ifs(filename, ios::binary);
  for(int z = 0; z < depth;z++) {
    klabels[z] = new sidType[sliceSize];
    ifs.read((char*)klabels[z],sliceSize*sizeof(sidType));
  }
  ifs.close();

  createIndexingStructures(klabels);

#ifndef USE_REVERSE_INDEXING
  for(int z = 0; z < depth;z++) {
    delete[] klabels[z];
  }
  delete[] klabels;
#endif
}

void Slice3d::importSupervoxelsFromBuffer(const uint* buffer, int _width, int _height, int _depth)
{
#ifndef USE_REVERSE_INDEXING
  sidType** klabels = 0;
#endif

  printf("[Slice3d] Importing supervoxel labels from buffer. size = (%d,%d,%d) =? (%d,%d,%d), supernode_step=%d\n",
         width, height, depth, _width, _height, _depth, supernode_step);

  assert(_width == width);
  assert(_height == height);
  assert(_depth == depth);

  if(klabels != 0) {
    printf("[Slice3d] Error in importSupervoxels : supervoxels have already been generated\n");
    return;
  }

  klabels = new sidType*[depth];
  ulong sliceSize = width*height;

  ulong idx = 0;
  for(int z=0; z<depth; z++) {
    klabels[z] = new sidType[sliceSize];
    for(ulong s = 0; s < sliceSize; ++s) {
      klabels[z][s] = buffer[idx];
      ++idx;
    }
  }

  createIndexingStructures(klabels);

#ifndef USE_REVERSE_INDEXING
  for(int z = 0; z < depth;z++)
    delete[] klabels[z];
  delete[] klabels;
#endif
}

void Slice3d::importSupervoxels(const char* filename)
{
#ifndef USE_REVERSE_INDEXING
  sidType** klabels = 0;
#endif

  if(klabels != 0) {
    printf("[Slice3d] Error in importSupervoxels : supervoxels have already been generated\n");
    return;
  }

  const int MAX_SIZE = 100;
  char line[MAX_SIZE];

  ifstream ifs(filename);
  ifs.getline(line,MAX_SIZE);
  sscanf(line,"%d %d %d %d",&depth,&height,&width,&supernode_step);

  PRINT_MESSAGE("[Slice3d] Importing supervoxel labels from %s. depth=%d, height=%d, width=%d, supernode_step=%d\n",
                filename,depth,height,width,supernode_step);

  klabels = new sidType*[depth];
  ulong slice_size = width*height;
  for(int z=0;z<depth;z++)
    klabels[z] = new sidType[slice_size];

  for(int z = 0;z<depth;z++) {
    for(int y = 0;y<height;y++) {
      for(int x = 0;x<width;x++) {
        ifs.getline(line,MAX_SIZE);
        klabels[z][y*width+x] = (sidType)atoi(line);
      }
    }
  }
  ifs.close();

  createIndexingStructures(klabels);

#ifndef USE_REVERSE_INDEXING
  for(int z = 0;z < depth;z++)
    delete[] klabels[z];
  delete[] klabels;
#endif
}

void Slice3d::createReverseIndexing(sidType**& _klabels)
{
  _klabels = new sidType*[depth];
  ulong slice_size = width*height;
  for(int z=0;z<depth;z++) {
    _klabels[z] = new sidType[slice_size];
  }

  supernode* s;
  int sid = 0;
  node n;
  for(map<sidType, supernode* >::iterator it = mSupervoxels->begin();
      it != mSupervoxels->end(); it++) {
    sid = it->first;
    s = it->second;
    nodeIterator ni = s->getIterator();
    ni.goToBegin();
    while(!ni.isAtEnd())
      {
        ni.get(n);
        ni.next();
        _klabels[n.z][n.y*width + n.x] = sid;
      }
  }
}

void Slice3d::exportSupervoxels(const char* filename)
{
  sidType** _klabels;
#ifndef USE_REVERSE_INDEXING
  createReverseIndexing(_klabels);
#else
  _klabels = this->klabels;
#endif

  if(_klabels == 0) {
    printf("[Slice3d] Error in exportSupervoxels : supervoxels have been generated yet\n");
    return;
  }

  ofstream ofs(filename, ios::binary);
  ofs << depth << " " << height << " " << width << " " << supernode_step << endl;

  // TODO : use write function to avoid looping
  // export data
  for(int z=0;z<depth;z++) {
    for(int y=0;y<height;y++) {
      for(int x=0;x<width;x++) {
        ofs << _klabels[z][y*width+x] << endl;
      }
    }
  }
  ofs.close();

  // NFO file used by VIVA
  stringstream snfo;
  snfo << filename << ".nfo";
  ofstream nfo(snfo.str().c_str());
  nfo << "cubeDepth " << depth << endl;
  nfo << "cubeHeight " << height << endl;
  nfo << "cubeWidth " << width << endl;
  nfo << "cubeFile " << getNameFromPath(filename) << endl;
  nfo << "type uchar" << endl;
}

void Slice3d::exportSupervoxelsToBinaryFile(const char* filename)
{
  sidType** _klabels;
#ifndef USE_REVERSE_INDEXING
  PRINT_MESSAGE("[Slice3d] Creating reverse index\n");
  createReverseIndexing(_klabels);
#else
  _klabels = this->klabels;
#endif

  if(_klabels == 0) {
    printf("[Slice3d] Error in exportSupervoxels : supervoxels have been generated yet\n");
    return;
  }

  ofstream ofs(filename, ios::binary);
  ulong sliceSize = width*height;
  for(int z=0;z<depth;z++) {
    ofs.write((char*)_klabels[z],sliceSize*sizeof(sidType));
  }
  ofs.close();

  // NFO file used by VIVA
  stringstream snfo;
  snfo << filename << ".nfo";
  ofstream nfo(snfo.str().c_str());
  nfo << "cubeDepth " << depth << endl;
  nfo << "cubeHeight " << height << endl;
  nfo << "cubeWidth " << width << endl;
  nfo << "cubeFile " << getNameFromPath(filename) << endl;
  nfo << "type int" << endl;
}


#ifdef USE_REVERSE_INDEXING
sidType Slice3d::getSID(uint x,uint y,uint z)
{
  return klabels[z][y*width+x];
}
#endif

void Slice3d::loadFromDir(const char* dir, int nImgs)
{
  node end;
  end.x = start_x + width;
  end.y = start_y + height;
  end.z = nImgs;
  node start;
  start.x = start_x;
  start.y = start_y;
  start.z = start_z;
  loadFromDir(dir, start, end);
}

void Slice3d::loadFromDir(const char* dir, const node& start, const node& end)
{
  const int bytes_per_pixel = 1;
  IplImage* img;
  IplImage* gray_img;
  int nImgs = end.z-start.z;

  // load files
  vector<string> files;
  getFilesInDir(dir, files, start.z, "png", true);
  if(files.size() == 0) {
    getFilesInDir(dir, files, start.z, "tif", true);
  }

  int nValidImgs = 0;
  for(vector<string>::iterator itFile = files.begin();
      itFile != files.end(); itFile++) {
    if(width == UNITIALIZED_SIZE) {
      IplImage* img_slice = cvLoadImage(itFile->c_str(),0);
      if(!img_slice)
        continue;

      width = img_slice->width;
      height = img_slice->height;

      cvReleaseImage(&img_slice);
    }

    nValidImgs++;
  }

  if(nImgs != -1) {
    if(nImgs > nValidImgs) {
      printf("[Slice3d] Warning : nImgs=%d > nValidImgs=%d\n", nImgs, nValidImgs);
      nImgs = nValidImgs;
    }
  }
  else {
    nImgs = nValidImgs;
  }

  // ask for enough memory for the texels and make sure we got it before proceeding
  depth = nImgs; //files.size();

  if(!isDirectory(dir)) {
    PRINT_MESSAGE("[Slice3d] Loading data from file %s\n", dir);
    importData(dir);
    return;
  }

  stringstream sVolDataFile;
  sVolDataFile << dir;
  sVolDataFile << "/volumedata";

  if(fileExists(sVolDataFile.str().c_str())) {
    PRINT_MESSAGE("[Slice3d] Loading data from %s\n", sVolDataFile.str().c_str());
    importData(sVolDataFile.str().c_str());
    return;
  }

  ulong n = width*height*sizeof(char);
  raw_data = new uchar[n*depth];

  PRINT_MESSAGE("[Slice3d] Loading %d %dx%d images (%ld pixels) from directory %s\n",
                nImgs, width, height, n, dir);
  int d = 0;
  for(int iImage = 0; iImage < nImgs; iImage++)
    {
      string image = files[iImage];

      // Load image in black and white
      // Do no handle 3d cubes in color for now !
      IplImage* img_slice = cvLoadImage(image.c_str(),0);

      if(!img_slice) {
	PRINT_MESSAGE("[Slice3d] Warning : image %s not loaded properly\n", image.c_str());
        continue;
      }

      if(img_slice->width != width || img_slice->height != height)
        {
          PRINT_MESSAGE("[Slice3d] Warning : (img_slice->width != width || img_slice->height != height)\n");
          if(img_slice->nChannels != bytes_per_pixel)
            {
              gray_img = cvCreateImage(cvSize(img_slice->width,img_slice->height),IPL_DEPTH_8U,bytes_per_pixel);
              cvCvtColor(img_slice,gray_img,CV_RGB2GRAY);
              img = cvCreateImage(cvSize(width,height),IPL_DEPTH_8U,bytes_per_pixel);
              cvResize(gray_img,img);
              cvReleaseImage(&gray_img);
            }
          else
            {
              img = cvCreateImage(cvSize(width,height),IPL_DEPTH_8U,bytes_per_pixel);
              cvResize(img_slice,img);
            }

          memcpy(raw_data+(d*n),img->imageData,n);
          cvReleaseImage(&img);
          cvReleaseImage(&img_slice);
        }
      else
        {
          img = img_slice;

          if(img_slice->nChannels != bytes_per_pixel)
            {
              // FIXME
              printf("[Slice3d] img_slice->nChannels %d\n", img_slice->nChannels);
              exit(-1);
            }

          if(img->widthStep != img->width)
            {
              int i = 0;
              for(int y = 0; y < img->height ; y++)
                for(int x = 0; x < img->width ; x++)
                  {
                    raw_data[d*n+i] = ((uchar*)(img->imageData + img->widthStep*y))[x*img->nChannels];
                    i++;
                  }
            }
          else
            memcpy(raw_data+(d*n),img->imageData,n);

          cvReleaseImage(&img);
        }
      ++d;
    }
}

void Slice3d::loadFromDir(const char* dir, uchar*& raw_data,
                          int& width, int& height, int* nImgs)
{
  const int bytes_per_pixel = 1;
  IplImage* img;
  IplImage* gray_img;

  // load files
  vector<string> files;
  getFilesInDir(dir, files,"png", true);
  if(files.size() == 0) {
    getFilesInDir(dir, files,"tif", true);
  }

  width = -1;
  int nValidImgs = 0;
  for(vector<string>::iterator itFile = files.begin();
      itFile != files.end(); itFile++)
    {
      if((itFile->c_str()[0] != '.') && ( (getExtension(*itFile) == "png") || (getExtension(*itFile) == "tif")) )
        {
          if(width == -1)
            {
              IplImage* img_slice = cvLoadImage(itFile->c_str(),0);
              if(!img_slice)
                continue;

              width = img_slice->width;
              height = img_slice->height;

              cvReleaseImage(&img_slice);
            }

          nValidImgs++;
        }
    }

  if(*nImgs != -1)
    {
      if(*nImgs > nValidImgs)
        {
          printf("[PixelData] Warning : nImgs=%d > nValidImgs=%d\n", *nImgs, nValidImgs);
          *nImgs = nValidImgs;
        }
    }
  else
    *nImgs = nValidImgs;

  // ask for enough memory for the texels and make sure we got it before proceeding
  uint n = width*height*sizeof(char);
  raw_data = new uchar[n*(*nImgs)];

  printf("[PixelData] Loading %d images from directory %s, width=%d, height=%d\n", *nImgs, dir, width, height);
  int d = 0;
  for(int iImage = 0; d < *nImgs; iImage++)
    {
      string image = files[iImage];

      // Load image in black and white
      // Do no handle 3d cubes in color for now !
      IplImage* img_slice = cvLoadImage(image.c_str(),0);

      if(!img_slice)
        continue;

      if(img_slice->width != width || img_slice->height != height)
        {
          if(img_slice->nChannels != bytes_per_pixel)
            {
              gray_img = cvCreateImage(cvSize(img_slice->width,img_slice->height),IPL_DEPTH_8U,bytes_per_pixel);
              cvCvtColor(img_slice,gray_img,CV_RGB2GRAY);
              img = cvCreateImage(cvSize(width,height),IPL_DEPTH_8U,bytes_per_pixel);
              cvResize(gray_img,img);
              cvReleaseImage(&gray_img);
            }
          else
            {
              img = cvCreateImage(cvSize(width,height),IPL_DEPTH_8U,bytes_per_pixel);
              cvResize(img_slice,img);
            }

          memcpy(raw_data+(d*n),img->imageData,n);
          cvReleaseImage(&img);
          cvReleaseImage(&img_slice);
        }
      else
        {
          img = img_slice;

          if(img_slice->nChannels != bytes_per_pixel)
            {
              // FIXME
              printf("[Slice3d] img_slice->nChannels %d\n", img_slice->nChannels);
              exit(-1);
            }

          if(img->widthStep != img->width)
            {
              int i = 0;
              for(int y = 0; y < img->height ; y++)
                for(int x = 0; x < img->width ; x++)
                  {
                    raw_data[d*n+i] = ((uchar*)(img->imageData + img->widthStep*y))[x*img->nChannels];
                    i++;
                  }
            }
          else
            memcpy(raw_data+(d*n),img->imageData,n);

          cvReleaseImage(&img);
        }

      d++;
    }
}


bool Slice3d::importData(const char* filename,
                          const int iwidth,
                          const int iheight,
                          const int idepth)
{
  string ext = getExtension(filename);
  //transform(ext.begin(),ext.end(), ext.begin(), tolower);
  //if(ext == "tif" || ext == "TIF")
  if(ext == "nrrd" || ext == "mha")
    {
      printf("[Slice3d] Import NRRD/MHA Cube\n");
#ifdef USE_ITK
      importCube(filename,
                 raw_data,
                 width,
                 height,
                 depth);
#else
      printf("[Slice3d] Error can not import NRRD Cube without ITK\n");
      return false;
#endif
    }
  else
    {
      if(ext == "tif" || ext == "TIF")
	{
	  printf("[Slice3d] Import tif Cube\n");
#ifdef USE_ITK
	  importTIFCube(filename,
			raw_data,
			width,
			height,
			depth);
#else
	  printf("[Slice3d] Error can not import TIF Cube without ITK\n");
          return false;
#endif
	}
      else
	{

	  if(iwidth == -1)
	    {
	      stringstream snfo;
	      snfo << filename << ".nfo";
	      ifstream info(snfo.str().c_str());
	      if(!info.good()) {
                printf("[Slice3d] Could not load %s\n", snfo.str().c_str());
                return false;
              }

	      string name;
	      string attribute;
	      while(info.good())
		{
		  info >> name;
		  info >> attribute;
		  if(!strcmp(name.c_str(), "cubeDepth"))
		    depth = atoi(attribute.c_str());
		  else if(!strcmp(name.c_str(), "cubeHeight"))
		    height = atoi(attribute.c_str());
		  else if(!strcmp(name.c_str(), "cubeWidth"))
		    width = atoi(attribute.c_str());
		  else if(!strcmp(name.c_str(), "type"))
		    {
		      if(strcmp(attribute.c_str(), "uchar"))
			{
			  printf("[Slice3d] Only uchar types are allowed\n");
                          return false;
			}
		    }
		}
	      info.close();
	    }
	  else
	    {
	      width = iwidth;
	      height = iheight;
	      depth = idepth;
	    }

	  ifstream ifs(filename,ios::binary);
	  ulong n = ((ulong)width)*height*depth*sizeof(char);

	  raw_data = new uchar[n];
	  ifs.read((char*)raw_data,n);
	  ifs.close();
	}
    }
  return true;
}

void Slice3d::exportData(const char* filename)
{
  ulong n = width*height*depth*sizeof(char);
  ofstream ofs(filename);
  ofs.write((char*)raw_data,n);
  ofs.close();

  // NFO file used by VIVA
  stringstream snfo;
  snfo << filename << ".nfo";
  ofstream nfo(snfo.str().c_str());
  nfo << "cubeDepth " << depth << endl;
  nfo << "cubeHeight " << height << endl;
  nfo << "cubeWidth " << width << endl;
  nfo << "cubeFile " << getNameFromPath(filename) << endl;
  nfo << "type uchar" << endl;

  nfo.close();
}

void Slice3d::rescaleRawData()
{
  // bin intensities
  ulong n = width*height*depth;
  const int nBuckets = 256;
  ulong buckets[nBuckets];
  for(int i = 0; i < nBuckets; ++i) {
    buckets[i] = 0;
  }
  for(ulong i = 0; i < n; ++i) {
    ++buckets[raw_data[i]];
  }

  // compute min and max
  int min_intensity_count = n*0.01;
  int cum_count = 0;
  double min_intensity = 0;
  for(int i = 0; i < nBuckets; ++i) {
    cum_count += buckets[i];

    if(cum_count >= min_intensity_count) {
      min_intensity = i;
      break;
    }
  }

  double max_intensity = min_intensity;
  int max_intensity_count = n*0.99;
  for(int i = min_intensity; i < nBuckets; ++i) {
    cum_count += buckets[i];

    if(cum_count >= max_intensity_count) {
      max_intensity = i;
      break;
    }
  }

  // rescale
  PRINT_MESSAGE("[Slice3d] Rescaling data. min=%g, max=%g\n", min_intensity,
                max_intensity);
  double new_intensity = 0;
  for(ulong i = 0; i < n; ++i) {
    new_intensity = (raw_data[i] - min_intensity) / (double)(max_intensity - min_intensity);
    new_intensity *= 255;
    if(new_intensity > 255) {
      new_intensity = 255;
    }
    if(new_intensity < 0) {
      new_intensity = 0;
    }
    raw_data[i] = new_intensity;
  }
}

supernode* Slice3d::getSupernode(sidType sid)
{
  if(sid < (sidType)mSupervoxels->size())
    return (*mSupervoxels)[sid];
  else
    return 0;
}

/*
 * This is different from computeSupernodeLabel(labels) that takes the labels
 * as input instead of the mask images
 */
inline labelType Slice3d::computeSupernodeLabel(sidType sid, uchar* mask_data)
{
  labelType label;
  int countObject = 0;
  int countBackground = 0;
  node n;
  supernode* s = getSupernode(sid);
  nodeIterator ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd()) {
    ni.get(n);
    ni.next();
    if(mask_data[n.z*sliceSize+n.y*width+n.x] == BACKGROUND_MASKVALUE) {
      countBackground++;
    } else {
      countObject++;
    }
  }

  if(includeOtherLabel) {
    int total = ni.size();
    if(countObject > (minPercentToAssignLabel*total)) {
        label = FOREGROUND;
    } else {
      if(countBackground > (minPercentToAssignLabel*total))
        label = BACKGROUND;
      else {
        label = OTHER_LABEL;
      }
    }
  } else {
    if( (countObject/(double)(countBackground + countObject)) > minPercentToAssignLabel) {
      label = FOREGROUND;
    } else {
      label = BACKGROUND;
    }
  }

  return label;
}

labelType Slice3d::computeSupernodeLabel_advanced(sidType sid, uchar* mask_data)
{
  labelType label;
  int countObject = 0;
  int countBackground = 0;
  supernode* s = getSupernode(sid);
  node n;
  nodeIterator ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd())
    {
      ni.get(n);
      ni.next();
      if(mask_data[n.z*sliceSize+n.y*width+n.x] == BACKGROUND_ADVANCED_MASKVALUE)
        countBackground++;
      else
	{
	  if(mask_data[n.z*sliceSize+n.y*width+n.x] == FOREGROUND_ADVANCED_MASKVALUE)
	    countObject++;
	}
    }

  if(countObject>(minPercentToAssignLabel*s->size()))
    label = FOREGROUND;
  else {
    if(countBackground>(minPercentToAssignLabel*s->size()))
      label = BACKGROUND;
    else
      label = OTHER_LABEL;
  }

  return label;
}

void Slice3d::generateSupernodeLabels(const char* fn_annotation,
                                      bool includeBoundaryLabels,
                                      bool useColorImages)
{
  generateSupernodeLabelFromMaskDirectory(fn_annotation,
                                          includeBoundaryLabels,
                                          useColorImages);
}

void Slice3d::generateSupernodeLabelFromMaskDirectory(const char* mask_dir,
                                                      bool includeBoundaryLabels,
                                                      bool useColorImages)
{
  if(!isDirectory(mask_dir)) {
    printf("[Slice3d] Error while loading ground-truth labels: %s is not a valid directory\n", mask_dir);
    return;
  }

  // load all the mask images in a cube of data
  uchar* mask_data;
  int _width;
  int _height;
  int _depth = -1;
  loadFromDir(mask_dir, mask_data,
              _width, _height, &_depth);

  // check if ground-truth cube has to be resized
  if(width != _width || height != _height || depth != _depth) {
    printf("[Slice3d] Resizing mask data (%d,%d,%d) -> (%d,%d,%d)\n",
           _width, _height, _depth, width, height, depth);
    ulong image_size = (ulong)_width*(ulong)_height;
    ulong new_image_size = (ulong)width*(ulong)height;
    ulong new_cube_size = new_image_size*(ulong)depth;
    uchar* new_mask_data = new uchar[new_cube_size];
    ulong cubeIdx = 0;
    for(int z = 0; z < depth; ++z) {
      for(int y = 0; y < height; ++y) {
        for(int x = 0; x < width; ++x) {
          new_mask_data[cubeIdx] = mask_data[z*image_size + y*_width + x];
          ++cubeIdx;
        }
      }
    }
    delete[] mask_data;
    mask_data = new_mask_data;
  }

  generateSupernodeLabelFromMaskImages(mask_data,
                                       includeBoundaryLabels,
                                       useColorImages);

  delete[] mask_data;
}

void Slice3d::generateSupernodeLabelFromMaskImages(uchar* mask_data,
                                                   bool includeBoundaryLabels,
                                                   bool useColorImages)
{
  if(supernodeLabelsLoaded) {      
    printf("[Slice3d] Warning : Supernode labels have already been loaded\n");
    return;
  }

  supernodeLabelsLoaded = true;

  supernode* s;
  int nr_classes = 2;
  if(includeBoundaryLabels) {
    nr_classes = 3;
  }

  // first pass to compute labels
  if(useColorImages) {
    for(map<sidType, supernode* >::iterator it = mSupervoxels->begin();
        it != mSupervoxels->end(); it++) {
      s = it->second;
      s->setData(computeSupernodeLabel_advanced(s->id, mask_data),nr_classes);
    }
  } else {
    for(map<sidType, supernode* >::iterator it = mSupervoxels->begin();
        it != mSupervoxels->end(); it++) {
      s = it->second;
      s->setData(computeSupernodeLabel(s->id, mask_data),nr_classes);
    }
  }

  int slabel;
  if(includeBoundaryLabels) {
    // second pass to change foreground labels to boundary labels if they are
    // touching background labels
    for(map<sidType, supernode* >::iterator it = mSupervoxels->begin();
        it != mSupervoxels->end(); it++) {
      s = it->second;
      slabel = s->getLabel();
      if(slabel == FOREGROUND) {
        for(vector < supernode* >::iterator itN = s->neighbors.begin();
            itN != s->neighbors.end();itN++) {
          supernode* ns = *itN;
          if(s->id == ns->id) {
            printf("[Slice3d] Error : supernode %d has a neighbor with the same sid\n", s->id);
            continue;
          }

          if(ns->getLabel() == BACKGROUND) {
            // BACKGROUND label was found among the neighbors
            // switch label to BOUNDARY
            s->setLabel(BOUNDARY);
            break;
          }
        }
      }

      if(s->getLabel() == FOREGROUND) {
        int countObject = 0;
        int countBackground = 0;
        node n;
        nodeIterator ni = s->getIterator();
        ni.goToBegin();
        while(!ni.isAtEnd()) {
          ni.get(n);
          ni.next();
          if(mask_data[n.z*sliceSize+n.y*width+n.x] == BACKGROUND_MASKVALUE) {
            countBackground++;
          } else {
            countObject++;
          }
        }
        //int minCount = minPercentToAssignLabel*ni.size();
        //if(countObject > minCount && countBackground > minCount) {
        if(countBackground > 0) {
          s->setLabel(BOUNDARY);
        }
      }

    }
  }
}

labelType* Slice3d::getSupernodeLabels()
{
  ulong cubeSize = sliceSize*depth;
  labelType* labelCube = new labelType[cubeSize];

  // Make prediction for every superpixel
  ulong idx = 0;
  labelType label;
  supernode* s;
  node n;
  for(map<sidType, supernode* >::iterator it = mSupervoxels->begin();
      it != mSupervoxels->end(); it++) {
    s = it->second;
    label = s->getLabel();

    nodeIterator ni = s->getIterator();
    ni.goToBegin();
    while(!ni.isAtEnd()) {
      ni.get(n);
      ni.next();
      idx = ((ulong)n.z*sliceSize) + n.y*width + n.x;
      //labelCube[idx] = label;
      labelCube[idx] = (label/2.0)*255; // temporary hack
    }
  }

  return labelCube;
}

void Slice3d::createSupernodeLabels(const uchar* nodeLabels,
                                    labelType*& labelCube,
                                    int nClasses)
{
  ulong cubeSize = sliceSize*(ulong)depth;
  labelCube = new labelType[cubeSize];
  nClasses = max(1,nClasses-1);

  // Create cube containing the labels of each supervoxel.
  ulong idx = 0;
  labelType label;
  supernode* s;
  node n;
  for(map<sidType, supernode* >::iterator it = mSupervoxels->begin();
      it != mSupervoxels->end(); it++) {
    s = it->second;
    label = nodeLabels[it->first];

    nodeIterator ni = s->getIterator();
    ni.goToBegin();
    while(!ni.isAtEnd())
      {
        ni.get(n);
        ni.next();
        idx = ((ulong)n.z*sliceSize) + n.y*width + n.x;
        //labelCube[idx] = label;
        labelCube[idx] = (label/(float)nClasses)*255;
      }
  }
}

void Slice3d::createOverlayAnnotationImage(const char* filename, int imageId)
{
  labelType label;
  const int nChannels = 3;
  uchar col[nChannels];
  col[0] = 255; col[1] = 0; col[2] = 0;

  IplImage* img = cvCreateImage(cvSize(width,height),IPL_DEPTH_8U, nChannels);
  uchar* pImg;
  ulong imageSize = width*height;
  ulong cubeIdx = imageId*imageSize;
  for(int y=0; y<height; y++) {
    for(int x=0; x<width; x++) {
      pImg = &((uchar*)(img->imageData + img->widthStep*y))[x*img->nChannels];
      label = getSupernodeLabel(getSid(x,y,imageId));

      for(int i=0; i < nChannels; i++) {
        if(col[i] > 0 && label != T_BACKGROUND)
          pImg[i] = col[i];
        else
          pImg[i] = raw_data[cubeIdx];
      }

      cubeIdx++;
    }
  }

  cvSaveImage(filename, img);
  cvReleaseImage(&img);
}

void Slice3d::exportOverlay(const char* filename)
{
#ifndef USE_ITK
  printf("[Slice3d] Warning : you have to set USE_ITK to true to export overlays\n");
#else

  ulong cubeSize = sliceSize*(ulong)depth;
  labelType* outputCube = new labelType[cubeSize*3];
  const int nChannels = 3;
  uchar col[nChannels];
  col[0] = 255; col[1] = 0; col[2] = 0;
  ulong outputCubeIdx = 0;
  ulong idx = 0;
  labelType label;
  supernode* s;
  node n;
  for(map<sidType, supernode* >::iterator it = mSupervoxels->begin();
      it != mSupervoxels->end(); it++) {
    s = it->second;
    label = s->getLabel();

    nodeIterator ni = s->getIterator();
    ni.goToBegin();
    while(!ni.isAtEnd()) {
      ni.get(n);
      ni.next();
      idx = ((ulong)n.z*sliceSize) + n.y*width + n.x;
      outputCubeIdx = idx*nChannels;
      
      for(int i=0; i < nChannels; i++) {
        if(col[i] > 0 && label != T_BACKGROUND)
          outputCube[outputCubeIdx] = col[i];
        else
          outputCube[outputCubeIdx] = raw_data[idx];
        
        outputCubeIdx++;
      }
    }
  }

  exportColorTIFCube(outputCube, filename,
                     depth, height, width);

  //const int firstImageToExtract = getDepth()/2;
  const int firstImageToExtract = 1;
  const int nImagesToExtract = 3;
  stringstream soutOverlayImage;
  soutOverlayImage << getDirectoryFromPath(filename);
  soutOverlayImage << "/snapshot_" << getNameFromPath(filename);
  PRINT_MESSAGE("[Slice3d] Exporting %d-th image from cube %s\n",
                firstImageToExtract, soutOverlayImage.str().c_str());
  exportImageFromColorCube(soutOverlayImage.str().c_str(),
                           outputCube,
                           getWidth(), getHeight(), getDepth(),
                           firstImageToExtract, nImagesToExtract);

  delete[] outputCube;
#endif
}

void Slice3d::exportOverlay(const char* filename, labelType* labels)
{
#ifndef USE_ITK
  printf("[Slice3d] Warning : you have to set USE_ITK to true to export overlays\n");
#else

  ulong cubeSize = sliceSize*(ulong)depth;
  labelType* outputCube = new labelType[cubeSize*3];
  const int nChannels = 3;
  uchar col[nChannels];
  col[0] = 255; col[1] = 0; col[2] = 0;
  ulong outputCubeIdx = 0;
  ulong idx = 0;
  labelType label;
  supernode* s;
  node n;
  for(map<sidType, supernode* >::iterator it = mSupervoxels->begin();
      it != mSupervoxels->end(); it++) {
    s = it->second;
    label = labels[it->first];

    nodeIterator ni = s->getIterator();
    ni.goToBegin();
    while(!ni.isAtEnd()) {
      ni.get(n);
      ni.next();
      idx = ((ulong)n.z*sliceSize) + n.y*width + n.x;
      outputCubeIdx = idx*nChannels;
      
      for(int i=0; i < nChannels; i++) {
        if(col[i] > 0 && label != T_BACKGROUND)
          outputCube[outputCubeIdx] = col[i];
        else
          outputCube[outputCubeIdx] = raw_data[idx];
        
        outputCubeIdx++;
      }
    }
  }

  exportColorTIFCube(outputCube, filename,
                     depth, height, width);

  //const int firstImageToExtract = getDepth()/2;
  const int firstImageToExtract = 1;
  const int nImagesToExtract = 3;
  stringstream soutOverlayImage;
  soutOverlayImage << getDirectoryFromPath(filename);
  soutOverlayImage << "/snapshot_" << getNameFromPath(filename);
  PRINT_MESSAGE("[Slice3d] Exporting %d-th image from cube %s\n",
                firstImageToExtract, soutOverlayImage.str().c_str());
  exportImageFromColorCube(soutOverlayImage.str().c_str(),
                           outputCube,
                           getWidth(), getHeight(), getDepth(),
                           firstImageToExtract, nImagesToExtract);

  delete[] outputCube;
#endif
}

void Slice3d::exportSupernodeLabels(const char* filename, int nClasses)
{
  ulong cubeSize = sliceSize*(ulong)depth;
  labelType* labelCube = new labelType[cubeSize];
  nClasses = max(1,nClasses-1);

  // Create cube containing the labels of each supervoxel.
  ulong idx = 0;
  labelType label;
  supernode* s;
  node n;
  for(map<sidType, supernode* >::iterator it = mSupervoxels->begin();
      it != mSupervoxels->end(); it++) {
    s = it->second;
    label = s->getLabel();

    nodeIterator ni = s->getIterator();
    ni.goToBegin();
    while(!ni.isAtEnd()) {
      ni.get(n);
      ni.next();
      idx = ((ulong)n.z*sliceSize) + n.y*width + n.x;
      labelCube[idx] = (label/(float)nClasses)*255;
    }
  }

  exportCube(labelCube, filename, depth, height, width);

  delete[] labelCube;
}

void Slice3d::exportSupervoxelsToNRRD(const char* filename)
{
#ifndef USE_ITK
  PRINT_MESSAGE("[Slice3d] Set USE_ITK to true to import NRRD cubes\n");
#else
  ulong cubeSize = sliceSize*(ulong)depth;
  uint* labelCube = new uint[cubeSize];

  // Create cube containing the labels of each supervoxel.
  ulong idx = 0;
  supernode* s;
  node n;
  for(map<sidType, supernode* >::iterator it = mSupervoxels->begin();
      it != mSupervoxels->end(); it++) {
    s = it->second;

    nodeIterator ni = s->getIterator();
    ni.goToBegin();
    while(!ni.isAtEnd()) {
      ni.get(n);
      ni.next();
      idx = ((ulong)n.z*sliceSize) + n.y*width + n.x;
      labelCube[idx] = it->first;
    }
  }

  exportNRRDCube(labelCube, filename, depth, height, width);

  delete[] labelCube;
#endif
}

void Slice3d::exportSupernodeLabels(const char* filename, int nClasses,
                                    vector<labelType>& labels,
                                    int _nLabels)
{
  ulong cubeSize = sliceSize*(ulong)depth;
  labelType* labelCube = new labelType[cubeSize];
  memset(labelCube,0,cubeSize*sizeof(labelType));
  nClasses = max(1,nClasses-1);

  // Create cube containing the labels of each supervoxel.
  ulong idx = 0;
  labelType label;
  supernode* s;
  node n;

  if(_nLabels == -1) {
    _nLabels = cubeSize;
  }

  for(int sid=0; sid < _nLabels; sid++)
    {
      s = getSupernode(sid);
      label = labels[sid];
      if(label>nClasses) {
        printf("[Slice3d] Error in exportSupernodeLabels : label=%d > nClasses=%d\n",label,nClasses+1);
        exit(-1);
      }

      nodeIterator ni = s->getIterator();
      ni.goToBegin();
      while(!ni.isAtEnd())
        {
          ni.get(n);
          ni.next();
          idx = ((ulong)n.z*sliceSize) + n.y*width + n.x;
          labelCube[idx] = (label/(float)nClasses)*255;
        }
    }

  exportCube(labelCube,
             filename,
             depth,
             height,
             width);

  delete[] labelCube;
}

void Slice3d::exportSupernodeLabels(const char* filename, int nClasses,
                                    labelType* labels,
                                    int _nLabels,
				    const map<labelType, ulong>* labelToClassIdx)
{
  ulong cubeSize = sliceSize*(ulong)depth;
  labelType* labelCube = new labelType[cubeSize];
  memset(labelCube,0,cubeSize*sizeof(labelType));
  nClasses = max(1,nClasses-1);

  // Create cube containing the labels of each supervoxel.
  ulong idx = 0;
  labelType label;
  supernode* s;
  node n;

  if(_nLabels == -1) {
    _nLabels = cubeSize;
  }

  for(int sid = 0; sid < _nLabels; sid++)
    {
      s = getSupernode(sid);
      label = labels[sid];
      if(label>nClasses) {
        printf("[Slice3d] Error in exportSupernodeLabels : label=%d > nClasses=%d\n",label,nClasses+1);
        exit(-1);
      }

      nodeIterator ni = s->getIterator();
      ni.goToBegin();
      while(!ni.isAtEnd())
        {
          ni.get(n);
          ni.next();
          idx = ((ulong)n.z*sliceSize) + n.y*width + n.x;
          labelCube[idx] = (label/(float)nClasses)*255;
        }
    }

  exportCube(labelCube,
             filename,
             depth,
             height,
             width);

  int firstImageToExtract = getDepth()/2;
  stringstream sBaseName;
  sBaseName << getDirectoryFromPath(filename) << "/";
  sBaseName << getNameFromPathWithoutExtension(filename);
  sBaseName << "_" << firstImageToExtract;
  sBaseName << ".png";

  PRINT_MESSAGE("[Slice3d] Exporting %d-th image from cube %s\n",
                firstImageToExtract, sBaseName.str().c_str());
  exportImageFromCube(sBaseName.str().c_str(),
                      labelCube,
                      getWidth(), getHeight(),
                      firstImageToExtract, 1);

  delete[] labelCube;
}

void Slice3d::resize(sizeSliceType w, sizeSliceType h, sizeSliceType d, map<sidType, sidType>* sid_mapping)
{
  assert(w <= width);
  assert(h <= height);
  assert(d <= depth);

  // copy data
  sizeSliceType new_sliceSize = w*h;
  sizeSliceType cubeSize = new_sliceSize*d;
  uchar* new_raw_data = new uchar[cubeSize];
  ulong new_idx = 0;
  for(int z = 0; z < d; ++z) {
    for(int y = 0; y < h; ++y) {
      for(int x = 0; x < w; ++x) {
        new_raw_data[new_idx] = raw_data[getIndex(x,y,z)];
        ++new_idx;
      }
    }
  }

  // iterate over all the nodes and generate new klabels array
  sidType** new_klabels = new sidType*[d];
  for(int z = 0;z < d;z++) {
    new_klabels[z] = new sidType[new_sliceSize];
  }

  supernode* s;
  node n;
  int sid = 0;
  const map<sidType, supernode* >& _supernodes = getSupernodes();
  for(map<sidType, supernode* >::const_iterator it = _supernodes.begin();
      it != _supernodes.end(); it++) {
    s = it->second;
    nodeIterator ni = s->getIterator();
    ni.goToBegin();
    bool used_supernode = false;
    while(!ni.isAtEnd()) {
      ni.get(n);
      ni.next();

      if(n.z < d && n.y < h && n.x < w) {
        new_klabels[n.z][n.y*w + n.x] = sid;
        used_supernode = true;
      }
    }
    if(used_supernode) {
      if(sid_mapping) {
        (*sid_mapping)[sid] = it->first;
      }
      ++sid;
    }
  }

  width = w;
  height = h;
  depth = d;
  sliceSize = new_sliceSize;

  delete[] raw_data;
  raw_data = new_raw_data;

  createIndexingStructures(new_klabels, true);

#ifndef USE_REVERSE_INDEXING
  for(int z = 0; z < depth;z++) {
    delete[] new_klabels[z];
  }
  delete[] new_klabels;
#else
  klabels = new_klabels;
#endif
}

void Slice3d::exportProbabilities(const char* filename, int nClasses,
                                  float* pbs)
{
  sizeSliceType cubeSize = sliceSize*depth;
  labelType* labelCube = new labelType[cubeSize];
  memset(labelCube,0,cubeSize*sizeof(labelType));
  nClasses = max(1, nClasses-1);

  // Make prediction for every superpixel
  sizeSliceType idx = 0;
  char pb;
  supernode* s;
  node n;

  for(uint sid = 0; sid < getNbSupernodes(); sid++) {
    s = getSupernode(sid);
    pb = pbs[sid]*255; 

    nodeIterator ni = s->getIterator();
    ni.goToBegin();
    while(!ni.isAtEnd()) {
      ni.get(n);
      ni.next();
      idx = n.z*sliceSize + n.y*width + n.x;
      labelCube[idx] = pb;
    }
  }

  exportCube(labelCube, filename, depth, height, width);
  delete[] labelCube;
}
