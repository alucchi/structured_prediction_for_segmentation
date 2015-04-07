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

#include "energyParam.h"
#include "inference_globals.h"
#include "graphInference.h"

#include <fstream>
#include <stdlib.h>
#include <map>

using namespace std;

//------------------------------------------------------------------------------

EnergyParam::EnergyParam()
{
  sizePsi = 0;
  nUnaryWeights = 0;
  nClasses = 0;
  nGradientLevels = 0;
  nOrientations = 0;
  nScalingCoefficients = 0;
  includeLocalEdges = false;
  weights = 0;
  nDistances = 1;
}

EnergyParam::EnergyParam(const EnergyParam& _param)
{
  sizePsi = _param.sizePsi;
  nUnaryWeights = _param.nUnaryWeights;
  nClasses = _param.nClasses;
  nGradientLevels = _param.nGradientLevels;
  nOrientations = _param.nOrientations;
  nScalingCoefficients = _param.nScalingCoefficients;
  includeLocalEdges = _param.includeLocalEdges;
  weights = _param.weights;
  nDistances = _param.nDistances;
}

EnergyParam& EnergyParam::operator=(EnergyParam const& _param)
{
  sizePsi = _param.sizePsi;
  nUnaryWeights = _param.nUnaryWeights;
  nClasses = _param.nClasses;
  nGradientLevels = _param.nGradientLevels;
  nOrientations = _param.nOrientations;
  nScalingCoefficients = _param.nScalingCoefficients;
  includeLocalEdges = _param.includeLocalEdges;
  weights = _param.weights;
  nDistances = _param.nDistances;

  return *this;
}

EnergyParam::EnergyParam(const char* filename)
{
  load(filename);
}

EnergyParam::~EnergyParam()
{
  if(weights) {
    delete[] weights;
  }
}

void EnergyParam::save(const char* filename)
{
  ofstream ofs(filename);
  ofs << sizePsi << endl;
  ofs << nClasses << endl;
  ofs << nGradientLevels << endl;
  ofs << nOrientations << endl;
  ofs << nScales << endl;
  ofs << nLocalScales << endl;
  ofs << nScalingCoefficients << endl;

  ofs.precision(9);

  for(int i = 0; i < sizePsi; i++) {
    ofs << weights[i] << endl;
  }

  ofs << endl;
  ofs.close();
}

void EnergyParam::load(const char* filename)
{
  // Load global stats
  ifstream ifsStats(filename);
  if(ifsStats.fail()) {
    printf("[EnergyParam] Error while loading %s\n", filename);
    return;
  }

  ifsStats >> sizePsi;
  ifsStats >> nClasses;
  ifsStats >> nGradientLevels;
  ifsStats >> nOrientations;
  ifsStats >> nScales;
  ifsStats >> nLocalScales;
  ifsStats >> nScalingCoefficients;

  // TODO
  //ifsStats >> nDistances;

  // TODO : might want to save those flags in a file too
  nUnaryWeights = (nClasses==2)?1:nClasses;
  useGlobalClassifier = (nLocalScales!=nScales);
  includeLocalEdges = (nGradientLevels!=0);

  printf("[EnergyParam] sizePsi=%d, nClasses=%d, nGradientLevels=%d, nOrientations=%d, nScales=%d, nLocalScales=%d, nScalingCoefficients=%d\n",
         sizePsi, nClasses, nGradientLevels, nOrientations, nScales, nLocalScales, nScalingCoefficients);

  INFERENCE_PRINT("[EnergyParam] sizePsi=%d, nClasses=%d, nGradientLevels=%d, nOrientations=%d, nScales=%d, nLocalScales=%d, nScalingCoefficients=%d\n",
         sizePsi, nClasses, nGradientLevels, nOrientations, nScales, nLocalScales, nScalingCoefficients);

  // allocate memory
  weights = new double[sizePsi];

  // fill table
  for(int i = 0; i < sizePsi; i++) {
    ifsStats >> weights[i];
  }

  ifsStats.close();
}

void EnergyParam::printMetaData()
{
  printf("[EnergyParam] sizePsi = %d, nClasses = %d, nGradientsLevels = %d, nOrientations = %d\n",
         sizePsi, nClasses, nGradientLevels, nOrientations);
}

void EnergyParam::print() const
{
  if(!weights) {
    return;
  }

  INFERENCE_PRINT("[EnergyParam] w=[");

  //unary terms
  for(int c=0; c < nUnaryWeights; c++) {
    INFERENCE_PRINT("%d:%.2g ", c, weights[c]);
  }

  /*
  for(int s=0; s < nScalingCoefficients; s++)
    INFERENCE_PRINT("* %d:%.2g",SVM_FEAT_INDEX0+s,weights[SVM_FEAT_INDEX0+s]);
  INFERENCE_PRINT("]\n");
  */

  if(nGradientLevels == 0)
  {
      INFERENCE_PRINT("%d:%.2g",nUnaryWeights, weights[nUnaryWeights]);
      INFERENCE_PRINT("\n");
  }
  else
  {
      // pairwise terms
      INFERENCE_PRINT("\t");
      for(int c=0; c < nClasses; c++)
        INFERENCE_PRINT("%d\t", c);
      INFERENCE_PRINT("\n");

      int idx = nUnaryWeights;

      for(int g=0; g < nGradientLevels; g++)
        {

    #if ORIENTATION_SENSITIVE
          for(int o=0; o < nOrientations; o++)
    #endif

          for(int c=0; c < nClasses; c++)
            {
              INFERENCE_PRINT("%d", c);
              for(int c2=0; c2 < nClasses; c2++)
                {
                  INFERENCE_PRINT("\t%.2g", weights[idx]);
                  idx++;
                }
              INFERENCE_PRINT("\n");
            }
        }
    }
}


void EnergyParam::loadLabelNames(const char* filename)
{
  ifstream ifsCol(filename);
  string line;
  int label = 0;
  int r,g,b;
  ulong classIdx;
  char labelName[50];
  while(getline(ifsCol, line))
    {
      sscanf(line.c_str(),"%d %lu %d %d %d %s", &label, &classIdx,&r,&g,&b,labelName);
      INFERENCE_PRINT("label=%d, classIdx=%ld, RGB=(%d,%d,%d) name=%s\n",
             label, classIdx, (int)r,(int)g,(int)b,labelName);
      labelNames.push_back((string)labelName);
    }
  ifsCol.close();
}
