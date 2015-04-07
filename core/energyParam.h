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

#ifndef ENERGY_PARAM_H
#define ENERGY_PARAM_H

#include <fstream>
#include <vector>
#include <string>

class EnergyParam
{
 public:

  EnergyParam();

   /**
     * Load global stats from a file
     */
  EnergyParam(const char* filename);

  EnergyParam(const EnergyParam& _param);

  EnergyParam& operator=(EnergyParam const& _param);

  ~EnergyParam();

  void load(const char* filename);
  void loadLabelNames(const char* filename);
  void print() const;
  void printMetaData();
  void save(const char* filename);

  int sizePsi;
  int nUnaryWeights;
  int nClasses;
  int nDistances;
  int nGradientLevels;
  int nOrientations;
  int nScales;
  int nLocalScales;
  int nScalingCoefficients;
  bool includeLocalEdges;
  bool useGlobalClassifier;

  // weight vector
  double* weights;

  std::vector<std::string> labelNames;
};

#endif //ENERGY_PARAM_H
