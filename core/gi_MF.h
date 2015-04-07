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

#ifndef GI_MF_H
#define GI_MF_H

// SliceMe
#include "Slice.h"

#include "graphInference.h"
#include "energyParam.h"

#include <map>
#include <vector>

//------------------------------------------------------------------------------

class GI_MF : public GraphInference
{
 public:
  /**
   * Constructor for SSVM framework
   */
  GI_MF(Slice_P* _slice,
        const EnergyParam* _param,
         double* _smw,
         labelType* _groundTruthLabels,
         double* _lossPerLabel,
         Feature* _feature,
         std::map<sidType, nodeCoeffType>* _nodeCoeffs);

  ~GI_MF();

  double run(labelType* inferredLabels,
             int id,
             size_t maxiter,
             labelType* nodeLabelsGroundTruth = 0,
             bool computeEnergyAtEachIteration = false,
             double* _loss = 0);

  void setBelieves(double** _b) { believes = _b; ownBelievesBuffer = false; }

 private:
  void exportBelieves(const char* filename);

  bool ownBelievesBuffer;
  double** believes;
};

#endif // GI_MF_H
