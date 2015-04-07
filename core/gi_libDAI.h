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

#ifndef GI_LIBDAI_H
#define GI_LIBDAI_H

// SliceMe
#include "Feature.h"
#include "Slice_P.h"

#include "graphInference.h"
#include "energyParam.h"

#include <dai/alldai.h>  // Include main libDAI header file
#include <dai/bp.h>

// standard libraries
#include <map>
#include <vector>

//------------------------------------------------------------------------------

class GI_libDAI : public GraphInference
{
 public:


  /**
   * Constructor for SSVM framework
   */
  GI_libDAI(Slice_P* _slice, 
            const EnergyParam* _param,
            double* _smw,
            labelType* _groundTruthLabels,
            double* _lossPerLabel,
    	    Feature* _feature,
            std::map<sidType, nodeCoeffType>* _nodeCoeffs,
            std::map<sidType, edgeCoeffType>* _edgeCoeffs
            );

  void addLocalEdges();

  void addNodes();

  void addUnaryNodes();

  void allocateGraph();

  void buildPottsModel(double gamma);

  void buildSSVMGraph();

  double computeEnergy(labelType* nodeLabels);

  void computeEdgePotentials();
  void computeNodePotentials();

  /**
   * Copy labels to nodeLabelsBP
   * Also computes the energy and the loss
   * @param nodeLabelsBP labels computed by using BP::findMaximum()
   * @param nodeLabelsGroundTruth if not null, the loss is also computed
   */
  void copyLabels_MSRC(std::vector<std::size_t>& labels,
                       labelType* nodeLabels,
                       dai::BP& bp,
                       dai::FactorGraph& fg);

  void createColoredAnnotationImageMax(std::vector<std::size_t>& labels,
                                       const char* outputFilename,
                                       map<labelType, ulong>* labelToClassIdx);

  void getLabels(dai::BP& bp,
                 dai::FactorGraph& fg,
                 std::vector<std::size_t>& labels);

  void getMarginals(float* marginals,
                    const int label);

  dai::Real** getUnaryPotentials() { return unaryPotentials; }

  void precomputePotentials();

  double run(labelType* inferredLabels,
             int id,
             size_t maxiter,
             labelType* nodeLabelsGroundTruth = 0,
             bool computeEnergyAtEachIteration = false,
             double* _loss = 0);

 private:
  std::vector<dai::Factor> factors;
  std::vector<dai::Var> vars;
  ulong nVars;
  ulong nFactors;
  // The factors are normalized so that the maximum potential value is equal to MAX_POTENTIAL
  dai::Real maxPotential;

  dai::Real** unaryPotentials;
  uint nUnaryPotentials;
  map<uint, dai::Real*> edgePotentials;
  uint nEdgePotentials;
};

#endif //GI_LIBDAI_H
