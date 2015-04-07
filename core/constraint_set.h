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

#ifndef CONSTRAINT_SET_H
#define CONSTRAINT_SET_H

#include "svm_struct_api_types.h"

// maximum number of constraints to be stored
#define CONSTRAINT_SET_DEFAULT_SIZE 100

struct c_item {
  SWORD* w; // psi
  double loss;
  //int* indices;
  int id;
};

typedef std::pair<c_item*, double> constraint;

template<template <typename> class P = std::less >
struct compare_pair_second {
    template<class T1, class T2> bool operator()(const std::pair<T1, T2>& left, const std::pair<T1, T2>& right) {
        return P<T2>()(left.second, right.second);
    }
};

typedef int cs_id_type;

enum eSortingType
  {
    CS_DISTANCE = 0,
    CS_SQUARE_DISTANCE,
    CS_MIN_DISTANCE,
    CS_MARGIN,
    CS_MARGIN_DISTANCE,
    CS_USE_MVC
  };

//------------------------------------------------------------------------------

class ConstraintSet
{
 public:
  static ConstraintSet* pInstance;

  static ConstraintSet* Instance() 
  {    
    if (pInstance == 0)  // is it the first call?
      {
        pInstance = new ConstraintSet; // create unique instance
      }
    return pInstance; // address of unique instance
  }

  static void setInstance(ConstraintSet* aInstance) {
    pInstance = aInstance;
  }

  ConstraintSet();

  ~ConstraintSet();

  /**
   * Returns true if constraint was added to the set.
   */
  bool add(cs_id_type id, SWORD* w, double loss, int sizePsi);

  bool add(cs_id_type id, SWORD* w, double loss, int sizePsi, double sorting_value);

  void clear();

  double computeDistance(cs_id_type id, SWORD* w);

  double computeLoss(cs_id_type id);

  double computeMargin(cs_id_type id, double* w);

  double computeScore(cs_id_type id, double* w);

  /**
   * Compute score = w*\psi
   */
  inline double computeScore(SWORD* psi, double* w)
  {
    double score = 0;
    SWORD* p = psi;
    while (p->wnum) {
      score += w[p->wnum]*(p->weight);
      ++p;
    }
    return score;
  }

  bool contains(cs_id_type id, SWORD* w);

  bool contains(vector< constraint >* _cs, SWORD* w);

  int count(cs_id_type id);

  void create_item(SWORD* w, int sizePsi, c_item* _item);

  bool isFull() { return constraints.size() == max_number_constraints; }

  /**
   * Get all constraints
   */
  void getConstraints(vector< constraint >& all_cs);

  /**
   * Returns constraints for a specific example
   */
  const vector< constraint >* getConstraints(cs_id_type id);

  vector< constraint >* getMutableConstraints(cs_id_type id);

  /**
   * Returns most violated constraints
   */
  const constraint* getMostViolatedConstraint(cs_id_type id, double* w, int* max_index = 0);

  ulong getCapacity() { return max_number_constraints; }

  ulong getSize();

  void getSortingValue(double& sorting_value, cs_id_type id, SWORD* w);

  void printStats(cs_id_type id, double* w);

  void save(const char* filename);

  void saveMargins(double* w, const char* filename);

  void setSortingAlgorithm(eSortingType _type) { sortingType = _type; }

 private:
  /**
   * Maximum number of constraints per example
   */
  ulong max_number_constraints;

  map<cs_id_type, vector< constraint >* > constraints;

  eSortingType sortingType;

  inline double computeDistance(SWORD* wa, SWORD* wb)
  {
    double d = 0;
    SWORD* pa = wa;
    SWORD* pb = wb;
    while (pa->wnum) {
      double t = wa->weight - wb->weight;
      d += t*t;
      ++pa;
      ++pb;
    }
    return sqrt(d);
  }

  inline double computeSquareDistance(SWORD* wa, SWORD* wb)
  {
    double d = 0;
    SWORD* pa = wa;
    SWORD* pb = wb;
    while (pa->wnum) {
      double t = wa->weight - wb->weight;
      d += t*t;
      ++pa;
      ++pb;
    }
    return d;
  }
};

#endif // CONSTRAINT_SET_H
