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

#include "constraint_set.h"
#include "Config.h"
#include "svm_struct_globals.h" // for SSVM_PRINT

ConstraintSet* ConstraintSet::pInstance = 0; // initialize pointer

ConstraintSet::ConstraintSet()
{
  max_number_constraints = CONSTRAINT_SET_DEFAULT_SIZE;
  sortingType = CS_DISTANCE;
  string config_tmp;
  if(Config::Instance()->getParameter("cs_max_number_constraints", config_tmp)) {
    max_number_constraints = atoi(config_tmp.c_str());
    printf("[ConstraintSet] max_number_constraints = %ld\n", max_number_constraints);
  }
}

ConstraintSet::~ConstraintSet()
{
  clear();
}

void ConstraintSet::clear()
{
  for(map<cs_id_type, vector< constraint >* >::iterator itC = constraints.begin();
      itC != constraints.end(); ++itC) {
    for(vector<constraint>::iterator it = itC->second->begin();
        it != itC->second->end(); ++it) {
      delete[] it->first->w;
      delete[] it->first;
    }
    delete itC->second;
  }
  constraints.clear();
}

void ConstraintSet::getConstraints(vector< constraint >& all_cs)
{
  //todo
  for(map<cs_id_type, vector< constraint >* >::iterator itC = constraints.begin();
      itC != constraints.end(); ++itC) {
    for(vector<constraint>::iterator it = itC->second->begin();
        it != itC->second->end(); ++it) {
      all_cs.push_back(*it);
    }
  }
}

const vector< constraint >* ConstraintSet::getConstraints(cs_id_type id)
{
  if(constraints.find(id) != constraints.end()) {
    return constraints[id];
  } else {
    return 0;
  }
}

vector< constraint >* ConstraintSet::getMutableConstraints(cs_id_type id)
{
  assert(constraints.find(id) != constraints.end());
  return constraints[id];
}

const constraint* ConstraintSet::getMostViolatedConstraint(cs_id_type id, double* w, int* max_index)
{
  constraint* c = 0;
  if(constraints.find(id) != constraints.end()) {
    vector< constraint >* _cs = constraints[id];
    double max_margin = 0;
    bool initialized = false;
    int i = 0;
    int i_max = 0;
    for(vector<constraint>::iterator it = _cs->begin();
        it != _cs->end(); ++it) {
      double margin = computeScore(it->first->w, w) + it->first->loss;
      if(!initialized || (margin > max_margin)) {
        max_margin = margin;
        c = &(*it);
        i_max = i;
        initialized = true;
      }
      ++i;
    }
    printf("[ConstraintSet] getMostViolatedConstraint(%d): id %d/%ld max_margin %g\n", id, i_max, _cs->size(), max_margin);
    if(max_index) {
      *max_index = i_max;
    }
  }
  return c;
}

double ConstraintSet::computeLoss(cs_id_type id)
{
  double loss = 0;
  if(constraints.find(id) != constraints.end()) {
    vector< constraint >* _cs = constraints[id];
    for(vector<constraint>::iterator it = _cs->begin();
        it != _cs->end(); ++it) {
      loss += it->first->loss;
    }
  }
  return loss;
}

double ConstraintSet::computeMargin(cs_id_type id, double* w)
{
  return computeScore(id, w) + computeLoss(id);
}

double ConstraintSet::computeScore(cs_id_type id, double* w)
{
  double score = 0;
  if(constraints.find(id) != constraints.end()) {
    vector< constraint >* _cs = constraints[id];
    for(vector<constraint>::iterator it = _cs->begin();
        it != _cs->end(); ++it) {
      score += computeScore(it->first->w, w);
    }
  }
  return score;
}

int ConstraintSet::count(cs_id_type id)
{
  if(constraints.find(id) == constraints.end()) {
    return 0;
  } else {
    return constraints[id]->size();
  }
}

double ConstraintSet::computeDistance(cs_id_type id, SWORD* w)
{
  double distance = 0;
  const vector< constraint >* _cs = getConstraints(id);

  switch(sortingType)
    {
    case CS_DISTANCE:
    case CS_MARGIN_DISTANCE:
    case CS_USE_MVC:
      for(vector<constraint>::const_iterator it = _cs->begin();
          it != _cs->end(); ++it) {
        distance += computeDistance(w, it->first->w);
      }
      break;
    case CS_SQUARE_DISTANCE:
      for(vector<constraint>::const_iterator it = _cs->begin();
          it != _cs->end(); ++it) {
        distance += computeSquareDistance(w, it->first->w);
      }
      break;
    case CS_MIN_DISTANCE:
      {
      vector<constraint>::const_iterator it = _cs->begin();
      distance = computeDistance(w, it->first->w);
      ++it;
      for(;it != _cs->end(); ++it) {
        distance = min(distance, computeDistance(w, it->first->w));
      }
      break;
      }
    default:
      printf("[constraint_set] Unknown sorting type %d\n", (int)sortingType);
      exit(-1);
      break;
    }

  return distance/(double)_cs->size();
}

bool ConstraintSet::add(cs_id_type id, SWORD* w, double loss, int sizePsi)
{
  return add(id, w, loss, sizePsi, 0);
}

void ConstraintSet::create_item(SWORD* w, int sizePsi, c_item* _item)
{
#if USE_SPARSE_VECTORS
  // count number of non zeros
  int n_non_zeros = 0;
  int i = 0;
  while(w[i].wnum) {
    if(w[i].weight != 0) {
      ++n_non_zeros;
    }
    ++i;
  }

  // add +1 for last element whose index is -1
  _item->w = new SWORD[n_non_zeros + 1];

  int si = 0; // index for sparse vector
  i = 0;
  while(w[i].wnum) {
    if(w[i].weight != 0) {
      _item->w[si] = w[i];
      ++si;
    }
    ++i;
  }
  _item->w[si].wnum = 0;
#else
  _item->w = new SWORD[sizePsi];
  int i = 0;
  while(w[i].wnum) {
    _item->w[i] = w[i];
    ++i;
  }
  _item->w[i].wnum = 0;
  _item->w[i].weight = 0;
#endif
  _item->loss = 0;
}

void ConstraintSet::getSortingValue(double& sorting_value, cs_id_type id, SWORD* w)
{
  switch(sortingType)
  {
  case CS_DISTANCE:
  case CS_SQUARE_DISTANCE:
  case CS_MIN_DISTANCE:
  case CS_USE_MVC:
    sorting_value = computeDistance(id, w);
    break;
  case CS_MARGIN_DISTANCE:
    {
      double distance = computeDistance(id, w);
      if(sorting_value == 0) {
        sorting_value = 1e-10;
      }
      sorting_value *= distance;
    }
    break;
  case CS_MARGIN:
    // do nothing. sorting_value already contains margin value
    break;
  default:
    printf("[ConstraintSet] Unknown sorting type %d\n", (int)sortingType);
    exit(-1);
    break;
  }
}

bool ConstraintSet::contains(cs_id_type id, SWORD* w)
{
  ulong n_constraints = count(id);
  vector< constraint >* _cs = 0;

  if(n_constraints == 0) {
    return false;
  } else {
    _cs = constraints[id];
  }

  return contains(_cs, w);
}

bool ConstraintSet::contains(vector< constraint >* _cs, SWORD* w)
{
  bool new_constraint = true;
  for(vector<constraint>::iterator it = _cs->begin(); it != _cs->end(); ++it) {
    SWORD* _w = it->first->w;
    bool is_diff = false;
    int i = 0;
    while(w[i].wnum) {
      if(w[i].weight != _w[i].weight) {
        is_diff = true;
        break;
      }
      ++i;
    }
    if(!is_diff) {
      new_constraint = false;
      break;
    }
  }
  return !new_constraint;
}

bool ConstraintSet::add(cs_id_type id, SWORD* w, double loss, int sizePsi, double sorting_value)
{
  double added = false;
  ulong n_constraints = count(id);
  vector< constraint >* _cs = 0;

  if(n_constraints == 0) {
    _cs = new vector< constraint >;
    constraints[id] = _cs;
    printf("[ConstraintSet] Creating new set for id %d\n", id);
  } else {
    _cs = constraints[id];
  }

  if(n_constraints < max_number_constraints) {

    // check if constraint is different from all the known constraints
    bool existing_constraint = contains(_cs, w);

    if(!existing_constraint) {
      c_item* _item = new c_item;
      create_item(w, sizePsi, _item);
      _item->loss = loss;
      _item->id = _cs->size();
      _cs->push_back(make_pair(_item, sorting_value)); // use distance of 0 for now
      ++n_constraints;
      added = true;
      SSVM_PRINT("[ConstraintSet] cs_id %d: Added new constraint with value %g. Set contains %ld constraints\n",
                 id, sorting_value, _cs->size());
      
      // sort if set is full
      if(_cs->size() == max_number_constraints) {

        printf("[ConstraintSet] cs_id %d: Constraint set is full. Computing distances...\n", id);
        int i = 0;
        for(vector<constraint>::iterator it = _cs->begin(); it != _cs->end(); ++it) {
          getSortingValue(it->second, id, it->first->w);
          printf("[ConstraintSet] cs_id %d: Sorting value %d = %g\n", id, i, it->second);
          ++i;
        }

        std::sort(_cs->begin(), _cs->end(), compare_pair_second<>());
        save("constraints_0_sorted.txt");
      }
    } else {
      SSVM_PRINT("[ConstraintSet] cs_id %d: Already existing constraint in set\n", id);
    }

  } else {

    bool existing_constraint = contains(_cs, w);

    if(!existing_constraint) {
      // add new constraint
      c_item* _item = new c_item;
      create_item(w, sizePsi, _item);
      _item->loss = loss;
      _item->id = _cs->size();
      _cs->push_back(make_pair(_item, sorting_value)); // use distance of 0 for now
      added = true;

      // remove constraint with smallest score
      for(vector<constraint>::iterator it = _cs->begin(); it != _cs->end(); ++it) {
        getSortingValue(it->second, id, it->first->w);
      }
      std::sort(_cs->begin(), _cs->end(), compare_pair_second<>());
      _cs->erase(_cs->begin());
    } else {
      SSVM_PRINT("[ConstraintSet] cs_id %d: Already existing constraint in set\n", id);
    }
  }


  return added;
}

void ConstraintSet::save(const char* filename)
{
  ofstream ofs(filename, ios::out);
  for(map<cs_id_type, vector< constraint >* >::iterator itC = constraints.begin();
      itC != constraints.end(); ++itC) {
    for(vector<constraint>::iterator it = itC->second->begin();
      it != itC->second->end(); ++it) {
      // output id + sorting value + loss
      ofs << itC->first << " " << it->second << " " << it->first->loss;
      SWORD* _w = it->first->w;
      while (_w->wnum) {
        ofs << " " << _w->weight;
        ++_w;
      }
      ofs << endl;
    }
  }
  ofs.close();
}

void ConstraintSet::saveMargins(double* w, const char* filename)
{
  ofstream ofs(filename, ios::out);
  for(map<cs_id_type, vector< constraint >* >::iterator itC = constraints.begin();
      itC != constraints.end(); ++itC) {
    for(vector<constraint>::iterator it = itC->second->begin();
      it != itC->second->end(); ++it) {
      // output id + margin
      double margin = computeScore(it->first->w, w) + it->first->loss;
      ofs << itC->first << " " << margin << endl;
    }
  }
  ofs.close();
}

ulong ConstraintSet::getSize()
{
  ulong size = 0;
  for(map<cs_id_type, vector< constraint >* >::iterator itC = constraints.begin();
      itC != constraints.end(); ++itC) {
    size += itC->second->size();
  }
  return size;
}

void ConstraintSet::printStats(cs_id_type id, double* w)
{
  if(constraints.find(id) != constraints.end()) {
    vector< constraint >* _cs = constraints[id];
    for(vector<constraint>::iterator it = _cs->begin();
        it != _cs->end(); ++it) {
      double score = computeScore(it->first->w, w);
      double loss = it->first->loss;
      double margin = loss + score;
      printf("[ConstraintSet] Score = %g, Loss = %g, Margin = %g\n", score, loss, margin);
    }
  }
}
