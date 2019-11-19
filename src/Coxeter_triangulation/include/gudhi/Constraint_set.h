/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Siargey Kachanovich
 *
 *    Copyright (C) 2019 Inria
 *
 *    Modification(s):
 *      - YYYY/MM Author: Description of the modification
 */

#ifndef CONSTRAINT_SET_H_
#define CONSTRAINT_SET_H_

class Constraint_set {

public:

  typedef std::set<std::size_t>::const_iterator iterator;
  
  void insert(std::size_t i) {
    set_.insert(i);
  }

  void erase(std::size_t i) {
    set_.erase(i);
  }
  
  Constraint_set with(std::size_t i) const {
    Constraint_set result = *this;
    result.insert(i);
    return result;
  }

  Constraint_set without(std::size_t i) const {
    Constraint_set result = *this;
    result.erase(i);
    return result;
  }

  iterator begin() const {
    return set_.begin();
  }

  iterator end() const {
    return set_.end();
  }

  std::size_t size() const {
    return set_.size();
  }
  
  Constraint_set() {}
  
  template<class Range>
  Constraint_set(const Range& range) {
    for (auto i: range)
      set_.insert(i);
  }
private:
  std::set<std::size_t> set_;
};

#endif
