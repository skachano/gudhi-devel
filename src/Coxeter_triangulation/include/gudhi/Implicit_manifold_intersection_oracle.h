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

#ifndef IMPLICIT_MANIFOLD_INTERSECTION_ORACLE_H_
#define IMPLICIT_MANIFOLD_INTERSECTION_ORACLE_H_

#include <Eigen/Dense>

#include <Gudhi/Query_result.h>

namespace Gudhi {

namespace coxeter_triangulation {

/** \class Implicit_manifold_intersection_oracle
 *  \brief An oracle that supports the intersection query on an implicit manifold.
 *
 *  \tparam Function_ The function template parameter. Should be a model of 
 *   the concept FunctionForImplicitManifold.
 *  \tparam Domain_ The domain template parameter. Should be a model of
 *   the concept DomainForManifoldTracing.
 *
 *  \ingroup coxeter_triangulation
 */
template<class Function_,
	 class Domain_ = Infinite_domain>
class Implicit_manifold_intersection_oracle {
  
public:

  template <class Simplex_handle,
	    class Triangulation>
  Query_result<Simplex_handle> intersects(const Simplex_handle& simplex,
					  const Triangulation& triangulation) {
    return Query_result();
  }

  template <class Simplex_handle,
	    class Triangulation>
  Query_result<Simplex_handle> intersects_boundary(const Simplex_handle& simplex,
						   const Triangulation& triangulation) {
    return Query_result();
  }

  bool lies_in_domain(const Eigen::VectorXd& p) {
    return true;
  }
  
  /** \brief Constructs an intersection oracle for an implicit manifold potentially 
   *   with boundary from given function and domain.
   *
   *  @param function The input function that represents the implicit manifold
   *   before the restriction with the domain.
   *  @param domain The input domain that can be used to define an implicit
   *   manifold with boundary.
   *  @param threshold The input parameter that defines the distance in terms
   *   of affine coordinates to a lower-dimensional face for the intersection 
   *   point to be considered lying at the lower-dimensional face.
   */
  Implicit_manifold_intersection_oracle(const Function& function,
					const Domain& domain,
					double threshold = 0)
    : fun_(function), domain_(domain), threshold_(threshold) {}

  /** \brief Constructs an intersection oracle for an implicit manifold 
   *   without boundary from a given function.
   *
   *  @param function The input function that represents the implicit manifold
   *   without boundary.
   *  @param threshold The input parameter that defines the distance in terms
   *   of affine coordinates to a lower-dimensional face for the intersection 
   *   point to be considered lying at the lower-dimensional face.
   */
  Implicit_manifold_intersection_oracle(const Function& function, double threshold = 0)
    : fun_(function), domain_(fun.amb_d()), threshold_(threshold) {}
  
private:
  Function fun_;
  Domain domain_;
  double threshold_ = 0;
};

template<class Function,
	 class Domain>
Implicit_manifold_intersection_oracle<Function, Domain> make_oracle(const Function& f, const Domain& dom, double threshold = 0){
  return Implicit_manifold_intersection_oracle<Function, Domain>(f, dom, threshold);
}


template<class Function>
Implicit_manifold_intersection_oracle<Function> make_oracle(const Function& f){
  return Implicit_manifold_intersection_oracle<Function>(f);
}

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif
