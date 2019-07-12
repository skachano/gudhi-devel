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

#include <gudhi/Functions/Domain_from_function.h>
#include <gudhi/Functions/Constant_function.h>
#include <gudhi/Query_result.h>

namespace Gudhi {

namespace coxeter_triangulation {

using Infinite_domain = Domain_from_function<Constant_function>;

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

  /* Computes the affine coordinates of the intersection point of the implicit manifold
   * and the affine hull of the simplex. */
  template <class Simplex_handle, 
	    class Triangulation>
  Eigen::VectorXd compute_lambda(const Simplex_handle& simplex,
				 const Triangulation& triangulation) const {
    std::size_t cod_d = this->cod_d();
    Eigen::MatrixXd matrix(cod_d + 1, cod_d + 1);
    for (std::size_t i = 0; i < cod_d + 1; ++i)
      matrix(0, i) = 1;
    std::size_t j = 0;
    for (auto v: simplex.vertex_range()) {
      Eigen::VectorXd v_coords = fun_(triangulation.cartesian_coordinates(v));
      for (std::size_t i = 1; i < cod_d + 1; ++i)
	matrix(i, j) = v_coords(i-1);
      j++;
    }
    Eigen::VectorXd z(cod_d + 1);
    z(0) = 1;
    for (std::size_t i = 1; i < cod_d + 1; ++i)
      z(i) = 0;
    return matrix.colPivHouseholderQr().solve(z);
  }

  /* Computes the affine coordinates of the intersection point of the boundary
   * of the implicit manifold and the affine hull of the simplex. */
  template <class Simplex_handle, 
	    class Triangulation>
  Eigen::VectorXd compute_boundary_lambda(const Simplex_handle& simplex,
					  const Triangulation& triangulation) const {
    std::size_t cod_d = this->cod_d();
    Eigen::MatrixXd matrix(cod_d + 2, cod_d + 2);
    for (std::size_t i = 0; i < cod_d + 2; ++i)
      matrix(0, i) = 1;
    std::size_t j = 0;
    for (auto v: simplex.vertex_range()) {
      Eigen::VectorXd v_coords = fun_(triangulation.cartesian_coordinates(v));
      for (std::size_t i = 1; i < cod_d + 1; ++i)
	matrix(i, j) = v_coords(i-1);
      matrix(cod_d + 1, j) = domain_.function_(triangulation.cartesian_coordinates(v))(0);
      j++;
    }
    Eigen::VectorXd z(cod_d + 2);
    z(0) = 1;
    for (std::size_t i = 1; i < cod_d + 2; ++i)
      z(i) = 0;
    return matrix.colPivHouseholderQr().solve(z);
  }

  
public:

  template <class Simplex_handle,
	    class Triangulation>
  Query_result<Simplex_handle> intersects(const Simplex_handle& simplex,
					  const Triangulation& triangulation) {
    return Query_result<Simplex_handle>();
  }

  template <class Simplex_handle,
	    class Triangulation>
  Query_result<Simplex_handle> intersects_boundary(const Simplex_handle& simplex,
						   const Triangulation& triangulation) {
    return Query_result<Simplex_handle>();
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
  Implicit_manifold_intersection_oracle(const Function_& function,
					const Domain_& domain,
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
  Implicit_manifold_intersection_oracle(const Function_& function, double threshold = 0)
    : fun_(function),
      domain_(Constant_function(function.amb_d(), 1, Eigen::VectorXd::Constant(1,-1))),
      threshold_(threshold) {}
  
private:
  Function_ fun_;
  Domain_ domain_;
  double threshold_ = 0;
};

template<class Function_,
	 class Domain_>
Implicit_manifold_intersection_oracle<Function_, Domain_> make_oracle(const Function_& f,
								     const Domain_& dom,
								     double threshold = 0){
  return Implicit_manifold_intersection_oracle<Function_, Domain_>(f, dom, threshold);
}


template<class Function_>
Implicit_manifold_intersection_oracle<Function_> make_oracle(const Function_& f){
  return Implicit_manifold_intersection_oracle<Function_>(f);
}

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif