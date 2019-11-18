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

#include <gudhi/Permutahedral_representation/face_from_indices.h>
#include <gudhi/Functions/Constant_function.h>
#include <gudhi/Functions/PL_approximation.h>
#include <gudhi/Query_result.h>

#include <vector>
#include <tuple>

namespace Gudhi {

namespace coxeter_triangulation {

/** \class Implicit_manifold_intersection_oracle
 *  \brief An oracle that supports the intersection query on an implicit manifold.
 *
 *  \tparam Function_ The function template parameter. Should be a model of 
 *   the concept FunctionForImplicitManifold.
 *  \tparam Domain_function_ The domain function template parameter. Should be a model of
 *   the concept FunctionForImplicitManifold.
 *
 *  \ingroup coxeter_triangulation
 */
template<class Function_,
	 class... Domain_functions_>
class Implicit_manifold_intersection_oracle {

  typedef std::tuple<Domain_functions_...> Domain_function_tuple;

public:
  /* Computes the affine coordinates of the intersection point of the implicit manifold
   * and the affine hull of the simplex. */
  template <class Simplex_handle,
	    class Constraint_set,
	    class Triangulation>
  Eigen::VectorXd compute_lambda(const Simplex_handle& simplex,
				 const Constraint_set& constraint_set,
				 const Triangulation& triangulation) const {
    std::size_t cod_d = this->cod_d();
    std::size_t matrix_dim = cod_d + 1 + constraint_set.size();
    if (simplex.dimension() != matrix_dim - 1) {
      std::cerr << "Implicit_manifold_intersection_oracle::compute lambda error: simplex of wrong dimension.\n";
      return Eigen::VectorXd();
    }
    Eigen::MatrixXd matrix(matrix_dim, matrix_dim);
    for (std::size_t i = 0; i < matrix_dim; ++i)
      matrix(0, i) = 1;
    std::size_t j = 0;
    for (auto v: simplex.vertex_range()) {
      Eigen::VectorXd v_coords = function_(triangulation.cartesian_coordinates(v));
      std::size_t i = 1;
      for (; i < cod_d + 1; ++i)
	matrix(i, j) = v_coords(i-1);
      for (std::size_t& I: constraint_set) {
	Eigen::VectorXd bv_coords = domain_function<I>()(triangulation.cartesian_coordinates(v));
	matrix(i++, j) = bv_coords(0);
      }
      j++;
    }
    Eigen::VectorXd z(matrix_dim);
    z(0) = 1;
    for (std::size_t i = 1; i < matrix_dim; ++i)
      z(i) = 0;
    Eigen::VectorXd lambda = matrix.colPivHouseholderQr().solve(z);
    return lambda;
  }

  /* Computes the intersection result for a given simplex in a triangulation. */
  template <class Simplex_handle,
	    class Triangulation>
  Query_result<Simplex_handle> intersection_result(const Eigen::VectorXd& lambda,
						   const Simplex_handle& simplex,
						   const Triangulation& triangulation) const {
    using QR = Query_result<Simplex_handle>;
    std::size_t amb_d = triangulation.dimension();
    std::size_t cod_d = simplex.dimension();

    for (std::size_t i = 0; i < (std::size_t)lambda.size(); ++i)
      if (lambda(i) < 0 || lambda(i) > 1)
	return QR({Eigen::VectorXd(), false});

    Eigen::MatrixXd vertex_matrix(cod_d + 1, amb_d);
    auto v_range = simplex.vertex_range();
    auto v_it = v_range.begin();
    for (std::size_t i = 0; i < cod_d + 1 && v_it != v_range.end(); ++v_it, ++i) {
      Eigen::VectorXd v_coords = triangulation.cartesian_coordinates(*v_it);
      for (std::size_t j = 0; j < amb_d; ++j)
	vertex_matrix(i, j) = v_coords(j);
    }
    Eigen::VectorXd intersection = lambda.transpose()*vertex_matrix;
    return QR({intersection, true});
  }
  
public:

  /** \brief Ambient dimension of the implicit manifold. */
  std::size_t amb_d() const {
    return function_.amb_d();
  }
  
  /** \brief Codimension of the implicit manifold. */
  std::size_t cod_d() const {
    return function_.cod_d();
  }

  /** \brief Intersection query with the relative interior of the manifold.
   *  
   *  \details The returned structure Query_result contains the boolean value
   *   that is true only if the intersection point of the query simplex and
   *   the relative interior of the manifold exists, the intersection point
   *   and the face of the query simplex that contains 
   *   the intersection point.
   *   
   *  \tparam Simplex_handle The class of the query simplex.
   *   Needs to be a model of the concept SimplexInCoxeterTriangulation.
   *  \tparam Triangulation The class of the triangulation.
   *   Needs to be a model of the concept TriangulationForManifoldTracing.
   *
   *  @param[in] simplex The query simplex. The dimension of the simplex
   *   should be the same as the codimension of the manifold 
   *   (the codomain dimension of the function).
   *  @param[in] triangulation The ambient triangulation. The dimension of 
   *   the triangulation should be the same as the ambient dimension of the manifold 
   *   (the domain dimension of the function).
   */
  template <class Simplex_handle,
	    class Constraint_set,
	    class Triangulation>
  Query_result<Simplex_handle> intersects(const Simplex_handle& simplex,
					  const Constraint_set& constraint_set,
					  const Triangulation& triangulation) const {
    Eigen::VectorXd lambda = compute_lambda(simplex, constraint_set, triangulation);
    return intersection_result(lambda, simplex, triangulation);
  }
  
  /** \brief Returns true if the input point lies inside the piecewise-linear
   *   domain induced by the given ambient triangulation that defines the relative
   *   interior of the piecewise-linear approximation of the manifold.
   *
   * @param p The input point. Needs to have the same dimension as the ambient
   *  dimension of the manifold (the domain dimension of the function).
   * @param triangulation The ambient triangulation. Needs to have the same
   *  dimension as the ambient dimension of the manifold 
   *  (the domain dimension of the function).
   */
  template <class Triangulation,
	    class Constraint_set>
  bool lies_in_domain(const Eigen::VectorXd& p,
		      const Constraint_set& constraint_set,
		      const Triangulation& triangulation) const {
    for (const std::size_t& I: constraint_set) {
      Eigen::VectorXd pl_p = make_pl_approximation(domain_function<I>(), triangulation)(p);
      if (pl_p(0) > 0) 
	return false;
    }
    return true;
  }

  /** \brief Returns the function that defines the interior of the manifold. */
  const Function_& function() const {
    return function_;
  }

  /** \brief Returns the I-th domain function for a given I. 
   *  @param[in] I Template parameter that represents the number of the domain function
   *   starting from 0.
   */
  template <std::size_t I>
  const typename std::tuple_element<I, Domain_function_tuple>::type& domain_function() const {
    return std::get<I>(domain_function_tuple_);
  }

  /** \brief Constructs an intersection oracle for an implicit manifold potentially 
   *   with boundary from given function and domain.
   *
   *  @param function The input function that represents the implicit manifold
   *   before the restriction with the domain.
   *  @param domain_functions The input domain functions that can be used to define an implicit
   *   manifold with boundary and corners.
   */
  Implicit_manifold_intersection_oracle(const Function_& function,
					const Domain_functions_&... domain_functions)
    : function_(function), domain_function_tuple_(std::make_tuple(domain_functions...)) {}
  
private:
  Function_ function_;
  Domain_function_tuple domain_function_tuple_;
};

/** \brief Static constructor of an intersection oracle from a function with a domain.
 *
 *  @param function The input function that represents the implicit manifold
 *   before the restriction with the domain.
 *  @param domain_functions The input domain functions that can be used to define an implicit
 *   manifold with boundary and corners.
 *
 *  \ingroup coxeter_triangulation
 */
template<class Function_,
	 class... Domain_functions_>
Implicit_manifold_intersection_oracle<Function_, Domain_functions_...>
make_oracle(const Function_& function,
	    const Domain_functions_&... domain_functions){
  return Implicit_manifold_intersection_oracle<Function_, Domain_functions_...>(function,
										domain_functions...);
}

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif
