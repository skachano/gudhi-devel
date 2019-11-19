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

#include <gudhi/Functions/Function.h>
#include <gudhi/Permutahedral_representation/face_from_indices.h>
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
 *  \tparam Constraint_functions_ The pack template parameter for constraint functions. All functions
 *   should be models of the concept FunctionForImplicitManifold.
 *
 *  \ingroup coxeter_triangulation
 */
class Implicit_manifold_intersection_oracle {

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
      Eigen::VectorXd v_coords = (*function_)(triangulation.cartesian_coordinates(v));
      std::size_t i = 1;
      for (; i < cod_d + 1; ++i)
	matrix(i, j) = v_coords(i-1);
      for (std::size_t I: constraint_set) {
	Eigen::VectorXd v_cart = triangulation.cartesian_coordinates(v);
	Eigen::VectorXd bv_coords = (*constraint_functions_.at(I))(v_cart);
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
    return function_->amb_d();
  }
  
  /** \brief Codimension of the implicit manifold. */
  std::size_t cod_d() const {
    return function_->cod_d();
  }

  /** \brief Intersection query with the manifold.
   *  
   *  \details The returned structure Query_result contains the boolean value
   *   that is true only if the intersection point of the query simplex and
   *   the piece of the manifold that saturates a given set of constraints exists
   *   and the intersection point.
   *   
   *  \tparam Simplex_handle The class of the query simplex.
   *   Needs to be a model of the concept SimplexInCoxeterTriangulation.
   *  \tparam Simplex_handle The class of the constraint set.
   *   Needs to be a model of the concept ConstraintSetForManifoldTracing.
   *  \tparam Triangulation The class of the triangulation.
   *   Needs to be a model of the concept TriangulationForManifoldTracing.
   *
   *  @param[in] simplex The query simplex. The dimension of the simplex
   *   should be the codimension of the manifold (the codomain dimension of the function)
   *   + the number of constraints in constraint_set.
   *  @param[in] constraint_set The set of constraints that defines a stratum on
   *   the manifold.
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
  
  /** \brief Returns true if the input point lies within the piecewise-linear
   *   domains defined by a given set of constraint functions.
   *
   * @param[in] p The input point. Needs to have the same dimension as the ambient
   *  dimension of the manifold (the domain dimension of the function).
   * @param[in] constraint_set The set of constraints that defines a stratum on
   *   the manifold.
   * @param[in] triangulation The ambient triangulation. Needs to have the same
   *  dimension as the ambient dimension of the manifold 
   *  (the domain dimension of the function).
   */
  template <class Triangulation,
	    class Constraint_set>
  bool lies_in_domain(const Eigen::VectorXd& p,
		      const Constraint_set& constraint_set,
		      const Triangulation& triangulation) const {
    for (const std::size_t& I: constraint_set) {
      Eigen::VectorXd pl_p = make_pl_approximation(*constraint_functions_.at(I), triangulation)(p);
      if (pl_p(0) > 0) 
	return false;
    }
    return true;
  }

  /** \brief Returns the function that defines the interior of the manifold. */
  Function* function() const {
    return function_;
  }

  /** \brief Returns the I-th constraint function for a given I. 
   *  @param[in] I Template parameter that represents the index of the constraint function
   *   starting from 0.
   */
  const std::vector<Function*>& constraint_functions() const {
    return constraint_functions_;
  }

  /** \brief Constructs an intersection oracle for an implicit manifold potentially 
   *   with boundary and corners from given function and constraints.
   *
   *  @param[in] function The input function that represents the implicit manifold
   *   before the applying the constraints.
   *  @param[in] constraint_functions The input constraint functions that can be used to define an implicit
   *   manifold with boundary and corners.
   */
  template <class Function_,
	    class Function_ptr_range>
  Implicit_manifold_intersection_oracle(Function_* function,
				        Function_ptr_range& constr_function_range)
    : function_(function) {
    for (auto f_ptr: constr_function_range)
      constraint_functions_.push_back(f_ptr);
  }

private:
  Function* function_;
  std::vector<Function*> constraint_functions_;
};

/** \brief Static constructor of an intersection oracle from a function and constraints.
 *
 *  @param[in] function The input function that represents the implicit manifold
 *   before the applying the constraints.
 *  @param[in] constraint_functions The input constraint functions that can be used to define an implicit
 *   manifold with boundary and corners.
 *
 *  \ingroup coxeter_triangulation
 */
template <class Function_,
	  class Function_ptr_range>
Implicit_manifold_intersection_oracle
make_oracle(Function_* function,
	    Function_ptr_range& constr_function_range) {
  return Implicit_manifold_intersection_oracle(function, constr_function_range);
}

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif
