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

#ifndef DOMAIN_FROM_FUNCTION_H_
#define DOMAIN_FROM_FUNCTION_H_

#include <cstdlib>

#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/** \class Domain_from_function
 * \brief Builds a domain from a given function.
 * The codomain of the function must be one-dimensional.
 * The resulting domain is a model of the concept DomainForManifoldTracing.
 *
 * \tparam Function The function template parameter. Should be a model of 
 * the concept FunctionForImplicitManifold.
 *
 * \ingroup coxeter_triangulation
 */
template <class Function>
struct Domain_from_function {
  bool operator()(const Eigen::VectorXd& p) const {
    return (function_(p)[0] < 0);
  }

  Domain_from_function(const Function& function)
    : function_(function) {
  }

  Function function_;
};


/** 
 * \brief Static constructor of a domain from a function.
 *
 * @param[in] function The input function.
 *
 * \tparam Function The function template parameter. Should be a model of 
 * the concept FunctionForImplicitManifold.
 */
template <typename Function>
Domain_from_function<Function> make_domain_from_function(const Function& function) {
  return Domain_from_function<Function>(function);
}

} // namespace coxeter_triangulation

} // namespace Gudhi

#endif
