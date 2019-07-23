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

#ifndef FUNCTIONS_TRANSLATE_H_
#define FUNCTIONS_TRANSLATE_H_

#include <cstdlib>

#include <gudhi/Functions/Function.h>
#include <Eigen/Dense>

namespace Gudhi {

namespace coxeter_triangulation {

/* \class Translate
 * \brief Translates the zero-set of the function by a vector.
 * The underlying function corresponds to f(x-off), where off is the offset vector.
 *
 * \tparam Function The function template parameter. Should be a model of 
 * the concept FunctionForImplicitManifold.
 *
 * \ingroup coxeter_triangulation
 */
template <class Function>
struct Translate : public Function {
  
  /** 
   * \brief Value of the function at a specified point.
   * @param[in] p The input point. The dimension needs to coincide with the ambient dimension.
   */
  void evaluate(const Eigen::VectorXd& p, Eigen::VectorXd& result) const {
    fun_.evaluate(p - off_, result);
  }

  /** \brief Returns the domain (ambient) dimension. */
  std::size_t amb_d() const {return fun_.amb_d();}

  /** \brief Returns the codomain dimension. */
  std::size_t cod_d() const {return fun_.cod_d();}

  /** \brief Returns a point on the zero-set. */
  void seed(Eigen::VectorXd& result) const {
    fun_.seed(result);
    result += off_;
  }

  /** 
   * \brief Constructor of the translated function.
   *
   * @param[in] function The function to be translated.
   * @param[in] off The offset vector. The dimension should correspond to the 
   * domain (ambient) dimension of 'function'.
   */
  Translate(const Function& function, const Eigen::VectorXd& off) :
    fun_(function), off_(off) {
  }
  Function fun_;
  Eigen::VectorXd off_;
};


/** 
 * \brief Static constructor of a translated function.
 *
 * @param[in] function The function to be translated.
 * @param[in] off The offset vector. The dimension should correspond to the 
 * domain (ambient) dimension of 'function'.
 *
 * \tparam Function The function template parameter. Should be a model of 
 * the concept FunctionForImplicitManifold.
 */
template <class Function>
Translate<Function> translate(const Function& function, Eigen::VectorXd off) {
  return Translate<Function>(function, off);
}

} // namespace coxeter_triangulation

} // namespace Gudhi

#endif
