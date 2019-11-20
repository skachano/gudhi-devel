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

#ifndef MANIFOLD_TRACING_H_
#define MANIFOLD_TRACING_H_

#include <gudhi/Query_result.h>

#include <boost/functional/hash.hpp>

#include <Eigen/Dense>

#include <queue>
#include <unordered_map>

namespace Gudhi {

namespace coxeter_triangulation {

/**
 *  \ingroup coxeter_triangulation
 */

/** \class Manifold_tracing 
 *  \brief A class that assembles methods for manifold tracing algorithm.
 *
 *  \tparam Triangulation_ The type of the ambient triangulation.
 *   Needs to be a model of the concept TriangulationForManifoldTracing.
 */
template <class Triangulation_>
class Manifold_tracing {
  
  typedef typename Triangulation_::Simplex_handle Simplex_handle;  
  struct Simplex_hash {
    typedef Simplex_handle argument_type;
    typedef std::size_t result_type;
    result_type operator()(const argument_type& s) const noexcept {
      return boost::hash<typename Simplex_handle::Vertex>()(s.vertex());
    }
  };
  
public:

  typedef std::unordered_set<std::size_t> Constraint_set;
  /** \brief Type of the output simplex map with keys of type Triangulation_::Simplex_handle
   *   and values of type Constraint_set_intersection_pair.
   *   This type should be used for the output in the method manifold_tracing_algorithm.
   */
  typedef std::unordered_map<Simplex_handle,
			     std::pair<Constraint_set, Eigen::VectorXd>,
			     Simplex_hash>
  Out_simplex_map;

  /**
   * \brief Computes the set of k-simplices that intersect
   * a boundaryless implicit manifold given by an intersection oracle, where k
   * is the codimension of the manifold.
   * The computation is based on the seed propagation --- it starts at the 
   * given seed points and then propagates along the manifold.
   *
   * \tparam Point_range Range of points of type Eigen::VectorXd.
   * \tparam Intersection_oracle Intersection oracle that represents the manifold.
   *  Needs to be a model of the concept IntersectionOracle.
   *
   * \param[in] seed_points The range of points on the manifold from which 
   * the computation begins.
   * \param[in] triangulation The ambient triangulation.
   * \param[in] oracle The intersection oracle for the manifold.
   * The ambient dimension needs to match the dimension of the
   * triangulation.
   * \param[out] out_simplex_map The output map, where the keys are k-simplices in
   * the input triangulation that intersect the input manifold and the mapped values 
   * are the intersection points.
   */
  template <class Intersection_oracle>
  void manifold_tracing_algorithm(const Eigen::VectorXd& seed_point,
				  Triangulation_& triangulation,
				  const Intersection_oracle& oracle,
				  Out_simplex_map& out_simplex_map) {
    std::size_t amb_d = oracle.amb_d();
    std::size_t cod_d = oracle.cod_d();
    std::size_t num_constraints = oracle.constraint_functions().size();
    std::queue<std::pair<Simplex_handle, Constraint_set> > queue;

    Constraint_set init_constraint_set;
    for (std::size_t I = 0; I < oracle.constraint_functions().size(); ++I)
      if ((*oracle.constraint_functions().at(I))(seed_point)(0) == 0)
	init_constraint_set.insert(I);
    typename Simplex_handle::Vertex y(amb_d, 0);
    typename Simplex_handle::OrderedSetPartition omega(cod_d + 1 + init_constraint_set.size());
    for (std::size_t i = 0; i < omega.size(); ++i) {
      omega[i].push_back(i);
    }
    for (std::size_t i = omega.size(); i < amb_d + 1; ++i) {
      omega[omega.size() - 1].push_back(i);
    }
    Simplex_handle init_s(y, omega);
    Eigen::VectorXd barycenter = triangulation.barycenter(init_s);
    out_simplex_map.emplace(std::make_pair(init_s,
  					   std::make_pair(init_constraint_set, seed_point)));
    queue.emplace(std::make_pair(init_s, init_constraint_set));
    Eigen::VectorXd p_shift = seed_point - barycenter;
    triangulation.change_offset(triangulation.offset() + p_shift);
    
    while (!queue.empty()) {
      Simplex_handle s;
      Constraint_set constr;
      std::tie(s, constr) = queue.front();
      queue.pop();
      if (s.dimension() != amb_d)
	for (auto cof: s.cofacet_range())
	  for (auto face: cof.facet_range()) {
	    Query_result<Simplex_handle> qr = oracle.intersects(face, constr, triangulation);
	    if (qr.success) {
	      Constraint_set face_constr(constr);
	      std::size_t excess_constraints = 0;
	      for (std::size_t I = 0; I < num_constraints; ++I)
		if (constr.find(I) == constr.end() &&
		    (*oracle.constraint_functions().at(I))(qr.intersection)(0) >= 0) {
		  face_constr.insert(I);
		  if (++excess_constraints > 1)
		    break;
		}
	      if (excess_constraints == 0 &&
		  out_simplex_map.emplace(std::make_pair(face,
							 std::make_pair(face_constr,
									qr.intersection))).second)
		queue.emplace(std::make_pair(face, constr));
	      if (excess_constraints == 1) {
		assert(face_constr.size() == constr.size() + 1);
		Query_result<Simplex_handle> qrb = oracle.intersects(cof,
								     face_constr,
								     triangulation);
		if (qrb.success &&
		    out_simplex_map.emplace(std::make_pair(cof,
							   std::make_pair(face_constr,
									  qrb.intersection))).second)
		  queue.emplace(std::make_pair(cof, face_constr));
	      }
	    }
	  }
      if (s.dimension() != cod_d) 
	for (auto facet: s.facet_range())
	  for (std::size_t I: constr) {
	    Constraint_set facet_constr(constr);
	    facet_constr.erase(I);
	    Query_result<Simplex_handle> qr = oracle.intersects(facet,
								facet_constr,
								triangulation);
	    if (qr.success &&
		(*oracle.constraint_functions().at(I))(qr.intersection)(0) < 0 &&
		out_simplex_map.emplace(std::make_pair(facet,
						       std::make_pair(facet_constr,
								      qr.intersection))).second)
	      queue.emplace(std::make_pair(facet, facet_constr));
	  }
    }
    
  }

  /** \brief Empty constructor */
  Manifold_tracing() {}
  
};

/**
 * \brief Static method for Manifold_tracing<Triangulation_>::manifold_tracing_algorithm
 * that computes the set of k-simplices that intersect
 * a boundaryless implicit manifold given by an intersection oracle, where k
 * is the codimension of the manifold.
 * The computation is based on the seed propagation --- it starts at the 
 * given seed points and then propagates along the manifold.
 *
 * \tparam Point_range Range of points of type Eigen::VectorXd.
 *  \tparam Triangulation_ The type of the ambient triangulation.
 *   Needs to be a model of the concept TriangulationForManifoldTracing.
 * \tparam Intersection_oracle Intersection oracle that represents the manifold.
 *  Needs to be a model of the concept IntersectionOracle.
 * \tparam Out_simplex_map Needs to be Manifold_tracing<Triangulation_>::Out_simplex_map.
 *
 * \param[in] seed_points The range of points on the manifold from which 
 * the computation begins.
 * \param[in] triangulation The ambient triangulation.
 * \param[in] oracle The intersection oracle for the manifold.
 * The ambient dimension needs to match the dimension of the
 * triangulation.
 * \param[out] out_simplex_map The output map, where the keys are k-simplices in
 * the input triangulation that intersect the input manifold and the mapped values 
 * are the intersection points.
 *
 * \ingroup coxeter_triangulation
 */
template <class Point_range,
	  class Triangulation,
	  class Intersection_oracle,
	  class Out_simplex_map>
void manifold_tracing_algorithm(const Point_range& seed_points,
			        Triangulation& triangulation,
				const Intersection_oracle& oracle,
				Out_simplex_map& out_simplex_map) {
  Manifold_tracing<Triangulation> mt;
  mt.manifold_tracing_algorithm(seed_points,
				triangulation,
				oracle,
				out_simplex_map);
}

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif
