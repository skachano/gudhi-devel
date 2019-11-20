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

#ifndef CELL_COMPLEX_H_
#define CELL_COMPLEX_H_

#include <Eigen/Dense> 

#include <algorithm>
#include <gudhi/Permutahedral_representation/Simplex_comparator.h>
#include <fstream> // for Hasse_diagram_persistence.h

#include <gudhi/Cell_complex/Hasse_diagram_cell.h> // for Hasse_cell

namespace Gudhi {

namespace coxeter_triangulation {

/**
 *  \ingroup coxeter_triangulation
 */

/** \class Cell_complex 
 *  \brief A class that constructs the cell complex from the output provided by the class 
 *   Gudhi::Manifold_tracing<Triangulation_>.
 *
 *  \tparam Out_simplex_map_ The type of a map from a simplex type that is a
 *   model of SimplexInCoxeterTriangulation to Eigen::VectorXd.
 *
 *  \ingroup coxeter_triangulation
 */
template <class Out_simplex_map_>
class Cell_complex {
public:
  
  /** \brief Type of a simplex in the ambient triangulation.
   *  Is a model of the concept SimplexInCoxeterTriangulation.
   */
  typedef typename Out_simplex_map_::key_type Simplex_handle;

  typedef typename Out_simplex_map_::mapped_type::first_type Constraint_set;

private:
  typedef typename Out_simplex_map_::hasher Simplex_hash;
  struct Pair_hash {
    typedef std::pair<Simplex_handle, Constraint_set> argument_type;
    typedef std::size_t result_type;
    result_type operator()(const argument_type& s) const noexcept {
      return Simplex_hash()(s.first);
    }
  };
  
public:
  /** \brief Type of a cell in the cell complex. 
   *  Always is Gudhi::Hasse_cell from the Hasse diagram module.
   *  The additional information is the boolean that is true if and only if the cell lies
   *  on the boundary.
   */
  typedef Gudhi::Hasse_diagram::Hasse_diagram_cell<int, double, int> Hasse_cell;
  /** \brief Type of a map from permutahedral representations of simplices in the
   *  ambient triangulation to the corresponding cells in the cell complex of some
   *  specific dimension.
   */
  typedef std::unordered_map<std::pair<Simplex_handle, Constraint_set>,
			     Hasse_cell*,
			     Pair_hash> Simplex_cell_map;
  /** \brief Type of a vector of maps from permutahedral representations of simplices in the
   *  ambient triangulation to the corresponding cells in the cell complex of various dimensions.
   */
  typedef std::vector<Simplex_cell_map> Simplex_cell_maps;
  
  /** \brief Type of a map from cells in the cell complex to the permutahedral representations
   *  of the corresponding simplices in the ambient triangulation.
   */
  typedef std::map<Hasse_cell*, std::pair<Simplex_handle, Constraint_set> > Cell_simplex_map;

  /** \brief Type of a map from vertex cells in the cell complex to the permutahedral representations
   *  of their Cartesian coordinates.
   */
  typedef std::map<Hasse_cell*, Eigen::VectorXd> Cell_point_map;

private:
    
  Hasse_cell* insert_cell(const Simplex_handle& simplex,
			  const Constraint_set& constr_set,
			  std::size_t cell_d) {
    Simplex_cell_map& simplex_cell_map = simplex_cell_maps_[cell_d];
    auto pair = std::make_pair(simplex, constr_set);
    auto map_it = simplex_cell_map.find(pair);
    if (map_it == simplex_cell_map.end()) {
      hasse_cells_.push_back(new Hasse_cell(cell_d));
      Hasse_cell* new_cell = hasse_cells_.back();    
      simplex_cell_map.emplace(std::make_pair(pair, new_cell));
      return new_cell;
    }
    return map_it->second;
  }

  void expand_level(std::size_t cell_d) {
    for (auto& sc_pair: simplex_cell_maps_[cell_d - 1]) {
      const Simplex_handle& simplex = sc_pair.first.first;
      const Constraint_set& constr_set = sc_pair.first.second;
      Hasse_cell* cell = sc_pair.second;
      std::size_t amb_d = simplex.vertex().size();
      if (simplex.dimension() < amb_d)
	for (Simplex_handle cofacet: simplex.cofacet_range()) {
	  Hasse_cell* new_cell = insert_cell(cofacet, constr_set, cell_d);
	  new_cell->get_boundary().emplace_back(std::make_pair(cell, 1));
	}
      for (std::size_t I: constr_set) {
	Constraint_set new_constr_set(constr_set);
        new_constr_set.erase(I);
	Hasse_cell* new_cell = insert_cell(simplex, new_constr_set, cell_d);
	new_cell->get_boundary().emplace_back(std::make_pair(cell, 1));
      }
    }
  }
  
  void construct_complex_(const Out_simplex_map_& out_simplex_map) {
    if (!out_simplex_map.empty()) {
      const Simplex_handle& simplex = out_simplex_map.begin()->first;
      const Constraint_set& constr_set = out_simplex_map.begin()->second.first;
      cod_d_ = simplex.dimension() - constr_set.size();
    }
    for (auto& os_pair: out_simplex_map) {
      const Simplex_handle& simplex = os_pair.first;
      const Constraint_set& constr_set = os_pair.second.first;
      const Eigen::VectorXd& point = os_pair.second.second;
      Hasse_cell* new_cell = insert_cell(simplex, constr_set, 0);
      cell_point_map_.emplace(std::make_pair(new_cell, point));
    }
    for (std::size_t cell_d = 1; cell_d < simplex_cell_maps_.size(); ++cell_d)
      expand_level(cell_d);
  }
  
public:  
  
  /**
   * \brief Constructs the the cell complex that approximates an \f$m\f$-dimensional manifold
   *  without boundary embedded in the \f$ d \f$-dimensional Euclidean space
   *  from the output of the class Gudhi::Manifold_tracing.
   *
   * \param[in] out_simplex_map A map from simplices of dimension \f$(d-m)\f$ 
   * in the ambient triangulation that intersect the relative interior of the manifold 
   * to the intersection points.
   */
  void construct_complex(const Out_simplex_map_& out_simplex_map) {
    simplex_cell_maps_.resize(intr_d_ + 1);
    construct_complex_(out_simplex_map);
  }
  
  /**
   * \brief Constructs the skeleton of the cell complex that approximates
   *  an \f$m\f$-dimensional manifold without boundary embedded 
   *  in the \f$d\f$-dimensional Euclidean space
   *  up to a limit dimension from the output of the class Gudhi::Manifold_tracing.
   *
   * \param[in] out_simplex_map A map from simplices of dimension \f$(d-m)\f$ 
   * in the ambient triangulation that intersect the relative interior of the manifold 
   * to the intersection points.
   * \param[in] limit_dimension The dimension of the constructed skeleton.
   */
  void construct_complex(const Out_simplex_map_& out_simplex_map,
			 std::size_t limit_dimension) {
    simplex_cell_maps_.resize(limit_dimension + 1);
    construct_complex_(out_simplex_map);
  }

  /**
   * \brief Returns the dimension of the cell complex.
   */
  std::size_t intrinsic_dimension() const {
    return intr_d_;
  }

  /**
   * \brief Returns a vector of maps from the cells of various dimensions in the interior
   *  of the cell complex of type Gudhi::Hasse_cell to the permutahedral representations 
   *  of the corresponding simplices in the ambient triangulation.
   */
  const Simplex_cell_maps& simplex_cell_maps() const {
    return simplex_cell_maps_;
  }
  
  /**
   * \brief Returns a map from the cells of a given dimension in the interior 
   *  of the cell complex of type Gudhi::Hasse_cell to the permutahedral representations 
   *  of the corresponding simplices in the ambient triangulation.
   *
   * \param[in] cell_d The dimension of the cells.
   */
  const Simplex_cell_map& simplex_cell_map(std::size_t cell_d) const {
    return simplex_cell_maps_[cell_d];
  }

  /**
   * \brief Returns a map from the vertex cells in the cell complex of type Gudhi::Hasse_cell
   *  to their Cartesian coordinates.
   */
  const Cell_point_map& cell_point_map() const {
    return cell_point_map_;
  }

  /**
   * \brief Conxtructor for the class Cell_complex.
   *
   * \param[in] intrinsic_dimension The dimension of the cell complex.
   */
  Cell_complex(std::size_t intrinsic_dimension)
    : intr_d_(intrinsic_dimension) {}

  ~Cell_complex() {
    for (Hasse_cell* hs_ptr: hasse_cells_)
      delete hs_ptr;
  }
  
private:
  std::size_t intr_d_, cod_d_;
  Simplex_cell_maps simplex_cell_maps_;
  Cell_point_map cell_point_map_;
  std::vector<Hasse_cell*> hasse_cells_;
};

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif
