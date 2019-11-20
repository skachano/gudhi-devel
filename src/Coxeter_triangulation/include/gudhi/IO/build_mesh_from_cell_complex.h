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

#ifndef IO_BUILD_MESH_FROM_CELL_COMPLEX_H_
#define IO_BUILD_MESH_FROM_CELL_COMPLEX_H_

#include <gudhi/IO/Mesh_medit.h>

namespace Gudhi {

namespace coxeter_triangulation {

template <class Hasse_cell,
	  class Simplex_cell_map,
	  class Toggle_range,
	  class Ref_range>
void populate_mesh(Mesh_medit& output,
		   const Simplex_cell_map& sc_map,
		   const Toggle_range& toggle_range,
		   const Ref_range& ref_range,
		   std::size_t amb_d,
		   std::map<Hasse_cell*, std::size_t> vi_map) {
  using Mesh_element_vertices = Mesh_medit::Mesh_elements::value_type::first_type;
  // std::map<Hasse_cell*, std::size_t> ci_map;
  // std::size_t index = vi_map.size() + 1; // current size of output.vertex_points
  // if (sc_map.size() >= 3) 
  //   for (const auto& sc_pair: sc_map[2]) {
  //     Eigen::VectorXd barycenter = Eigen::VectorXd::Zero(amb_d);
  //     std::set<std::size_t> vertex_indices;
  //     Hasse_cell* cell = sc_pair.second;
  //     for (const auto& ei_pair: cell->get_boundary())
  // 	for (const auto& vi_pair: ei_pair.first->get_boundary())
  // 	  vertex_indices.emplace(vi_map[vi_pair.first]);
  //     for (const std::size_t& v: vertex_indices)
  // 	barycenter += output.vertex_points[v-1];
  //     ci_map.emplace(std::make_pair(cell, index++));
  //     output.vertex_points.emplace_back((1./vertex_indices.size()) * barycenter);    
  //   }

  // if (configuration.toggle_edges && sc_map.size() >= 2)
  //   for (const auto& sc_map: sc_map[1]) {
  //     Hasse_cell* edge_cell = sc_map.second;
  //     Mesh_element_vertices edge;
  //     for (const auto& vi_pair: edge_cell->get_boundary())
  // 	edge.push_back(vi_map[vi_pair.first]);
  //     output.edges.emplace_back(std::make_pair(edge, configuration.ref_edges));
  //   }
  
  // if (configuration.toggle_triangles && sc_map.size() >= 3)
  //   for (const auto& sc_pair: sc_map[2]) {
  //     for (const auto& ei_pair: sc_pair.second->get_boundary()) {
  // 	Mesh_element_vertices triangle(1, ci_map[sc_pair.second]);
  // 	for (const auto& vi_pair: ei_pair.first->get_boundary())
  // 	  triangle.push_back(vi_map[vi_pair.first]);
  // 	output.triangles.emplace_back(std::make_pair(triangle, configuration.ref_triangles));
  //     }
  //   }
  
  // if (configuration.toggle_tetrahedra && sc_map.size() >= 4)
  //   for (const auto& sc_pair: sc_map[3]) {
  //     Eigen::VectorXd barycenter = Eigen::VectorXd::Zero(amb_d);
  //     std::set<std::size_t> vertex_indices;
  //     Hasse_cell* cell = sc_pair.second;
  //     for (const auto& ci_pair: cell->get_boundary())
  // 	for (const auto& ei_pair: ci_pair.first->get_boundary())
  // 	  for (const auto& vi_pair: ei_pair.first->get_boundary())
  // 	    vertex_indices.emplace(vi_map[vi_pair.first]);
  //     for (const std::size_t& v: vertex_indices)
  // 	barycenter += output.vertex_points[v-1];
  //     output.vertex_points.emplace_back((1./vertex_indices.size()) * barycenter);

  //     for (const auto& ci_pair: cell->get_boundary())
  // 	for (const auto& ei_pair: ci_pair.first->get_boundary()) {
  // 	  Mesh_element_vertices tetrahedron = {index, ci_map[sc_pair.second]};
  // 	  for (const auto& vi_pair: ei_pair.first->get_boundary())
  // 	    tetrahedron.push_back(vi_map[vi_pair.first]);
  // 	  output.tetrahedra.emplace_back(std::make_pair(tetrahedron, configuration.ref_tetrahedra));
  // 	}
  //     index++;
  //   }  
}

template <class Cell_complex,
	  class Toggle_range,
	  class Ref_range>
Mesh_medit build_mesh_from_cell_complex(const Cell_complex& cell_complex,
					const Toggle_range& toggle_vector,
					const Ref_range& ref_vector) {
  using Hasse_cell = typename Cell_complex::Hasse_cell;
  Mesh_medit output;
  std::map<Hasse_cell*, std::size_t> vi_map;
  std::size_t index = 1; // current size of output.vertex_points

  if (cell_complex.cell_point_map().empty())
    return output;
  std::size_t amb_d = std::min((int) cell_complex.cell_point_map().begin()->second.size(), 3);
  
  for (const auto& cp_pair: cell_complex.cell_point_map()) {
    vi_map.emplace(std::make_pair(cp_pair.first, index++));
    output.vertex_points.push_back(cp_pair.second);
    output.vertex_points.back().conservativeResize(amb_d);
  }  

  populate_mesh(output, cell_complex.simplex_cell_maps(), toggle_vector, ref_vector, amb_d, vi_map);  
  return output;
}

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif
