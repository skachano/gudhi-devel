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

#ifndef IO_OUTPUT_MESHES_TO_OFF_H_
#define IO_OUTPUT_MESHES_TO_OFF_H_

#include <gudhi/IO/Mesh_medit.h>
#include <Eigen/Dense>
#include <fstream>

namespace Gudhi {

namespace coxeter_triangulation {

using Vertex_points = Mesh_medit::Vertex_points;
using Mesh_elements = Mesh_medit::Mesh_elements;
using Scalar_field_range = Mesh_medit::Scalar_field_range;

template <std::size_t I = 0,
	  typename... Meshes>
typename std::enable_if<I == sizeof... (Meshes), void>::type
fill_meshes (Vertex_points& vertex_points,
	     Mesh_elements& edges,
	     Mesh_elements& triangles,
	     Mesh_elements& tetrahedra,
	     Scalar_field_range& triangles_scalar_range,
	     Scalar_field_range& tetrahedra_scalar_range,
	     std::size_t index,
	     const Meshes&... meshes) {
}

template <std::size_t I = 0,
	  typename... Meshes>
typename std::enable_if<I != sizeof... (Meshes), void>::type
fill_meshes (Vertex_points& vertex_points,
	     Mesh_elements& edges,
	     Mesh_elements& triangles,
	     Mesh_elements& tetrahedra,
	     Scalar_field_range& triangles_scalar_range,
	     Scalar_field_range& tetrahedra_scalar_range,
	     std::size_t index,
	     const Meshes&... meshes) {
  auto mesh = std::get<I>(std::forward_as_tuple(meshes...));
  for (const auto& v: mesh.vertex_points)
    vertex_points.push_back(v);
  for (const auto& e: mesh.edges) {
    std::vector<std::size_t> edge;
    for (const auto& v_i: e.first)
      edge.push_back(v_i + index);
    edges.emplace_back(std::make_pair(edge, e.second));
  }
  for (const auto& t: mesh.triangles) {
    std::vector<std::size_t> triangle;
    for (const auto& v_i: t.first)
      triangle.push_back(v_i + index);
    triangles.emplace_back(std::make_pair(triangle, t.second));
  }
  for (const auto& t: mesh.tetrahedra) {
    std::vector<std::size_t> tetrahedron;
    for (const auto& v_i: t.first)
      tetrahedron.push_back(v_i + index);
    tetrahedra.emplace_back(std::make_pair(tetrahedron, t.second));
  }
  for (const auto& b: mesh.triangles_scalar_range)
    triangles_scalar_range.push_back(b);
  for (const auto& b: mesh.tetrahedra_scalar_range)
    tetrahedra_scalar_range.push_back(b);
  fill_meshes<I+1, Meshes...>(vertex_points,
			      edges,
			      triangles,
			      tetrahedra,
			      triangles_scalar_range,
			      tetrahedra_scalar_range,
			      index + mesh.vertex_points.size(),
			      meshes...);
}

/** \brief Outputs meshes in the Object File Format (OFF).
 *  
 *  @param[in] n_format If true, then the output is multidimensional (nOFF). If false, then 
 *   the output is considered to be three-dimensional, which is necessary for example in MeshLab.
 *  @param[in] file_name The name of the output file.
 *  @param[in] meshes A pack of meshes to be specified separated by commas.
 */
template <typename... Meshes>
void output_meshes_to_off(bool n_format, std::string file_name, const Meshes&... meshes) {
  Vertex_points vertex_points;
  Mesh_elements edges, triangles, tetrahedra;
  Scalar_field_range triangles_scalar_range, tetrahedra_scalar_range;
  fill_meshes(vertex_points,
	      edges,
	      triangles,
	      tetrahedra,
	      triangles_scalar_range,
	      tetrahedra_scalar_range,
	      -1,
	      meshes...);
  
  std::ofstream ofs (file_name + ".off", std::ofstream::out);

  if (n_format)
    ofs << "n";
  ofs << "OFF\n";
  if (vertex_points.empty()) {
    ofs << "0 0 0\n";
    return;
  }
  if (n_format)
    ofs << vertex_points.begin()->size() << "\n";
  ofs << vertex_points.size() << " " << triangles.size() << " 0\n";

  // vertices
  for (auto p: vertex_points) {
    if (n_format) {
      ofs << p[0];
      for (int i = 1; i < p.size(); i++)
	ofs << " " << p[i];
    }
    else
      if (p.size() < 3) {
	ofs << p[0];
	for (int i = 1; i < p.size(); i++)
	  ofs << " " << p[i];
	for (int i = p.size(); i < 3; i++)
	  ofs << " 0.0";
      }
      else
	ofs << p[0] << " " << p[1] << " " << p[2];
    ofs << "\n";
  }
  
  // // edges
  // for (auto e: edges) {
  //   ofs << "2 ";
  //   for (auto v: e.first)      
  //     ofs << v << " ";
  //   ofs << "\n";
  //   // ofs << e.second << std::endl;
  // }
  // triangles
  for (auto t: triangles) {
    ofs << "3 ";
    for (auto v: t.first)      
      ofs << v << " ";
    ofs << "\n";
    // ofs << t.second << std::endl;
  }
  
  ofs.close();
}

} // namespace coxeter_triangulation 

} // namespace Gudhi

#endif
