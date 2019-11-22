#include <iostream>
#include <complex>
#include <unordered_set>

#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Functions/Constant_function.h>
#include <gudhi/Functions/Function_affine_plane_in_Rd.h>
#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Constraint_set.h>
#include <gudhi/Implicit_manifold_intersection_oracle.h>
#include <gudhi/Manifold_tracing.h>
#include <gudhi/Cell_complex.h>
#include <gudhi/Functions/random_orthogonal_matrix.h>
#include <gudhi/Functions/Linear_transformation.h>
#include <gudhi/Functions/Cartesian_product.h>
#include <gudhi/Clock.h>

#include <gudhi/IO/build_mesh_from_cell_complex.h>
#include <gudhi/IO/output_meshes_to_medit.h>

using namespace Gudhi::coxeter_triangulation;

int main(int argc, char** argv) {
  // Creating a circle S1 in R2 of specified radius
  std::size_t d = 2;
  std::size_t k = 0;
  Constant_function fun(d, k, Eigen::VectorXd::Zero(0));  
  Eigen::VectorXd seed = Eigen::VectorXd::Zero(d);
  Function_Sm_in_Rd fun_bound(1.0, 1);
  Eigen::MatrixXd normal_matrix_2 = Eigen::MatrixXd::Zero(2, 1);
  normal_matrix_2(1, 0) = 1.;
  Function_affine_plane_in_Rd fun_bound2(normal_matrix_2, Eigen::Vector2d(0, 0.5));

  Eigen::MatrixXd normal_matrix_3 = Eigen::MatrixXd::Zero(2, 1);
  normal_matrix_3(0, 0) = 1.;
  Function_affine_plane_in_Rd fun_bound3(normal_matrix_3, Eigen::Vector2d(0.5, 0));

  std::vector<Function*> constraint_functions_ = {&fun_bound, &fun_bound2, &fun_bound3};
  std::cout << "fun1(seed) =\n" << fun_bound(seed) << "\n"; 
  std::cout << "fun2(seed) =\n" << fun_bound2(seed) << "\n"; 
  
  // // Creating a flat torus S1xS1 in R4 from two circle functions
  // auto fun_flat_torus = make_product_function(fun_circle, fun_circle);

  // // Apply a random rotation in R4
  // auto matrix = random_orthogonal_matrix(4);
  // auto fun_flat_torus_rotated = make_linear_transformation(fun_flat_torus, matrix);

  // // Computing the seed of the function fun_flat_torus
  // Eigen::VectorXd seed = fun_flat_torus_rotated.seed();    

  // // Defining a domain function that defines the boundary, which is a hyperplane passing by the origin and orthogonal to x.
  // Eigen::MatrixXd normal_matrix = Eigen::MatrixXd::Zero(4, 1);
  // for (std::size_t i = 0; i < 4; i++)
  //   normal_matrix(i,0) = -seed(i);
  // Function_affine_plane_in_Rd fun_bound(normal_matrix, -seed/2);

  
  // Defining the intersection oracle
  auto oracle = make_oracle(&fun, constraint_functions_);
  // Function* fun_ptr = oracle.function();
  // std::cout << fun_ptr->amb_d() << " " << fun_ptr->cod_d() << "\n" << fun_ptr->seed() << "\n";
  // fun = oracle.constraint_functions().at(0);
  // std::cout << fun_ptr->amb_d() << " " << fun_ptr->cod_d() << "\n" << fun_ptr->seed() << "\n";
  
  // Define a Coxeter triangulation scaled by a factor lambda.
  // The triangulation is translated by a random vector to avoid violating the genericity hypothesis.
  double lambda = 0.6;
  Coxeter_triangulation<> cox_tr(oracle.amb_d());
  cox_tr.change_offset(Eigen::VectorXd::Random(oracle.amb_d()));
  cox_tr.change_matrix(lambda * cox_tr.matrix());

  // Manifold tracing algorithm
  using MT = Manifold_tracing<Coxeter_triangulation<> >;
  using Out_simplex_map = typename MT::Out_simplex_map;
  Out_simplex_map out_simplex_map;
  manifold_tracing_algorithm(seed, cox_tr, oracle, out_simplex_map);

  std::cout << "Output size = " << out_simplex_map.size() << "\n";
  // for (const auto& o_pair: out_simplex_map) {
  //   std::cout << o_pair.first << ": ";
  //   std::cout << "{";
  //   auto o_it = o_pair.second.first.begin();
  //   if (!o_pair.second.first.empty())
  //     std::cout << *o_it++;
  //   while (o_it != o_pair.second.first.end())
  //     std::cout << ", " << *o_it++;
  //   std::cout << "},\n" << o_pair.second.second << "\n";
  // }
  std::vector<std::size_t> dim_lists(fun.amb_d() + 1);
  for (const auto& o_pair: out_simplex_map) {
    dim_lists[o_pair.first.first.dimension()]++;
  }
  for (std::size_t i = 0; i < dim_lists.size(); ++i)
    std::cout << i << ": " << dim_lists[i] << "\n";
  Cell_complex<Out_simplex_map> cell_complex(fun.amb_d() - fun.cod_d());
  cell_complex.construct_complex(out_simplex_map);

  int chi = 0;
  for (std::size_t i = 0; i < cell_complex.simplex_cell_maps().size(); ++i) {
    std::cout << " Cells of dim " << i << ": " << cell_complex.simplex_cell_map(i).size() << "\n";
    chi += (2*((i+1)%2)-1)*cell_complex.simplex_cell_map(i).size();
  }
  std::cout << "Euler characteristic: " << chi << "\n";

  std::vector<std::vector<bool> >
    toggle_vectors(4, std::vector<bool>(std::pow(2, constraint_functions_.size()), true));
  std::vector<std::vector<std::size_t> >
    ref_vectors(4, std::vector<std::size_t>(std::pow(2, constraint_functions_.size()), 1));
  output_meshes_to_medit(2,
  			 "corner_manifold",
  			 build_mesh_from_cell_complex(cell_complex,
  						      toggle_vectors,
  						      ref_vectors));

}
