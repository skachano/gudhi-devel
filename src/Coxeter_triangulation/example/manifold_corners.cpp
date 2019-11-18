#include <iostream>
#include <complex>

#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Functions/Function_affine_plane_in_Rd.h>
#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Constraint_set.h>
#include <gudhi/Implicit_manifold_intersection_oracle.h>
#include <gudhi/Manifold_tracing.h>
#include <gudhi/Cell_complex_corners.h>
#include <gudhi/Functions/random_orthogonal_matrix.h>
#include <gudhi/Functions/Linear_transformation.h>
#include <gudhi/Functions/Cartesian_product.h>
#include <gudhi/Clock.h>

#include <gudhi/IO/build_mesh_from_cell_complex_corners.h>
// #include <gudhi/IO/output_meshes_to_medit.h>
#include <gudhi/IO/output_meshes_to_off.h>

using namespace Gudhi::coxeter_triangulation;

int main(int argc, char** argv) {
  // Creating a circle S1 in R2 of specified radius
  double radius = 1.0;
  Function_Sm_in_Rd fun_circle(radius, 1);

  // Creating a flat torus S1xS1 in R4 from two circle functions
  auto fun_flat_torus = make_product_function(fun_circle, fun_circle);

  // Apply a random rotation in R4
  auto matrix = random_orthogonal_matrix(4);
  auto fun_flat_torus_rotated = make_linear_transformation(fun_flat_torus, matrix);

  // Computing the seed of the function fun_flat_torus
  Eigen::VectorXd seed = fun_flat_torus_rotated.seed();    
  
  // Defining a domain function that defines the boundary, which is a hyperplane passing by the origin and orthogonal to x.
  Eigen::MatrixXd normal_matrix = Eigen::MatrixXd::Zero(4, 1);
  for (std::size_t i = 0; i < 4; i++)
    normal_matrix(i,0) = -seed(i);
  Function_affine_plane_in_Rd fun_bound(normal_matrix, -seed/2);
    
  // Defining the intersection oracle
  auto oracle = make_oracle(fun_flat_torus_rotated, fun_bound);

  auto oracle_function = oracle.function();
  std::cout << "oracle_function.amb_d() = " << oracle_function.amb_d() << "\n";
  std::cout << "oracle_function.cod_d() = " << oracle_function.cod_d() << "\n";
  std::cout << "oracle_function.seed() =\n" << oracle_function.seed() << "\n";
  auto oracle_domain_function = oracle.domain_function<0>();
  std::cout << "oracle_domain_function.amb_d() = " << oracle_domain_function.amb_d() << "\n";
  std::cout << "oracle_domain_function.cod_d() = " << oracle_domain_function.cod_d() << "\n";
  std::cout << "oracle_domain_function.seed() =\n" << oracle_domain_function.seed() << "\n";
}
