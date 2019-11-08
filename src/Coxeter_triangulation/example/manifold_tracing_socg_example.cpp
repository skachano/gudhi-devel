// workaround for the annoying boost message in boost 1.69
#define BOOST_PENDING_INTEGER_LOG2_HPP
#include <boost/integer/integer_log2.hpp>
// end workaround 

#include <iostream>

#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Functions/Function_affine_plane_in_Rd.h>
#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Functions/Function_iron_in_R3.h>
#include <gudhi/Functions/Cartesian_product.h>
#include <gudhi/Functions/Linear_transformation.h>
#include <gudhi/Functions/Embed_in_Rd.h>
#include <gudhi/Implicit_manifold_intersection_oracle.h>
#include <gudhi/Manifold_tracing.h>
#include <gudhi/Cell_complex.h>
#include <gudhi/Functions/random_orthogonal_matrix.h>

#include <gudhi/Clock.h>

#include <gudhi/IO/build_mesh_from_cell_complex.h>
#include <gudhi/IO/output_meshes_to_medit.h>

using namespace Gudhi::coxeter_triangulation;

int main(int argc, char** argv) {

  std::size_t d = 6;
  if (argc > 1)
    d = atoi(argv[1]);
  double lambda = 0.2;
  if (argc > 2)
    lambda = atof(argv[2]);
  
  // Creating the iron surface
  auto fun_iron = make_embedding(Function_iron_in_R3(), d);

  // Apply a random rotation in R4
  auto matrix = random_orthogonal_matrix(d);
  auto fun_iron_rotated = make_linear_transformation(fun_iron, matrix);

  // Computing the seed of the function fun_iron_rotated
  Eigen::VectorXd seed = fun_iron_rotated.seed();    
  
  // Defining a domain function that defines the boundary, which is a hyperplane passing by the origin and orthogonal to x.
  Eigen::MatrixXd normal_matrix = Eigen::MatrixXd::Zero(4, 1);
  for (std::size_t i = 0; i < 4; i++)
    normal_matrix(i,0) = -seed(i);
  Function_affine_plane_in_Rd fun_bound(normal_matrix, -seed/2);
    
  // Defining the intersection oracle
  auto oracle = make_oracle(fun_iron_rotated, fun_bound);

  // Define a Coxeter triangulation scaled by a factor lambda.
  // The triangulation is translated by a random vector to avoid violating the genericity hypothesis.
  Coxeter_triangulation<> cox_tr(oracle.amb_d());
  cox_tr.change_offset(Eigen::VectorXd::Random(oracle.amb_d()));
  cox_tr.change_matrix(lambda * cox_tr.matrix());

  // Manifold tracing algorithm
  using MT = Manifold_tracing<Coxeter_triangulation<> >;
  using Out_simplex_map = typename MT::Out_simplex_map;
  std::vector<Eigen::VectorXd> seed_points(1, seed);
  Out_simplex_map interior_simplex_map, boundary_simplex_map;

  Gudhi::Clock mta_clock("Manifold tracing algorithm");
  manifold_tracing_algorithm(seed_points, cox_tr, oracle, interior_simplex_map, boundary_simplex_map);
  std::cout << mta_clock << std::endl;
  std::cout << "Interior vertices: " << interior_simplex_map.size() << "\n";
  std::cout << "Boundary vertices: " << boundary_simplex_map.size() << "\n";
  
  // Constructing the cell complex
  std::size_t intr_d = oracle.amb_d() - oracle.cod_d();
  Cell_complex<Out_simplex_map> cell_complex(intr_d);
  cell_complex.construct_complex(interior_simplex_map, boundary_simplex_map);

  // Output the cell complex to a file readable by medit
  output_meshes_to_medit(3, 
			 "iron_with_boundary",
			 build_mesh_from_cell_complex(cell_complex,
						      Configuration(false, true, true, 1, 1, 2),
						      Configuration(true, true, true, 6, 13, 14)));

  return 0;
}
