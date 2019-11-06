#include <iostream>

#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Functions/Function_affine_plane_in_Rd.h>
#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Implicit_manifold_intersection_oracle_corners.h>
#include <gudhi/Manifold_tracing.h>
#include <gudhi/Cell_complex_corners.h>
#include <gudhi/Functions/random_orthogonal_matrix.h>
#include <gudhi/Functions/Linear_transformation.h>

#include <gudhi/IO/build_mesh_from_cell_complex_corners.h>
#include <gudhi/IO/output_meshes_to_medit.h>

using namespace Gudhi::coxeter_triangulation;

/* A definition of a function that defines a 2d surface embedded in R^4, but that normally
 * lives on a complex projective plane.
 * In terms of harmonic coordinates [x:y:z] of points on the complex projective plane,
 * the equation of the manifold is x^3*y + y^3*z + z^3*x = 0.
 * The embedding consists of restricting the manifold to the affine subspace z = 1.
 */
// struct Function_custom_function : public Function {

//   Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
//     // The real and imaginary parts of the variables x and y
//     double x = p(0), y = p(1), z = p(2), w = p(3);
//     Eigen::VectorXd result(cod_d());

//     result(0) = std::pow(std::sqrt(x*x + y*y) - R, 2) + z*z + w*w - r*r;
//     // result(1) = y/x - 2*(w/z)/(1+(w/z)*(w/z));
//     result(1) = y - w*std::sin(w);
//     return result;
//   }

//   std::size_t amb_d() const {return 4;};
//   std::size_t cod_d() const {return 2;};

//   Eigen::VectorXd seed() const {
//     Eigen::VectorXd result = Eigen::VectorXd::Zero(amb_d());
//     result(0) = R;
//     result(2) = r;
//     return result;
//   }

//   Function_custom_function() {}
// private:
//   const double R = 1;
//   const double r = 1;
// };

struct Function_x_abs : public Function {

  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    double x = p(0);
    Eigen::VectorXd result(cod_d());
    result(0) = std::abs(x) - 1;
    return result;
  }

  std::size_t amb_d() const {return 3;};
  std::size_t cod_d() const {return 1;};

  Eigen::VectorXd seed() const {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(amb_d());
    result(0) = 1;
    return result;
  }

  Function_x_abs() {}
};

struct Function_y_abs : public Function {

  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    double y = p(1);
    Eigen::VectorXd result(cod_d());
    result(0) = std::abs(y) - 1;
    return result;
  }

  std::size_t amb_d() const {return 3;};
  std::size_t cod_d() const {return 1;};

  Eigen::VectorXd seed() const {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(amb_d());
    result(0) = 1;
    return result;
  }

  Function_y_abs() {}
};


int main(int argc, char** argv) {

  // The function for the (non-compact) manifold
  // Function_custom_function fun_custom;
  // Eigen::MatrixXd matrix = random_orthogonal_matrix(fun_custom.amb_d());
  // auto fun = make_linear_transformation(fun_custom, matrix);

  // // Seed of the function
  // Eigen::VectorXd seed = fun.seed();
  
  // // Defining the intersection oracle
  // auto oracle = make_oracle(fun);

  // Test for the constant function and the MTA
  std::size_t d = 3;
  std::size_t k = 1;
  Eigen::MatrixXd normal_matrix = Eigen::MatrixXd::Zero(d,k);
  normal_matrix << 0, 0, 1;
  Function_affine_plane_in_Rd fun(normal_matrix);
  Eigen::VectorXd seed = Eigen::VectorXd::Zero(d);
  Function_x_abs fun_x;
  Function_y_abs fun_y;
  auto oracle = make_oracle(fun, fun_x, fun_y);

  // Define a Coxeter triangulation scaled by a factor lambda.
  // The triangulation is translated by a random vector to avoid violating the genericity hypothesis.
  double lambda = 0.1;
  if (argc > 1)
    lambda = atof(argv[1]);
  Eigen::MatrixXd rot_matrix = random_orthogonal_matrix(d);
  Coxeter_triangulation<> cox_tr(oracle.amb_d());
  cox_tr.change_offset(Eigen::VectorXd::Random(oracle.amb_d()));
  cox_tr.change_matrix(lambda * rot_matrix * cox_tr.matrix());
  
  // Manifold tracing algorithm
  using MT = Manifold_tracing<Coxeter_triangulation<> >;
  using Out_simplex_map = typename MT::Out_simplex_map;
  std::vector<Eigen::VectorXd> seed_points(1, seed);
  Out_simplex_map
    interior_simplex_map,
    boundary1_simplex_map,
    boundary2_simplex_map,
    corner_simplex_map;
  manifold_tracing_algorithm(seed_points, cox_tr, oracle,
			     interior_simplex_map,
			     boundary1_simplex_map,
			     boundary2_simplex_map,
			     corner_simplex_map);
  // std::cout << "Output size (interior) = " << interior_simplex_map.size() << "\n";
  // for (auto& m_pair: interior_simplex_map) 
  //   std::cout << " " << m_pair.first << "\n";
  // std::cout << "Output size (boundary1) = " << boundary1_simplex_map.size() << "\n";
  // for (auto& m_pair: boundary1_simplex_map)
  //   std::cout << " " << m_pair.first << "\n";
  // std::cout << "Output size (boundary2) = " << boundary2_simplex_map.size() << "\n";
  // for (auto& m_pair: boundary2_simplex_map)
  //   std::cout << " " << m_pair.first << "\n";
  // std::cout << "Output size (corner) = " << corner_simplex_map.size() << "\n";
  // for (auto& m_pair: corner_simplex_map)
  //   std::cout << " " << m_pair.first << "\n";
  
  // Constructing the cell complex
  std::size_t intr_d = oracle.amb_d() - oracle.cod_d();
  Cell_complex<Out_simplex_map> cell_complex(intr_d);
  cell_complex.construct_complex(interior_simplex_map,
				 boundary1_simplex_map,
				 boundary2_simplex_map,
				 corner_simplex_map);
  for (std::size_t i = 0; i < cell_complex.interior_simplex_cell_maps().size(); i++)
    std::cout << "Cells (interior[" << i << "]) = "
	      << cell_complex.interior_simplex_cell_map(i).size() << "\n";
  for (std::size_t i = 0; i < cell_complex.boundary1_simplex_cell_maps().size(); i++)
    std::cout << "Cells (boundary1[" << i << "]) = "
	      << cell_complex.boundary1_simplex_cell_map(i).size() << "\n";
  for (std::size_t i = 0; i < cell_complex.boundary2_simplex_cell_maps().size(); i++)
    std::cout << "Cells (boundary2[" << i << "]) = "
	      << cell_complex.boundary2_simplex_cell_map(i).size() << "\n";
  for (std::size_t i = 0; i < cell_complex.corner_simplex_cell_maps().size(); i++)
    std::cout << "Cells (corner[" << i << "]) = "
	      << cell_complex.corner_simplex_cell_map(i).size() << "\n";
  std::cout << "Euler characteristic = " <<
    cell_complex.interior_simplex_cell_map(0).size()
    + cell_complex.boundary1_simplex_cell_map(0).size()
    + cell_complex.boundary2_simplex_cell_map(0).size()
    + cell_complex.corner_simplex_cell_map(0).size()
    - cell_complex.interior_simplex_cell_map(1).size()
    - cell_complex.boundary1_simplex_cell_map(1).size()
    - cell_complex.boundary2_simplex_cell_map(1).size()
    + cell_complex.interior_simplex_cell_map(2).size() << "\n";

  // Output the cell complex to a file readable by medit
  output_meshes_to_medit(3,
  			 "custom_manifold",
  			 build_mesh_from_cell_complex(cell_complex,
  						      Configuration(true, true, true, 1, 5, 3),
  						      Configuration(true, true, true, 2, 13, 14)));
  return 0;
}
