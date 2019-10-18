#include <iostream>

#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Implicit_manifold_intersection_oracle.h>
#include <gudhi/Manifold_tracing.h>
#include <gudhi/Cell_complex.h>
#include <gudhi/Functions/random_orthogonal_matrix.h>
#include <gudhi/Functions/Linear_transformation.h>

#include <gudhi/IO/build_mesh_from_cell_complex.h>
#include <gudhi/IO/output_meshes_to_medit.h>

using namespace Gudhi::coxeter_triangulation;

/* A definition of a function that defines a 2d surface embedded in R^4, but that normally
 * lives on a complex projective plane.
 * In terms of harmonic coordinates [x:y:z] of points on the complex projective plane,
 * the equation of the manifold is x^3*y + y^3*z + z^3*x = 0.
 * The embedding consists of restricting the manifold to the affine subspace z = 1.
 */
struct Function_custom_function : public Function {

  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    // The real and imaginary parts of the variables x and y
    double x = p(0), y = p(1), z = p(2), w = p(3);
    Eigen::VectorXd result(cod_d());

    result(0) = std::pow(std::sqrt(x*x + y*y) - R, 2) + z*z + w*w - r*r;
    result(1) = y/x - 2*(w/z)/(1+(w/z)*(w/z));
    return result;
  }

  std::size_t amb_d() const {return 4;};
  std::size_t cod_d() const {return 2;};

  Eigen::VectorXd seed() const {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(amb_d());
    result(0) = R;
    result(2) = r;
    return result;
  }

  Function_custom_function() {}
private:
  const double R = 1;
  const double r = 1;
};

int main(int argc, char** argv) {

  // The function for the (non-compact) manifold
  Function_custom_function fun_custom;
  Eigen::MatrixXd matrix = random_orthogonal_matrix(fun_custom.amb_d());
  auto fun = make_linear_transformation(fun_custom, matrix);

  // Seed of the function
  Eigen::VectorXd seed = fun.seed();

  // Defining the intersection oracle
  auto oracle = make_oracle(fun);

  // Define a Coxeter triangulation scaled by a factor lambda.
  // The triangulation is translated by a random vector to avoid violating the genericity hypothesis.
  double lambda = 0.1;
  if (argc > 1)
    lambda = atof(argv[1]);
  Coxeter_triangulation<> cox_tr(oracle.amb_d());
  cox_tr.change_offset(Eigen::VectorXd::Random(oracle.amb_d()));
  cox_tr.change_matrix(lambda * cox_tr.matrix());
  
  // Manifold tracing algorithm
  using MT = Manifold_tracing<Coxeter_triangulation<> >;
  using Out_simplex_map = typename MT::Out_simplex_map;
  std::vector<Eigen::VectorXd> seed_points(1, seed);
  Out_simplex_map interior_simplex_map;
  manifold_tracing_algorithm(seed_points, cox_tr, oracle, interior_simplex_map);
  std::cout << "Output size = " << interior_simplex_map.size() << "\n";
  
  // Constructing the cell complex
  std::size_t intr_d = oracle.amb_d() - oracle.cod_d();
  Cell_complex<Out_simplex_map> cell_complex(intr_d);
  cell_complex.construct_complex(interior_simplex_map);

  // Output the cell complex to a file readable by medit
  output_meshes_to_medit(3,
			 "custom_manifold",
			 build_mesh_from_cell_complex(cell_complex,
						      Configuration(true, true, true, 1, 5, 3),
						      Configuration(true, true, true, 2, 13, 14)));
  return 0;
}
