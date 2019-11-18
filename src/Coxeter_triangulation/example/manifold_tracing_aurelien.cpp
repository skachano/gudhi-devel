#include <iostream>
#include <complex>

#include <gudhi/Coxeter_triangulation.h>
#include <gudhi/Functions/Function_affine_plane_in_Rd.h>
#include <gudhi/Functions/Function_Sm_in_Rd.h>
#include <gudhi/Implicit_manifold_intersection_oracle_corners.h>
#include <gudhi/Manifold_tracing.h>
#include <gudhi/Cell_complex_corners.h>
#include <gudhi/Functions/random_orthogonal_matrix.h>
#include <gudhi/Functions/Linear_transformation.h>
#include <gudhi/Clock.h>

#include <gudhi/IO/build_mesh_from_cell_complex_corners.h>
// #include <gudhi/IO/output_meshes_to_medit.h>
#include <gudhi/IO/output_meshes_to_off.h>

using namespace Gudhi::coxeter_triangulation;

int binomialCoeff(int n, int k) 
{ 
  int res = 1; 
  
    // Since C(n, k) = C(n, n-k) 
  if ( k > n - k ) 
    k = n - k; 
  
  // Calculate value of [n * (n-1) *---* (n-k+1)] / [k * (k-1) *----* 1] 
  for (int i = 0; i < k; ++i) { 
    res *= (n - i); 
    res /= (i + 1); 
  }   
  return res; 
} 

/* A definition of a function that defines a 2d surface embedded in R^4, but that normally
 * lives on a complex projective plane.
 * In terms of harmonic coordinates [x:y:z] of points on the complex projective plane,
 * the equation of the manifold is x^3*y + y^3*z + z^3*x = 0.
 * The embedding consists of restricting the manifold to the affine subspace z = 1.
 */
struct Function_aurelien_x3y_y3z_z3x : public Function {

  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    // The real and imaginary parts of the variables x and y
    double xr = p(0), xi = p(1), yr = p(2), yi = p(3);
    Eigen::VectorXd result(cod_d());
    
    // Squares and cubes of real and imaginary parts used in the computations
    double
      xr2 = xr*xr, xi2 = xi*xi, yr2 = yr*yr, yi2 = yi*yi,
      xr3 = xr2*xr, xi3 = xi2*xi, yr3 = yr2*yr, yi3 = yi2*yi;

    // The first coordinate of the output is Re(x^3*y + y^3 + x)
    result(0) = 
      xr3*yr - 3*xr*xi2*yr - 3*xr2*xi*yi + xi3*yi
      + yr3 - 3*yr*yi2 + xr;
    // The second coordinate of the output is Im(x^3*y + y^3 + x)
    result(1) =
      3*xr2*xi*yr + xr3*yi - 3*xr*xi2*yi - xi3*yr
      + 3*yr2*yi - yi3 + xi;
    return result;
  }
  
  std::size_t amb_d() const {return 4;};
  std::size_t cod_d() const {return 2;};

  Eigen::VectorXd seed() const {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(4);
    // Real solutions to x^3*y + y^3 + x = 0
    result(0) = -0.124045635085; result(2) = 0.5;
    return result;
  }

  Function_aurelien_x3y_y3z_z3x() {}
private:
  const double R = 1;
  const double r = 1;
};

struct Function_aurelien_x2Y_y2Z_z2X : public Function {
// Cells (interior[0]) = 6935
// Cells (interior[1]) = 16791
// Cells (interior[2]) = 9857
// Cells (boundary1[0]) = 377
// Cells (boundary1[1]) = 384
// Cells (boundary2[0]) = 423
// Cells (boundary2[1]) = 430
// Cells (corner[0]) = 14
// Euler characteristic (interior) = 1
// Euler characteristic (boundary) = -14
// Euler characteristic (corners) = 14
// Euler characteristic (whole surface) = -4
  
  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    // The real and imaginary parts of the variables x and y
    double xr = p(0), xi = p(1), yr = p(2), yi = p(3);
    Eigen::VectorXd result(cod_d());
    
    // Squares and cubes of real and imaginary parts used in the computations
    double
      xr2 = xr*xr, xi2 = xi*xi, yr2 = yr*yr, yi2 = yi*yi;

    // The first coordinate of the output is Re(x^2*'y + y^2 + 'x)
    result(0) = 
      xr2*yr - xi2*yr + 2*xr*xi*yi
      + yr2 - yi2 + xr;
    // The second coordinate of the output is Im(x^2*'y + y^2 + 'x)
    result(1) =
      2*xr*xi*yr - xr2*yi + xi2*yi
      + 2*yr*yi - xi;
    return result;
  }
  
  std::size_t amb_d() const {return 4;};
  std::size_t cod_d() const {return 2;};

  Eigen::VectorXd seed() const {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(4);
    // Real solutions to x^2*y + y^2 + x = 0
    // 0.5*x^2 + x + 0.25 = 0
    // x^2 + 2*x + 0.5 = 0
    // D = 1 - 0.5 = 0.5
    // x = -1 + sqrt(0.5)
    // 1 + 0.5 - 2*sqrt(0.5) - 2 + 2*sqrt(0.5) + 0.5 = 0
    result(0) = -1 + std::sqrt(0.5); result(2) = 0.5;
    return result;
  }

  Function_aurelien_x2Y_y2Z_z2X() {}
private:
  const double R = 1;
  const double r = 1;
};

struct Function_aurelien_xY_yZ_zX : public Function {
// Cells (interior[0]) = 4112
// Cells (interior[1]) = 9937
// Cells (interior[2]) = 5826
// Cells (boundary1[0]) = 222
// Cells (boundary1[1]) = 225
// Cells (boundary2[0]) = 222
// Cells (boundary2[1]) = 225
// Cells (corner[0]) = 6
// Euler characteristic (interior) = 1
// Euler characteristic (boundary) = -6
// Euler characteristic (corners) = 6
// Euler characteristic (whole surface) = 0
  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    // The real and imaginary parts of the variables x and y
    double xr = p(0), xi = p(1), yr = p(2), yi = p(3);
    Eigen::VectorXd result(cod_d());
    
    // The first coordinate of the output is Re(x*'y + y + 'x)
    result(0) = 
      xr*yr + xi*yi
      + yr + xr;
    // The second coordinate of the output is Im(x*'y + y + 'x)
    result(1) =
      xi*yr - xr*yi
      + yi - xi;
    return result;
  }
  
  std::size_t amb_d() const {return 4;};
  std::size_t cod_d() const {return 2;};

  Eigen::VectorXd seed() const {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(4);
    // Real solutions to x*y + y + x = 0
    // 1.5*x + 0.5 = 0
    // x = -1/3
    result(0) = -1./3; result(2) = 0.5;
    return result;
  }

  Function_aurelien_xY_yZ_zX() {}
private:
  const double R = 1;
  const double r = 1;
};

struct Function_aurelien_xYd_yZd_zXd : public Function {
  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    // The real and imaginary parts of the variables x and y
    std::complex<double> x(p(0), p(1)), y(p(2), p(3));
    std::complex<double> res_complex = x * std::pow(std::conj(y), d_) + y + std::pow(std::conj(x), d_);
    Eigen::VectorXd result(cod_d());
    result << std::real(res_complex), std::imag(res_complex);
    // The first coordinate of the output is Re(x*'y^d + y + 'x^d)
    // for (std::size_t i = 0, sign = 1; i < d_; i += 2, sign = -sign)
    //   result(0) += sign * xr * binomialCoeff(d_, i) * std::pow(-yi, i) * std::pow(yr, d_ - i);
    // for (std::size_t i = 1, sign = -1; i < d_; i += 2, sign = -sign)
    //   result(0) += sign * xi * binomialCoeff(d_, i) * std::pow(-yi, i) * std::pow(yr, d_ - i);
    // result(0) += yr;
    // for (std::size_t i = 0, sign = 1; i < d_; i += 2, sign = -sign)
    //   result(0) += sign * binomialCoeff(d_, i) * std::pow(-xi, i) * std::pow(xr, d_ - i);

    // The second coordinate of the output is Im(x*'y^d + y + 'x^d)
    // for (std::size_t i = 1, sign = 1; i < d_; i += 2, sign = -sign)
    //   result(1) += sign * xr * binomialCoeff(d_, i) * std::pow(-yi, i) * std::pow(yr, d_ - i);
    // for (std::size_t i = 0, sign = 1; i < d_; i += 2, sign = -sign)
    //   result(1) += sign * xi * binomialCoeff(d_, i) * std::pow(-yi, i) * std::pow(yr, d_ - i);
    // result(1) += yi;
    // for (std::size_t i = 1, sign = 1; i < d_; i += 2, sign = -sign)
    //   result(1) += sign * binomialCoeff(d_, i) * std::pow(-xi, i) * std::pow(xr, d_ - i);
    return result;
  }
  
  std::size_t amb_d() const {return 4;};
  std::size_t cod_d() const {return 2;};

  Eigen::VectorXd seed() const {
    return seed_;
  }

  Function_aurelien_xYd_yZd_zXd(std::size_t d, Eigen::VectorXd seed)
    : d_(d), seed_(seed) {}
private:
  std::size_t d_;
  Eigen::VectorXd seed_;  
};

struct Function_x_y : public Function {

  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    double xr = p(0), xi = p(1), yr = p(2), yi = p(3);
    Eigen::VectorXd result(cod_d());
    // result(0) = xr*xr + xi*xi - yr*yr - yi*yi;
    result(0) = xr*xr + xi*xi - 1;
    return result;
  }

  std::size_t amb_d() const {return 4;};
  std::size_t cod_d() const {return 1;};

  Eigen::VectorXd seed() const {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(amb_d());
    return result;
  }

  Function_x_y() {}
};

struct Function_y_abs : public Function {

  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    double yr = p(2), yi = p(3);
    Eigen::VectorXd result(cod_d());
    result(0) = yr*yr + yi*yi - 1;
    return result;
  }

  std::size_t amb_d() const {return 4;};
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
  // Function_aurelien_x3y_y3z_z3x fun_custom;
  // Function_aurelien_x2Y_y2Z_z2X fun_custom;
  // Function_aurelien_xY_yZ_zX fun_custom;
  // Eigen::Vector4d seed0(-1./3, 0, 0.5, 0);
  // Function_aurelien_xYd_yZd_zXd fun_custom(1, seed0);
  Eigen::Vector4d seed0(-0.90388511704530858879416532, 0, 0.5, 0);
  Function_aurelien_xYd_yZd_zXd fun_custom(7, seed0);

  Function_x_y fun_xy;
  Function_y_abs fun_y;
  Eigen::MatrixXd matrix = random_orthogonal_matrix(fun_custom.amb_d());
  auto fun = make_linear_transformation(fun_custom, matrix);
  auto fun_boundary1 = make_linear_transformation(fun_xy, matrix);
  auto fun_boundary2 = make_linear_transformation(fun_y, matrix);

  // // Seed of the function
  // Eigen::VectorXd seed = fun.seed();
  
  // // Defining the intersection oracle
  // auto oracle = make_oracle(fun);

  // Test for the constant function and the MTA
  std::size_t d = 4;
  std::size_t k = 2;
  Eigen::MatrixXd normal_matrix = Eigen::MatrixXd::Zero(d,k);
  normal_matrix << 0, 0, 1;
  // Function_aurelien fun;
  Eigen::VectorXd seed = fun.seed();
  auto oracle = make_oracle(fun, fun_boundary1, fun_boundary2);

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
  Gudhi::Clock mta_clock("Manifold tracing algorithm");
  manifold_tracing_algorithm(seed_points, cox_tr, oracle,
			     interior_simplex_map,
			     boundary1_simplex_map,
			     boundary2_simplex_map,
			     corner_simplex_map);
  std::cout << mta_clock << std::endl;

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

  int i0 = cell_complex.interior_simplex_cell_map(0).size();
  int i1 = cell_complex.interior_simplex_cell_map(1).size();
  int i2 = cell_complex.interior_simplex_cell_map(2).size();
  int a0 = cell_complex.boundary1_simplex_cell_map(0).size();
  int a1 = cell_complex.boundary1_simplex_cell_map(1).size();
  int b0 = cell_complex.boundary2_simplex_cell_map(0).size();
  int b1 = cell_complex.boundary2_simplex_cell_map(1).size();
  int c0 = cell_complex.corner_simplex_cell_map(0).size();

  int ei = i0 - i1 + i2;
  int eb = a0 - a1 + b0 - b1;
  int ec = c0;
  
  std::cout << "Euler characteristic (interior) = " << ei << "\n";
  std::cout << "Euler characteristic (boundary) = " << eb << "\n";
  std::cout << "Euler characteristic (corners) = " << ec << "\n";
  std::cout << "Euler characteristic (whole surface) = " << 3*ei + 3*eb/2 + ec << "\n";

  // std::cout << "Euler characteristic (total) = " <<
  //   i0 + a0 + b0 + c0 - i1 - a1 - b1 + i2 << "\n";

  // Output the cell complex to a file readable by medit
  // output_meshes_to_medit(3,
  // 			 "custom_manifold",
  // 			 build_mesh_from_cell_complex(cell_complex,
  // 						      Configuration(true, true, true, 1, 5, 3),
  // 						      Configuration(true, true, true, 2, 13, 14)));
  output_meshes_to_off(false,
		       "custom_manifold",
		       build_mesh_from_cell_complex(cell_complex,
						    Configuration(true, true, true, 1, 5, 3),
						    Configuration(true, true, true, 2, 13, 14)));
  return 0;
}
