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
    result(1) = 1;
    return result;
  }

  Function_y_abs() {}
};

struct Function_z_abs : public Function {

  Eigen::VectorXd operator()(const Eigen::VectorXd& p) const {
    double z = p(2);
    Eigen::VectorXd result(cod_d());
    result(0) = std::abs(z) - 1;
    return result;
  }

  std::size_t amb_d() const {return 3;};
  std::size_t cod_d() const {return 1;};

  Eigen::VectorXd seed() const {
    Eigen::VectorXd result = Eigen::VectorXd::Zero(amb_d());
    result(2) = 1;
    return result;
  }

  Function_z_abs() {}
};


int main(int argc, char** argv) {
  // Creating a circle S1 in R2 of specified radius
  std::size_t d = 3;
  std::size_t k = 0;
  Constant_function fun(d, k, Eigen::VectorXd::Zero(0));  
  Eigen::VectorXd seed = Eigen::VectorXd::Zero(d);
  Function_x_abs fun_bound_x;
  Function_y_abs fun_bound_y;
  Function_z_abs fun_bound_z;

  std::vector<Function*> constraint_functions_ = {&fun_bound_x, &fun_bound_y, &fun_bound_z};
  std::cout << "fun_x(seed) =\n" << fun_bound_x(seed) << "\n"; 
  std::cout << "fun_y(seed) =\n" << fun_bound_y(seed) << "\n"; 
  std::cout << "fun_z(seed) =\n" << fun_bound_z(seed) << "\n"; 
  
  // Defining the intersection oracle
  auto oracle = make_oracle(&fun, constraint_functions_);
  
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
  ref_vectors[1][0] = 2;
  ref_vectors[1][1] = 3;
  ref_vectors[1][2] = 4;
  ref_vectors[1][3] = 5;
  ref_vectors[1][4] = 6;
  ref_vectors[1][5] = 7;
  ref_vectors[1][6] = 8;
  ref_vectors[1][7] = 9;
  output_meshes_to_medit(3,
  			 "corner_manifold",
  			 build_mesh_from_cell_complex(cell_complex,
  						      toggle_vectors,
  						      ref_vectors));

}
