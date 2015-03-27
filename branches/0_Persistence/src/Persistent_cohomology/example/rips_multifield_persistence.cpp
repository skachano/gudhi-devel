 /*    This file is part of the Gudhi Library. The Gudhi library 
  *    (Geometric Understanding in Higher Dimensions) is a generic C++ 
  *    library for computational topology.
  *
  *    Author(s):       Clément Maria
  *
  *    Copyright (C) 2014  INRIA Sophia Antipolis-Méditerranée (France)
  *
  *    This program is free software: you can redistribute it and/or modify
  *    it under the terms of the GNU General Public License as published by
  *    the Free Software Foundation, either version 3 of the License, or
  *    (at your option) any later version.
  *
  *    This program is distributed in the hope that it will be useful,
  *    but WITHOUT ANY WARRANTY; without even the implied warranty of
  *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *    GNU General Public License for more details.
  *
  *    You should have received a copy of the GNU General Public License
  *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
  */

#include "gudhi/reader_utils.h"
#include "gudhi/graph_simplicial_complex.h"
#include "gudhi/distance_functions.h"
#include "gudhi/Simplex_tree.h"
#include "gudhi/Persistent_cohomology.h"
#include "gudhi/Persistent_cohomology/Multi_field.h"

#include <boost/program_options.hpp>

using namespace Gudhi;
using namespace Gudhi::persistent_cohomology;

typedef int    Vertex_handle;
typedef double Filtration_value;

void program_options( int argc, char *   argv[] 
                    , std::string      & filepoints
                    , std::string      & filediag
                    , Filtration_value & threshold
                    , int              & dim_max
                    , int              & min_p
                    , int              & max_p
                    , Filtration_value & min_persistence );

int main (int argc, char * argv[]) 
{
  std::string      filepoints;
  std::string      filediag;
  Filtration_value threshold;
  int              dim_max;
  int              min_p;
  int              max_p;
  Filtration_value min_persistence;

  program_options(argc,argv,filepoints,filediag,threshold,dim_max,min_p,max_p,min_persistence);

// Extract the points from the file filepoints
  typedef std::vector<double> Point_t;
  std::vector< Point_t > points;
  read_points( filepoints, points );

// Compute the proximity graph of the points
  Graph_t prox_graph = compute_proximity_graph( points, threshold
                                              , euclidean_distance<Point_t> );

// Construct the Rips complex in a Simplex Tree
  Simplex_tree<> st;        
  st.insert_graph(prox_graph); // insert the proximity graph in the simplex tree
  st.expansion( dim_max ); // expand the graph until dimension dim_max

  std::cout << "st=" << st << std::endl;

// Sort the simplices in the order of the filtration
  st.initialize_filtration();

// Compute the persistence diagram of the complex
  Persistent_cohomology< Simplex_tree<>, Multi_field > pcoh( st );
  pcoh.init_coefficients( min_p, max_p ); //initilizes the coefficient field for homology
  // compute persistent homology, disgarding persistent features of life shorter than min_persistence
  pcoh.compute_persistent_cohomology( min_persistence );

// Output the diagram in filediag
  if(filediag.empty()) {  pcoh.output_diagram(); }
  else {
   std::ofstream out(filediag);
   pcoh.output_diagram(out);
   out.close(); }

  return 0;
}



void program_options( int argc, char *   argv[] 
                    , std::string      & filepoints
                    , std::string      & filediag
                    , Filtration_value & threshold
                    , int              & dim_max
                    , int              & min_p
                    , int              & max_p
                    , Filtration_value & min_persistence )
{
  namespace po = boost::program_options;
  po::options_description     hidden("Hidden options");
  hidden.add_options()
        ("input-file",  po::value<std::string>(&filepoints), 
          "Name of file containing a point set. Format is one point per line:   X1 ... Xd \n");
    
  po::options_description visible("Allowed options");
  visible.add_options()
        ("help,h",          "produce help message")
        ("output-file,o", po::value<std::string>(&filediag)->default_value(std::string()), 
                   "Name of file in which the persistence diagram is written. Default print in std::cout")
        ("max-edge-length,r", po::value<Filtration_value>(&threshold)->default_value(0), 
                   "Maximal length of an edge for the Rips complex construction.")
        ("cpx-dimension,d", po::value<int>(&dim_max)->default_value(1),        
                  "Maximal dimension of the Rips complex we want to compute.")
        ("min-field-charac,p", po::value<int>(&min_p)->default_value(2),                  
                  "Minimal characteristic p of the coefficient field Z/pZ.")
        ("max-field-charac,q", po::value<int>(&max_p)->default_value(1223),                  
                  "Minimial characteristic q of the coefficient field Z/pZ.")
        ("min-persistence,m", po::value<Filtration_value>(&min_persistence), 
                  "Minimal lifetime of homology feature to be recorded. Default is 0");

  po::positional_options_description pos;
  pos.add("input-file", 1);
  
  po::options_description all; all.add(visible).add(hidden);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
                options(all).positional(pos).run(), vm);
  po::notify(vm);

  if (vm.count("help") || !vm.count("input-file"))
  { 
      std::cout << std::endl;
      std::cout << "Compute the persistent homology with various coefficient fields \n";
      std::cout << "of a Rips complex defined on a set of input points. The coefficient \n";
      std::cout << "fields are all the Z/rZ for a prime number r contained in the \n";
      std::cout << "specified range [p,q]\n \n";
      std::cout << "The output diagram contains one bar per line, written with the convention: \n";
      std::cout << "   p1*...*pr   dim b d \n";
      std::cout << "where dim is the dimension of the homological feature,\n";
      std::cout << "b and d are respectively the birth and death of the feature and \n";
      std::cout << "p1*...*pr is the product of prime numbers pi such that the homology \n";
      std::cout << "feature exists in homology with Z/piZ coefficients."<< std::endl << std::endl;

      std::cout << "Usage: " << argv[0] << " [options] input-file" << std::endl << std::endl;
      std::cout << visible << std::endl; 
      std::abort();
  }
}
