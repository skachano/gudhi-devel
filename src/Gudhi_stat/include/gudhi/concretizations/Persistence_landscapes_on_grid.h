/*    This file is part of the Gudhi Library. The Gudhi library
 *    (Geometric Understanding in Higher Dimensions) is a generic C++
 *    library for computational topology.
 *
 *    Author(s):       Pawel Dlotko
 *
 *    Copyright (C) 2015  INRIA Sophia-Saclay (France)
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

#pragma once
#ifndef Persistence_landscape_on_grid_on_grid_H
#define Persistence_landscape_on_grid_on_grid_H

//standard include
#include <iostream>
#include <vector>
#include <limits>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <unistd.h>
#include <cmath>


//gudhi include
#include <gudhi/abstract_classes/Abs_Vectorized_topological_data.h>
#include <gudhi/abstract_classes/Abs_Topological_data_with_averages.h>
#include <gudhi/abstract_classes/Abs_Topological_data_with_distances.h>
#include <gudhi/abstract_classes/Abs_Real_valued_topological_data.h>
#include <gudhi/abstract_classes/Abs_Topological_data_with_scalar_product.h>
#include <gudhi/concretizations/read_persitence_from_file.h>



using namespace std;


namespace Gudhi
{
namespace Gudhi_stat
{

/**
 * Given two points in R^2, the procedure compute the parameters A and B of the line y = Ax + B that crosses those two points.
**/
std::pair<double,double> compute_parameters_of_a_line( std::pair<double,double> p1 , std::pair<double,double> p2 )
{
    double a = (p2.second-p1.second)/( p2.first - p1.first );
    double b = p1.second - a*p1.first;
    return std::make_pair(a,b);
}

struct greater
{
    template<class T>
    bool operator()(T const &a, T const &b) const { return a > b; }
};

//double epsi = std::numeric_limits<double>::epsilon();
double epsi = 0.000005;


/**
 *  A procedure used to compare doubles. Typically gien two doubles A and B, comparing A == B is not good idea. In this case, we use the procedure almostEqual with the epsi defined at
 *  the top of the file. Setting up the epsi give the user a tolerance on what should be consider equal.
**/
inline bool almost_equal( double a , double b )
{
    if ( fabs(a-b) < epsi )
        return true;
    return false;
}




class Persistence_landscape_on_grid  :
								public Abs_Vectorized_topological_data ,
								public Abs_Topological_data_with_distances,
								public Abs_Real_valued_topological_data,
								public Abs_Topological_data_with_averages,
								public Abs_Topological_data_with_scalar_product
{
public:
	/**
	 * Default constructor.
	**/
    Persistence_landscape_on_grid()
    {
		this->set_up_numbers_of_functions_for_vectorization_and_projections_to_reals();
		this->grid_min = this->grid_max = 0;
	}

    /**
	 * Constructor that takes as an input a vector of birth-death pairs.
	**/
    Persistence_landscape_on_grid( const std::vector< std::pair< double , double > >& p , double grid_min_ , double grid_max_ , size_t number_of_points_ );

    /**
	 * Assignement operator.
	**/
    Persistence_landscape_on_grid& operator=( const Persistence_landscape_on_grid& org );

    /**
	 * Copy constructor.
	**/
    Persistence_landscape_on_grid(const Persistence_landscape_on_grid&);

    /**
	 * Constructor that reads persistence intervals from file and creates persistence landscape. The format of the input file is the following: in each line we put birth-death pair. Last line is assumed
	 * to be empty. Even if the points within a line are not ordered, they will be ordered while the input is read. The additional parameters of this procedure are: ranges of grid, resoltion of a grid
	 * and the dimension of intervals that are need to be read from a file (in case of Gudhi format files). 
	**/
    Persistence_landscape_on_grid(const char* filename , double grid_min_, double grid_max_ , size_t number_of_points_ , size_t dimension_ = 0 );
        
    /**
	 * Constructor that reads persistence intervals from file and creates persistence landscape. The format of the input file is the following: in each line we put birth-death pair. Last line is assumed
	 * to be empty. Even if the points within a line are not ordered, they will be ordered while the input is read. The additional parameter is the resoution of a grid. The remaning parameters are 
	 * calculated based on data. 
	**/
    Persistence_landscape_on_grid(const char* filename , size_t number_of_points );


    /**
     * This procedure loads a landscape from file. It erase all the data that was previously stored in this landscape.
    **/
    void load_landscape_from_file( const char* filename );


    /**
     * The procedure stores a landscape to a file. The file can be later used by a procedure load_landscape_from_file.
    **/
    void print_to_file( const char* filename )const;



	/**
	 * This function compute integral of the landscape (defined formally as sum of integrals on R of all landscape functions)
	**/
    double compute_integral_of_landscape()const
    {
		size_t maximal_level = this->number_of_nonzero_levels();
		double result = 0;
		for ( size_t i = 0 ; i != maximal_level ; ++i )
		{
			result += this->compute_integral_of_landscape(i);
		}
		return result;
	}


    /**
	 * This function compute integral of the 'level'-level of a landscape.
	**/
    double compute_integral_of_landscape( size_t level )const
    {
		bool dbg = false;
		double result = 0;		
		double dx = (this->grid_max - this->grid_min)/(double)(this->values_of_landscapes.size()-1);
		
		if ( dbg )
		{			
			cerr << "this->grid_max : " << this->grid_max << endl;
			cerr << "this->grid_min : " << this->grid_min << endl;
			cerr << "this->values_of_landscapes.size() : " << this->values_of_landscapes.size() << endl;
			getchar();
		}
		
		
		double previous_x = this->grid_min-dx;
		double previous_y = 0;
		for ( size_t i = 0 ; i != this->values_of_landscapes.size() ; ++i )
		{
			double current_x = previous_x + dx;
			double current_y = 0;
			if ( this->values_of_landscapes[i].size() > level )current_y = this->values_of_landscapes[i][level];
		
			if ( dbg )
			{
				cerr << "this->values_of_landscapes[i].size() : " << this->values_of_landscapes[i].size() << " , level : " << level << endl;	
				if ( this->values_of_landscapes[i].size() > level )cerr << "this->values_of_landscapes[i][level] : " << this->values_of_landscapes[i][level] << endl;
				cerr << "previous_y : " << previous_y << endl;
				cerr << "current_y : " << current_y << endl;
				cerr << "dx : " << dx << endl;
				cerr << "0.5*dx*( previous_y + current_y ); " << 0.5*dx*( previous_y + current_y ) << endl;
			}
			
			result += 0.5*dx*( previous_y + current_y );			
			previous_x = current_x;
			previous_y = current_y;
		}
		return result;
	}

	/**
	 * This function compute integral of the landscape p-th power of a landscape (defined formally as sum of integrals on R of p-th powers of all landscape functions)
	**/
	double compute_integral_of_landscape( double p )const
	{
		size_t maximal_level = this->number_of_nonzero_levels();
		double result = 0;
		for ( size_t i = 0 ; i != maximal_level ; ++i )
		{
			result += this->compute_integral_of_landscape(p,i);
		}
		return result;
	}

    /**
	 * This function compute integral of the landscape p-th power of a level of a landscape (defined formally as sum of integrals on R of p-th powers of all landscape functions)
	**/	
	double compute_integral_of_landscape( double p , size_t level )const
	{
		bool dbg = false;
		
		double result = 0;				
		double dx = (this->grid_max - this->grid_min)/(double)(this->values_of_landscapes.size()-1);		
		double previous_x = this->grid_min;
		double previous_y = 0;
		if ( this->values_of_landscapes[0].size() > level )previous_y = this->values_of_landscapes[0][level];
		
		if ( dbg )
		{
			cerr << "dx : " << dx << endl;
			cerr << "previous_x : " << previous_x << endl;
			cerr << "previous_y : " << previous_y << endl;
			cerr << "power : " << p << endl;
			getchar();
		}
		
		for ( size_t i = 0 ; i != this->values_of_landscapes.size() ; ++i )
		{
			double current_x = previous_x + dx;
			double current_y = 0;
			if ( this->values_of_landscapes[i].size() > level )current_y = this->values_of_landscapes[i][level];
			
			if ( dbg )cerr << "current_y : " << current_y << endl;
			
			if ( current_y == previous_y )continue; 
			
			std::pair<double,double> coef = compute_parameters_of_a_line( std::make_pair( previous_x , previous_y ) , std::make_pair( current_x , current_y ) );
			double a = coef.first;
			double b = coef.second;
			
			if ( dbg )
			{
				cerr << "A line passing through points : (" << previous_x << "," << previous_y << ") and (" << current_x << "," << current_y << ") is : " << a << "x+" << b << endl;
			}
			
			//In this interval, the landscape has a form f(x) = ax+b. We want to compute integral of (ax+b)^p = 1/a * (ax+b)^{p+1}/(p+1)
			double value_to_add = 0;
			if ( a != 0 )
			{			
				value_to_add = 1/(a*(p+1)) * ( pow((a*current_x+b),p+1) - pow((a*previous_x+b),p+1));				
			}
			else
			{
				value_to_add = ( current_x - previous_x )*( pow(b,p) );
			}				
			result += value_to_add;
			if ( dbg )
			{
				cerr << "Increasing result by : " << value_to_add << endl;
				cerr << "restult : " << result << endl;
				getchar();
			}		
			previous_x = current_x;
			previous_y = current_y;
		}		
		if ( dbg )cerr << "The total result is : " << result << endl;
		return result;
	}
	
	/**
     * Writing landscape into a stream. A i-th level landscape starts with a string "lambda_i". Then the discontinuity points of the landscapes follows.
     * Shall those points be joined with lines, we will obtain the i-th landscape function.
    **/
    friend std::ostream& operator<<(std::ostream& out, const Persistence_landscape_on_grid& land )
	{
		double dx = (land.grid_max - land.grid_min)/(double)(land.values_of_landscapes.size()-1);
		double x = land.grid_min;
		for ( size_t i = 0 ; i != land.values_of_landscapes.size() ; ++i )
		{
			out << x << " : ";
			for ( size_t j = 0 ; j != land.values_of_landscapes[i].size() ; ++j )
			{
				out << land.values_of_landscapes[i][j] << " ";
			}
			out << endl;
			x += dx;
		}
		return out;
	}


    /**
     * A function that computes the value of a landscape at a given point. The parameters of the function are: unsigned level and double x.
     * The procedure will compute the value of the level-landscape at the point x.
    **/
    double compute_value_at_a_given_point( unsigned level , double x )const
    {
		bool dbg = false;
		if ( (x < this->grid_min) || (x > this->grid_max) )return 0;
		
		//find a position of a vector closest to x:
		double dx = (this->grid_max - this->grid_min)/(double)(this->values_of_landscapes.size()-1);
		size_t position = size_t((x-this->grid_min)/dx);
		
		if ( dbg )
		{
			std::cerr << "This is a procedure compute_value_at_a_given_point \n";
			std::cerr << "level : " << level << endl;
			std::cerr << "x : " << x << endl;
			std::cerr << "psoition : " << position << endl;
		}
		//check if we are not exacly in the grid point:
		if ( almost_equal( position*dx+ this->grid_min ,  x) )
		{
			if ( this->values_of_landscapes[position].size() < level )
			{
				return this->values_of_landscapes[position][level];
			}
			else
			{
				return 0;
			}
		}
		//in the other case, approximate with a line:
		std::pair<double,double> line;
		if ( (this->values_of_landscapes[position].size() > level) && (this->values_of_landscapes[position+1].size() > level) )
		{
			line = compute_parameters_of_a_line( std::make_pair( position*dx+ this->grid_min , this->values_of_landscapes[position][level] ) , std::make_pair( (position+1)*dx+ this->grid_min , this->values_of_landscapes[position+1][level] ) );
		}
		else
		{
			if ( (this->values_of_landscapes[position].size() > level) || (this->values_of_landscapes[position+1].size() > level) )
			{
				if ( (this->values_of_landscapes[position].size() > level) )
				{
						line = compute_parameters_of_a_line( std::make_pair( position*dx+ this->grid_min , this->values_of_landscapes[position][level] ) , std::make_pair( (position+1)*dx+ this->grid_min , 0 ) );
				}
				else
				{
					//(this->values_of_landscapes[position+1].size() > level)
					line = compute_parameters_of_a_line( std::make_pair( position*dx+ this->grid_min , 0 ) , std::make_pair( (position+1)*dx+ this->grid_min , this->values_of_landscapes[position+1][level] ) );	
				}
			}
			else
			{
				return 0;
			}
		}
		//compute the value of the linear function parametrized by line on a point x:
		return line.first*x+line.second;
	}
    

	/**
	 * A function that compute sum of two landscapes.
	**/
    friend Persistence_landscape_on_grid add_two_landscapes ( const Persistence_landscape_on_grid& land1 ,  const Persistence_landscape_on_grid& land2 )
    {
        return operation_on_pair_of_landscapes_on_grid< std::plus<double> >(land1,land2);
    }

    /**
	 * A function that compute difference of two landscapes.
	**/
    friend Persistence_landscape_on_grid subtract_two_landscapes ( const Persistence_landscape_on_grid& land1 ,  const Persistence_landscape_on_grid& land2 )
    {
        return operation_on_pair_of_landscapes_on_grid< std::minus<double> >(land1,land2);
    }

	/**
	 * An operator +, that compute sum of two landscapes.
	**/
    friend Persistence_landscape_on_grid operator+( const Persistence_landscape_on_grid& first , const Persistence_landscape_on_grid& second )
    {
        return add_two_landscapes( first,second );
    }

	/**
	 * An operator -, that compute difference of two landscapes.
	**/
    friend Persistence_landscape_on_grid operator-( const Persistence_landscape_on_grid& first , const Persistence_landscape_on_grid& second )
    {
        return subtract_two_landscapes( first,second );
    }

	/**
	 * An operator * that allows multipilication of a landscape by a real number.
	**/
    friend Persistence_landscape_on_grid operator*( const Persistence_landscape_on_grid& first , double con )
    {
        return first.multiply_lanscape_by_real_number_not_overwrite(con);
    }

	/**
	 * An operator * that allows multipilication of a landscape by a real number (order of parameters swapped).
	**/
    friend Persistence_landscape_on_grid operator*( double con , const Persistence_landscape_on_grid& first  )
    {
        return first.multiply_lanscape_by_real_number_not_overwrite(con);
    }

	friend bool check_if_defined_on_the_same_domain( const Persistence_landscape_on_grid& land1, const Persistence_landscape_on_grid& land2 )
	{
		if ( land1.values_of_landscapes.size() != land2.values_of_landscapes.size() )return false;
		if ( land1.grid_min != land2.grid_min )return false;
		if ( land1.grid_max != land2.grid_max )return false;
		return true;
	}

	/**
	 * Operator +=. The second parameter is persistnece landwscape.
	**/
    Persistence_landscape_on_grid operator += ( const Persistence_landscape_on_grid& rhs )
    {
        *this = *this + rhs;
        return *this;
    }

	/**
	 * Operator -=. The second parameter is persistnece landwscape.
	**/
    Persistence_landscape_on_grid operator -= ( const Persistence_landscape_on_grid& rhs )
    {
        *this = *this - rhs;
        return *this;
    }


	/**
	 * Operator *=. The second parameter is a real number by which the y values of all landscape functions are multiplied. The x-values remain unchanged. 
	**/
    Persistence_landscape_on_grid operator *= ( double x )
    {
        *this = *this*x;
        return *this;
    }

	/**
	 * Operator /=. The second parameter is a real number.
	**/
    Persistence_landscape_on_grid operator /= ( double x )
    {
        if ( x == 0 )throw( "In operator /=, division by 0. Program terminated." );
        *this = *this * (1/x);
        return *this;
    }

	/**
	 * An operator to compare two persistence landscapes.
	**/
    bool operator == ( const Persistence_landscape_on_grid& rhs  )const
    {
		bool dbg = true;
		if ( ! this->values_of_landscapes.size() == rhs.values_of_landscapes.size() )
		{
			if (dbg) cerr << "values_of_landscapes of incompatable sizes\n";
			return false;
		}
		if ( !almost_equal( this->grid_min , rhs.grid_min ) )
		{
			if (dbg) cerr << "grid_min not equal\n";
			return false; 
		}
		if ( !almost_equal(this->grid_max,rhs.grid_max ) )
		{	
			if (dbg) cerr << "grid_max not equal\n";
			return false;		
		}
		for ( size_t i = 0 ; i != this->values_of_landscapes.size() ; ++i )
		{					
			for ( size_t aa = 0 ; aa != this->values_of_landscapes[i].size() ; ++aa )
			{
				if ( !almost_equal( this->values_of_landscapes[i][aa]  , rhs.values_of_landscapes[i][aa]  ) )
				{					
					if (dbg) 
					{
						cerr << "Problem in the position : " << i << " of values_of_landscapes. \n";
						cerr << this->values_of_landscapes[i][aa]   << " " <<  rhs.values_of_landscapes[i][aa] << endl;
					}
					return false;
				}
			}				
		}	
		return true;
	}


    /**
	 * An operator to compare two persistence landscapes.
	**/
    bool operator != ( const Persistence_landscape_on_grid& rhs  )const
    {
		return !((*this) == rhs);
	}


	/**
	 * Computations of maximum (y) value of landscape.
	**/
    double compute_maximum()const
    {
       //since the function can only be entirely positive or negative, the maximal value will be an extremal value in the arrays:
		double max_value = -std::numeric_limits<double>::max();
		for ( size_t i = 0 ; i != this->values_of_landscapes.size() ; ++i )
		{
			if ( this->values_of_landscapes[i].size() )
			{
				if ( this->values_of_landscapes[i][0] > max_value )max_value = this->values_of_landscapes[i][0];
				if ( this->values_of_landscapes[i][ this->values_of_landscapes[i].size()-1 ] > max_value )max_value = this->values_of_landscapes[i][ this->values_of_landscapes[i].size()-1 ];
			}
		}
		return max_value;        
    }
    
    /**
	 * Computations of minimum and maximum value of landscape.
	**/
    std::pair<double,double> compute_minimum_maximum()const
    {        
       //since the function can only be entirely positive or negative, the maximal value will be an extremal value in the arrays:
		double max_value = -std::numeric_limits<double>::max();
		double min_value = 0;
		for ( size_t i = 0 ; i != this->values_of_landscapes.size() ; ++i )
		{
			if ( this->values_of_landscapes[i].size() )
			{
				if ( this->values_of_landscapes[i][0] > max_value )max_value = this->values_of_landscapes[i][0];
				if ( this->values_of_landscapes[i][ this->values_of_landscapes[i].size()-1 ] > max_value )max_value = this->values_of_landscapes[i][ this->values_of_landscapes[i].size()-1 ];
				
				if ( this->values_of_landscapes[i][0] < min_value )min_value = this->values_of_landscapes[i][0];
				if ( this->values_of_landscapes[i][ this->values_of_landscapes[i].size()-1 ] < min_value )min_value = this->values_of_landscapes[i][ this->values_of_landscapes[i].size()-1 ];
			}
		}
		return std::make_pair(min_value , max_value);
    }
    
    /**
     * This function computes maximal lambda for which lambda-level landscape is nonzero.
    **/
    size_t number_of_nonzero_levels()const
    {
		size_t result = 0;
		for ( size_t i = 0 ; i != this->values_of_landscapes.size() ; ++i )
		{
			if ( this->values_of_landscapes[i].size() > result )result = this->values_of_landscapes[i].size();
		}
		return result;
	}

	/**
	 * Computations of a L^i norm of landscape, where i is the input parameter.
	**/
    double compute_norm_of_landscape( double i )
    {
		std::vector< std::pair< double , double > > p;
        Persistence_landscape_on_grid l(p,this->grid_min,this->grid_max,this->values_of_landscapes.size()-1);       
        
        if ( i != -1 )
        {
            return compute_discance_of_landscapes_on_grid(*this,l,i);
        }
        else
        {
            return compute_max_norm_discance_of_landscapes(*this,l);
        }
    }

	/**
 	 * An operator to compute the value of a landscape in the level 'level' at the argument 'x'.
	**/
    double operator()(unsigned level,double x)const{return this->compute_value_at_a_given_point(level,x);}

	/**
	 * Computations of L^{\infty} distance between two landscapes.
	**/
    friend double compute_max_norm_discance_of_landscapes( const Persistence_landscape_on_grid& first, const Persistence_landscape_on_grid& second );
    //friend double compute_max_norm_discance_of_landscapes( const Persistence_landscape_on_grid& first, const Persistence_landscape_on_grid& second , unsigned& nrOfLand , double&x , double& y1, double& y2 );




	/**
	 * Function to compute absolute value of a PL function. The representation of persistence landscapes allow to store general PL-function. When computing distance betwen two landscapes, we compute difference between
	 * them. In this case, a general PL-function with negative value can appear as a result. Then in order to compute distance, we need to take its absolute value. This is the purpose of this procedure.
	**/
    void abs()
    {
		//Be careful here. We assume that the functions are either entirely positive or negative. They do not change signs. That is why I can implemnt abs in the way presented below:
		for ( size_t i = 0 ; i != this->values_of_landscapes.size() ; ++i )
		{
			for ( size_t j = 0 ; j != this->values_of_landscapes[i].size() ; ++j )
			{
				this->values_of_landscapes[i][j] = std::abs( this->values_of_landscapes[i][j] );
			}
		}
	}

	/**
	 * Computes the number of landscape functions.
	**/
    size_t size()const{return this->number_of_nonzero_levels(); }

	/**
	 *  Computate maximal value of lambda-level landscape.
	**/
    double find_max( unsigned lambda )const
    {		
		double max_value = -std::numeric_limits<double>::max();
		for ( size_t i = 0 ; i != this->values_of_landscapes.size() ; ++i )
		{
			if ( this->values_of_landscapes[i].size() > lambda )
			{
				if ( this->values_of_landscapes[i][lambda] > max_value )max_value = this->values_of_landscapes[i][lambda];				
			}
		}
		return max_value;
	}

	/**
	 * Function to compute inner (scalar) product of two landscapes.
	**/
    friend double compute_inner_product( const Persistence_landscape_on_grid& l1 , const Persistence_landscape_on_grid& l2 )
    {
		if ( !check_if_defined_on_the_same_domain(l1,l2) )throw "Landscapes are not defined on the same grid, the program will now terminate";
		size_t maximal_level = l1.number_of_nonzero_levels();
		double result = 0;
		for ( size_t i = 0 ; i != maximal_level ; ++i )
		{
			result += compute_inner_product(l1,l2,i);
		}
		return result;
	}

	
	
//SOMETHIND IS WRONG OVER HERE!!! THE SCALAR PRODUC DEPENDS A LOT ON THE RESOLUTION OF A GRID> THIS SHOULD NOT BE LIKE THIS!!!
	/**
	 * Function to compute inner (scalar) product of given levels of two landscapes.
	**/
    friend double compute_inner_product( const Persistence_landscape_on_grid& l1 , const Persistence_landscape_on_grid& l2 , size_t level )
    {
		if ( !check_if_defined_on_the_same_domain(l1,l2) )throw "Landscapes are not defined on the same grid, the program will now terminate";
		double result = 0;
		
		double dx = (l1.grid_max - l1.grid_min)/(double)(l1.values_of_landscapes.size()-1);
		
		double previous_x = l1.grid_min-dx;
		double previous_y_l1 = 0;
		double previous_y_l2 = 0;
		for ( size_t i = 0 ; i != l1.values_of_landscapes.size() ; ++i )
		{
			double current_x = previous_x + dx;
			double current_y_l1 = 0;
			if ( l1.values_of_landscapes[i].size() > level )current_y_l1 = l1.values_of_landscapes[i][level];
			
			double current_y_l2 = 0;
			if ( l2.values_of_landscapes[i].size() > level )current_y_l2 = l2.values_of_landscapes[i][level];
			
			std::pair<double,double> l1_coords = compute_parameters_of_a_line( std::make_pair( previous_x , previous_y_l1 ) , std::make_pair( current_x , current_y_l1 ) );			
			std::pair<double,double> l2_coords = compute_parameters_of_a_line( std::make_pair( previous_x , previous_y_l2 ) , std::make_pair( current_x , current_y_l2 ) );
			
			//let us assume that the first line is of a form y = ax+b, and the second one is of a form y = cx + d. Then here are a,b,c,d:
			double a = l1_coords.first;
			double b = l1_coords.second;
			
			double c = l2_coords.first;
			double d = l2_coords.second;
			
			//now, to compute the inner product in this interval we need to compute the integral of (ax+b)(cx+d) = acx^2 + (ad+bc)x + bd in the interval from previous_x to current_x:
			//The integal is ac/3*x^3 + (ac+bd)/2*x^2 + bd*x
		
			result += a*c/3*current_x*current_x*current_x + (a*c+b*d)/2*current_x*current_x + b*d*current_x-
			          a*c/3*previous_x*previous_x*previous_x + (a*c+b*d)/2*previous_x*previous_x + b*d*previous_x;
			
			previous_x = current_x;
			previous_y_l1 = current_y_l1;
			previous_y_l2 = current_y_l2;
			
		}		
		return result;
	}




	//concretization of abstract functions:

	/**
	 * The number of projections to R is defined to the number of nonzero landscape functions. I-th projection is an integral of i-th landscape function over whole R.
	**/
    double project_to_R( int number_of_function )
    {
		return this->compute_integral_of_landscape( (size_t)number_of_function );
	}

    std::vector<double> vectorize( int number_of_function )
    {
		//TODO, think of something smarter over here
		if ( ( number_of_function < 0 ) || ( (size_t)number_of_function >= this->values_of_landscapes.size() ) )
		{
			throw "Wrong number of function\n";
		}
		std::vector<double> v = this->values_of_landscapes[ number_of_function ];
		return v;
	}
    void compute_average( std::vector< Abs_Topological_data_with_averages* > to_average )
    {
		
		bool dbg = false;
		//After execution of this procedure, the average is supposed to be in the current object. To make sure that this is the case, we need to do some cleaning first.
		this->values_of_landscapes.clear();
		this->grid_min = this->grid_max = 0;
		
		//if there is nothing to averate, then the average is a zero landscape. 
		if ( to_average.size() == 0 )return;
		
		//now we need to check if the grids in all objects of to_average are the same:
		for ( size_t i = 0 ; i != to_average.size() ; ++i )
		{	
			if ( !check_if_defined_on_the_same_domain(*((Persistence_landscape_on_grid*)(to_average[0])),*((Persistence_landscape_on_grid*)(to_average[i]))) )throw "Two grids are not compatible"; 				
		}
		
		this->values_of_landscapes = std::vector< std::vector<double> >( ((Persistence_landscape_on_grid*)(to_average[0]))->values_of_landscapes.size() );
		this->grid_min = ((Persistence_landscape_on_grid*)(to_average[0]))->grid_min;
		this->grid_max = ((Persistence_landscape_on_grid*)(to_average[0]))->grid_max;
		
		if ( dbg )
		{
			cerr << "Computations of average. The data from the current landscape have been cleared. We are ready to do the computations. \n";
		}
						
		//for every point in the grid:
		for ( size_t grid_point = 0 ; grid_point != ((Persistence_landscape_on_grid*)(to_average[0]))->values_of_landscapes.size() ; ++grid_point )
		{
			
			//set up a vector of the correct size:
			size_t maximal_size_of_vector = 0;
			for ( size_t land_no = 0 ; land_no != to_average.size() ; ++land_no )
			{
				if ( ((Persistence_landscape_on_grid*)(to_average[land_no]))->values_of_landscapes[grid_point].size() > maximal_size_of_vector )
				maximal_size_of_vector = ((Persistence_landscape_on_grid*)(to_average[land_no]))->values_of_landscapes[grid_point].size();
			}
			this->values_of_landscapes[grid_point] = std::vector<double>( maximal_size_of_vector );
			
			if ( dbg )
			{
				cerr << "We are considering the point : " << grid_point << " of the grid. In this point, there are at most : " << maximal_size_of_vector << " nonzero landscape functions \n";
			}
		
			//and compute an arythmetic average:
			for ( size_t land_no = 0 ; land_no != to_average.size() ; ++land_no )
			{
				//summing:				
				for ( size_t i = 0 ; i != ((Persistence_landscape_on_grid*)(to_average[land_no]))->values_of_landscapes[grid_point].size() ; ++i )
				{
					//compute the average in a smarter way.
					this->values_of_landscapes[grid_point][i] += ((Persistence_landscape_on_grid*)(to_average[land_no]))->values_of_landscapes[grid_point][i];
				}								
			}
			//normalizing:
			for ( size_t i = 0 ; i != this->values_of_landscapes[grid_point].size() ; ++i )
			{
				this->values_of_landscapes[grid_point][i] /= (double)to_average.size();
			}			 
		}
	}//compute_average
	
	
	/**
	 * Computations of L^{p} distance between two landscapes on a grid. p is the parameter of the procedure.
	**/
	friend double compute_discance_of_landscapes_on_grid( const Persistence_landscape_on_grid& first, const Persistence_landscape_on_grid& second , int p )
	{		
		bool dbg = false;
		//This is what we want to compute: (\int_{- \infty}^{+\infty}| first-second |^p)^(1/p). We will do it one step at a time:

		if ( dbg )
		{
			cerr << "first : " << first << endl;
			cerr << "second : " << second << endl;
			getchar();
		}

		//first-second :
		Persistence_landscape_on_grid lan = first-second;
		
		if ( dbg )
		{
			cerr << "Difference : " << lan << endl;
		}

		//| first-second |:
		lan.abs();
		
		if ( dbg )
		{
			cerr << "Abs : " << lan << endl;
		}
		
		if ( p != -1 )
		{
			//\int_{- \infty}^{+\infty}| first-second |^p
			double result;
			if ( p != 1 )
			{
				if (dbg){cerr << "p : " << p << endl; getchar();}
				result = lan.compute_integral_of_landscape( (double)p );				
				if (dbg){cerr << "integral : " << result << endl;getchar();}
			}
			else
			{
				result = lan.compute_integral_of_landscape();			
				if (dbg){cerr << "integral, wihtout power : " << result << endl;getchar();}
			}
			//(\int_{- \infty}^{+\infty}| first-second |^p)^(1/p)
			return pow( result , 1/(double)p );
		}
		else
		{
			//p == -1
			return lan.compute_maximum();
		}
	}


    double distance( const Abs_Topological_data_with_distances* second , double power = 1 )
    {
		if ( power != -1 )
		{
			return compute_discance_of_landscapes_on_grid( *this , *((Persistence_landscape_on_grid*)second) , power );
		}
		else
		{
			return compute_max_norm_discance_of_landscapes( *this , *((Persistence_landscape_on_grid*)second) );
		}
	}


	double compute_scalar_product( const Abs_Topological_data_with_scalar_product* second )
	{
		return compute_inner_product( (*this) , *((Persistence_landscape_on_grid*)second) );
	}


	std::vector< std::vector< double > > output_for_visualization()
	{
		return this->values_of_landscapes;
	}
	
	
	//a function used to create a gnuplot script for visualization of landscapes
	void plot( const char* filename , size_t from_ = std::numeric_limits<size_t>::max(), size_t to_ = std::numeric_limits<size_t>::max() );


private:
	double grid_min;
	double grid_max;
	std::vector< std::vector< double > > values_of_landscapes;
	
	void set_up_numbers_of_functions_for_vectorization_and_projections_to_reals()
	{
		//warning, this function can be only called after filling in the values_of_landscapes vector.
		this->number_of_functions_for_vectorization = this->values_of_landscapes.size();
		this->number_of_functions_for_projections_to_reals = this->values_of_landscapes.size();
	}
	
	void set_up_values_of_landscapes( const std::vector< std::pair< double , double > >& p , double grid_min_ , double grid_max_ , size_t number_of_points_ );
	template < typename oper > friend Persistence_landscape_on_grid operation_on_pair_of_landscapes_on_grid( const Persistence_landscape_on_grid& land1 ,  const Persistence_landscape_on_grid& land2 );
	Persistence_landscape_on_grid multiply_lanscape_by_real_number_not_overwrite( double x )const;
};


void Persistence_landscape_on_grid::set_up_values_of_landscapes( const std::vector< std::pair< double , double > >& p , double grid_min_ , double grid_max_ , size_t number_of_points_ )
{
	bool dbg = false;
	if ( dbg )
	{
		std::cerr << "Here is the procedure : set_up_values_of_landscapes. The parameters are : grid_min_ : " << grid_min_ << ", grid_max_ : " << grid_max_ << ", number_of_points_ : " << number_of_points_ << endl;
		//getchar();
		std::cerr << "Here are the intervals at our disposal : \n";
		for ( size_t i = 0 ; i != p.size() ; ++i )
		{
			std::cerr << p[i].first << " , " << p[i].second << endl;
		} 
	}
	
	this->values_of_landscapes = std::vector< std::vector< double > >( number_of_points_+1 );
	this->grid_min = grid_min_;
	this->grid_max = grid_max_;
	
	if ( grid_max_ <= grid_min_ )
	{
		throw "Wrong parameters of grid_min and grid_max given to the procedure. THe grid have negative, or zero size. The program will now terminate.\n";
	}
	
	double dx = ( grid_max_ - grid_min_ )/(double)(number_of_points_);		
	//for every interval in the diagram:
	for ( size_t int_no = 0 ; int_no != p.size() ; ++int_no )
	{		
		size_t grid_interval_begin = (p[int_no].first-grid_min_)/dx;
		size_t grid_interval_end = (p[int_no].second-grid_min_)/dx;
		size_t grid_interval_midpoint = (size_t)(0.5*(grid_interval_begin+grid_interval_end));
		
		if ( dbg )
		{
			cerr << "Considering an interval : " << p[int_no].first << "," << p[int_no].second << endl;
		
			std::cerr << "grid_interval_begin : " << grid_interval_begin << std::endl;
			std::cerr << "grid_interval_end : " << grid_interval_end << std::endl;
			std::cerr << "grid_interval_midpoint : " << grid_interval_midpoint << std::endl;	
		}
		
		double landscape_value = dx;
		for ( size_t i = grid_interval_begin+1 ; i < grid_interval_midpoint ; ++i )
		{
			if ( dbg )
			{
				std::cerr << "Adding landscape value (going up) for a point : " << i << " equal : " << landscape_value << std::endl;
			}
			this->values_of_landscapes[i].push_back( landscape_value );
			landscape_value += dx;
		}
		for ( size_t i = grid_interval_midpoint ; i <= grid_interval_end ; ++i )
		{						
			if ( landscape_value > 0 ) 
			{ 
				this->values_of_landscapes[i].push_back( landscape_value );	
				if ( dbg )
				{				
					std::cerr << "AAdding landscape value (going down) for a point : " << i << " equal : " << landscape_value << std::endl;
				}
			}
			landscape_value -= dx;	
		}  
	}
	
	//and now we need to sort the valuesL
	for ( size_t pt = 0 ; pt != this->values_of_landscapes.size() ; ++pt )
	{
		std::sort( this->values_of_landscapes[pt].begin() , this->values_of_landscapes[pt].end() , greater() );	
	}
}//set_up_values_of_landscapes

Persistence_landscape_on_grid::Persistence_landscape_on_grid( const std::vector< std::pair< double , double > >& p , double grid_min_ , double grid_max_ , size_t number_of_points_ )
{
	this->set_up_values_of_landscapes( p , grid_min_ , grid_max_ , number_of_points_ );
}//Persistence_landscape_on_grid

Persistence_landscape_on_grid& Persistence_landscape_on_grid::operator=( const Persistence_landscape_on_grid& org )
{
	this->grid_min = org.grid_min;
	this->grid_max = org.grid_max;
	this->values_of_landscapes = org.values_of_landscapes;
	return (*this);
}//operator=

Persistence_landscape_on_grid::Persistence_landscape_on_grid(const Persistence_landscape_on_grid& org)
{
	this->grid_min = org.grid_min;
	this->grid_max = org.grid_max;
	this->values_of_landscapes = org.values_of_landscapes;
}//copy constructor

Persistence_landscape_on_grid::Persistence_landscape_on_grid(const char* filename , double grid_min_, double grid_max_ , size_t number_of_points_ , size_t dimension )
{
	//standard file with barcode
    std::vector< std::pair< double , double > > p = read_standard_file( filename );    
    //gudhi file with barcode
    //std::vector< std::pair< double , double > > p = read_gudhi_file( filename , dimension );        
	
	this->set_up_values_of_landscapes( p , grid_min_ , grid_max_ , number_of_points_ );
}

Persistence_landscape_on_grid::Persistence_landscape_on_grid(const char* filename , size_t number_of_points_ )
{
	//standard file with barcode
    std::vector< std::pair< double , double > > p = read_standard_file( filename );    
    //gudhi file with barcode
    //std::vector< std::pair< double , double > > p = read_gudhi_file( filename , dimension );     
    
    double grid_min_ = std::numeric_limits<double>::max();
    double grid_max_ = -std::numeric_limits<double>::max();
    for ( size_t i = 0 ; i != p.size() ; ++i )
    {
		if ( p[i].first < grid_min_ )grid_min_ = p[i].first;
		if ( p[i].second > grid_max_ )grid_max_ = p[i].second;
	}	
	this->set_up_values_of_landscapes( p , grid_min_ , grid_max_ , number_of_points_ );
}

void Persistence_landscape_on_grid::load_landscape_from_file( const char* filename )
{
	//check if the file exist.
	if ( !( access( filename, F_OK ) != -1 ) )
	{
		cerr << "The file : " << filename << " do not exist. The program will now terminate \n";
		throw "The file from which you are trying to read the persistence landscape do not exist. The program will now terminate \n";
	}	
	std::ifstream in;
	in.open( filename );
	
	size_t number_of_points_in_the_grid = 0;
	in >> this->grid_min >> this->grid_max >> number_of_points_in_the_grid;
	
	std::vector< std::vector< double > > v(number_of_points_in_the_grid);
	std::string line;
	std::getline(in, line);
	double number;
	for ( size_t i = 0 ; i != number_of_points_in_the_grid ; ++i )
	{
		//read a line of a file and convert it to a vector.
		std::vector< double > vv;		
		std::getline(in, line);		
		//cerr << "Reading line : " << line << endl;getchar();		
		std::istringstream stream(line);
		while (stream >> number)
		{
			vv.push_back(number);
		}		
		v[i] = vv;
	}
	this->values_of_landscapes = v;
	in.close();
}

void Persistence_landscape_on_grid::print_to_file( const char* filename )const
{
	std::ofstream out;
	out.open( filename );
	
	//first we store the parameters of the grid:
	out << grid_min << std::endl << grid_max << std::endl << this->values_of_landscapes.size() << std::endl;
	
	//and now in the following lines, the values of this->values_of_landscapes for the following arguments:
	for ( size_t i = 0 ; i != this->values_of_landscapes.size() ; ++i )
	{
		for ( size_t j = 0 ; j != this->values_of_landscapes[i].size() ; ++j )
		{
			out << this->values_of_landscapes[i][j] << " ";
		}
		out << std::endl;
	}
	
	out.close();
}

void Persistence_landscape_on_grid::plot( const char* filename , size_t from_ , size_t to_ )
{		
    //this program create a gnuplot script file that allows to plot persistence diagram.
    ofstream out;

    std::ostringstream nameSS;
    nameSS << filename << "_GnuplotScript";
    std::string nameStr = nameSS.str();
    out.open( (char*)nameStr.c_str() );

    std::pair<double,double> min_max = compute_minimum_maximum();
    out << "set xrange [" << this->grid_min << " : " << this->grid_max << "]" << endl;
    out << "set yrange [" << min_max.first << " : " << min_max.second << "]" << endl;
    
    size_t number_of_nonzero_levels = this->number_of_nonzero_levels();
    double dx = ( this->grid_max - this->grid_min )/((double)this->values_of_landscapes.size()-1);
    
    
    size_t from = 0;
    if ( from_ != std::numeric_limits<size_t>::max() )
    {
		if ( from_ < number_of_nonzero_levels )
		{
			from = from_;
		}
		else
		{
			return;
		}
	}	
	size_t to = number_of_nonzero_levels;
	if ( to_ != std::numeric_limits<size_t>::max() )
    {
		if ( to_ < number_of_nonzero_levels )
		{
			to = to_;
		}
	}
    
    
    out << "plot ";    
    for ( size_t lambda= from ; lambda != to ; ++lambda )
    {
        //out << "     '-' using 1:2 title 'l" << lambda << "' with lp";
        out << "     '-' using 1:2 notitle with lp";
        if ( lambda+1 != to )
        {
            out << ", \\";
        }
        out << endl;
    }
    
    for ( size_t lambda = from ; lambda != to ; ++lambda )
    {
		double point = this->grid_min;
        for ( size_t i = 0 ; i != this->values_of_landscapes.size() ; ++i )
        {
			double value = 0;
			if ( this->values_of_landscapes[i].size() > lambda )			
			{
				value = this->values_of_landscapes[i][lambda];
			}
            out << point << " " << value << endl;
            point += dx;
        }
        out << "EOF" << endl;
    }
    cout << "Gnuplot script to visualize persistence diagram written to the file: " << nameStr << ". Type load '" << nameStr << "' in gnuplot to visualize." << endl;
}

template < typename T > 
Persistence_landscape_on_grid operation_on_pair_of_landscapes_on_grid ( const Persistence_landscape_on_grid& land1 ,  const Persistence_landscape_on_grid& land2 )
{
	//first we need to check if the domains are the same:
	if ( !check_if_defined_on_the_same_domain(land1,land2) )throw "Two grids are not compatible"; 	
	
	T oper;
	Persistence_landscape_on_grid result;
	result.values_of_landscapes = std::vector< std::vector< double > >( land1.values_of_landscapes.size()  );
	result.grid_min = land1.grid_min;
	result.grid_max = land1.grid_max;
	
	//now we perorm the operations:
	for ( size_t grid_point = 0 ; grid_point != land1.values_of_landscapes.size() ; ++grid_point )
	{
		result.values_of_landscapes[grid_point] = std::vector< double >( std::max( land1.values_of_landscapes[grid_point].size() , land2.values_of_landscapes[grid_point].size() ) );
		for ( size_t lambda = 0 ; lambda != std::max( land1.values_of_landscapes[grid_point].size() , land2.values_of_landscapes[grid_point].size() ) ; ++lambda )
		{
			double value1 = 0;
			double value2 = 0;
			if ( lambda < land1.values_of_landscapes[grid_point].size() )value1 = land1.values_of_landscapes[grid_point][lambda];
			if ( lambda < land2.values_of_landscapes[grid_point].size() )value2 = land2.values_of_landscapes[grid_point][lambda];
			result.values_of_landscapes[grid_point][lambda] = oper( value1 , value2 );
		}
	}
	
	return result;
}

Persistence_landscape_on_grid Persistence_landscape_on_grid::multiply_lanscape_by_real_number_not_overwrite( double x )const
{
	Persistence_landscape_on_grid result;
	result.values_of_landscapes = std::vector< std::vector< double > >( this->values_of_landscapes.size()  );
	result.grid_min = this->grid_min;
	result.grid_max = this->grid_max;
	
	for ( size_t grid_point = 0 ; grid_point != this->values_of_landscapes.size() ; ++grid_point )
	{
		result.values_of_landscapes[grid_point] = std::vector< double >( this->values_of_landscapes[grid_point].size() );
		for ( size_t i = 0 ; i != this->values_of_landscapes[grid_point].size() ; ++i )
		{
			result.values_of_landscapes[grid_point][i] = x*this->values_of_landscapes[grid_point][i];
		}
	}
	
	return result;
}

double compute_max_norm_discance_of_landscapes( const Persistence_landscape_on_grid& first, const Persistence_landscape_on_grid& second )
{
	double result = 0;
	
	//first we need to check if first and second is defined on the same domain"
	if ( !check_if_defined_on_the_same_domain(first, second) )throw "Two grids are not compatible"; 
	
	for ( size_t i = 0 ; i != first.values_of_landscapes.size() ; ++i )
	{
		for ( size_t j = 0 ; j != std::min( first.values_of_landscapes[i].size() , second.values_of_landscapes[i].size() ) ; ++j )
		{
			if ( result < abs( first.values_of_landscapes[i][j] - second.values_of_landscapes[i][j] ) )
			{
				result = abs( first.values_of_landscapes[i][j] - second.values_of_landscapes[i][j] );
			}
		}
		if ( first.values_of_landscapes[i].size() == std::min( first.values_of_landscapes[i].size() , second.values_of_landscapes[i].size() ) )
		{
			for ( size_t j = first.values_of_landscapes[i].size() ; j != second.values_of_landscapes[i].size() ; ++j )
			{
				if ( result < second.values_of_landscapes[i][j] )result = second.values_of_landscapes[i][j];
			}
		}
		if ( second.values_of_landscapes[i].size() == std::min( first.values_of_landscapes[i].size() , second.values_of_landscapes[i].size() ) )
		{
			for ( size_t j = second.values_of_landscapes[i].size() ; j != first.values_of_landscapes[i].size() ; ++j )
			{
				if ( result < first.values_of_landscapes[i][j] )result = first.values_of_landscapes[i][j];
			}
		}
	}  
	return result;
}



}//namespace Gudhi_stat
}//namespace Gudhi

#endif