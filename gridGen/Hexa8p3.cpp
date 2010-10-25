/*
    8 noded Hexahedral element 
    shape function, Jacobian function, 
    Gaussian quadrature formulas for 
    volume integration

    Copyright (C) 2010  Pavanakumar Mohanamuraly

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.s

*/

# include "Hexa8p3.h"

Hexa8p3::Hexa8p3(){
	spMat.init();
}

Hexa8p3::Hexa8p3(tvmet::Matrix<double,8,3>& x_i){
 spMat.init();
 setHexNodes(x_i);
}

void Hexa8p3::setHexNodes(tvmet::Matrix<double,8,3>& x_i){
 J_i = 0.0;
 x_i_j_quad = 0.0;
 x_i_j = 0.0;
 N = 0.0;
   
 // Copy these values
	x_i_j = x_i; 
 tvmet::Vector<double,8> N_ret;
 tvmet::Vector<double,3> temp_point;
 // Step 1
 // Calculate shape functions 
	for(int i = 0 ; i < 27 ; ++i ){
  temp_point = tvmet::row( spMat.xi_i_j , i );
  N_ret = shapeFun( temp_point );
		for(int j = 0 ; j < 8 ; ++j ){
   N(i,j) = N_ret(j);
  }
 }
 // Find out quadrature points in physical space using shape funtion
	x_i_j_quad = N * x_i_j;
	// Step 2 
	// Calculate Jacobian
	tvmet::Matrix<double,3,3> jacobian(0.0);
	tvmet::Matrix<double,8,3> shape_der(0.0);
	for(int i = 0 ; i < 27 ; ++i ){
		// Calculate the derivative of N for each quadrature point and store it in a temp variable
		temp_point = tvmet::row( spMat.xi_i_j , i );
		shape_der = shapeFunDer( temp_point );
  jacobian = tvmet::trans(shape_der) * x_i_j;
		// Store the jacobian determinant in J_i
		J_i(i) = determinant(jacobian);
		jacobian = 0.0;
	}
}

tvmet::Vector<double,8> Hexa8p3::shapeFun(tvmet::Vector<double,3> &xi_i){
	tvmet::Vector<double,8> N_ret(1.0);
	for(int i = 0 ; i < 8 ; ++i )
		for(int j = 0 ; j < 3 ; ++j )
			N_ret(i) *= 0.5 * ( 1.0 + spMat.permuteShapeFun(i,j) * xi_i(j) );
	return N_ret;
}

tvmet::Matrix<double,8,3> Hexa8p3::shapeFunDer(tvmet::Vector<double,3> &xi_i){
	tvmet::Matrix<double,8,3> dN_ret(1.0);
	for(int j = 0 ; j < 8 ; ++j ){
		for(int k = 0 ; k < 3 ; ++k ){
			for(int i = 0 ; i < 3 ; ++i ){
				dN_ret(j,k) *= 0.5 * ( 1 + spMat.permute(k,i) * xi_i(i) * spMat.permuteShapeFun(j,i) );
   }
			dN_ret(j,k) *= spMat.permuteShapeFun(j,k);
  }
  //////////////////////////////////////////
  //          DEBUG                    ////
  //////////////////////////////////////////
  /*std::cout << "dN = " << "1/8 * ";
		for(int k = 0 ; k < 3 ; ++k ){
			for(int i = 0 ; i < 3 ; ++i ){
				std::cout << "( 1 + " << spMat.permute(k,i) << "*" << xi_i(i) << "*" << spMat.permuteShapeFun(j,i) << ")";
   }
			std::cout << "*" << spMat.permuteShapeFun(j,k) << std::endl;
  }*/
   //////////////////////////////////////////
  //          DEBUG                    ////
  //////////////////////////////////////////
 }
	return dN_ret;
}

std::string Hexa8p3::tecplotZone(){
 std::stringstream ss;
 ss << "VARIABLES = \"X\", \"Y\", \"Z\"\nZONE NODES=8, ELEMENTS=1, DATAPACKING=POINT, ZONETYPE=FEBRICK";
	for(int i = 0 ; i < 8 ; ++i){
		ss << "\n";
		for(int j = 0 ; j < 3 ; ++j){
			ss << x_i_j(i,j) << " ";
		}
	}
	for(int i = 1 ; i < 9 ; ++i)
		ss << i << " ";
 // Quadrature nodes in physical space
	ss << "\nZONE";
	for(int i = 0 ; i < 27 ; ++i){
		ss << "\n";
		for(int j = 0 ; j < 3 ; ++j){
			ss << x_i_j_quad(i,j) << " ";
		}
	}
 // Quadrature nodes in mapped space
	ss << "\nZONE";
	for(int i = 0 ; i < 27 ; ++i){
		ss << "\n";
		for(int j = 0 ; j < 3 ; ++j){
			ss << spMat.xi_i_j(i,j) << " ";
		}
	}
 // Quadrature nodes in mapped space
	ss << "\nZONE NODES=8, ELEMENTS=1, DATAPACKING=POINT, ZONETYPE=FEBRICK";
	for(int i = 0 ; i < 8 ; ++i){
		ss << "\n";
		for(int j = 0 ; j < 3 ; ++j){
			ss << spMat.permuteShapeFun(i,j) << " ";
		}
	}
 ss << "\n" << "1 2 3 4 5 6 7 8";
	return ss.str();
}

std::string Hexa8p3::asciiOutput(){
	std::stringstream ss;
	for(int i = 0 ; i < 27 ; ++i){
  ss << x_i_j_quad(i,0) << " " << x_i_j_quad(i,1) << " " << x_i_j_quad(i,2) << " " << J_i(i) << "\n";
	}
	return ss.str();
}

std::string Hexa8p3::colocationPointOutput(){
	std::stringstream ss;
	for(int i = 0 ; i < 27 ; ++i){
  ss << x_i_j_quad(i,0) << " " << x_i_j_quad(i,1) << " " << x_i_j_quad(i,2) << "\n";
	}
	return ss.str();
}

double Hexa8p3::determinant( tvmet::Matrix<double,3,3> &A ) {
  tvmet::Vector<double,3> _1(row(A,0));
  tvmet::Vector<double,3> _2(row(A,1));
  tvmet::Vector<double,3> _3(row(A,2));
  return ( _1[0] * (_2[1] * _3[2] - _2[2] * _3[1]) - 
				_1[1] * (_2[0] * _3[2] - _2[2] * _3[0]) + 
				_1[2] * (_2[0] * _3[1] - _2[1] * _3[0]));
}

double Hexa8p3::vtp( tvmet::Vector<double,3> &_1 , tvmet::Vector<double,3> &_2 , tvmet::Vector<double,3> &_3 ) {
  return ( _1[0] * (_2[1] * _3[2] - _2[2] * _3[1]) - 
				_1[1] * (_2[0] * _3[2] - _2[2] * _3[0]) + 
				_1[2] * (_2[0] * _3[1] - _2[1] * _3[0]));
}


tvmet::Matrix<double,27,3> *Hexa8p3::getGaussPoints(){
return &x_i_j_quad;
}

double Hexa8p3::volIntegral( tvmet::Vector<double,27> &funVal){
 double integral_val = 0.0;
 for(int i = 0 ; i < 27 ; ++i )
  integral_val += spMat.w_i(i) * J_i(i) * funVal(i);
 return integral_val;
}
	
tvmet::Vector<double,3> Hexa8p3::volIntegral( tvmet::Matrix<double,27,3> &funVal){
 tvmet::Vector<double,3> integral_val(0.0);
	for(int dimension = 0 ; dimension < 3 ; ++dimension )
		for(int i = 0 ; i < 27 ; ++i )
		integral_val(dimension) += spMat.w_i(i) * J_i(i) * funVal(i,dimension);
 integral_val /= tvmet::sum(J_i);
 return integral_val;
}

void Hexa8p3::setJacobian( tvmet::Vector<double,27> &jack ){
 J_i = jack; 
}

double Hexa8p3::volumeOfHex(){
	tvmet::Vector<double,3> y0, y1, y2, y3, y4, y5, y6;
	y0 = tvmet::row(x_i_j,6) - tvmet::row(x_i_j,0);
	y1 = tvmet::row(x_i_j,1) - tvmet::row(x_i_j,0);
	y2 = tvmet::row(x_i_j,2) - tvmet::row(x_i_j,5);
	y3 = tvmet::row(x_i_j,4) - tvmet::row(x_i_j,0);
	y4 = tvmet::row(x_i_j,5) - tvmet::row(x_i_j,7);
	y5 = tvmet::row(x_i_j,3) - tvmet::row(x_i_j,0);
	y6 = tvmet::row(x_i_j,7) - tvmet::row(x_i_j,2);
	return 1.0 / 6.0 * ( vtp( y0, y1, y2 ) +	vtp( y0, y3, y4 ) + vtp( y0, y5, y6 ) );
}

