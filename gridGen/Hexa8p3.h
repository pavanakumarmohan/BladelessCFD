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

# include <tvmet/Vector.h>
# include <tvmet/Matrix.h>
# include <iostream>
# include <list>
# include <string>
# include <sstream>

typedef struct {
 tvmet::Matrix<int,8,3> permuteShapeFun;
 tvmet::Matrix<double,27,3> xi_i_j;
 tvmet::Vector<double,27> w_i;
 tvmet::Matrix<int,3,3> permute;
	tvmet::Matrix<int,3,3> unit3by3;
 
void init(){
 permuteShapeFun = -1, -1 , -1 ,	
                    1, -1 , -1 ,	
                    1,  1 , -1 ,	
                   -1,  1 , -1 ,	
																			-1, -1 ,  1 ,	
                    1, -1 ,  1 ,	
                    1,  1 ,  1 ,	
                   -1,  1 ,  1 ;

 permute = 0 , 1 , 1 ,
           1 , 0 , 1 ,
           1 , 1 , 0 ;

 unit3by3 = 1 , 0 , 0 ,
            0 , 1 , 0 ,
            0 , 0 , 1 ;


	xi_i_j =          -.7745966692414835  ,    -.7745966692414835  ,    -.7745966692414835,    
		-.7745966692414835  ,    -.7745966692414835   ,    .0000000000000000    ,
		-.7745966692414835  ,    -.7745966692414835   ,    .7745966692414835    ,
		-.7745966692414835  ,     .0000000000000000   ,   -.7745966692414835    ,
		-.7745966692414835  ,     .0000000000000000   ,    .0000000000000000    ,
		-.7745966692414835  ,     .0000000000000000   ,    .7745966692414835    ,
		-.7745966692414835  ,     .7745966692414835   ,   -.7745966692414835    ,
		-.7745966692414835  ,     .7745966692414835   ,    .0000000000000000    ,
		-.7745966692414835  ,     .7745966692414835   ,    .7745966692414835    ,
		.0000000000000000  ,    -.7745966692414835   ,   -.7745966692414835    ,
		.0000000000000000  ,    -.7745966692414835   ,    .0000000000000000    ,
		.0000000000000000  ,    -.7745966692414835   ,    .7745966692414835    ,
		.0000000000000000  ,     .0000000000000000   ,   -.7745966692414835    ,
		.0000000000000000  ,     .0000000000000000   ,    .0000000000000000    ,
		.0000000000000000  ,     .0000000000000000   ,    .7745966692414835    ,
		.0000000000000000  ,     .7745966692414835   ,   -.7745966692414835    ,
		.0000000000000000  ,     .7745966692414835   ,    .0000000000000000    ,
		.0000000000000000  ,     .7745966692414835   ,    .7745966692414835    ,
		.7745966692414835  ,    -.7745966692414835   ,   -.7745966692414835    ,
		.7745966692414835  ,    -.7745966692414835   ,    .0000000000000000    ,
		.7745966692414835  ,    -.7745966692414835   ,    .7745966692414835    ,
		.7745966692414835  ,     .0000000000000000   ,   -.7745966692414835    ,
		.7745966692414835  ,     .0000000000000000   ,    .0000000000000000    ,
		.7745966692414835  ,     .0000000000000000   ,    .7745966692414835    ,
		.7745966692414835  ,     .7745966692414835   ,   -.7745966692414835    ,
		.7745966692414835  ,     .7745966692414835   ,    .0000000000000000    ,
		.7745966692414835  ,     .7745966692414835   ,    .7745966692414835; 


	w_i =.1714677640603567,    
					.2743484224965707,    
					.1714677640603567,    
					.2743484224965707,    
					.4389574759945131,    
					.2743484224965707,    
					.1714677640603567,    
					.2743484224965707,    
					.1714677640603567,    
					.2743484224965707,    
					.4389574759945131,    
					.2743484224965707,    
					.4389574759945131,    
					.7023319615912209,    
					.4389574759945131,    
					.2743484224965707,    
					.4389574759945131,    
					.2743484224965707,    
					.1714677640603567,    
					.2743484224965707,    
					.1714677640603567,    
					.2743484224965707,    
					.4389574759945131,    
					.2743484224965707,    
					.1714677640603567,    
					.2743484224965707,    
					.1714677640603567;
}

}specialMatrices;

class Hexa8p3{
	public:
		// Empty Constructor
  Hexa8p3(); 
		// Constructor takes the 8 vertices as arguements
		Hexa8p3(tvmet::Matrix<double,8,3>&); 
		// shape function evaluation 
		tvmet::Vector<double,8> shapeFun(tvmet::Vector<double,3> &); 
		// shape function derivative
		tvmet::Matrix<double,8,3> shapeFunDer(tvmet::Vector<double,3> &); 
  std::string tecplotZone();
  std::string asciiOutput();
  std::string colocationPointOutput();
		std::string binaryOutput();
		double determinant( tvmet::Matrix<double,3,3> & );
		double vtp( tvmet::Vector<double,3> &, tvmet::Vector<double,3> &, tvmet::Vector<double,3> &); 
  double volIntegral( tvmet::Vector<double,27> &funVal);
  double volumeOfHex();
 	tvmet::Vector<double,3> volIntegral( tvmet::Matrix<double,27,3> &funVal);
		tvmet::Matrix<double,27,3> *getGaussPoints();
  void setJacobian(tvmet::Vector<double,27> &);
  void setHexNodes(tvmet::Matrix<double,8,3>&);
  // Access to special matrices like the permutation matrix, element sign matrix, etc.	
		specialMatrices spMat;

	private:
		// Shape function storage for the 27 points
		tvmet::Matrix<double,27,8> N; 
		// The coordinates of the 27 quadrature nodes in physical space (for fun eval)
		tvmet::Matrix<double,27,3> x_i_j_quad; 
		// determiant of Jacobian of the 27 nodes
		tvmet::Vector<double,27> J_i; 
		// The coordinates of the 8 nodes forming the hexa element
		tvmet::Matrix<double,8,3> x_i_j; 
};

