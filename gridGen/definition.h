#ifndef __DEFINITION_3DCYL_EULER_H__

#define __DEFINITION_3DCYL_EULER_H__

#include<fstream>
#include<iostream>
#include<blitz/array.h>
#include<sstream>
#include<cstdlib>
#include<tvmet/Vector.h>
#include<tvmet/Matrix.h>
#include<string>
#include<vector>
#include<cmath>
#include<list>
#include<map>
#include <gsl/gsl_linalg.h>


#ifdef __GNUC__
#include<stdint.h>
#define UINTPTR_MAX UINT32_MAX
#define UINT32_MAX 4294967295U
#endif

#ifdef USE_MPI
typedef struct{
public:
	void sleepPIDForGDBAttach(){
		int i = 0;
		char hostname[256];
		gethostname(hostname, sizeof(hostname));
		printf("PID %d on %s ready for attach\n", getpid(), hostname);
		fflush(stdout);
		while (0 == i)
			sleep(5);
	}
}MPIDebugHelper;
#endif

namespace IDX {
 // for state varialbes
 static int rho = 0;
 static int u = 1;
 static int v = 2;
 static int w = 3;
 static int p = 4;
 // for normals and points
 static int x = 0;
 static int y = 1;
 static int z = 2;

}

enum {
 RHO = 0, U, V, W, P
};
namespace airStandardAtmosphere {
 static const double gamma = 1.4;
 static const double gm1 = 0.4;
 static const double gbygm1 = 1.4 / 0.4;
 static const double R = 287.058; // J/kg/K
 static const double pressure = 101300; // Pa
 static const double temperature = 293.15; // K
 static const double density = pressure / (R * temperature);
 static const double speedOfSound = sqrt(gamma * pressure / density); // m/s
}

 //\param A 5x3 matrix
typedef tvmet::Matrix<double,5,3> _5x3_;
// \param flux This is the basic data type that is an array of 5 real values
typedef tvmet::Vector<double, 5 > flux;
// \param fluxArray3 is a blitz++ 3D array of Flux data
typedef blitz::Array<flux, 3 > fluxArray3;
// \param fluxArray1 is a blitz++ 1D array of Flux data
typedef blitz::Array<flux, 1 > fluxArray1;
// \param point This is the basic data type that is an array of 3 real values used for representing coordinates, normals, etc.
typedef tvmet::Vector<double, 3 > point;
// \param pointArray a blitz++ array of point data
typedef blitz::Array<point, 3 > pointArray;
// \param cell This is the basic data type to hold the pointer reference to a node of a Hexahedral cell
typedef tvmet::Vector<uintptr_t, 8 > cell;
// \param cellArray a blitz++ array of cell data
typedef blitz::Array<cell, 3 > cellArray;
// \param cellFace This is the basic data type to hold the pointer reference to a face of a Hexahedral cell
typedef tvmet::Vector<uintptr_t, 6 > cellFace;
// \param polyCell data structure storing the polygonal cell
// \param cellFaceArray a blitz++ array of cellFace constant data pointers
typedef blitz::Array<cellFace, 3 > cellFaceArray_ptr;

typedef struct polyCell {
 int nfaces; // number of faces of the polygon
 uintptr_t *faces; // memory addr pointer array - allocated based on the the # of faces
 uintptr_t *cells; // memory addr pointer array - the neighbouring cell (cells[x]) sharing the face (faces[x])
 uintptr_t *nodes; // memory addr pointer array - the nodes that form this cell
} polyCell;
// \param polyCellArray is a blitz++ one dimensional array of polyCell data structure (along the z-direction only)
typedef blitz::Array<polyCell, 1 > polyCellArray;

// \param Face is the primitive data structue of a face

typedef struct Face {
 uintptr_t Left; // pointer reference to Left neighbour
 uintptr_t Right; // pointer reference to Right neighbour
 point Normal; // nx ny and nz normal components and magintude is the area of the face
 // To add Centroid
 point Centroid;
} Face;

// \param int3 This is a basic data type used to represent index (i,j,k)
typedef tvmet::Vector<int, 3 > int3;
// \param Range This is to avid long variable names
typedef blitz::Range Range;
// \param solutionRotationMatrix This is the rotation matrix used to rotate the flux vectors as the grid rotates
typedef tvmet::Matrix<double, 5, 5 > solutionRotationMatrix;
// \param gridRotationMatrix This is the rotation matrix used to rotate the grid
typedef tvmet::Matrix<double, 3, 3 > gridRotationMatrix;


// \brief namespace numerals contain math helper functions like vector triple product and 3x3 determinant
namespace utility {
 // \brief vtp function returns the vector triple product of three vectors _1, _2 and _3 (1x3)

 inline double vtp(point &_1, point &_2, point &_3) {
  return ( _1[0] * (_2[1] * _3[2] - _2[2] * _3[1]) - _1[1] * (_2[0] * _3[2] - _2[2] * _3[0]) + _1[2] * (_2[0] * _3[1] - _2[1] * _3[0]));
 }
 // \brief determinant function returns the determinan of 3x3 matrix formed by row vectors _1, _2 and _3

 inline double determinant(point &_1, point &_2, point &_3) {
  return ( _1[0] * (_2[1] * _3[2] - _2[2] * _3[1]) - _1[1] * (_2[0] * _3[2] - _2[2] * _3[0]) + _1[2] * (_2[0] * _3[1] - _2[1] * _3[0]));
 }
 // \brief cross function returns the corss produict of two vectors formed by the clc cyclic points _1,_2 and _3

 inline point cross(point &_1, point &_2, point &_3) {
  point a = point(_2 - _1), b = point(_3 - _1);
  return point(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
 }
 // \brife mag function returns the magnitude of vector p

 inline double mag(point &p) {
  return sqrt(tvmet::dot(p, p));
 }

 inline double mag(flux &f) {
  return sqrt(tvmet::dot(f, f));
 }

 // \brief normalise function normalises a vector to unit magnitude

 inline point normalize(point vec) {
  double mag = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
  return point(vec[0] / mag, vec[1] / mag, vec[2] / mag);
 }

 inline double flux_point_dot(flux &q, point &normal) {
  return double( (q[1] * normal[0] + q[2] * normal[1] + q[3] * normal[2]));
 }

 inline double flux_point_dot(point &normal, flux &q) {
  return double( (q[1] * normal[0] + q[2] * normal[1] + q[3] * normal[2]));
 }

 inline std::map<std::string,std::string> getVariablesInNamelist(std::string filename) {
 std::ifstream namelist(filename.c_str());
 std::string testStr;
 std::string parsedStr;
 int section = 0;
 int testSection = 0;
 size_t pos;
 std::map<std::string,std::string> mapStr;

 while (!namelist.eof()) {
  // Get One line from file
  std::getline(namelist, testStr);
  // Strip all spaces
  // testStr.erase(std::remove(testStr.begin(), testStr.end(), ' '), testStr.end());
  // Check
  // For Section name
  if (*(testStr.begin()) == '&') {
   section++;
   testSection++;
  }// For Section end
  else if (*(testStr.begin()) == '/') {
   testSection--;
  }//For comment or new line at beginning of line
  else if ( testStr[size_t(0)] == '!' || testStr[size_t(0)] == '\n') {
   //std::cout << "Skipping" << std::endl;
  }// Parse the string
   else {
    pos = testStr.find('!');
    if (pos != std::string::npos) { // remove comment
     testStr = testStr.substr(0, pos);
    }
    parsedStr.append(testStr);
    parsedStr.append("\n");
   }
  }
  parsedStr.erase(parsedStr.length() - 1, 1);
  std::cout << "Processed String :" << parsedStr;

  std::istringstream parser(parsedStr);
  std::string first, second;

  while (!parser.eof()) {
   getline(parser, testStr);
   pos = testStr.find('=');
   first = testStr.substr(0, pos);
   second = testStr.substr(pos + 1, testStr.length() - 1);
   mapStr[first] = second;
  }
  return mapStr;
 }
 
template<class T> class Interpolator{
  public:
/*       Linear interpolation
        target  - the target point, 0.0 - 1.0
        v       - a pointer to an array of size 2 containg the two values
*/
  inline static T Linear(double target, T *v){
   return (T) ( target * (v[1]) + ( T(1.0f) - target ) * (v[0]) );
  }
/* 
        BiLinear interpolation, linear interpolation in 2D
        target  - a 2D point (X,Y)
        v       - an array of size 4 containg values cockwise around the square starting from bottom left

        cost: performs 3 linear interpolations
    */
  inline static T Bilinear(double *target, T *v){
   T v_prime[2] = { Linear(target[1], &(v[0])), Linear(target[1], &(v[2])) };      
   return Linear( target[0] , v_prime );
  }
  /* 
     TriLinear interpolation, linear interpolation in 2D
     target  - a 3D point (X,Y,Z)
     v       - an array of size 8 containg the values of the 8 corners 
               of a cube defined as two faces: 0-3 face one (front face) 
                                               4-7 face two (back face)

    cost: 7 linear interpolations
   */
  inline static T Trilinear(double *target, T *v){
   T v_prime[2] = { Bilinear(&(target[0]), &(v[0])) , Bilinear(&(target[1]), &(v[4])) };
   return Linear( target[2] , v_prime );
  }
 };

 // LU Decompostion to solve the linear algebra problem Ax = b for a 4x4 system
 // b is a matrix of size (4,A.extent(0)) 4 is the # of flux variables
 inline void ludecomp( blitz::Array<double,2> &A, blitz::Array<double,2> &b ){
  if(A.extent(0) != A.extent(1)){
   std::cout << "Error: Matrix not square" << "\n";
  }
  if( b.extent(1) != A.extent(0) || b.extent(0) != 5 ){
   std::cout << "Error: Matrix not square" << "\n";
  }
  gsl_matrix_view m = gsl_matrix_view_array ( A.data(), A.extent(0), A.extent(1) );
  gsl_vector_view B = gsl_vector_view_array ( b.data(), A.extent(0) );
  gsl_linalg_cholesky_decomp(&m.matrix);
		//std::cout << "LS Matrix = " << A << "\n";
  // b1
  gsl_linalg_cholesky_svx(&m.matrix, &B.vector);
  // b2
  B = gsl_vector_view_array ( b.data() + A.extent(0) , A.extent(0) );
  gsl_linalg_cholesky_svx(&m.matrix, &B.vector);
  //b3
  B = gsl_vector_view_array ( b.data() + 2 * A.extent(0) , A.extent(0) );
  gsl_linalg_cholesky_svx(&m.matrix, &B.vector);
  //b4
  B = gsl_vector_view_array ( b.data() + 3 * A.extent(0) , A.extent(0) );
  gsl_linalg_cholesky_svx(&m.matrix, &B.vector);
  //b5
  B = gsl_vector_view_array ( b.data() + 4 * A.extent(0) , A.extent(0) );
  gsl_linalg_cholesky_svx(&m.matrix, &B.vector);
 }
}

// \brief namespace numerals contain numeric constants which are regularly used in the code
namespace numerals {
 // \param one_by_2 contains the value 1/2
 static const double one_by_2 = 1.0 / 2.0;
 // \param one_by_6 contains the value 1/6
 static const double one_by_6 = 1.0 / 6.0;
 // \param kap2 second order Jameson type dissipation coefficient
 static const double kap2 = 1.0 / 4.0;
 // \param kap2 fourth order Jameson type dissipation coefficient
 static const double kap4 = 1.0 / 256.0;
 // \param f_q These are the Runge-Kutta 4th order scheme parameters
 static const tvmet::Vector<double, 5 > f_q(0.5, 0.5, 1.0, 1.0);
 static const tvmet::Vector< double , 2 > rk2_f_q( 1.0 , 1.0 );
 // \param fac These are the Runge-Kutta 4th order scheme parameters
 static const tvmet::Vector<double, 5 > fac(1.0e0 / 6.0e0, 1.0e0 / 3.0e0, 1.0e0 / 3.0e0, 1.0e0 / 6.0e0);
 static const tvmet::Vector<double, 2 > rk2_fac( 0.5 , 0.5 );
 static const uintptr_t zero_ptr = UINTPTR_MAX;
 static const double _2_MPI = 2.0 * M_PI;
}

// \param These are solver specific so declared them as seperate typedef's

typedef struct cellFaceConnectivity {
 std::list<int> faceConnectivity;
 std::list<int> cellConnectivity;
 double volume;
} cellFaceConnectivity;
typedef blitz::Array<cellFaceConnectivity, 1 > cellFaceConnectivityArray;

typedef struct faceData {
 point normal;
 point centroid;
 int32_t leftCell;
 int32_t rightCell;
} faceData;
typedef blitz::Array<faceData, 1 > faceDataArray;

namespace source {

 typedef struct cell {
  std::vector<int32_t> nodeConnectivity;
  point centroid;
  double weight;
 } cell;
 // \param polyCellArray is a blitz++ one dimensional array of polyCell data structure (along the z-direction only)
 typedef blitz::Array<cell, 1 > cellArray;
 typedef blitz::Array<point, 1 > pointArray;
}

namespace helperSolver {

 inline double enthalpy(flux &F) {
  return airStandardAtmosphere::gbygm1 * F[IDX::p] / F[IDX::rho];
 }

 inline double energy(flux & F) {
  return F[IDX::p] / (F[IDX::rho] * airStandardAtmosphere::gm1);
 }

 inline double temperature(flux & F) {
  return double( F[IDX::p] / (F[IDX::rho] * airStandardAtmosphere::R));
 }

 inline flux normalFlux(flux &q, point &normal) {
  flux temp;
  point vel = point(q[IDX::u], q[IDX::v], q[IDX::w]);
  double h0 = enthalpy(q) + 0.5 * tvmet::dot(vel, vel);
  double vn = tvmet::dot(vel, normal);
  temp[0] = q[IDX::rho] * vn;
  temp[1] = q[IDX::rho] * q[IDX::u] * vn + q[IDX::p] * normal[0];
  temp[2] = q[IDX::rho] * q[IDX::v] * vn + q[IDX::p] * normal[1];
  temp[3] = q[IDX::rho] * q[IDX::w] * vn + q[IDX::p] * normal[2];
  temp[4] = q[IDX::rho] * h0 * vn;
  //  temp[0] = q[IDX::rho] * vn;
  //  temp[1] = q[IDX::rho] * q[IDX::u] * vn;
  //  temp[2] = q[IDX::rho] * q[IDX::v] * vn;
  //  temp[3] = q[IDX::rho] * q[IDX::w] * vn;
  //  temp[4] = q[IDX::rho] * h0 * vn;
  return temp;
 }

 inline double velocityDot(flux &q1, flux &q2) {
  return double(q1[IDX::u] * q2[IDX::u] + q1[IDX::v] * q2[IDX::v] + q1[IDX::w] * q2[IDX::w]);
 }

 inline double normalVelocity(flux &q, point & n) {
  return double(q[IDX::u] * n[IDX::x] + q[IDX::v] * n[IDX::y] + q[IDX::w] * n[IDX::z]);
 }

}

namespace bladeObject {
 typedef tvmet::Vector<double, 2 > pair;
 typedef std::pair< double, blitz::Array<point, 1 > > loadPair;

 typedef struct {

  void addLoad(blitz::Array<point, 1 > &load, double t) {
   if (load.size() != nSpan) {
    std::cerr << "Incomaptible load array size = " << load.size() << std::endl;
   }
   loads.insert(loadPair(t, load));
  }

  void setSpan(int n) {
   nSpan = n;
  }

  point find(int loc, double t) {
   //std::cout << "Finding Time Key :" << t << std::endl;
   std::map< double, blitz::Array<point, 1 > >::iterator it = loads.find(t);
   if (it == loads.end()) {
    std::cerr << "Unknown time key " << std::endl;
    exit(1);
   }
   //std::cout << "Found Load for the key " << std::endl;
   //std::cout << "data:" << it->second << std::endl;
   return (it->second)(loc);
  }

 std::map< double, blitz::Array<point, 1 > >::iterator begin(){
  return loads.begin();
 }

 std::map< double, blitz::Array<point, 1 > >::iterator end(){
  return loads.end();
 }
 
 private:
  int nSpan;
  std::map< double, blitz::Array<point, 1 > > loads;
 } wopwopLoads;

 typedef struct {
  // intrisic types
  bool tipLossFactorFlag;
  bool tipJetFlag;
  bool bladeDefFlag;
  char InflowFileName[256];
  char BEMTFileName[256];
  int nBlade;
  int nLE;
  int nTE;
  int nSpan;
  double psiZero;
  double axisValue;
  double CT;
  double dPsi;
  double Solidity;
  double omega;
  double rootChord;
  double rootCutOut;
  double rotorRadius;
  double CollPitch;
  double FWHM_Z;
  double FWHM_Y;
  double area;
  // Objects
  point location;
  point axis;
  std::vector<pair> taperRatio;
  std::vector<pair> twist;
  std::vector<pair> sweep;
  std::vector<pair> flex;
  std::vector<pair> offset;
  std::vector<pair> dx;
  // Stores Time history of loads
  wopwopLoads *load;
  // Stores Current time blade loads
  blitz::Array<point, 1 > *interpLoad;
 } BEMT_Rotor_Blade;


}

#endif

