#ifndef __USCYL_GRID_H__
#define __USCYL_GRID_H__

#ifdef __GNUC__
#include <ext/hash_map>
namespace std {
 using namespace __gnu_cxx;
}
#else
#include <hash_map>
#endif

#include "definition.h"

class StretchedCylindricalGrid {
public:
 // Constructor for equally spaced cylindrical grids
 StretchedCylindricalGrid(int jmax, blitz::Array<double,1> rdist, blitz::Array<double,1> zdist);
 //Empty Constructor to create the cylindrical mesh 
 StretchedCylindricalGrid();
 
 // Functions returning value
 int getHexCellStart();
 int getPolyCellStart();
 int getTotalFaces();
 int getTotalHexCells();
 int getTotalNodes();
 int getTotalPolyCells();
 void setZOffset(double off);
 // Functions returning void
 void debugFaceList(int i, int j, int k);
 void debugCell(int i, int j, int k);
 void Constructor(int jmax, blitz::Array<double,1> rdist, blitz::Array<double,1> zdist);
 void writeGridFile(std::string filename);
 void writeNodesAndCells(std::string filename);
 void writeUnstructuredData(std::string filename);
 void writeUnstructuredDataGzipped(std::string filename);
 void writeUnstructredGridForSolver(std::string filename);
 void writeDebugCentroid(std::string filename);

 // Depreciated Functions
 void getCellMapFromGrid(cellFaceConnectivityArray &_cellFaceConnectivityArray); // Depreacted due to logic change but just left it for usefulness sake
 void getFaceDataFromGrid(faceDataArray &_faceDataArray); // Depreacted due to logic change but just left it for usefulness sake
 void getVolumeFromGrid(blitz::Array<double, 1 > &); // Depreacted due to logic change but just left it for usefulness sake
 void outputHigherOrderNodes(std::string filename);
 ~StretchedCylindricalGrid();

private:

 // _param _imax contains the maximum # of nodes along r for the grid
 int _imax;
 // _param _jmax contains the maximum # of nodes along theta for the grid
 int _jmax;
 // _param _kmax contains the maximum # of nodes along z for the grid
 int _kmax;
 // _param _dr contains the radial spacing
 blitz::Array<double,1> _rdist;
 // _param _dtheta contains azimuthal spacing
 double _dtheta;
 // _param _dz contains the spacing along z
 blitz::Array<double,1> _zdist;
 // _param _zoffset contains the offset along z
 double _zoffset;
 // \param _cell contains pointers to the 8 nodes connected to a cell located at (i,j,k) where i,j,k are cell center locations
 cellArray _cell;
 // \param _mapFace is a blitz::array containing pointer to face and left/right cell for i,j,k coordinate of the cell
 cellFaceArray_ptr _mapFace;
 // \param _node contains the x,y,z coordinates of the node at (i,j,k) where i,j,k are nodal location
 pointArray _node;
 // \param _pcell contains pointers to the _jmax * 2 nodes connected to a cell located at (i,j,k) where i,j,k are cell center locations
 polyCellArray _pcell;
 //\param _volume contains the cell volume
 blitz::Array<double, 3 > _volume;
 // To add _cell_centroid
 //\param _cell_centroid contains the cell centroid
 blitz::Array<point, 3 > _cell_centroid;
 // \param _Face is a hash_map containing struct Face data for a given face number
 std::hash_map<int, Face> _Face;
 // \param _Faceit is a _Face iterator
 std::hash_map<int, Face>::iterator _Faceit;
 /////////////////////////////////////////////////////////////////////////
 std::hash_map<uintptr_t, int> _hashMapFace;
 std::hash_map<uintptr_t, int>::iterator _hashMapFaceit;
 ////////////////////////////////////////////////////////////////////////
 // \param _mapCells is a stl::list containing the map of memory address of face to face #
 std::hash_map<uintptr_t, int> _mapCells;
 // \param _mapCellsit is a _mapNodes iterator
 std::hash_map<uintptr_t, int>::iterator _mapCellsit;
 // \param _mapNodes is a stl::list containing the node number to the i,j,k coordinate of the node
 std::hash_map<uintptr_t, int> _mapNodes;
 // \param _mapNodesit is a _mapNodes iterator
 std::hash_map<uintptr_t, int>::iterator _mapNodesit;
 // Private Functions returing value
 int findNode(int i, int j, int k);
 int findNode(uintptr_t addr);
 int findCell(int i, int j, int k);
 int findCell(uintptr_t addr);
 int findFace(int i, int j, int k);
 int findFace(uintptr_t addr);
 double calculateVolume(point &_0, point &_1, point &_2, point &_3, point &_4, point &_5, point &_6, point &_7);
 uintptr_t correctNgbr(Face *F, uintptr_t comp);
 // Private Functions returing void
 void assignCells2Faces();
 void assignCells2Nodes();
 void calculateNodes();
 void createNodeList();
 void calculateVolume();
 void calculateFaceNormal();
 void formPolygonalCells();
 void mapCells();
 //void calculateGauss
};

#endif
