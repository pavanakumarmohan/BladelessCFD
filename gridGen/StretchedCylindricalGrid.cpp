#include "StretchedCylindricalGrid.h"
#include "Hexa8p3.h"

#include<fstream>
#include <sstream>

/*
    \brief Empty Constructor to create the cylindrical mesh 
*/ 
StretchedCylindricalGrid::StretchedCylindricalGrid() { }
 
/*
    \brief Constructor to create the cylindrical mesh 
    \parm imax contains the maximum # of nodes along r direction
    \parm jmax contains the maximum # of nodes along \theta direction
    \parm kmax contains the maximum # of nodes along z direction
    \parm rmax is the maximum distance of the boundary location in r direction
    \parm zmax is the maximum distance of the boundary location in z direction

    Calculates the following :
    1) Node coordinates -> _node
    2) A cell data structure that contains the link to the 6 nodes connected to it -> cell
    3) cell volume of each cell -> cellVolume

    for a grid of imax,jmax and kmax there are 
    1) imax ( +1 for ghost point ) cell centres along r (there is no rmin=0)
    2) jmax ( no need for periodic point ) cell centres along theta ( preiodic ) 
    3) kmax + 1 ( +2 for ghost point ) cell centres along z ( two boundaries zmax and zmin )

    for a grid of imax,jmax and kmax there are 
    1) imax + 1 ( +1 for ghost node and  +1 for boundary node ) nodes along r (there is no rmin=0)
    2) jmax + 1 ( +1 for periodic node ) nodes along theta ( periodic ) 
    3) kmax + 2 ( +2 for ghost node ) nodes along z ( two boundaries zmax and zmin )

 */
StretchedCylindricalGrid::StretchedCylindricalGrid(int jmax, blitz::Array<double,1> rdist, blitz::Array<double,1> zdist) {
 Constructor(jmax, rdist, zdist);
}

/*
    \brief Constructor to create the cylindrical mesh 
    \parm imax contains the maximum # of nodes along r direction
    \parm jmax contains the maximum # of nodes along \theta direction
    \parm kmax contains the maximum # of nodes along z direction
    \parm rmax is the maximum distance of the boundary location in r direction
    \parm zmax is the maximum distance of the boundary location in z direction

    Calculates the following :
    1) Node coordinates -> _node
    2) A cell data structure that contains the link to the 6 nodes connected to it -> cell
    3) cell volume of each cell -> cellVolume

    for a grid of imax,jmax and kmax there are 
    1) imax ( +1 for ghost point ) cell centres along r (there is no rmin=0)
    2) jmax ( no need for periodic point ) cell centres along theta ( preiodic ) 
    3) kmax + 1 ( +2 for ghost point ) cell centres along z ( two boundaries zmax and zmin )

    for a grid of imax,jmax and kmax there are 
    1) imax + 1 ( +1 for ghost node and  +1 for boundary node ) nodes along r (there is no rmin=0)
    2) jmax + 1 ( +1 for periodic node ) nodes along theta ( periodic ) 
    3) kmax + 2 ( +2 for ghost node ) nodes along z ( two boundaries zmax and zmin )

 */
void StretchedCylindricalGrid::Constructor(int jmax, blitz::Array<double,1> rdist, blitz::Array<double,1> zdist) {
 // Assign the imax,jmax and kmax information to the class
 _imax = rdist.size() - 1;
 _jmax = jmax;
 _kmax = zdist.size() - 2;
 _dtheta = 2.0 * M_PI / (jmax);
 _rdist.resize(rdist.size());
 _zdist.resize(zdist.size());
 _zdist = zdist;
 _rdist = rdist;
 //_zoffset = 0.0;

 /* Memory Allocation for the Nodes and Cell data */
 _node.resize(Range(1, _imax + 1), Range(1, _jmax + 1), Range(1, _kmax + 2));
 _cell.resize(Range(1, _imax), Range(1, _jmax), Range(1, _kmax + 1));
 _pcell.resize(Range(1, _kmax + 1));
 _volume.resize(Range(1, _imax), Range(1, _jmax), Range(1, _kmax + 1));
 _cell_centroid.resize(Range(1, _imax), Range(1, _jmax), Range(1, _kmax + 1));
 _mapFace.resize(Range(1, _imax), Range(1, _jmax), Range(1, _kmax + 1));

 _mapFace = cellFace(UINTPTR_MAX);

 mapCells(); // Map memroy address of cells to cell numbers assigned
 calculateNodes(); // Using the dr, dtheta and dz information we form the nodes at i,j,k
 assignCells2Nodes(); // The formed nodes are then assigned to cells
 createNodeList(); // A hasp_map of Nodes for each i,j,k pointing to a node number is created
 calculateVolume(); // The volume of the hexahedral elements are calculated
 assignCells2Faces(); // Unique faces are formed and they are assiged to the respective cells
 formPolygonalCells(); // The polygonal centre volume is created and the faces and area are also calculated

 // Constructor Stat's (for debug purpose can be removed if too trouble some)
 std::cout << "\nCreated Unstructured Cylindrical grid\n";
 std::cout << "-------------------------------------\n";
 std::cout << "Total Nodes         = " << (_imax + 1)*(_jmax)*(_kmax + 2) << std::endl;
 std::cout << "Total Elements      = " << _imax * _jmax * (_kmax + 1) << std::endl;
 std::cout << "Total Faces Created = " << _Face.size() << std::endl;
 std::cout << "dtheta      = " << _dtheta << std::endl;
 std::cout << "imax, jmax, kmax    = " << _imax << " " << _jmax << " " << _kmax << std::endl;
 std::cout << "---------------- END ---------------------" << std::endl;
}

/*
\brief This function sets the offset along Z direction
 */
void StretchedCylindricalGrid::setZOffset(double off) {
 _zoffset = off;
}

/*
\brief This function calculates the nodes and stores it in the _node array of the class.
       Note that the loop is over the nodes and not cell centers
 */
void StretchedCylindricalGrid::calculateNodes() {
 for (int j = 1; j <= _jmax; ++j) {
  for (int i = 1; i <= _imax + 1; ++i) {
   for (int k = 1; k <= _kmax + 2; ++k) {
    _node(i, j, k) = point( _rdist(i - 1) * cos(j * _dtheta), _rdist(i - 1) * sin(j * _dtheta), _zdist(k - 1));
   }
  }
  // Take care of periodic node
  for (int i = 1; i <= _imax + 1; ++i) {
   for (int k = 1; k <= _kmax + 2; ++k) {
    _node(i, _jmax + 1, k) = _node(i, 1, k);
   }
  }
 }
}

/*
\brief This function assigns the nodes calculated from calculateNodes() to the respective cells.
       Note that the loop is over the cell centers and not nodes
*/
void StretchedCylindricalGrid::assignCells2Nodes() {
 // Loop through all cell centres (alos includes the ghost and centre-line degenerate triangles)
 for (int i = 1; i <= _imax; ++i) {
  for (int j = 1; j <= _jmax; ++j) {
   for (int k = 1; k <= _kmax + 1; ++k) {
    _cell(i, j, k)[0] = (uintptr_t)&(_node(i, j, k)); // Node 1
    _cell(i, j, k)[1] = (uintptr_t)&(_node(i + 1, j, k)); // Node 2
    _cell(i, j, k)[2] = (uintptr_t)&(_node(i + 1, j + 1, k)); // Node 3
    _cell(i, j, k)[3] = (uintptr_t)&(_node(i, j + 1, k)); // Node 4
    _cell(i, j, k)[4] = (uintptr_t)&(_node(i, j, k + 1)); // Node 5
    _cell(i, j, k)[5] = (uintptr_t)&(_node(i + 1, j, k + 1)); // Node 6
    _cell(i, j, k)[6] = (uintptr_t)&(_node(i + 1, j + 1, k + 1)); // Node 7
    _cell(i, j, k)[7] = (uintptr_t)&(_node(i, j + 1, k + 1)); // Node 8
   }
  }
 }
}

int StretchedCylindricalGrid::getTotalHexCells() {
 return (_imax - 1) * _jmax * (_kmax + 1);
}

int StretchedCylindricalGrid::getTotalPolyCells() {
 return _kmax + 1;
}

int StretchedCylindricalGrid::getTotalNodes() {
 return (_imax + 1) * _jmax * (_kmax + 2);
}

int StretchedCylindricalGrid::getTotalFaces() {
 return _Face.size() - 1;
}

int StretchedCylindricalGrid::getHexCellStart() {
 return 1;
}

int StretchedCylindricalGrid::getPolyCellStart() {
 return _mapCells.size() - _kmax;
}

/*
\brief Write the unstructured grid to file
\param filename type(std::string) contains the filename string to write to
 */
void StretchedCylindricalGrid::writeUnstructuredData(std::string filename) {
 std::ofstream fout(filename.c_str());
 uintptr_t addr;
 if (fout.fail()) {
  std::cerr << "Error : Unable to open file " << filename << std::endl;
  exit(-1);
 }
 std::cout << "\nSuccessfully wrote Cylindrical grid to disk\n";
 std::cout << "-------------------------------------------\n";
 std::cout << "Total Nodes    = " << (_imax + 1)*(_jmax)*(_kmax + 2) << std::endl;
 std::cout << "Total Elements = " << (_imax - 1) * _jmax * (_kmax + 1) << std::endl;
 std::cout << "File Name      = " << filename << std::endl;
 std::cout << "File Format    = Tecplot Unstructured (ASCII)" << std::endl;

 fout << "TITLE = \"FE-Volume Grid\"\nFILETYPE=GRID\nVARIABLES = \"X\", \"Y\", \"Z\"\n";
 fout << "ZONE NODES=" << (_imax + 1)*(_jmax)*(_kmax + 2) << ", ELEMENTS=" << (_imax - 1) * _jmax * (_kmax + 1) << ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK\n";
 // Write all nodes
 for (int i = 1; i <= _imax + 1; ++i) {
  for (int j = 1; j <= _jmax; ++j) {
   for (int k = 1; k <= _kmax + 2; ++k) {
    fout << _node(i, j, k)[0] << "\t";
   }
   fout << "\n";
  }
 }
 for (int i = 1; i <= _imax + 1; ++i) {
  for (int j = 1; j <= _jmax; ++j) {
   for (int k = 1; k <= _kmax + 2; ++k) {
    fout << _node(i, j, k)[1] << "\t";
   }
   fout << "\n";
  }
 }
 for (int i = 1; i <= _imax + 1; ++i) {
  for (int j = 1; j <= _jmax; ++j) {
   for (int k = 1; k <= _kmax + 2; ++k) {
    fout << _node(i, j, k)[2] + _zoffset << "\n";
   }
   fout << "\n";
  }
 }
 // Write all cells
 for (int i = 2; i <= _imax; ++i) {
  for (int j = 1; j <= _jmax; ++j) {
   for (int k = 1; k <= _kmax + 1; ++k) {
    for (int l = 0; l < 8; ++l) {
     addr = (uintptr_t) _cell(i, j, k)[l];
     fout << findNode(addr) << "\t";
    }
    fout << "\n";
   }
  }
 }
 // Write Centerline solution as separate zone
	fout << "ZONE I=" << _kmax + 1 << "\n";
	for (int k = 1; k <= _kmax + 1; ++k) 
  fout << k << " " << k << " " << k << "\n";
}


/*
\brief Generates a stl::list of node numbers that maps to the node at (i,j,k)
 */
void StretchedCylindricalGrid::createNodeList() {
 int dummy = 1;
 uintptr_t addr;
 for (int i = 1; i <= _imax + 1; ++i) {
  for (int j = 1; j <= _jmax; ++j) {
   for (int k = 1; k <= _kmax + 2; ++k) {
    addr = (uintptr_t) & _node(i, j, k);
    _mapNodes[addr] = dummy;
    dummy++;
   }
  }
 }
 for (int i = 1; i <= _imax + 1; ++i) {
  for (int k = 1; k <= _kmax + 2; ++k) {
   addr = (uintptr_t) & _node(i, _jmax + 1, k); //
   _mapNodes[addr] = findNode(i, 1, k);
  }
 }
}

/*
\brief Searches the stl::hasp_map of {memory address <=> node numbers} and return the node number
\param i nodal coordianate r direction
\param j nodal coordianate theta direction
\param k nodal coordianate z direction
 */
int StretchedCylindricalGrid::findNode(int i, int j, int k) {
 uintptr_t addr = (uintptr_t) & _node(i, j, k);
 _mapNodesit = _mapNodes.find(addr);
 return _mapNodesit->second;
}

/*
\brief Overloaded Function findNode Searches the stl::hash_map of {memory address <=> node numbers} and return the node number
\param addr memory address of the node as uintptr_t
 */
int StretchedCylindricalGrid::findNode(uintptr_t addr) {
 _mapNodesit = _mapNodes.find(addr);
 return _mapNodesit->second;
}

/*
\brief Searches the stl::hasp_map of {memory address <=> node numbers} and return the cell number
\param i nodal coordianate r direction
\param j nodal coordianate theta direction
\param k nodal coordianate z direction
 */
int StretchedCylindricalGrid::findCell(int i, int j, int k) {
 uintptr_t addr = (uintptr_t) & _cell(i, j, k);
 _mapCellsit = _mapCells.find(addr);
 return _mapCellsit->second;
}

/*
\brief Overloaded Function findCell Searches the stl::hash_map of {memory address <=> node numbers} and return the cell number
\param addr memory address of the node as uintptr_t
 */
int StretchedCylindricalGrid::findCell(uintptr_t addr) {
 _mapCellsit = _mapCells.find(addr);
 return _mapCellsit->second;
}

/*
\brief Overloaded Function findFace Searches the stl::hash_map of {memory address <=> node numbers} and return the node number
\param addr memory address of the node as uintptr_t
 */
int StretchedCylindricalGrid::findFace(uintptr_t addr) {
 _hashMapFaceit = _hashMapFace.find(addr);
 return _hashMapFaceit->second;
}

/*
\brief Overloaded Function findNode Searches the stl::hash_map of {memory address <=> node numbers} and return the node number
\param addr memory address of the node as uintptr_t
 */
uintptr_t StretchedCylindricalGrid::correctNgbr(Face *F, uintptr_t comp) {
 if (F->Right == comp) return F->Left;
 else return F->Right;
}

/*
\brief Function calculateVolume calculates the volume of a general Hexahedron using vector triple product for all cells and stores it in _volume
 */
void StretchedCylindricalGrid::calculateVolume() {
 // Loop through all cell centres (also includes the ghost)
 point y0, y1, y2, y3, y4, y5, y6;
 for (int i = 1; i <= _imax; ++i) {
  for (int j = 1; j <= _jmax; ++j) {
   for (int k = 1; k <= _kmax + 1; ++k) {
    // Get access to the cell nodes
    y0 = *((point *) _cell(i, j, k)[6]) - *((point *) _cell(i, j, k)[0]);
    y1 = *((point *) _cell(i, j, k)[1]) - *((point *) _cell(i, j, k)[0]);
    y2 = *((point *) _cell(i, j, k)[2]) - *((point *) _cell(i, j, k)[5]);
    y3 = *((point *) _cell(i, j, k)[4]) - *((point *) _cell(i, j, k)[0]);
    y4 = *((point *) _cell(i, j, k)[5]) - *((point *) _cell(i, j, k)[7]);
    y5 = *((point *) _cell(i, j, k)[3]) - *((point *) _cell(i, j, k)[0]);
    y6 = *((point *) _cell(i, j, k)[7]) - *((point *) _cell(i, j, k)[2]);
    _volume(i, j, k) = numerals::one_by_6 * (utility::vtp(y0, y1, y2) + utility::vtp(y0, y3, y4) + utility::vtp(y0, y5, y6));
    /******************************************************************** */
    // To add -- calculate the cell centroids
    _cell_centroid(i,j,k) = 0.0;
    for(int node_count = 0 ; node_count < 8 ; ++node_count )
	    _cell_centroid(i,j,k) += *((point *) _cell(i, j, k)[node_count]);
    _cell_centroid(i,j,k) *= 0.125; // There are 8 nodes
    /******************************************************************** */
   }
  }
 }
}

/*
\brief Overloaded Function calculateVolume calculates the volume of a general Hexahedron using vector triple product for given points _1 through _8
 */
double StretchedCylindricalGrid::calculateVolume(point &_0, point &_1, point &_2, point &_3, point &_4, point &_5, point &_6, point &_7) {
 // Loop through all cell centres (also includes the ghost)
 point y0, y1, y2, y3, y4, y5, y6;
 // Get access to the cell nodes
 y0 = _6 - _0;
 y1 = _1 - _0;
 y2 = _2 - _5;
 y3 = _4 - _0;
 y4 = _5 - _7;
 y5 = _3 - _0;
 y6 = _7 - _2;
 return double(numerals::one_by_6 * (utility::vtp(y0, y1, y2) + utility::vtp(y0, y3, y4) + utility::vtp(y0, y5, y6)));
}

/*
\brief assignCells2Faces creates the face heap and assigns the face from that heap to the cells ( long and dirty so could have lots of bugs! needs big change )
 * Notes:
 * 1) 0 -> i+1 face
 * 2) 1 -> i-1 face
 * 3) 2 -> j+1 face
 * 4) 3 -> j-1 face
 * 5) 4 -> k+1 face
 * 6) 5 -> k-1 face
 */
// TO DO: Correct the kmax + 1 and i = 1 cell connectivity to pcell instead of cell

void StretchedCylindricalGrid::assignCells2Faces() {
 // loop through cells to find the faces and assign the cells to the left and right of the face
 Face tempFace;
 int counter = 1;
 double ds;
 point _0, _1, _2, _3, _4, _5, _6, _7, unit_vector = point(1.0); // _0 through _7 are the 8 nodes of a hex cell (temp storage)

 // Periodic boundary cells (special case)
 for (int i = 2; i <= _imax - 1; ++i) {
  for (int k = 2; k <= _kmax; ++k) {
   // Assign the appropriate nodes
   _0 = *((point *) _cell(i, _jmax, k)[0]);
   _1 = *((point *) _cell(i, _jmax, k)[1]);
   _2 = *((point *) _cell(i, _jmax, k)[2]);
   _3 = *((point *) _cell(i, _jmax, k)[3]);
   _4 = *((point *) _cell(i, _jmax, k)[4]);
   _5 = *((point *) _cell(i, _jmax, k)[5]);
   _6 = *((point *) _cell(i, _jmax, k)[6]);
   _7 = *((point *) _cell(i, _jmax, k)[7]);
   // Push the face to face heap and then assign it to the left and right cell
   // the face node numbering sequence is 2-3-7-6
   tempFace.Left = findCell((uintptr_t) & _cell(i, _jmax, k));
   tempFace.Right = findCell((uintptr_t) & _cell(i, 1, k));
   point X(_2[0], _3[0], _7[0]), Y(_2[1], _3[1], _7[1]), Z(_2[2], _3[2], _7[2]);
   ds = 0.5 * sqrt(pow(utility::determinant(X, Y, unit_vector), 2) +
           pow(utility::determinant(Y, Z, unit_vector), 2) +
           pow(utility::determinant(Z, X, unit_vector), 2));
   X = point(_2[0], _7[0], _6[0]);
   Y = point(_2[1], _7[1], _6[1]);
   Z = point(_2[2], _7[2], _6[2]);
   ds += 0.5 * sqrt(pow(utility::determinant(X, Y, unit_vector), 2) +
           pow(utility::determinant(Y, Z, unit_vector), 2) +
           pow(utility::determinant(Z, X, unit_vector), 2));
   tempFace.Normal = ds * utility::normalize(utility::cross(_2, _3, _7));
   /******************************************************************** */
   // To add -- Face centroid calculation
   tempFace.Centroid = _2 + _7 + _6 + _3;
   tempFace.Centroid *= 0.25;
   /******************************************************************** */
   _Face[counter] = tempFace;
   _mapFace(i, _jmax, k)[2] = (uintptr_t) & _Face[counter];
   _mapFace(i, 1, k)[3] = (uintptr_t) & _Face[counter];
   counter++;
  }
 }
 // Interior cells
 for (int i = 2; i <= _imax - 1; ++i) {
  for (int j = 1; j <= _jmax; ++j) {
   for (int k = 2; k <= _kmax; ++k) {
    // Assign the appropriate nodes
    _0 = *((point *) _cell(i, j, k)[0]);
    _1 = *((point *) _cell(i, j, k)[1]);
    _2 = *((point *) _cell(i, j, k)[2]);
    _3 = *((point *) _cell(i, j, k)[3]);
    _4 = *((point *) _cell(i, j, k)[4]);
    _5 = *((point *) _cell(i, j, k)[5]);
    _6 = *((point *) _cell(i, j, k)[6]);
    _7 = *((point *) _cell(i, j, k)[7]);
    if (_mapFace(i, j, k)[0] == UINTPTR_MAX) { // No i+1 face found
     // Push the face to face heap and then assign it to the left and right cell
     // the face node numbering sequence is 1-2-6-5
     tempFace.Left = findCell((uintptr_t) & _cell(i, j, k));
     tempFace.Right = findCell((uintptr_t) &(_cell(i + 1, j, k)));
     point X(_1[0], _2[0], _6[0]), Y(_1[1], _2[1], _6[1]), Z(_1[2], _2[2], _6[2]);
     ds = 0.5 * sqrt(pow(utility::determinant(X, Y, unit_vector), 2) +
             pow(utility::determinant(Y, Z, unit_vector), 2) +
             pow(utility::determinant(Z, X, unit_vector), 2));
     X = point(_1[0], _6[0], _5[0]);
     Y = point(_1[1], _6[1], _5[1]);
     Z = point(_1[2], _6[2], _5[2]);
     ds += 0.5 * sqrt(pow(utility::determinant(X, Y, unit_vector), 2) +
             pow(utility::determinant(Y, Z, unit_vector), 2) +
             pow(utility::determinant(Z, X, unit_vector), 2));
     tempFace.Normal = ds * utility::normalize(utility::cross(_1, _2, _6)); // Calculate Normal
     /******************************************************************** */
     // To add -- Face centroid calculation
     tempFace.Centroid = _1 + _2 + _6 + _5;
     tempFace.Centroid *= 0.25;
     /******************************************************************** */
     _Face[counter] = tempFace;
     _mapFace(i, j, k)[0] = (uintptr_t) & _Face[counter];
     _mapFace(i + 1, j, k)[1] = (uintptr_t) & _Face[counter];
     counter++;
    }
    if (_mapFace(i, j, k)[1] == UINTPTR_MAX) {
     // Push the face to face heap and then assign it to the left and right cell
     // the face node numbering sequence is 0-3-7-4
     if (i - 1 != 1) {
      tempFace.Left = findCell((uintptr_t) &(_cell(i - 1, j, k)));
     } else {
      tempFace.Left = findCell((uintptr_t) &(_pcell(k)));
     }
     tempFace.Right = findCell((uintptr_t) &(_cell(i, j, k)));
     point X(_0[0], _3[0], _7[0]), Y(_0[1], _3[1], _7[1]), Z(_0[2], _3[2], _7[2]);
     ds = 0.5 * sqrt(pow(utility::determinant(X, Y, unit_vector), 2) +
             pow(utility::determinant(Y, Z, unit_vector), 2) +
             pow(utility::determinant(Z, X, unit_vector), 2));
     X = point(_0[0], _7[0], _4[0]);
     Y = point(_0[1], _7[1], _4[1]);
     Z = point(_0[2], _7[2], _4[2]);
     ds += 0.5 * sqrt(pow(utility::determinant(X, Y, unit_vector), 2) +
             pow(utility::determinant(Y, Z, unit_vector), 2) +
             pow(utility::determinant(Z, X, unit_vector), 2));
     tempFace.Normal = ds * utility::normalize(utility::cross(_0, _3, _7)); // Calculate Normal
     /******************************************************************** */
     // To add -- Face centroid calculation
     tempFace.Centroid = _0 + _3 + _7 + _4;
     tempFace.Centroid *= 0.25;
     /******************************************************************** */
     _Face[counter] = tempFace;
     _mapFace(i, j, k)[1] = (uintptr_t) & _Face[counter];
     _mapFace(i - 1, j, k)[0] = (uintptr_t) & _Face[counter];
     counter++;
    }
    if (_mapFace(i, j, k)[2] == UINTPTR_MAX) {
     // Push the face to face heap and then assign it to the left and right cell
     // the face node numbering sequence is 2-3-7-6
     tempFace.Left = findCell((uintptr_t) &(_cell(i, j, k)));
     tempFace.Right = findCell((uintptr_t) &(_cell(i, j + 1, k)));
     point X(_2[0], _3[0], _7[0]), Y(_2[1], _3[1], _7[1]), Z(_2[2], _3[2], _7[2]);
     ds = 0.5 * sqrt(pow(utility::determinant(X, Y, unit_vector), 2) +
             pow(utility::determinant(Y, Z, unit_vector), 2) +
             pow(utility::determinant(Z, X, unit_vector), 2));
     X = point(_2[0], _7[0], _6[0]);
     Y = point(_2[1], _7[1], _6[1]);
     Z = point(_2[2], _7[2], _6[2]);
     ds += 0.5 * sqrt(pow(utility::determinant(X, Y, unit_vector), 2) +
             pow(utility::determinant(Y, Z, unit_vector), 2) +
             pow(utility::determinant(Z, X, unit_vector), 2));
     tempFace.Normal = ds * utility::normalize(utility::cross(_2, _3, _7));
     /******************************************************************** */
     // To add -- Face centroid calculation
     tempFace.Centroid = _2 + _3 + _7 + _6;
     tempFace.Centroid *= 0.25;
     /******************************************************************** */
     _Face[counter] = tempFace;
     _mapFace(i, j, k)[2] = (uintptr_t) & _Face[counter];
     _mapFace(i, j + 1, k)[3] = (uintptr_t) & _Face[counter];
     counter++;

    }
    if (_mapFace(i, j, k)[3] == UINTPTR_MAX) {
     // Push the face to face heap and then assign it to the left and right cell
     // the face node numbering sequence is 1-5-4-0
     tempFace.Left = findCell((uintptr_t) &(_cell(i, j - 1, k)));
     tempFace.Right = findCell((uintptr_t) &(_cell(i, j, k)));
     point X(_1[0], _5[0], _4[0]), Y(_1[1], _5[1], _4[1]), Z(_1[2], _5[2], _4[2]);
     ds = 0.5 * sqrt(pow(utility::determinant(X, Y, unit_vector), 2) +
             pow(utility::determinant(Y, Z, unit_vector), 2) +
             pow(utility::determinant(Z, X, unit_vector), 2));
     X = point(_1[0], _4[0], _0[0]);
     Y = point(_1[1], _4[1], _0[1]);
     Z = point(_1[2], _4[2], _0[2]);
     ds += 0.5 * sqrt(pow(utility::determinant(X, Y, unit_vector), 2) +
             pow(utility::determinant(Y, Z, unit_vector), 2) +
             pow(utility::determinant(Z, X, unit_vector), 2));
     tempFace.Normal = ds * utility::normalize(utility::cross(_1, _4, _5)); // Assign the appropriate nodes
     /******************************************************************** */
     // To add -- Face centroid calculation
     tempFace.Centroid = _1 + _5 + _4 + _0;
     tempFace.Centroid *= 0.25; 
     /******************************************************************** */
     _Face[counter] = tempFace;
     _mapFace(i, j, k)[3] = (uintptr_t) & _Face[counter];
     _mapFace(i, j - 1, k)[2] = (uintptr_t) & _Face[counter];
     counter++;
    }
    if (_mapFace(i, j, k)[4] == UINTPTR_MAX) {
     // Push the face to face heap and then assign it to the left and right cell
     // the face node numbering sequence is 6-7-4-5
     tempFace.Left = findCell((uintptr_t) &(_cell(i, j, k)));

     // -----------------------------------------------------------------------------
     //                        Added recently (needs testing)
     // -----------------------------------------------------------------------------
     tempFace.Right = findCell((uintptr_t) &(_cell(i, j, k + 1)));
     // -----------------------------------------------------------------------------
     point X(_6[0], _7[0], _4[0]), Y(_6[1], _7[1], _4[1]), Z(_6[2], _7[2], _4[2]);
     ds = 0.5 * sqrt(pow(utility::determinant(X, Y, unit_vector), 2) +
             pow(utility::determinant(Y, Z, unit_vector), 2) +
             pow(utility::determinant(Z, X, unit_vector), 2));
     X = point(_6[0], _4[0], _5[0]);
     Y = point(_6[1], _4[1], _5[1]);
     Z = point(_6[2], _4[2], _5[2]);
     ds += 0.5 * sqrt(pow(utility::determinant(X, Y, unit_vector), 2) +
             pow(utility::determinant(Y, Z, unit_vector), 2) +
             pow(utility::determinant(Z, X, unit_vector), 2));
     tempFace.Normal = ds * utility::normalize(utility::cross(_6, _7, _4));
     /******************************************************************** */
     // To add -- Face centroid calculation
     tempFace.Centroid = _6 + _7 + _4 + _5;
     tempFace.Centroid *= 0.25; 
     /******************************************************************** */
     _Face[counter] = tempFace;
     _mapFace(i, j, k)[4] = (uintptr_t) & _Face[counter];
     _mapFace(i, j, k + 1)[5] = (uintptr_t) & _Face[counter];
     counter++;
    }
    if (_mapFace(i, j, k)[5] == UINTPTR_MAX) {
     // Push the face to face heap and then assign it to the left and right cell
     // the face node numbering sequence is 2-1-0-3
     tempFace.Left = findCell((uintptr_t) &(_cell(i, j, k - 1)));
     tempFace.Right = findCell((uintptr_t) &(_cell(i, j, k)));
     point X(_2[0], _1[0], _0[0]), Y(_2[1], _1[1], _0[1]), Z(_2[2], _1[2], _0[2]);
     ds = 0.5 * sqrt(pow(utility::determinant(X, Y, unit_vector), 2) +
             pow(utility::determinant(Y, Z, unit_vector), 2) +
             pow(utility::determinant(Z, X, unit_vector), 2));
     X = point(_2[0], _0[0], _3[0]);
     Y = point(_2[1], _0[1], _3[1]);
     Z = point(_2[2], _0[2], _3[2]);
     ds += 0.5 * sqrt(pow(utility::determinant(X, Y, unit_vector), 2) +
             pow(utility::determinant(Y, Z, unit_vector), 2) +
             pow(utility::determinant(Z, X, unit_vector), 2));
     tempFace.Normal = ds * utility::normalize(utility::cross(_2, _0, _1));
     /******************************************************************** */
     // To add -- Face centroid calculation
     tempFace.Centroid = _2 + _1 + _0 + _3;
     tempFace.Centroid *= 0.25; 
     /******************************************************************** */
     _Face[counter] = tempFace;
     _mapFace(i, j, k)[5] = (uintptr_t) & _Face[counter];
     _mapFace(i, j, k - 1)[4] = (uintptr_t) & _Face[counter];
     counter++;

    }
   }
  }
 }
 // std::cout << "Face (2,1,2)[0] : " << _mapFace(2, 1, 2)[0] << "Max uptr :" << UINTPTR_MAX << std::endl;
}

void StretchedCylindricalGrid::formPolygonalCells() {
 double ds; // To store area of the face
 Face tempFace; // A temp variable to store the face array
 int counter = _Face.size() + 1; // set counter to the last face in the face hash_map

 // Allocate the nfaces memory pointer array of faces and cells
 for (int k = 1; k <= _kmax + 1; ++k) {
  _pcell(k).nfaces = _jmax;
  _pcell(k).faces = new uintptr_t[_pcell(k).nfaces + 2];
  _pcell(k).cells = new uintptr_t[_pcell(k).nfaces];
  _pcell(k).nodes = new uintptr_t[_pcell(k).nfaces * 2];
 }

 // The following Task to be done :
 // 1) Link nodes to the cell correctly
 // 2) Link the side faces using existing face data
 // 3) Create the top and bottom face calculate area, normal and centroid and push it to face heap and link to cell
 // 4) calculate volume and use the existing volume array to store the volume data

 // Task 1:
 // Node numbering convention is same as that for the hex cells start from bottom and end in top
 for (int k = 1; k <= _kmax + 1; ++k) {
  for (int j = 1; j <= _jmax; ++j) {
   _pcell(k).nodes[j - 1] = (uintptr_t)&(_node(2, j, k));
   _pcell(k).nodes[_jmax + j - 1] = (uintptr_t)&(_node(2, j, k + 1));
  }
 }

 // Task 2:
 // Link using existing face data and put the reference to this cell
 for (int k = 1; k <= _kmax + 1; ++k) {
  for (int j = 1; j <= _jmax; ++j) {
   _pcell(k).faces[j - 1] = _mapFace(2, j, k)[1];
  }
  _pcell(k).faces[_jmax] = UINTPTR_MAX;
  _pcell(k).faces[_jmax + 1] = UINTPTR_MAX;
 }

 // Task 3 and 4:
 // Top and Bottom Face
 // Note: Top face is _jmax and bottom face is _jmax + 1 th index
 for (int k = 2; k <= _kmax; ++k) {
  // Top face not assigned
  if (_pcell(k).faces[_jmax] == UINTPTR_MAX) {
   ds = 0.5 * _jmax * _rdist(1) * _rdist(1) * sin(2.0 * M_PI / double(_jmax));
   tempFace.Normal = point(ds * point(0.0, 0.0, 1.0)); // Calculate Normal
   tempFace.Left = findCell((uintptr_t) &(_pcell(k)));
   tempFace.Right = findCell((uintptr_t) &(_pcell(k + 1)));
   /******************************************************************** */
   // To add -- Face centroid calculation
   tempFace.Centroid = point( 0.0 , 0.0 , _zdist(k) );
   /******************************************************************** */
   _Face[counter] = tempFace;
   // Assign the face to the cell
   _pcell(k).faces[_jmax] = (uintptr_t) &(_Face[counter]);
   _pcell(k + 1).faces[_jmax + 1] = (uintptr_t) &(_Face[counter]);
   _volume(1, 1, k) = ds * ( _zdist(k) - _zdist(k-1));
   /******************************************************************** */
   // To add cell centroid calculation
   _cell_centroid(1,1,k) = point(0.0,0.0, 0.5 * ( _zdist(k) + _zdist(k-1) ) );
   /******************************************************************** */
   counter++;
  }
  // Bottom face not assigned
  if (_pcell(k).faces[_jmax + 1] == UINTPTR_MAX) {
   ds = 0.5 * _jmax * _rdist(1) * _rdist(1) * sin(2 * M_PI / double(_jmax));
   tempFace.Normal = point(ds * point(0.0, 0.0, 1.0)); // Calculate Normal
   tempFace.Right = findCell((uintptr_t) &(_pcell(k)));
   tempFace.Left = findCell((uintptr_t) &(_pcell(k - 1)));
   /******************************************************************** */
   // To add -- Face centroid calculation
   tempFace.Centroid = point(0.0,0.0, 0.5 * ( _zdist(k) + _zdist(k-1) ) );
   /******************************************************************** */
   _Face[counter] = tempFace;
   // Assign the face to the cell
   _pcell(k).faces[_jmax + 1] = (uintptr_t) &(_Face[counter]);
   _pcell(k - 1).faces[_jmax] = (uintptr_t) &(_Face[counter]);
   _volume(1, 1, k) = ds * ( _zdist(k) - _zdist(k-1) );
   /******************************************************************** */
   // To add cell centroid calculation
   _cell_centroid(1, 1, k) = point(0.0,0.0, 0.5 * ( _zdist(k) + _zdist(k-1) ) );
   /******************************************************************** */
   counter++;
  }
 }
 _volume(1, 1, _kmax + 1) = _volume(1, 1, 2);
 _volume(1, 1, 1) = _volume(1, 1, 2);
 // To do -- periodic face centroid mapping
 _cell_centroid(1, 1, _kmax + 1) = point(0.0, 0.0, 0.5 * ( _zdist(_zdist.size() - 1) + _zdist(_zdist.size() - 2) ) );
 _cell_centroid(1, 1, 1) = point(0.0, 0.0, 0.5 * ( _zdist(0) + _zdist(1) ) );

 // Create the hash map of face node # to face memory address
 for (int i = 1; i < _Face.size(); ++i) {
  _hashMapFace[(uintptr_t) &(_Face[i])] = i;
 }
}

/*
\brief This function maps the cells using the memory address to a unique cell number
 */
void StretchedCylindricalGrid::mapCells() {
 uintptr_t addr;
 int counter = 1;
 // Map all Hex Cells
 for (int i = 2; i <= _imax; ++i) {
  for (int j = 1; j <= _jmax; ++j) {
   for (int k = 1; k <= _kmax + 1; ++k) {
    addr = (uintptr_t) &(_cell(i, j, k));
    _mapCells[addr] = counter;
    counter++;
   }
  }
 }
 // Map all Poly Cells
 for (int k = 1; k <= _kmax + 1; ++k) {
  addr = (uintptr_t) &(_pcell(k));
  _mapCells[addr] = counter;
  counter++;
 }
}

void StretchedCylindricalGrid::writeUnstructredGridForSolver(std::string filename) {
 std::ofstream fout(filename.c_str(), std::ios::binary);
 uintptr_t addr;
 int32_t dummy;
 double real_dummy;
 if (fout.fail()) {
  std::cerr << "Error : Unable to open file " << filename << std::endl;
  exit(-1);
 }
 std::cout << "Writing unstructued grid information to solver " << std::endl;
 std::cout << "---------------------------------------------- " << std::endl;
 // ----------------------------------------------------------------------------------------------------------------------------------------------
 // Write Header Information
 // ----------------------------------------------------------------------------------------------------------------------------------------------
 dummy = 7;
 fout.write((char *) & dummy, sizeof (int32_t)); // Magic cookie
 dummy = 0;
 fout.write((char *) & dummy, sizeof (int32_t)); // Major version number
 dummy = 1;
 fout.write((char *) & dummy, sizeof (int32_t)); // Minor version number
 dummy = 0;
 fout.write((char *) & dummy, sizeof (int32_t)); // Revision number
 dummy = getTotalNodes();
 fout.write((char *) & dummy, sizeof (int32_t)); // Total Number of nodes
 std::cout << "Total Nodes = " << dummy << std::endl;
 dummy = _mapCells.size();
 fout.write((char *) & dummy, sizeof (int32_t)); // Total Number of cells
 std::cout << "Total Elements = " << dummy << std::endl;
 dummy = _Face.size();
 fout.write((char *) & dummy, sizeof (int32_t)); // Total Number of faces
 std::cout << "Total Faces = " << dummy << std::endl;
 dummy = getTotalHexCells();
 fout.write((char *) & dummy, sizeof (int32_t)); // Total Number of hex cells
 std::cout << "Total Hex Cells = " << dummy << std::endl;
 dummy = _jmax + 2;
 fout.write((char *) & dummy, sizeof (int32_t)); // Faces contained by one poly cell
 std::cout << "One Poly Cell Faces = " << dummy << std::endl;
 // ----------------------------------------------------------------------------------------------------------------------------------------------
 // Face data
 // ----------------------------------------------------------------------------------------------------------------------------------------------
 std::cout << "First Face Data = [" << (_Face.begin())->second.Normal[0] << " " << (_Face.begin())->second.Normal[1] << " " << (_Face.begin())->second.Normal[2] << " " << (_Face.begin())->second.Left << " " << (_Face.begin())->second.Right << (_Face.begin())->second.Centroid << "]" << std::endl;
 for (std::hash_map<int, Face>::iterator _iFace = _Face.begin(); _iFace != _Face.end(); ++_iFace) {
  real_dummy = _iFace->second.Normal[0];
  fout.write((char *) & real_dummy, sizeof (double));
  real_dummy = _iFace->second.Normal[1];
  fout.write((char *) & real_dummy, sizeof (double));
  real_dummy = _iFace->second.Normal[2];
  fout.write((char *) & real_dummy, sizeof (double));
  dummy = _iFace->second.Left;
  fout.write((char *) & dummy, sizeof (int32_t));
  dummy = _iFace->second.Right;
  fout.write((char *) & dummy, sizeof (int32_t));
  /******************************************************************** */
  // To add -- write face centroid
  /******************************************************************** */
  real_dummy = _iFace->second.Centroid[0];
  fout.write((char *) & real_dummy, sizeof (double));
  real_dummy = _iFace->second.Centroid[1];
  fout.write((char *) & real_dummy, sizeof (double));
  real_dummy = _iFace->second.Centroid[2];
  fout.write((char *) & real_dummy, sizeof (double));
  /******************************************************************** */
 }
 std::cout << "Last Face Data = [" << _Face[_Face.size()].Normal[0] << " " << _Face[_Face.size()].Normal[1] << " " << _Face[_Face.size()].Normal[2] << " " << _Face[_Face.size()].Left << " " << _Face[_Face.size()].Right << _Face[_Face.size()].Centroid << "]" << std::endl;
 // ----------------------------------------------------------------------------------------------------------------------------------------------
 // Cell Volume and centroid data
 // ----------------------------------------------------------------------------------------------------------------------------------------------
 std::cout << "First Cell volume = " << _volume(2, 1, 1) << std::endl;
 for (int i = 2; i <= _imax; ++i) {
  for (int j = 1; j <= _jmax; ++j) {
   for (int k = 1; k <= _kmax + 1; ++k) {
    real_dummy = _volume(i, j, k);
    fout.write((char *) & real_dummy, sizeof (double));
   }
  }
 }
 for (int k = 1; k <= _kmax + 1; ++k) {
  real_dummy = _volume(1, 1, k);
  fout.write((char *) & real_dummy, sizeof (double));
 }
 std::cout << "Last Cell volume = " << _volume(1, 1, _kmax + 1) << std::endl;

 /******************************************************************** */
 // To add -- write cell centroid
 /******************************************************************** */
 std::cout << "First Cell centroid = " << _cell_centroid(2, 1, 1) << std::endl;
 for (int i = 2; i <= _imax; ++i) {
  for (int j = 1; j <= _jmax; ++j) {
   for (int k = 1; k <= _kmax + 1; ++k) {
    // X coordinate
    real_dummy = _cell_centroid(i, j, k)[0];
    fout.write((char *) & real_dummy, sizeof (double));
    // Y coordiante
    real_dummy = _cell_centroid(i, j, k)[1];
    fout.write((char *) & real_dummy, sizeof (double));
    // Z coordinate
    real_dummy = _cell_centroid(i, j, k)[2];
    fout.write((char *) & real_dummy, sizeof (double));
   }
  }
 }
 for (int k = 1; k <= _kmax + 1; ++k) {
   // X coordinate
   real_dummy = _cell_centroid(1, 1, k)[0];
   fout.write((char *) & real_dummy, sizeof (double));
   // Y coordiante
   real_dummy = _cell_centroid(1, 1, k)[1];
   fout.write((char *) & real_dummy, sizeof (double));
   // Z coordinate
   real_dummy = _cell_centroid(1, 1, k)[2];
   fout.write((char *) & real_dummy, sizeof (double));
 }
 std::cout << "Last Cell centroid = " << _cell_centroid(1, 1, _kmax + 1) << std::endl;
 /******************************************************************** */


 // ----------------------------------------------------------------------------------------------------------------------------------------------
 // write boundary cells and ghost cells list
 // ----------------------------------------------------------------------------------------------------------------------------------------------
 // ghost cell (sides)
 dummy = _jmax * (_kmax + 1) + 2 + 2 * _jmax * (_imax - 2); // Total Ghost cells
 fout.write((char *) & dummy, sizeof (int32_t)); // neighbour cell 4
 std::cout << "Total Ghost Cells in File = " << dummy << std::endl;
 for (int k = 1; k <= _kmax + 1; ++k) {
  for (int j = 1; j <= _jmax; ++j) {
   addr = (uintptr_t) &(_cell(_imax, j, k));
   dummy = findCell(addr);
   fout.write((char *) & dummy, sizeof (int32_t)); // neighbour cell 4
   //std::cout << dummy << std::endl;
  }
 }
 std::cout << "First Ghost Cell in File = " << findCell((uintptr_t) &(_cell(_imax, 1, 1))) << std::endl;
 // ghost cell (top and bottom)
 for (int i = 2; i < _imax; ++i) {
  for (int j = 1; j <= _jmax; ++j) {
   addr = (uintptr_t) &(_cell(i, j, _kmax + 1));
   dummy = findCell(addr);
   fout.write((char *) & dummy, sizeof (int32_t)); // neighbour cell 4
   //std::cout << dummy << std::endl;
   addr = (uintptr_t) &(_cell(i, j, 1));
   dummy = findCell(addr);
   fout.write((char *) & dummy, sizeof (int32_t)); // neighbour cell 4
   //std::cout << dummy << std::endl;
  }
 }
 // The two ghost polygonal cells
 addr = (uintptr_t) &(_pcell(_kmax + 1));
 dummy = findCell(addr);
 fout.write((char *) & dummy, sizeof (int32_t)); // neighbour cell 4
 //std::cout << dummy << std::endl;
 addr = (uintptr_t) &(_pcell(1));
 dummy = findCell(addr);
 fout.write((char *) & dummy, sizeof (int32_t)); // neighbour cell 4
 //std::cout << dummy << std::endl;
 fout.close();
 std::cout << "Last Ghost Cell in File = " << dummy << std::endl;
}

void StretchedCylindricalGrid::writeNodesAndCells(std::string filename) {
 std::ofstream fout(filename.c_str(), std::ios::binary);
 int32_t dummy;
 double real_dummy;
 if (fout.fail()) {
  std::cerr << "Error : Unable to open file " << filename << std::endl;
  exit(-1);
 }
 // ----------------------------------------------------------------------------------------------------------------------------------------------
 // Write Header Information
 // ----------------------------------------------------------------------------------------------------------------------------------------------
 dummy = 7;
 fout.write((char *) & dummy, sizeof (int32_t)); // Magic cookie
 dummy = 0;
 fout.write((char *) & dummy, sizeof (int32_t)); // Major version number
 dummy = 1;
 fout.write((char *) & dummy, sizeof (int32_t)); // Minor version number
 dummy = 0;
 fout.write((char *) & dummy, sizeof (int32_t)); // Revision number
 dummy = getTotalNodes();
 fout.write((char *) & dummy, sizeof (int32_t)); // Total Number of nodes
 dummy = _mapCells.size();
 fout.write((char *) & dummy, sizeof (int32_t)); // Total Number of cells
 // ----------------------------------------------------------------------------------------------------------------------------------------------
 // Node data
 // ----------------------------------------------------------------------------------------------------------------------------------------------
 for (int i = 1; i <= _imax + 1; ++i) {
  for (int j = 1; j <= _jmax; ++j) {
   for (int k = 1; k <= _kmax + 2; ++k) {
    real_dummy = _node(i, j, k)[0];
    fout.write((char *) & real_dummy, sizeof (real_dummy));
    real_dummy = _node(i, j, k)[1];
    fout.write((char *) & real_dummy, sizeof (real_dummy));
    real_dummy = _node(i, j, k)[2] + _zoffset;
    fout.write((char *) & real_dummy, sizeof (real_dummy));
   }
  }
 }
 // ----------------------------------------------------------------------------------------------------------------------------------------------
 // Cell data
 // ----------------------------------------------------------------------------------------------------------------------------------------------
 dummy = getTotalHexCells();
 fout.write((char *) & dummy, sizeof (int32_t)); // Total Hex Cells
 dummy = getTotalPolyCells();
 fout.write((char *) & dummy, sizeof (int32_t)); // Total Poly Cells
 dummy = _jmax * 2;
 fout.write((char *) & dummy, sizeof (int32_t)); // # of Nodes connecting to form Poly Cell
 dummy = _jmax + 2;
 fout.write((char *) & dummy, sizeof (int32_t)); // # of neighbour cells connecting to Poly Cell
 // ----------------------------------------------------------------------------------------------------------------------------------------------
 // Cell connectivity
 // ----------------------------------------------------------------------------------------------------------------------------------------------
 for (int i = 2; i <= _imax; ++i) {
  for (int j = 1; j <= _jmax; ++j) {
   for (int k = 1; k <= _kmax + 1; ++k) {
    // Node connectivity
    for (int l = 0; l < 8; ++l) {
     dummy = findNode((uintptr_t) _cell(i, j, k)[l]);
     fout.write((char *) & dummy, sizeof (int32_t)); // # of neighbour cells connecting to Poly Cell
    }
   }
  }
 }
 // ----------------------------------------------------------------------------------------------------------------------------------------------
}

void StretchedCylindricalGrid::outputHigherOrderNodes(std::string filename){
 bool debug(false);
 std::ofstream fout(filename.c_str());
 if (fout.fail()) {
  std::cerr << "Error : Unable to open file " << filename << std::endl;
  exit(-1);
 }
 // Map all Hex Cells
 Hexa8p3 hex;
 tvmet::Matrix<double,8,3> hex_nodes;
 int counter = 0;
 for (int i = 2; i <= _imax; ++i) {
  for (int j = 1; j <= _jmax; ++j) {
   for (int k = 1; k <= _kmax + 1; ++k) {
    for (int m = 0; m < 3; ++m){ 
     for (int n = 0; n < 8; ++n){
      hex_nodes(n,m) = (*((point *) _cell(i, j, k)[n]))(m) + hex.spMat.unit3by3(2,m) * _zoffset;
     }
    }
    // Form the Higher order element and print it to file
    hex.setHexNodes(hex_nodes);
    // output the calculated 27 Gaussian nodes, the determinant of Jacobian, 
    //fout << counter + 1 << " " << hex.asciiOutput();
    fout << hex.asciiOutput();
    counter++;	
   }
  }
 }
 std::cout << "A Total of " << counter << " cells quadrature node written to file " << filename << std::endl;
 fout.close();
 if(debug){
  std::ofstream fdebug_out("debugData");
  for (int i = 2; i <= _imax; ++i) {
   for (int j = 1; j <= _jmax; ++j) {
   int k = 10;
    for (int n = 0; n < 8; ++n){
     for (int m = 0; m < 3; ++m){ 
      fdebug_out << (*((point *) _cell(i, j, k)[n]))(m) + hex.spMat.unit3by3(2,m) * _zoffset << " ";	
     }
     fdebug_out << "\n";
    }
   }
  }
 }
}

void StretchedCylindricalGrid::writeDebugCentroid(std::string filename){
 std::ofstream fout(filename.c_str());
 if (fout.fail()) {
  std::cerr << "Error : Unable to open file " << filename << std::endl;
  exit(-1);
 }
 fout << "VARIABLES=X,Y,Z\nZONE\n";
 for (int i = 2; i <= _imax; ++i) {
  for (int j = 1; j <= _jmax; ++j) {
   for (int k = 1; k <= _kmax + 1; ++k) {
    fout << _cell_centroid(i,j,k)[0] << "\t" << _cell_centroid(i,j,k)[1] << "\t" << _cell_centroid(i,j,k)[2] << "\n";
   }
  }
 }
 fout << "ZONE\n";
 for (int k = 1; k <= _kmax + 1 ; ++k)
  fout << _cell_centroid(1,1,k)[0] << "\t" << _cell_centroid(1,1,k)[1] << "\t" << _cell_centroid(1,1,k)[2] << "\n";
 fout << "ZONE\n";
 for (std::hash_map<int, Face>::iterator _iFace = _Face.begin(); _iFace != _Face.end(); ++_iFace) {
  /******************************************************************** */
  // To add -- write face centroid
  /******************************************************************** */
  fout << _iFace->second.Centroid[0] << "\t" << _iFace->second.Centroid[1] << "\t" << _iFace->second.Centroid[2] << "\n";
  /******************************************************************** */
 }
}

StretchedCylindricalGrid::~StretchedCylindricalGrid() {
 ;
}
