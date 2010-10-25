#include <algorithm>
#include "StretchedStaticCylindricalGridSolver.h"
#include "Hexa8p3.h"


StretchedStaticCylindricalGridSolver::StretchedStaticCylindricalGridSolver(std::string filename)
{
  Constructor(filename);
  _outflowBC = false;
}

void StretchedStaticCylindricalGridSolver::Constructor(std::string filename)
{
// Constructor Stat's (for debug purpose can be removed if too trouble some)
//std::cout <<  "\nCreating Unstructured Cylindrical grid solver \n";
//std::cout <<  "-------------------------------------------------\n";

  std::ifstream fin(filename.c_str(), std::ios::binary);
  if (fin.fail()) {
    std::cout << "Error : Unable to open file " << filename << "\n";
    exit(-1);
  }
  int32_t dummy;
  double real_dummy;
  int nNodes;
// ---------------------------------------------------------------------------------------------------------------------------------------
// Read Header Information and verify identity
// ---------------------------------------------------------------------------------------------------------------------------------------
  fin.read((char *) & dummy, sizeof (int32_t)); // Magic cookie
  if (dummy != 7) {
    std::cout << "Magic cookie does not match in binary file :" << filename << "\n";
    exit(0);
  }
//std::cout <<  "Binary File Version : ";
  fin.read((char *) & dummy, sizeof (int32_t)); // Major version number
  if (dummy != 0) {
    std::cout << "Major Verison does not match in binary file :" << filename << "\n";
    std::cout << "Current Major Verison = " << 0 << " Major Version in file = " << dummy << "\n";
    exit(0);
  }
//std::cout <<  dummy << ".";
  fin.read((char *) & dummy, sizeof (int32_t)); // Minor version number
  if (dummy != 1) {
    std::cout << "warning : Minor Verison does not match that in binary file :" << filename << "\n";
    std::cout << "Current Major Verison = " << 0 << " Major Version in file = " << dummy << "\n";
  }
//std::cout <<  dummy << ".";
  fin.read((char *) & dummy, sizeof (int32_t)); // Revision version number
  if (dummy != 0) {
    std::cout << "warning : Binary file created at revision number = " << dummy << "the file might not be compatible " << "\n";
  }
//std::cout <<  dummy << "\n";
  fin.read((char *) & dummy, sizeof (int32_t)); // Total Nodes
//std::cout <<  "Total Nodes         = " << dummy << "\n";
  nNodes = dummy;
  _totalNodes = dummy;

  fin.read((char *) & dummy, sizeof (int32_t)); // Total Cells
//std::cout <<  "Total Elements      = " << dummy << "\n";
  _totalCells = dummy;

  fin.read((char *) & dummy, sizeof (int32_t)); // Total Faces
//std::cout <<  "Total Faces         = " << dummy << "\n";
  _totalFaces = dummy;

  fin.read((char *) & dummy, sizeof (int32_t)); // Total Hex Cells
//std::cout <<  "Total Hexa elements  = " << dummy << "\n";
  _totalHexCells = dummy;

  fin.read((char *) & dummy, sizeof (int32_t)); // one Poly Cell faces
//std::cout <<  "Faces conneciting Poly Cell = " << dummy << "\n";
  _onePolyCellFaces = dummy;

// Allocate data based on Info read so far
  _U.resize(Range(1, _totalCells));
  _Uavg.resize(Range(1, _totalCells));
  _URK.resize(Range(1, _totalCells));
  _rhsRK.resize(Range(1, _totalCells));
  _rhs.resize(Range(1, _totalCells));
  _volume.resize(Range(1, _totalCells));
  _cell_centroid.resize(Range(1, _totalCells));
  _faceDataArray.resize(Range(1, _totalFaces));
  _localDt.resize(Range(1, _totalCells));
  _charLength.resize(Range(1, _totalCells));
  _source.resize(Range(1, _totalCells));
// Fill Arrays with appropriate values
  _U = flux(0.0);
  _Uavg = flux(0.0);
  _URK = flux(0.0);
  _rhs = flux(0.0);
  _rhsRK = flux(0.0);
  _volume = double(0.0);

// ---------------------------------------------------------------------------------------------------------------------------------------
//  Face Data ( very very very very very very very very very very very very critical part of the code)
// ---------------------------------------------------------------------------------------------------------------------------------------
  for (int counter = 1; counter <= _totalFaces; ++counter) {

    /*********************
      Read Normals
    *********************/
    fin.read((char *) & real_dummy, sizeof (double));
    //std::cout <<  real_dummy << " ";
    _faceDataArray(counter).normal[0] = real_dummy;

    fin.read((char *) & real_dummy, sizeof (double));
    //std::cout <<  real_dummy << " ";
    _faceDataArray(counter).normal[1] = real_dummy;

    fin.read((char *) & real_dummy, sizeof (double));
    //std::cout <<  real_dummy << " ";
    _faceDataArray(counter).normal[2] = real_dummy;

    fin.read((char *) & dummy, sizeof (int32_t));
    //std::cout <<  dummy << " ";
    _faceDataArray(counter).leftCell = dummy;
    if (dummy <= 0) {
      std::cout <<  "Zero/-ve Index identified in input at Face # " << counter << " Left = " << dummy << "Right = " << _faceDataArray(counter).rightCell << " Normal = " << _faceDataArray(counter).normal << "\n";
      exit(0);
    }
    fin.read((char *) & dummy, sizeof (int32_t));
    _faceDataArray(counter).rightCell = dummy;
    if (dummy <= 0) {
      std::cout <<  "Zero/-ve Index identified in input at Face # " << counter << " Left = " << dummy << "Right = " << _faceDataArray(counter).rightCell << " Normal = " << _faceDataArray(counter).normal << "\n";
      exit(0);
    }

    /***********************
      To add - Read Centroids
     ***********************/
    fin.read((char *) & real_dummy, sizeof (double));
    //std::cout <<  real_dummy << " ";
    _faceDataArray(counter).centroid[0] = real_dummy;

    fin.read((char *) & real_dummy, sizeof (double));
    //std::cout <<  real_dummy << " ";
    _faceDataArray(counter).centroid[1] = real_dummy;

    fin.read((char *) & real_dummy, sizeof (double));
    //std::cout <<  real_dummy << " ";
    _faceDataArray(counter).centroid[2] = real_dummy;

  }
//std::cout <<  "First Face Data = [" << _faceDataArray(1).normal[0] << " " << _faceDataArray(1).normal[1] << " " << _faceDataArray(1).normal[2] << " " << _faceDataArray(1).leftCell << " " << _faceDataArray(1).rightCell << "]" << "\n";
//std::cout <<  "Last Face Data = [" << _faceDataArray(_totalFaces).normal[0] << " " << _faceDataArray(_totalFaces).normal[1] << " " << _faceDataArray(_totalFaces).normal[2] << " " << _faceDataArray(_totalFaces).leftCell << " " << _faceDataArray(_totalFaces).rightCell << "]" << "\n";
// ---------------------------------------------------------------------------------------------------------------------------------------
//  Volume data
// ---------------------------------------------------------------------------------------------------------------------------------------
  for (int i = 1; i <= _totalCells; ++i) {
    fin.read((char *) & real_dummy, sizeof (double));
    _volume(i) = real_dummy;
  }
//std::cout <<  "First cell Volume = " << _volume(1) << "\n";
//std::cout <<  "Last cell Volume = " << _volume(_totalCells) << "\n";
// ---------------------------------------------------------------------------------------------------------------------------------------
//  Centroid data
// ---------------------------------------------------------------------------------------------------------------------------------------
  for (int i = 1; i <= _totalCells; ++i) {
    fin.read((char *) & real_dummy, sizeof (double));
    _cell_centroid(i)[0] = real_dummy;
    fin.read((char *) & real_dummy, sizeof (double));
    _cell_centroid(i)[1] = real_dummy;
    fin.read((char *) & real_dummy, sizeof (double));
    _cell_centroid(i)[2] = real_dummy;
  }
//std::cout <<  "First cell centroid = " << _cell_centroid(1) << "\n";
//std::cout <<  "Last cell centroid = " << _cell_centroid(_totalCells) << "\n";

// Read in ghost cells
  fin.read((char *) & dummy, sizeof (int32_t));
  int nGhost = dummy;
//std::cout <<  "Total Ghost Cells in File = " << dummy << "\n";
  for (int i = 0; i < nGhost; ++i) {
    fin.read((char *) & dummy, sizeof (int32_t));
    _ghostCells.push_back(dummy);
    //std::cout <<  dummy << "\n";
  }
//std::cout <<  "First Ghost Cell in File  = " << *(_ghostCells.begin()) << "\n";
//std::cout <<  "Last  Ghost Cell in File  = " << (_ghostCells.back()) << "\n";
//std::cout <<  "---------------- END ---------------------" << "\n";

// ---------------------------------------------------------------------------------------------------------------------------------------
// Check if the grid has been read in correctly
// ---------------------------------------------------------------------------------------------------------------------------------------
//gridCheck();
// ---------------------------------------------------------------------------------------------------------------------------------------
// Mark all the faces that belong to or contain left/right neighbour as a Ghost cell
// Mark all interior cell faces and keep them
// ---------------------------------------------------------------------------------------------------------------------------------------
// Optimize to speed up the vector addition (pre allocate using a guess value )
  _ghostFaces.reserve(_ghostCells.size() * 4);
  _interiorFaces.reserve(_totalFaces);
//	std::cout <<  "Grouping Faces" << "\n";
// groupFaces();
  _CFL = 0.0;
  _currentTime = 0.0;
  init_flag = false;
  _free = flux(0.0);
  _ntime_samples = 1;

}

void StretchedStaticCylindricalGridSolver::calculateLocalDt()
{
  for (std::vector<faceData *>::iterator it = _interiorFaces.begin(); it != _interiorFaces.end(); ++it) {
    double dA = utility::mag((*it)->normal);
    _charLength((*it)->leftCell) += dA;
    _charLength((*it)->rightCell) += dA;
  }
  for (std::vector<faceMapGhost>::iterator it = _ghostFaces.begin(); it != _ghostFaces.end(); ++it) {
    double dA = utility::mag(it->fd_ptr->normal);
    _charLength(it->fd_ptr->leftCell) += dA;
    _charLength(it->fd_ptr->leftCell) += dA;
  }
  for (int i = 1; i <= _totalCells; ++i) {
    _charLength(i) /= _volume(i);
    _localDt(i) = _CFL * _charLength(i) / airStandardAtmosphere::speedOfSound;
    //std::cout <<  "Cell # " << i << " Local Dt = " << _localDt(i) << "\n";
    if (_localDt(i) < 0.0)
      std::cout <<  "Negative Time Step Detected dt = " << _localDt(i) << "\n";
  }
// Make the characteristic length of all ghost cells to 0.0
  for (std::vector<int>::iterator iter = _ghostCells.begin(); iter != _ghostCells.end(); ++iter)
    _localDt((*iter)) = 0.0;
}

double StretchedStaticCylindricalGridSolver::calculateDtCFL(double CFL, double vel)
{
  for (int i = 1; i <= _totalFaces; ++i) {
    double dA = utility::mag(this->_faceDataArray(i).normal);
    _charLength(_faceDataArray(i).leftCell) += dA;
    _charLength(_faceDataArray(i).rightCell) += dA;
    //std::cout <<  "Area = " << dA << "\n";
  }
  double smallestDt = 1e20;
  double tempDt;
  for (int i = 1; i <= _totalCells; ++i) {
    _charLength(i) = _volume(i) / _charLength(i);
  }
// Make all ghost RHS as zero
  for (std::vector<int>::iterator it = _ghostCells.begin(); it != _ghostCells.end(); ++it) {
    _charLength(*it) = 1e50;
  }

  for (int i = 1; i <= _totalCells; ++i) {
    tempDt = CFL * _charLength(i) / (airStandardAtmosphere::speedOfSound + vel);
    if (smallestDt > tempDt) {
      smallestDt = tempDt;
    }
    if (smallestDt <= 0.0)
      std::cout <<  "Negative Time Step Detected dt = " << smallestDt
                << " i = " << i << " charlength = " << _charLength(i) << "\n";
  }

// Make the characteristic length of all ghost cells to 0.0
  _localDt = smallestDt;
  return smallestDt;
}

void StretchedStaticCylindricalGridSolver::setCFL(double CFL)
{
  _CFL = CFL;
  calculateLocalDt();
}

void StretchedStaticCylindricalGridSolver::setDt(double dt)
{
  std::cout <<  "Local Dt = " << dt << "\n";
  _CFL = 0.7;
  _localDt = dt;
}

void StretchedStaticCylindricalGridSolver::setCurrentTime(double t)
{
  std::cout <<  "Current Time set to = " << t << "\n";
  _currentTime = t;
}

void StretchedStaticCylindricalGridSolver::Init(flux &f)
{
  _U = f;
  _URK = f;
  init_flag = true;
}

void StretchedStaticCylindricalGridSolver::Init(double rho,double u,double v,double w,double p)
{
  flux f = flux(rho,u,v,w,p);
  std::cout<<"Setting free stream values "<<f<<"\n";
  _U = f;
  _URK = f;
  init_flag = true;
}


void StretchedStaticCylindricalGridSolver::writeHexCellsSolution(fluxArray1 &RHS, std::string filename)
{
  std::ofstream fout(filename.c_str());
  if(fout.fail()) {
    std::cout<<"Error Opening file "<<filename<<"\n";
    exit(1);
  }
  fout << "variables = rho , u , v , w , p\n";
  fout << "FILETYPE = SOLUTION\nZONE NODES = " << _totalNodes << " , ELEMENTS = " << _totalHexCells << ", DATAPACKING = BLOCK , ZONETYPE = FEBRICK\n";
  fout << "SOLUTIONTIME = " << getTime() << "\n";
  fout << "VARLOCATION = ([1-5] = CELLCENTERED)\n";
  for (int i = 1; i <= _totalHexCells; ++i)
    fout << RHS(i)[IDX::rho] << "\n";
  for (int i = 1; i <= _totalHexCells; ++i)
    fout << RHS(i)[IDX::u] << "\n";
  for (int i = 1; i <= _totalHexCells; ++i)
    fout << RHS(i)[IDX::v] << "\n";
  for (int i = 1; i <= _totalHexCells; ++i)
    fout << RHS(i)[IDX::w] << "\n";
  for (int i = 1; i <= _totalHexCells; ++i)
    fout << RHS(i)[IDX::p] << "\n";
// Write Centerline solution as separate zone
  fout << "ZONE I=" << _totalCells - _totalHexCells << "\n";
  fout << "SOLUTIONTIME = " << getTime() << "\n";
  for (int i = _totalHexCells + 1; i <= _totalCells; ++i)
    fout << RHS(i)[IDX::rho] << " " << RHS(i)[IDX::u] << " "<< RHS(i)[IDX::v] << " "<< RHS(i)[IDX::w] << " "<< RHS(i)[IDX::p] << "\n";
  fout.close();
}

void StretchedStaticCylindricalGridSolver::writeSolution(std::string path)
{
  std::stringstream s;
  s.precision(2);
  s << std::scientific << path << "/solution" << _currentTime <<".dat" ;
  path = s.str();
  std::string plus("p");
  std::string minus("m");
  std::string point("pt");

  size_t j = path.find('+');
  if(j != std::string::npos)
    path.replace(j,1,plus);
  j = path.find('-');
  if(j != std::string::npos)
    path.replace(j,1,minus);
  j = path.find('.');
  if(j != std::string::npos)
    path.replace(j,1,point);
  writeHexCellsSolution(_U, path);
}

// Needs testing
// \brief roe is the Roe's approximate Riemann solver
// \param UL     --> Left Cell neighbour cell flux
// \param UR     --> Right Cell neighbour cell flux
// \param normal --> unit normal vector
// \param F      --> Flux to be evaluated at the interface

inline void StretchedStaticCylindricalGridSolver::roe(flux &qL, flux &qR, point &normal, flux & F)
{
  double h0L, h0R, Rfac, rhoRoe, aRoe, h0Roe, drho, dp, VnL, VnR, VnRoe, qSqRoe, dV, drhoTerm;
  flux fL, fR, df1, df2, df3;
  point dvel, velRoe;
  short i;
// ---------------------------------------------------------
//                Left and Right Flow values
// ---------------------------------------------------------
  h0L = helperSolver::enthalpy(qL) + 0.5 * helperSolver::velocityDot(qL, qL);
  h0R = helperSolver::enthalpy(qR) + 0.5 * helperSolver::velocityDot(qR, qR);
  VnL = helperSolver::normalVelocity(qL, normal);
  VnR = helperSolver::normalVelocity(qR, normal);
// ---------------------------------------------------------
//                Roe Averaged States
// ---------------------------------------------------------
  Rfac = sqrt(qR[IDX::rho] / qL[IDX::rho]);
  rhoRoe = Rfac * qL[IDX::rho];
  for (i = 0; i < 3; i++)
    velRoe[i] = (qL[i + 1] + Rfac * qR[i + 1]) / (1.0 + Rfac);
  h0Roe = (h0L + Rfac * h0R) / (1.0 + Rfac);
  qSqRoe = tvmet::dot(velRoe, velRoe);
  aRoe = sqrt(airStandardAtmosphere::gm1 * (h0Roe - 0.5 * qSqRoe));
  VnRoe = (VnL + Rfac * VnR) / (1.0 + Rfac);
// ---------------------------------------------------------
//         Differnetial values for characteristics
// ---------------------------------------------------------
  drho = qR[IDX::rho] - qL[IDX::rho];
  for (i = 0; i < 3; ++i)
    dvel[i] = qR[i + 1] - qL[i + 1];
  dp = qR[IDX::p] - qL[IDX::p];
  dV = tvmet::dot(dvel, normal);
  fL = helperSolver::normalFlux(qL, normal);
  fR = helperSolver::normalFlux(qR, normal);
// ---------------------------------------------------------
//         Form Inteface roe averaged fluxes df1 - df3
// ---------------------------------------------------------
// -----------
//     df1
// -----------
  drhoTerm = drho - dp / pow(aRoe, 2);
  df1[0] = fabs(VnRoe) * drhoTerm;
  for (i = 0; i < 3; i++) df1[i + 1] = fabs(VnRoe)*(velRoe[i] * drhoTerm + rhoRoe * (dvel[i] - normal[i] * dV));
  df1[4] = fabs(VnRoe)*(0.5 * qSqRoe * drhoTerm + rhoRoe * (tvmet::dot(velRoe, dvel) - VnRoe * dV));
// -----------
//     df2
// -----------
  df2[0] = 0.5 * fabs(VnRoe + aRoe)*(dp / pow(aRoe, 2) + rhoRoe * dV / aRoe);
  for (i = 0; i < 3; i++) df2[i + 1] = df2[0]*(velRoe[i] + normal[i] * aRoe);
  df2[4] = df2[0]*(h0Roe + VnRoe * aRoe);
// -----------
//     df3
// -----------
  df3[0] = 0.5 * fabs(VnRoe - aRoe)*(dp / pow(aRoe, 2) - rhoRoe * dV / aRoe);
  for (i = 0; i < 3; i++) df3[i + 1] = df3[0] * (velRoe[i] - normal[i] * aRoe);
  df3[4] = df3[0] * (h0Roe - VnRoe * aRoe);
// ----------------------------------------------------------------------------------------
//         Sum fluxes 0.5 * (fL + fR - df1 - df2 -df3 ) to yield flux at the interface
// ----------------------------------------------------------------------------------------
  F = 0.5 * (fL + fR - df1 - df2 - df3);
//F = 0.5 * (fL + fR);
}

void StretchedStaticCylindricalGridSolver::gridCheck()
{

  blitz::Array<int32_t, 1 > tempCell(Range(_start, _start + _totalCells - 1));
  blitz::Array<double, 1 > tempCellDouble(Range(_start, _start + _totalCells - 1));
  tempCell = 0;
  tempCellDouble = double(0.0);

  int32_t right, left;
  bool okFlag = true;
  std::ofstream fout("./data/debugFaceDouble.dat");
  fout << "variables=\"scalar\"\n";
  fout << "FILETYPE=SOLUTION\nZONE NODES=" << _totalNodes << ", ELEMENTS=" << _totalHexCells << ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK\nVARLOCATION = ([1] = CELLCENTERED)\n";
  std::cout << "\n\nChecking for Face consitency and possible errors" << "\n";
  std::cout << "-------------------------------------------------" << "\n";
  std::cout << "1) Checking for Double Counting of Face : ";
// ---------------------------------------------------------------------------------------------------------------------------------------
// Check if some face is double counted
// ---------------------------------------------------------------------------------------------------------------------------------------
  for (int i = 1; i <= _totalFaces; ++i) {
    right = _faceDataArray(i).rightCell;
    left = _faceDataArray(i).leftCell;

    if (right != -1) {
      tempCell(right) += 1;
    }
    if (left != -1) {
      tempCell(left) += 1;
    }
  }
  for (int i = _start; i < _start + _totalHexCells; ++i) {
    std::vector<int>::iterator result = std::find(_ghostCells.begin(), _ghostCells.end(), i);
    if (result != _ghostCells.end()) tempCell(i) = 6;
    fout << tempCell(i) << "\n";
    if (tempCell(i) != 6) {
      std::cout << "warning : Hexa Cell " << i << " value :" << tempCell(i) << "\n";
      okFlag == false;
    }
  }
  for (int i = _start + _totalHexCells; i < _start + _totalCells; ++i) {
    std::vector<int>::iterator result = std::find(_ghostCells.begin(), _ghostCells.end(), i);
    if (result != _ghostCells.end()) tempCell(i) = _onePolyCellFaces;
    fout << tempCell(i) << "\n";
    if (tempCell(i) != _onePolyCellFaces) {
      std::cout << "warning : PolyHedral Cell " << i << " value :" << tempCell(i) << "\n";
      okFlag == false;
    }
  }
  if (okFlag == false) {
    std::cout << "<failed>" << "\n";
    exit(0);
  } else
    std::cout << "<passed>" << "\n";
  tempCell.free();
  fout.close();
// ---------------------------------------------------------------------------------------------------------------------------------------
// Check for zero or negative volume
// ---------------------------------------------------------------------------------------------------------------------------------------
  okFlag = true;
  std::cout << "2) Checking for zero or negative volume : ";
  for (int i = _start; i < _start + _totalCells; ++i) {
    if (_volume(i) < 1e-10) {
      std::cout << " Volume is close to zero or negative in cell # : " << i << " volume = " << _volume(i) << "\n";
      okFlag = false;
    }
  }
  if (okFlag == false) {
    std::cout << "<failed>" << "\n";
    exit(0);
  } else
    std::cout << "<passed>" << "\n";
// ---------------------------------------------------------------------------------------------------------------------------------------
// Check if face area vectors added up to zero
// ---------------------------------------------------------------------------------------------------------------------------------------
  fout.open("./data/debugFaceAreaAdd.dat");
  fout << "variables=\"dotProduct\"\n";
  fout << "FILETYPE=SOLUTION\nZONE NODES=" << _totalNodes << ", ELEMENTS=" << _totalHexCells << ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK\nVARLOCATION = ([1] = CELLCENTERED)\n";

  okFlag = true;
  std::cout << "3) Checking if face area vectors add up to zero in a cell : ";
  for (int i = 1; i <= _totalFaces; ++i) {
    right = _faceDataArray(i).rightCell;
    left = _faceDataArray(i).leftCell;

    if (right != -1) {
      tempCellDouble(right) -= tvmet::dot(_faceDataArray(i).normal, point(1.0));
    }
    if (left != -1) {
      tempCellDouble(left) += tvmet::dot(_faceDataArray(i).normal, point(1.0));
    }
  }
// ---------------------------------------------------------------------------------------------------------------------------------------
// Check if dot product sum of area vector and some constant vector adds up to zero
// ---------------------------------------------------------------------------------------------------------------------------------------
  for (int i = _start; i < _start + _totalHexCells; ++i) {
    std::vector<int>::iterator result = std::find(_ghostCells.begin(), _ghostCells.end(), i);
    if (result != _ghostCells.end()) tempCellDouble(i) = (0.0);
    fout << tempCellDouble(i) << "\n";
    if (fabs(tempCellDouble(i)) > 1e-14) {
      std::cout << "warning : Hex Cell " << i << " dot product sum  :" << tempCellDouble(i) << "\n";
      okFlag = false;
    }
  }
  for (int i = _start + _totalHexCells; i < _start + _totalCells; ++i) {
    std::vector<int>::iterator result = std::find(_ghostCells.begin(), _ghostCells.end(), i);
    if (result != _ghostCells.end()) tempCellDouble(i) = (0.0);
    fout << tempCellDouble(i) << "\n";
    if (fabs(tempCellDouble(i)) > 1e-14) {
      std::cout << "warning : PolyHedral Cell " << i << " dot product sum  :" << tempCellDouble(i) << "\n";
      okFlag = false;
    }
  }
  if (okFlag == false) {
    std::cout << "<failed>" << "\n";
    //exit(0);
  } else
    std::cout << "<passed>" << "\n";

  fout.close();
  tempCellDouble.free();
  std::cout <<  "---------------- END ---------------------" << "\n";
// ---------------------------------------------------------------------------------------------------------------------------------------
}

void StretchedStaticCylindricalGridSolver::groupFaces()
{
  int right, left;
  std::vector<int>::iterator result;
  faceMapGhost temp;
// loop through face to find if any ghost cell is the left/right neighbour
  for (int i = 1; i <= _totalFaces; ++i) {
    right = _faceDataArray(i).rightCell;
    left = _faceDataArray(i).leftCell;
    // Fist check if the right cell neighbours is a ghost cell
    result = std::find(_ghostCells.begin(), _ghostCells.end(), right);
    if (result == _ghostCells.end()) {
      // If not check if the left cell neighbours is a ghost cell
      result = std::find(_ghostCells.begin(), _ghostCells.end(), left);
      if (result == _ghostCells.end()) {
        // The face is interior
        _interiorFaces.push_back(&_faceDataArray(i));
      } else {
        // The face has a left ghost cell
        temp.fd_ptr = &_faceDataArray(i);
        temp.isLeftGhost = true;
        _ghostFaces.push_back(temp);
      }
    } else {
      // The face has a right ghost cell
      temp.fd_ptr = &_faceDataArray(i);
      temp.isLeftGhost = false;
      _ghostFaces.push_back(temp);
    }
  }
}

void StretchedStaticCylindricalGridSolver::setFreeStream(flux free)
{
  _free = free;
}

void StretchedStaticCylindricalGridSolver::setFreeStream(double rho,double u,double v,double w,double p)
{
  _free = flux(rho,u,v,w,p);
}


StretchedStaticCylindricalGridSolver::~StretchedStaticCylindricalGridSolver()
{
  ;
}

bool StretchedStaticCylindricalGridSolver::takeOneRK4Substep()
{
  // if (_CFL < 1e-10) {
  //  std::cout << "Time Step not set or set to a very small value CFL = \n" << _CFL << "\n";
  //  std::cout << "Setting to defualt CFL value of 1 and calculating local time \n" << "\n";
  //  _CFL = 1.0;
  //  calculateLocalDt();
  // }
  if (init_flag == false) {
    std::cout << "Solver not initilized with values (CFL = " << _CFL << ")" << "\n";
    return true;
  }
  if (_free[0] < 1e-6 || _free[4] < 1e-6) {
    std::cout << "Free stream values not initilized or set to zero values _free = \n" << _free << "\n";
    return true;
  }
  flux q, src;
  point vel;
  bool done(false);

  // Flow field RK pass
  calculateRHS();
  std::cout <<  "pass = " << _pass << "\n";
  if (_pass <= 2) {
    for (int i = 1; i <= _totalCells; ++i) {
      vel = point( _U(i)[IDX::u] , _U(i)[IDX::v] , _U(i)[IDX::w] );
      src = flux(0.0, _source(i)[0], _source(i)[1], _source(i)[2],tvmet::dot(vel, _source(i)));
      q = -_localDt(i) * ( _rhs(i) + src );
      _rhsRK(i) += numerals::fac[_pass] * q;
      _U(i) = _URK(i) + numerals::f_q[_pass] * q;
    }
  } else {
    for (int i = 1; i <= _totalCells; ++i) {
      vel = point( _U(i)[IDX::u] , _U(i)[IDX::v] , _U(i)[IDX::w] );
      src = flux(0.0, _source(i)[0], _source(i)[1], _source(i)[2], tvmet::dot(vel, _source(i)));
      q = -_localDt(i) * ( _rhs(i) + src );
      _rhsRK(i) += numerals::fac[_pass] * q;
      _U(i) = _URK(i) + _rhsRK(i);
    }
    done = true;
    _rhs = flux(0.0);
    _rhsRK = _rhs;
    _URK = _U;
    _currentTime += _localDt(1);
  }
  ++_pass;
  return done;
}

bool StretchedStaticCylindricalGridSolver::takeOneRK2Substep()
{
// if (_CFL < 1e-10) {
//  std::cout << "Time Step not set or set to a very small value CFL = \n" << _CFL << "\n";
//  std::cout << "Setting to defualt CFL value of 1 and calculating local time \n" << "\n";
//  _CFL = 1.0;
//  calculateLocalDt();
// }
  if (init_flag == false) {
    std::cout << "Solver not initilized with values (CFL = " << _CFL << ")" << "\n";
    return true;
  }
  if (_free[0] < 1e-6 || _free[4] < 1e-6) {
    std::cout << "Free stream values not initilized or set to zero values _free = \n" << _free << "\n";
    return true;
  }
  flux q, src;
  point vel;
  bool done(false);
// Flow field RK pass
  calculateRHS();
  if (_pass < 1) {
    for (int i = 1; i <= _totalCells; ++i) {
      vel = point(_U(i)[IDX::u], _U(i)[IDX::v], _U(i)[IDX::w]);
      src = flux(0.0, _source(i)[0], _source(i)[1], _source(i)[2], tvmet::dot(vel, _source(i)));
      q = -_localDt(i) * (_rhs(i) + src);
      _rhsRK(i) += 0.5 * q;
      _U(i) = _URK(i) + q;
    }
  } else {
    for (int i = 1; i <= _totalCells; ++i) {
      vel = point(_U(i)[IDX::u], _U(i)[IDX::v], _U(i)[IDX::w]);
      src = flux(0.0, _source(i)[0], _source(i)[1], _source(i)[2], tvmet::dot(vel, _source(i)));
      q = -_localDt(i) * (_rhs(i) + src);
      _rhsRK(i) += 0.5 * q;
      _U(i) = _URK(i) + _rhsRK(i);
    }
    _rhs = flux(0.0);
    _rhsRK = _rhs;
    _URK = _U;
    _currentTime += _localDt(1);
    done = true;
  }
  ++_pass;
  return done;
}

void StretchedStaticCylindricalGridSolver::setupForTimeStep()
{
  _rhs = flux(0.0);
  _rhsRK = _rhs;
  _URK = _U;
  _pass = 0;
  for (int i = 1; i <= _totalCells; ++i) {
    for(int j = 0; j <= 4; ++j) {
      if(std::isnan(_U(i)[j]) != 0) {
        std::cout<< "Error Nan in finite volume grid\n";
        exit(1);
      }
    }
  }
}


void StretchedStaticCylindricalGridSolver::setCellValue(int cellNumber, flux &F)
{
///if (cellNumber > _totalCells || cellNumber < 1)
//std::cout << "Error : Bad Index or Index Exceeded the cell index list i =" << cellNumber << "\n";
//else
  _U(cellNumber) = F;
}

int StretchedStaticCylindricalGridSolver::getTotalCells()
{
  return _totalCells;
}

int StretchedStaticCylindricalGridSolver::getTotalNodes()
{
  int dummy = _totalNodes;
  return dummy;
}

point StretchedStaticCylindricalGridSolver::getVelocityAtCell(int cellNum)
{
  return point(_U(cellNum)[IDX::u], _U(cellNum)[IDX::v], _U(cellNum)[IDX::w]);
}

flux StretchedStaticCylindricalGridSolver::getFluxAtCell(int cellNum)
{
  return _U(cellNum);
}

double StretchedStaticCylindricalGridSolver::getCellVolume(int cellNum)
{
  return _volume(cellNum);
}


flux StretchedStaticCylindricalGridSolver::getFluxAtCellConservative(int cellNum)
{
  double ru = _U(cellNum)[IDX::u] * _U(cellNum)[IDX::rho];
  double rv = _U(cellNum)[IDX::v] * _U(cellNum)[IDX::rho];
  double rw = _U(cellNum)[IDX::w] * _U(cellNum)[IDX::rho];
  double et = _U(cellNum)[IDX::p] / ( airStandardAtmosphere::gm1 ) +
              0.5 *_U(cellNum)[IDX::rho] * ( pow(_U(cellNum)[IDX::u],2) +
                  pow(_U(cellNum)[IDX::v],2) + pow(_U(cellNum)[IDX::w],2) );
  return flux( _U(cellNum)[IDX::rho] , ru , rv , rw , et );
}

double StretchedStaticCylindricalGridSolver::getTime()
{
  return _currentTime;
}

void StretchedStaticCylindricalGridSolver::setSourceValue(int cellNum, point &src)
{
  _source(cellNum) += src;
}

point StretchedStaticCylindricalGridSolver::getSourceValue(int cellNum)
{
  return _source(cellNum);
}

void StretchedStaticCylindricalGridSolver::zeroSource()
{
  _source = point(0.0);
}

void StretchedStaticCylindricalGridSolver::readFromTecplotFile(std::string filename)
{
  std::ifstream fin(filename.c_str());
  char dummy[2000];
  std::string temp;
  int tempNodes, tempElements;
  double tempValue;
// File header
  fin.getline(dummy, 2000);
  temp = dummy;
  if (temp != "variables = rho , u , v , w , p") {
    std::cout << "Not an acceptable solution file as file header is corrupted (header : " << temp << ")" << "\n";
    return;
  }
// File type
  fin.getline(dummy, 2000);
  temp = dummy;
  if (temp != "FILETYPE = SOLUTION") {
    std::cout << "Not an acceptable solution file (filetype : " << temp << ")" << "\n";
    return;
  }
// Zone header
  fin.getline(dummy, 2000);
  std::cout <<  "Read zone header :" << dummy << "\n";
  sscanf(dummy, "ZONE NODES = %d , ELEMENTS = %d , DATAPACKING = BLOCK, ZONETYPE = FEBRICK", &tempNodes, &tempElements);
  std::cout <<  "Read from file : contains totally " << tempNodes << " nodes and " << tempElements << " elements data" << "\n";
  if (this->_totalHexCells != tempElements) {
    std::cout << "The total elements in solution file does not match that in the grid file" << "\n";
    return;
  }
  if (this->_totalNodes != tempNodes) {
    std::cout << "The total nodes in solution file does not match that in the grid file" << "\n";
    return;
  }
// solution time stamp
  fin.getline(dummy, 2000);
  sscanf(dummy, "SOLUTIONTIME = %lf", &tempValue);
  std::cout <<  "Read from file : solution file written at time " << tempValue << "\n";
  this->_currentTime = tempValue;
// solution varlocation
  fin.getline(dummy, 2000);
  temp = dummy;
  if (temp != "VARLOCATION = ([1-5] = CELLCENTERED)") {
    std::cout << "Not an acceptable solution file (varlocation : " << temp << ")" << "\n";
    return;
  }
  std::cout <<  "File looks ok ... reading data " << "\n";
  for (int i = 1; i <= _totalHexCells; ++i)
    fin >> this->_U(i)[IDX::rho];
  for (int i = 1; i <= _totalHexCells; ++i)
    fin >> this->_U(i)[IDX::u];
  for (int i = 1; i <= _totalHexCells; ++i)
    fin >> this->_U(i)[IDX::v];
  for (int i = 1; i <= _totalHexCells; ++i)
    fin >> this->_U(i)[IDX::w];
  for (int i = 1; i <= _totalHexCells; ++i)
    fin >> this->_U(i)[IDX::p];
  // Get Centerline Polyhedral cell values
  fin.getline(dummy, 2000);
  fin.getline(dummy, 2000);
  for (int i = _totalHexCells + 1; i <= _totalCells; ++i)
    fin >> _U(i)[IDX::rho] >> _U(i)[IDX::u] >> _U(i)[IDX::v] >> _U(i)[IDX::w] >> _U(i)[IDX::p];
  fin.close();
}

void StretchedStaticCylindricalGridSolver::setOutflowBC()
{
  _outflowBC = true;
}

void StretchedStaticCylindricalGridSolver::unsetOutflowBC()
{
  _outflowBC = false;
}

std::vector<int> & StretchedStaticCylindricalGridSolver::ghostList()
{
  return _ghostCells;
}

//    ---------------------------------------
//      Higher order solver implementation
//    ---------------------------------------
void StretchedStaticCylindricalGridSolver::gradientReconstruct()
{
  _DU = _5x3_(0.0);
// Green-Gauss reconstruction
  for (int i = 1 ; i <= _totalFaces ; ++i) {
    _5x3_ temp;
    point dA = _faceDataArray(i).normal;
    flux Uf = flux(_U(_faceDataArray(i).leftCell) + _U(_faceDataArray(i).rightCell));
    Uf *= 0.5;
    temp = Uf[0]*dA[0], Uf[0]*dA[1], Uf[0]*dA[2],
    Uf[1]*dA[0], Uf[1]*dA[1], Uf[1]*dA[2],
    Uf[2]*dA[0], Uf[2]*dA[1], Uf[2]*dA[2],
    Uf[3]*dA[0], Uf[3]*dA[1], Uf[3]*dA[2],
    Uf[4]*dA[0], Uf[4]*dA[1], Uf[4]*dA[2];
    _DU(_faceDataArray(i).leftCell)  += temp;
    _DU(_faceDataArray(i).rightCell) -= temp;
  }
}

void StretchedStaticCylindricalGridSolver::calculateRHS()
{
// To add -- Higher order reconstruction
  for (int i = 1 ; i <= _totalFaces ; ++i) {
    int leftCell = _faceDataArray(i).leftCell;
    int rightCell = _faceDataArray(i).rightCell;
    point normal = utility::normalize(_faceDataArray(i).normal);
    double dA = utility::mag(_faceDataArray(i).normal);
    point dr_L = point(_cell_centroid(leftCell) - _faceDataArray(i).centroid);
    point dr_R = point(_cell_centroid(rightCell) - _faceDataArray(i).centroid);
    flux leftState = flux(_U(leftCell) + tvmet::prod( _DU(leftCell) , dr_L ));
    flux rightState = flux(_U(rightCell) + tvmet::prod( _DU(leftCell) , dr_L ));
    flux fluxInterface;
    roe(leftState, rightState, normal, fluxInterface);
    _rhs(leftCell) += 0.5 * dA / _volume(leftCell) * fluxInterface;
    _rhs(rightCell) -= 0.5 * dA / _volume(rightCell) * fluxInterface;
  }
  for (std::vector<int>::iterator iter = _ghostCells.begin(); iter != _ghostCells.end(); ++iter)
    _rhs(*iter) = flux(0.0);
}

