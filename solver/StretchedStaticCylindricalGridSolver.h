#ifndef __STRETCHEDSTATICCYLINDRICALGRIDSOLVER_H__

#define __STRETCHEDSTATICCYLINDRICALGRIDSOLVER_H__

#include "definition.h"
#include<stdio.h>

class StretchedStaticCylindricalGridSolver
{
public:
// Empty constructr
  StretchedStaticCylindricalGridSolver() {
    ;
  }
  StretchedStaticCylindricalGridSolver(std::string);
  void Constructor(std::string);
// Functions returing value
  int getTotalCells();
  int getTotalNodes();
  double getTime();
  double getCellVolume(int);
  void takeTimeAverageData();
  void writeTimeAverageSolution(std::string path);
// Functions returing object
  point getSourceValue(int cellNum);
  point getVelocityAtCell(int cellNum);
  flux getFluxAtCell(int cellNum);
  flux getFluxAtCellConservative(int cellNum);
  std::vector<int> &ghostList();
// Functions returing void
  void takeOneRKTimeStep();
  void takeOneRKTimeStepSplit();
  bool takeOneRK4Substep();
  bool takeOneRK2Substep();
  void setOutflowBC();
  void unsetOutflowBC();
  void takeOneEulerTimeStep();
  void setCFL(double CFL);
  void Init(flux &f);
  void Init(double rho,double u,double v,double w,double p);
  void setFreeStream(flux free);
  void setFreeStream(double rho,double u,double v,double w,double p);
  void setCellValue(int cellNumber, flux &F);
  void setCurrentTime(double t);
  void setDt(double t);
  void setupForTimeStep();
  void setSourceValue(int cellNum, point &src);
  void zeroSource();
  void readFromTecplotFile(std::string filename);
  void writeRestart(std::string path);
  void readRestart(std::string path);
  void writeSolution(std::string file);
  void gradientReconstruct();
  double calculateDtCFL(double CFL, double vel);

//Depreciated Function
//StretchedStaticCylindricalGridSolver(UniformStaticCylindricalGrid &G, flux &FreeStream); // Depreacted just here for usefulness of the code
// Destructor RAII

  ~StretchedStaticCylindricalGridSolver();

private:
  int _ntime_samples;
// \param init_flag is used to check if the class is active
  bool init_flag;
// \param _pass keeps track of which RK pass we are currently residing
  int _pass;
// \param _totalCells stores the total number of cells in the computational domain
  int _totalCells;
// \param _totalFaces stores the total number of Faces in the computational domain
  int _totalFaces;
// \param _totalFaces stores the total number of Faces in the computational domain
  int _totalNodes;
// \param _start stores the starting index of the cell in case there is zero offset
  int _start;
// \param _hexStart stores the starting index of the hex cells
  int _hexStart;
// \param _hexStart stores the starting index of the hex cells
  int _totalHexCells;
// \param _onePolyCellFaces conatins the number of faces that form one single Poly Cell
  int _onePolyCellFaces;
// \param _onePolyCellNodes conatins the number of nodes that form one single Poly Cell
  int _onePolyCellNodes;
  double _CFL;
  double _currentTime;
// \param _q is the 1D blitz++ array of flux data type and stores state vector
  fluxArray1 _U,_Uavg;
// \param _qBase is the 1D blitz++ array of flux data type and stores state vector for RK time-stepping
  fluxArray1 _URK;
// \param _rhs is the 1D blitz++ array of flux data type to store the right hand side of the area integral over the control volume
  fluxArray1 _rhs;
// \param _rhs is the 1D blitz++ array of flux data type to store the right hand side of the area integral over the control volume for RK time-stepping
  fluxArray1 _rhsRK;
// \param _FreeStream stores the free stream flux used for the Riemann boundary condition
  flux _FreeStream;
// \param _CellConnectivity contains the cell connectivity from the Grid Object
  cellFaceConnectivityArray _cellFaceConnectivityArray;
// \param _CellConnectivity contains the cell connectivity from the Grid Object
  faceDataArray _faceDataArray;
// \param _volume contains the cell volume from the Grid Object
  blitz::Array<double, 1 > _volume;
// \param _5x3_ matrix containing the reconstruction matrix
  blitz::Array<_5x3_, 1 > _DU;
// \param _cell_centroid contains the cell centroid from the Grid Object
  blitz::Array<point, 1 > _cell_centroid;
//\param _boundaryCells list of boundary nodes
  std::vector<int> _boundaryCells;
//\param _ghostCells list of boundary nodes
  std::vector<int> _ghostCells;
//\param _free the free stream values
  flux _free;
//\param _outflowBC a boolean indicating if the outflow BC needs to be switched on/off
  bool _outflowBC;

  typedef struct {
    faceData *fd_ptr;
    bool isLeftGhost;
  } faceMapGhost;
//

  typedef struct {
    faceData *fd_ptr;
    bool isLeftBoundary;
  } faceMapBoundary;
//
  std::vector<faceData *> _interiorFaces;
//
  std::vector<faceMapGhost> _ghostFaces;
// Characteristic length for local time stepping
  blitz::Array<double, 1 > _charLength;
// local time step caluclated based on CFL criteria
  blitz::Array<double, 1 > _localDt;
// Source term
  blitz::Array<point, 1 > _source;
// Private Functions returning void
  void gridCheck();
  void calculateRHS_OBC();
  void calculateRHS();
  void roe(flux &UL, flux &UR, point &normal, flux &F);
  void writeHexCellsSolution(fluxArray1 &RHS, std::string);
  void groupFaces();
  void calculateLocalDt();

#ifdef USE_MPI
  // Cells below the top ghost cells
  fluxArray1 _topSendBuf , _topRecvBuf;
  // Cells above the bottom ghost cells
  fluxArray1 _botSendBuf , _botRecvBuf;
#endif

  inline double ramp(double t) {
    if (t < 0.1)
      return t / 0.1;
    else
      return 1.0;
  }
// Depreciated
  std::vector<faceMapBoundary> _boundaryFaces;
};

#endif
