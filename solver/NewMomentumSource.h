#ifndef _NEW_MOMENTUMSOURCE_H

#define	_NEW_MOMENTUMSOURCE_H

#include "definition.h"
#include "StretchedStaticCylindricalGridSolver.h"
#include <stdio.h>

typedef bladeObject::pair pair;

class NewMomentumSource
{
public:
  NewMomentumSource();
  NewMomentumSource(StretchedStaticCylindricalGridSolver *,std::string filename );
// Functions returning void
  void Constructor(StretchedStaticCylindricalGridSolver *,std::string filename);
  void setTipCorLen(double tipCorLen);
  void recalRevs();
  void writeSource(std::string filename);
  void writeLoad(std::string filename);
  void calculateHigherOrderDiscreteSourceIntegral();
  void readQuadratureNodes(std::string);
  void setBEMTBlade(bladeObject::BEMT_Rotor_Blade & );
  void parseGridFile(std::string filename);
// Destructor RAII
  virtual ~NewMomentumSource();

private:
// Private Data
  bool _bladeDefined;
  int _revs;
  int _imax;
  int _jmax;
  int _kmax;
  int _totalCells; // This not the total cells count only the Hex cell count
  int _rotorZIdx; // The Z index (of the array zdist) in which the rotor is placed
  int _rotorZIdxPad; // The padding above and below the Z index (of the array zdist)
                     // to limit the volume integration to only few cells along the 
                     // z direction (I want this to be done automatically but due to
                     // time constraints I am making this a manual input)
  int _rStart,_rEnd; // The starting and ending index of the rdist array - used to
		     // limit the number of source cells used for volume integration 
                     // along the r direction
  double _dtheta;
  double _tipCorLen;
  double _omegaWeight;
// Private Objects
  blitz::Array<double,1> _rdist;
  blitz::Array<double,1> _zdist;
  blitz::Array<point,2> _inflow;
  blitz::Array<int,2> _inflowCells;
  blitz::Array<point,1> _source;
  bladeObject::BEMT_Rotor_Blade _BEMT_Rotor_Blade;
  StretchedStaticCylindricalGridSolver *_S;
  blitz::Array<tvmet::Matrix<double,27,3>,1> _quadNodes;
  blitz::Array<tvmet::Vector<double,27>,1> _jacobian;
  blitz::Array<int,3> _ijk2cellNum;
  std::list<int> _sourceCells;
// Private Functions
// Function returing void
  void calculateBladeLoads(double *);
  void convertPeriodic(double *T);
  void countRevs(double theta);
  void interpolateInflow();
  void mapijk2CellNum();
// Function returning value
  double X(int i, int j);
  double Y(int i, int j);
// Function returning objects
  point getLoadAtBladeRadius(int blade, double r);
  std::string getStrippedInputFile(const char *filename);
  tvmet::Vector<double,3> getCellCentroidAtijk(int i,int j,int k);
};

#endif	/* _MOMENTUMSOURCE_H */

