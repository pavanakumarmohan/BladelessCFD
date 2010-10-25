#include "NewMomentumSource.h"
#include "Hexa8p3.h"

NewMomentumSource::NewMomentumSource() {}

NewMomentumSource::NewMomentumSource(StretchedStaticCylindricalGridSolver *S,std::string filename)
{
  Constructor(S,filename);
}

void NewMomentumSource::setTipCorLen(double tipCorLen)
{
  _tipCorLen = tipCorLen;
}

void NewMomentumSource::Constructor(StretchedStaticCylindricalGridSolver *S,std::string filename)
{
  _tipCorLen = 0.8;
  _S = S;
  // Parse the grid inoput file to get the necessary connectivity and node distribution information
  parseGridFile(filename);
  readQuadratureNodes("QuadratureNodes.dat");
  // Map all the i,j,k cell index to cell number (structured to unstructued data)
  _ijk2cellNum.resize(Range(1,_imax),Range(1,_jmax),Range(1,_kmax+1));
  mapijk2CellNum();
  // Calculate the _rStart and _rEnd index
  _rStart = -1;
  _rEnd = -1;
  for(int i = 0; i < _rdist.size() ; ++i){
   if( _BEMT_Rotor_Blade.rootCutOut * _BEMT_Rotor_Blade.rotorRadius >=  _rdist(i) && _rStart < 0 )
    _rStart = i;
   if( _BEMT_Rotor_Blade.rotorRadius >= _rdist(i) && _rEnd < 0 )
    _rEnd = i;
  }
  if( _rEnd < 0 || _rEnd == 1 ){
   std::cerr << "Error in defining blade radius.\nIt" << 
                  "is inconsistent with grid extent\n" << "_rEnd = " << _rEnd << std::endl;
   exit(2);
  }
  if( _rStart < 0 || _rStart == 1 ){
   std::cerr << "Error in defining blade root cut out.\nIt" << 
                  "is inconsistent with grid extent\n" << "_rStart = " << _rStart << std::endl;
   exit(2);
  }
  if( _rotorZIdx - _rotorZIdxPad < 0 ){
   std::cerr << "Error in defining the blade padding along Z.\nIt" << 
                  "exceeds grid index/extent\n" << "_rotorZIdx - _rotorZIdxPad = " << 
                  _rotorZIdx - _rotorZIdxPad << std::endl;
   exit(2);
  }
  if( _rotorZIdx + _rotorZIdxPad > _kmax ){
   std::cerr << "Error in defining the blade padding along Z.\nIt" << 
                  "exceeds grid index/extent\n" << "_rotorZIdx + _rotorZIdxPad = " <<
                  _rotorZIdx + _rotorZIdxPad << std::endl;
   exit(2);
  }
  if( _rotorZIdx < 1 || _rotorZIdx > _kmax ){
   std::cerr << "Error in defining the blade Z location index.\nIt" << 
                  "exceeds grid index/extent\n" << "_rotorZIdx = " <<
                  _rotorZIdx << std::endl;
   exit(2);
  }

  _BEMT_Rotor_Blade.interpLoad = new blitz::Array<point, 1 > [_BEMT_Rotor_Blade.nBlade];
  _BEMT_Rotor_Blade.FWHM_Y /= 2.35482005; // 2.0 ln(2)
  _BEMT_Rotor_Blade.FWHM_Z /= 2.35482005; // 2.0 ln(2)
  _BEMT_Rotor_Blade.FWHM_Y = 2.0 * pow(_BEMT_Rotor_Blade.FWHM_Y * _BEMT_Rotor_Blade.rootChord , 2);
  _BEMT_Rotor_Blade.FWHM_Z = 2.0 * pow(_BEMT_Rotor_Blade.FWHM_Y * _BEMT_Rotor_Blade.rootChord, 2);
 
  // Get the inflow cells and map the i,j --> cell number for the _inflowCells connectivity
  _inflowCells.resize(Range(2, _imax), Range(1, _jmax));
  _inflow.resize(Range(2, _imax), Range(1, _jmax));
  for (int i = 2; i <= _imax; ++i) {
   for (int j = 1; j <= _jmax; ++j) {
    int k = _rotorZIdx;
    _inflowCells(i,j) = _ijk2cellNum(i,j,k);
   }
  }
  // 4) Free the i,j,k --> cell number connectivity as it is no longer required
  _revs = 1;
}

void NewMomentumSource::recalRevs()
{
  double time = _S->getTime();
  _omegaWeight = _BEMT_Rotor_Blade.omega * time / ( 4.0 * M_PI ) ;
  if(_omegaWeight > 1.0) {
    _omegaWeight = 1.0;
    double t_in_ramp = 4.0 * M_PI / _BEMT_Rotor_Blade.omega;
    _revs = ceil( (time - t_in_ramp) * _BEMT_Rotor_Blade.omega / numerals::_2_MPI + 0.5 );
  } else
    _revs = ceil( pow( _BEMT_Rotor_Blade.omega*time , 2.0 ) / (8.0 * M_PI * M_PI) );
}

void NewMomentumSource::countRevs(double theta)
{
  if (theta > _revs * 2.0 * M_PI)
    ++_revs;
}


point NewMomentumSource::getLoadAtBladeRadius(int blade, double r)
{
  // To implement
  // 1) Find the i location corresponding to the rotor r location
  // 2) Do linear interpolation using the iLoc and iLoc + 1 radial values
  //    a) Check if the given location is less than radius of the blade
  //       i) If yes return the blade load or
  // 	 ii) Else return 0.0
  /*if( r < _BEMT_Rotor_Blade.rotorRadius ){
    point loadI = _BEMT_Rotor_Blade.interpLoad[blade](iLoc);
    point loadIp1 = _BEMT_Rotor_Blade.interpLoad[blade](iLoc + 1);
    return point(loadI + (loadI - loadIp1) * (r - (iLoc - 0.5) * _dr) / _dr);
    }else
    return point(0.0);*/
}

void NewMomentumSource::calculateBladeLoads(double *theta)
{
  // Get Inflow at the current blade position
  blitz::Array< point , 1 > tempLoad;
  tempLoad.resize(Range(_rStart, _rEnd));
  tempLoad = point(0.0);
  interpolateInflow();
  int jLoc;
  double F, f;
  for (int i = 0; i < _BEMT_Rotor_Blade.nBlade; ++i) {
    jLoc = int(( theta[i] / _dtheta )) + 1;
    if (jLoc > _jmax && jLoc <= _jmax + 1)
      jLoc = 1;
    for (int counter = _rStart; counter <= _rEnd ; ++counter) {
      // Calculate the inflow using linear interpolation
      // used to smooth the induced velocity variation when changing from one cell to another
      double linInflow;
      if (jLoc == 1) {
        linInflow = std::min(_inflow(counter, jLoc)[2], _inflow(counter, _jmax)[2]);
      } else {
        linInflow = std::min(_inflow(counter, jLoc)[2], _inflow(counter, jLoc - 1)[2]);
      }
      // To implement
      // 1) Find r based on the counter value
      // double r = ( counter - 0.5 ) * _dr;
      double r = _rdist(counter) ; /*Fixed */
      double vel = sqrt( pow(linInflow, 2) + pow(_BEMT_Rotor_Blade.omega * r * _omegaWeight, 2) );
      // phi = atan( v_i / (Omega * r) )
      double phi = atan2( linInflow , _BEMT_Rotor_Blade.omega * r * _omegaWeight );
      // alpha = theta - phi
      double alpha = _BEMT_Rotor_Blade.CollPitch + phi;
      double rho = _S->getFluxAtCell(_inflowCells(counter, jLoc))[IDX::rho];
      // dL = 0.5 * V^2 * c * (2 * pi * alpha) * dr
      double dL = rho * pow(vel, 2) * _BEMT_Rotor_Blade.rootChord * M_PI * alpha;
      F =  0.5 * ( 1.0 - cos( ( r - _tipCorLen *  _BEMT_Rotor_Blade.rotorRadius  ) /
                                (_BEMT_Rotor_Blade.rotorRadius * ( 1.0 - _tipCorLen ) ) * M_PI - M_PI ) );
      // Fx = F_r Cos(\theta) - F_{\theta} Sin(\theta)
      tempLoad(counter)[0] = dL * sin(phi) * sin(theta[i]);
      // Fy = F_r Sin(\theta) + F_{\theta} Cos(\theta)
      tempLoad(counter)[1] = -dL * sin(phi) * cos(theta[i]);
      // Fz
      tempLoad(counter)[2] = dL * cos(phi) * F;
      if(std::isnan(tempLoad(counter)[2]))
        std::cout << r << " " << vel << " " << phi << " " << alpha << " " << rho << " " << dL << "\n";
    }
    //_BEMT_Rotor_Blade.load[i].addLoad(tempLoad, _S->getTime());
    _BEMT_Rotor_Blade.interpLoad[i].resize(tempLoad.shape());
    _BEMT_Rotor_Blade.interpLoad[i] = tempLoad;
  }
}

void NewMomentumSource::convertPeriodic(double* T)
{
  for (int i = 0; i < _BEMT_Rotor_Blade.nBlade; ++i) {
    if (T[i] - 2.0 * M_PI >= 0.0) {

      int fac = int(T[i] / (2.0 * M_PI));
      T[i] -= fac * 2.0 * M_PI;
    }
  }
}

/*
void NewMomentumSource::writeSource(std::string filename) {
 std::stringstream s;
 s << filename << "/source.dat" ;
 std::ofstream fout(s.str().c_str());
 fout << "variables=Sx,Sy,Sz\n";
 fout << "FILETYPE=SOLUTION\nZONE NODES=" << _totalNodes <<
  ", ELEMENTS=" << _totalCells <<
  ", DATAPACKING=BLOCK, ZONETYPE=FEBRICK\nVARLOCATION = ([1-3] = CELLCENTERED)\n";

 for (int i = 1; i <= _totalCells; ++i)
  fout << _S->getSourceValue(i)[0] << "\n";
 for (int i = 1; i <= _totalCells; ++i)
  fout << _S->getSourceValue(i)[1] << "\n";
 for (int i = 1; i <= _totalCells; ++i)
  fout << _S->getSourceValue(i)[2] << "\n";

 fout << "ZONE I=" << _kmax + 1 << "\n";
 for (int i = _totalCells + 1; i <= _totalCells + _kmax + 1; ++i)
  fout << _S->getSourceValue(i)[0] << " " << _S->getSourceValue(i)[1] << " "<< _S->getSourceValue(i)[2] << "\n";
 fout.close();
}*/


/* Have to change this print the r values instead of r index values.
   This requires change in converter program.
*/
void NewMomentumSource::writeLoad(std::string filename)
{
  std::ofstream fout;
  fout.open(filename.c_str(), std::ios::app);
  if(fout.fail()) {
    std::cout << "Error Opening loads file " << filename.c_str() << std::endl;
    exit(1);
  }
  fout << "ZONE I=" << _BEMT_Rotor_Blade.nBlade << " J=" << _BEMT_Rotor_Blade.nSpan - 1 << "\n";
  fout << "SOLUTIONTIME=" << _S->getTime() << "\n";
  for (int i = 0; i < _BEMT_Rotor_Blade.nBlade; ++i) {
    for (int j = _rStart; j <= _rEnd ; ++j ) {
      fout << _rdist[j] << " " << _BEMT_Rotor_Blade.interpLoad[i](j)[0] << " " <<
           _BEMT_Rotor_Blade.interpLoad[i](j)[1] << " " << _BEMT_Rotor_Blade.interpLoad[i](j)[2] << "\n";
    }
  }
}

void NewMomentumSource::interpolateInflow()
{
  for (int i = 2; i <= _imax; ++i) {
    for (int j = 1; j <= _jmax; ++j) {
      _inflow(i, j) = _S->getVelocityAtCell(_inflowCells(i, j));
    }
  }
}

inline double NewMomentumSource::X(int i, int j)
{
  // To be implemented
}

inline double NewMomentumSource::Y(int i, int j)
{
  // To be implemented
}

NewMomentumSource::~NewMomentumSource() {}

void NewMomentumSource::readQuadratureNodes(std::string filename)
{
  bool debug = false;
  std::ifstream fin(filename.c_str());
  if (fin.fail()) {
    std::cout << "Error : Unable to open file " << filename << "\n";
    exit(-1);
  }
  //To implement
  // 1) Do something here to fix the _totalCells
  _quadNodes.resize(Range(1,_totalCells));
  _jacobian.resize(Range(1,_totalCells));
  // 2) Read quad nodes from file till end of file
  std::cout << "Reading In Quadrature nodes" << std::endl;
  for(int i = 1 ; i <= _totalCells ; ++i ){
   for(int j = 0 ; j < 27 ; ++j ){
    fin >> _quadNodes(i)(j,0) >> _quadNodes(i)(j,1) >> _quadNodes(i)(j,2) >> _jacobian(i)(j);
   }
  }
  fin.close();
}

void NewMomentumSource::calculateHigherOrderDiscreteSourceIntegral()
{
  if (_bladeDefined == false) {
    std::cout << "ERROR: Rotor and blade inputs not set" << "\n";
    exit(0);
  }
  ////////////////////////////////////////////////////////////////////////////
  // Some variables used in the Higer Order Volume integartion
  // source --> Used to store the blade load value at the 27 quadrature nodes
  // zeroSource --> This is just a dummy variable to store zero's
  // t --> Current time step value
  // theta --> The curent azimuthal postion of each individual blade
  // xy_p --> In plane (blade) coordinates of the quadrature nodes
  // tgt,norm --> The in plane (blade) tangential and normal direction of Actuator line (blade)
  // hex --> An instance of the Hexa8p3 class used for the volume integration
  // temp_point --> A dummy variable used to transfer values
  ////////////////////////////////////////////////////////////////////////////
  tvmet::Matrix<double,27,3> source,zeroSource(0.0);
  double t = _S->getTime(), theta[_BEMT_Rotor_Blade.nBlade],norm_prod,d,dz;
  tvmet::Vector<double, 2 > xy_p, tgt, norm;
  Hexa8p3 hex;
  point temp_point;
  /////////////////////////////////////////////////////////////////////////////
  // The weighting function _omegaWeight to ramp the source value slowly from
  // zero.....And also calculate the current rotor azimuthal position based on t.
  /////////////////////////////////////////////////////////////////////////////
  _omegaWeight = _BEMT_Rotor_Blade.omega * t / ( 4.0 * M_PI ) ;
  if(_omegaWeight > 1.0) {
    _omegaWeight = 1.0;
    double t_in_ramp = 4.0 * M_PI / _BEMT_Rotor_Blade.omega;
    for (int i = 0; i < _BEMT_Rotor_Blade.nBlade; ++i)
      theta[i] = _BEMT_Rotor_Blade.psiZero + i * (2.0 * M_PI / _BEMT_Rotor_Blade.nBlade) + (t - t_in_ramp) * _BEMT_Rotor_Blade.omega + M_PI;
  } else {
    for (int i = 0; i < _BEMT_Rotor_Blade.nBlade; ++i)
      theta[i] = _BEMT_Rotor_Blade.psiZero + i * (2.0 * M_PI / _BEMT_Rotor_Blade.nBlade) + pow( _BEMT_Rotor_Blade.omega * t , 2.0 ) / (4.0 * M_PI) ;
  }
  ///////////////////////////////////////////////////////////////////////////////
  // Keep track of number of revolutions
  // And make sure the current azimuthal postion of the rotor is between 0 - 2 \pi
  // Print some debug information
  // Calcualte the blade loads using the inflow values at the cell center
  ///////////////////////////////////////////////////////////////////////////////
  countRevs(theta[0]);
  convertPeriodic(theta);
  std::cout << "Blade Position : (" << _revs << " th rev)" << "\n";
  for (int i = 0; i < _BEMT_Rotor_Blade.nBlade; ++i)
    std::cout << "Blade : " << i << " Pos : " << theta[i] << "omegaWeight : " << _omegaWeight << "\n";
  std::cout << "Time : " << t << "\n";
  calculateBladeLoads(theta);
  ////////////////////////////////////////////////////////////////////////////////
  // Loop over all cells to find out the volume integral of source term
  ////////////////////////////////////////////////////////////////////////////////
  source = zeroSource;
  // To implement
  // for each cell in _inflowCells (identified at the constructor) do the source calculation
  
   for (int i = 1; i <= 10/* it should be the total source cells*/; ++i) {

    temp_point = point(0.0);
    for (int quad_count = 0; quad_count < 27; ++quad_count) {
      ////////////////////////////////////////////////////////////////
      // I am pre-computing the dz value here because there was some
      // anomaly when I tired computing the dz for the individual
      // blade as each blade gave a different dz value creating
      // asymmetry in the source terms. I have not figured out the
      // cause for the problem so I have change the way I am calculating
      // the dz value. But this is not a generalized way to calculate
      // the dz value because this assumes that the rotor has no
      // coning (a circular disk). If someone plans to incorporate
      // coning then they have to put this in the loop blade=1:nBlade
      ////////////////////////////////////////////////////////////////
      dz = _quadNodes(i)(quad_count,2) - _zdist(_rotorZIdx);
      /////////////////////////////////////////////////////////////////
      // Evaluate the source value at the 27 quadrature nodes for each
      // rotor and the cumulative value is stored in source(27,3)
      /////////////////////////////////////////////////////////////////
      for (int blade = 0; blade < _BEMT_Rotor_Blade.nBlade; ++blade) {
        xy_p[0] = _quadNodes(i)(quad_count,0);
        xy_p[1] = _quadNodes(i)(quad_count,1);
        norm[0] = cos(theta[blade]);
        norm[1] = sin(theta[blade]);
        norm_prod = tvmet::dot(norm, xy_p);
        if (norm_prod > 0.0 && norm_prod >= _BEMT_Rotor_Blade.rootCutOut * _BEMT_Rotor_Blade.rotorRadius ) {
          tgt[0] = -sin(theta[blade]);
          tgt[1] = cos(theta[blade]);
          d = tvmet::dot(tgt, xy_p);
          temp_point =	 getLoadAtBladeRadius(blade, norm_prod) *
                         exp( - dz * dz / _BEMT_Rotor_Blade.FWHM_Z ) *
                         exp( - d * d / _BEMT_Rotor_Blade.FWHM_Y ) /
                         ( M_PI * sqrt( _BEMT_Rotor_Blade.FWHM_Y * _BEMT_Rotor_Blade.FWHM_Z ) );
          // Transfer source value to fun array
          for (int dimension = 0; dimension < 3; ++dimension)
            source(quad_count,dimension) += temp_point(dimension);
        }
      }
      /////////////////////////////////////////////////////////////////
      // The Hexa8p3 class is used to do the volume integration.
      /////////////////////////////////////////////////////////////////
    }
    hex.setJacobian( _jacobian(i) );
    temp_point = hex.volIntegral(source);
    _S->setSourceValue(i, temp_point);
    source = zeroSource;
  }
  ////////////////////////////////////////////////////////////////
}


// Newly added
void NewMomentumSource::parseGridFile(std::string filename)
{
  // Parse the grid input file
  std::string dd=getStrippedInputFile(filename.c_str());
  std::istringstream inputFile(dd);
  double zoffset,trans_xyz;
  std::string pathPrefix;
  inputFile >> pathPrefix >> _jmax >> _imax >> _kmax;
  _rdist.resize(_imax);
  _zdist.resize(_kmax);
  for(int i = 0 ; i < _imax ; ++i ) {
    inputFile >> _rdist(i);
  }
  for(int i = 0 ; i < _kmax ; ++i ) {
    inputFile >> _zdist(i);
  }
  // Read translation and rotation info
  inputFile >> trans_xyz;
  inputFile >> trans_xyz;
  inputFile >> trans_xyz;

  inputFile >> trans_xyz;
  inputFile >> trans_xyz;
  inputFile >> trans_xyz;
  // Get grid inputs like rotor position along Z and padding along Z 
  inputFile >> _rotorZIdx;
  inputFile >> _rotorZIdxPad;
  // Get grid inputs like root cut out and end of blade  
  inputFile >> _BEMT_Rotor_Blade.rootCutOut;
  inputFile >> _BEMT_Rotor_Blade.rotorRadius;
  inputFile >> _BEMT_Rotor_Blade.nBlade;
  inputFile >> _BEMT_Rotor_Blade.psiZero;
  inputFile >> _BEMT_Rotor_Blade.omega;
  inputFile >> _BEMT_Rotor_Blade.FWHM_Z;
  inputFile >> _BEMT_Rotor_Blade.FWHM_Y;
  inputFile >> _BEMT_Rotor_Blade.CollPitch;
}

// Parser
std::string NewMomentumSource::getStrippedInputFile(const char *filename)
{
  std::ifstream namelist(filename);
  if(namelist.fail()) {
    std::cerr<<"Unable to open input file :"<<filename<<std::endl;
    exit(0);
  }
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
    testStr.erase(std::remove(testStr.begin(), testStr.end(), ' '), testStr.end());
    // Check
    //For comment or new line at beginning of line
    if ( testStr[size_t(0)] == '#' || testStr[size_t(0)] == '\n') {
      //std::cout << "Skipping" << std::endl;
    }// Parse the string
    else {
      pos = testStr.find('#');
      if (pos != std::string::npos) { // remove comment
        testStr = testStr.substr(0, pos);
      }
      parsedStr.append(testStr);
      parsedStr.append("\n");
    }
  }
  parsedStr.erase(parsedStr.length() - 1, 1);
  return parsedStr;
}

// Newly added
void NewMomentumSource::mapijk2CellNum(){
 int counter = 1;
 // Map all Hex Cells
 for (int i = 2; i <= _imax; ++i) {
  for (int j = 1; j <= _jmax; ++j) {
   for (int k = 1; k <= _kmax + 1; ++k) {
    _ijk2cellNum(i,j,k) = counter;
    counter++;
   }
  }
 }
 // Not mapping the Poly Cells
}

// Newly added
tvmet::Vector<double,3> NewMomentumSource::getCellCentroidAtijk(int i,int j,int k){
 return tvmet::Vector<double,3>(0.0);
}

