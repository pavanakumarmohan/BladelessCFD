#include "StretchedStaticCylindricalGridSolver.h"
#include "NewMomentumSource.h"
#include<iostream>
#include<sstream>
#include<string>

typedef struct {
// Solver Inputs common to all process
  double dt;
  int endIter;
  int rstFlag; // Restart flag
  char rstPrefix[2048]; // The iteration number for restart
  int timeIntegrationType; // Get the correct time integration type
  int nBlocks;
  int solOutIter,loadOutIter;
  char solOutPath[2048],loadOutPath[2048],filePrefix[2048];
  double tipCorLen;
} inputInfo;

std::string getStrippedInputFile(const char *filename);
typedef bladeObject::pair pair;

int main(int nargs, char **args)
{
  inputInfo ifo;
  StretchedStaticCylindricalGridSolver S;
  NewMomentumSource NM;
//NewMomentumSource MS;
  bladeObject::BEMT_Rotor_Blade blade;
  std::string dd;
  /****************************************************************
    Command line argument checking
   *****************************************************************/
  if(nargs < 2)
    dd=getStrippedInputFile("input.dat");
  else
    dd=getStrippedInputFile(args[1]);
  std::istringstream inputFile(dd);

  /****************************************************************
    Parse Input and read input
   *****************************************************************/
  inputFile >> ifo.nBlocks >> ifo.dt
            >> ifo.endIter >> ifo.rstFlag
            >> ifo.rstPrefix >> ifo.timeIntegrationType
            >> ifo.solOutIter >> ifo.solOutPath
            >> ifo.loadOutIter >> ifo.loadOutPath
            >> ifo.filePrefix;

  /****************************************************************
  Read Blade from input file
  *****************************************************************/
  inputFile >> blade.nBlade >> blade.omega >> blade.psiZero >> blade.rootCutOut >> blade.rotorRadius
            >> blade.location[2] >> blade.FWHM_Y >> blade.FWHM_Z >> blade.CollPitch >> blade.rootChord
            >> blade.CT >> ifo.tipCorLen ;

  /****************************************************************
    Print all the read input for verification
   *****************************************************************/
  std::cout << "*****************************************************************\n";
  std::cout << "\t\tnBlocks : " << ifo.nBlocks << "\n";
  std::cout << "\t\tDt : " << ifo.dt << "\n";
  std::cout << "\t\tTotal iterations to run : " << ifo.endIter << "\n";
  std::cout << "\t\tOutput Solution after iteration : " << ifo.solOutIter << "\n";
  std::cout << "\t\tSolution Output Path : " << ifo.solOutPath << "\n";
  std::cout << "\t\tOutput Load after iteration : " << ifo.loadOutIter << "\n";
  std::cout << "\t\tLoad Output Path : " << ifo.loadOutPath << "\n";
  if( nargs < 2 )
    std::cout << "\t\tInput File Argument : input.dat" << "\n";
  else
    std::cout << "\t\tInput File Argument : " << args[1] << "\n";
  std::cout << "\t\tInput File Path : " << ifo.filePrefix << "\n";
  std::cout << "*****************************************************************\n";
  std::cout << "\t\tNo:of blades : " << blade.nBlade << "\n";
  std::cout << "\t\tOmega : " << blade.omega << "\n";
  std::cout << "\t\tZ Location : " << blade.location[2] << "\n";
  std::cout << "\t\tPsi-Zero : " << blade.psiZero << "\n";
  std::cout << "\t\tFWHM Y : " << blade.FWHM_Y << "\n";
  std::cout << "\t\tFWHM Z : " << blade.FWHM_Z << "\n";
  std::cout << "\t\tRoot Chord : " << blade.rootChord << "\n";
  std::cout << "\t\tRoot Cut Out : " << blade.rootCutOut << "\n";
  std::cout << "\t\tRotor Radius : " << blade.rotorRadius << "\n";
  std::cout << "\t\tTip Correction Length : " << ifo.tipCorLen << "\n";
  std::cout << "*****************************************************************\n"<<std::endl;
  // File parsed success !!!
  NM.Constructor(&S,"inputGrid.dat");

  return 0;
}


std::string getStrippedInputFile(const char *filename)
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


