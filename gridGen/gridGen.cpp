#include "definition.h"
#include "StretchedCylindricalGrid.h"

std::string getStrippedInputFile(const char *filename); 

int main(int nargs, char **args) {
	std::string dd;
	if(nargs < 2)
		dd=getStrippedInputFile("input.dat");
	else
		dd=getStrippedInputFile(args[1]);

	std::istringstream inputFile(dd);
	int imax,jmax,kmax;
	double zoffset;
	std::string pathPrefix;
	std::stringstream dummy;
	tvmet::Vector<double,3> Translate;
	tvmet::Vector<double,3> Rotate;
	blitz::Array<double,1> rdist;
	blitz::Array<double,1> zdist;
	StretchedCylindricalGrid G;
	bool debug=true;

	inputFile >> pathPrefix >> jmax >> imax >> kmax;
	rdist.resize(imax);
	zdist.resize(kmax);
	for(int i = 0 ; i < imax ; ++i ){
		inputFile >> rdist(i);
	}
	for(int i = 0 ; i < kmax ; ++i ){
		inputFile >> zdist(i);
	}
	inputFile >> Translate[0] >> Translate[1] >> Translate[2];
	inputFile >> Rotate[0] >> Rotate[1] >> Rotate[2];

	std::cout << std::endl;
	std::cout << "Stripped input file " << std::endl;
	std::cout << inputFile.str() << std::endl;

	std::cout << "Read from input file" << std::endl;
	std::cout << "--------------------" << std::endl;
	std::cout << "imax = " << imax << std::endl;
	std::cout << "jmax = " << jmax << std::endl;
	std::cout << "kmax = " << kmax << std::endl;
	std::cout << "r distribution = ";
	for(int i = 0 ; i < imax ; ++i )
		std::cout << rdist(i) << "\\" ;
	std::cout << std::endl;
	std::cout << "z distribution = ";
	for(int i = 0 ; i < kmax ; ++i )
		std::cout << zdist(i) << "\\";
	std::cout << std::endl;
	G.Constructor(jmax,rdist,zdist);

	std::cout << "Writing block grid to file" << std::endl;
	std::cout << "--------------------------" << std::endl;

	pathPrefix = dummy.str() + "CylGrid.dat";
	G.writeUnstructredGridForSolver(pathPrefix);

	pathPrefix = dummy.str() + "TecplotGrid.dat";
	G.writeUnstructuredData(pathPrefix);

        if(debug){
		pathPrefix = dummy.str() + "Centorid.dat";
		G.writeDebugCentroid(pathPrefix);
	}
	
	// Higher Order Integration
	pathPrefix = dummy.str() + "QuadratureNodes.dat";
	G.outputHigherOrderNodes(pathPrefix);
	return 0;
}

std::string getStrippedInputFile(const char *filename) {
	std::ifstream namelist(filename);
	if(namelist.fail()){
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
