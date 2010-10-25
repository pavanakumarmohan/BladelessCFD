# -----------------------------------------------------------
# The folder for the blitz and tvmet libraries
# someday I want to remove this and wrap these two
# into the code
# -----------------------------------------------------------
BZ_DIR=/home/e495299/codes/
TVMET_DIR=/home/e495299/codes/
HEXA8P3=/home/e495299/codes/gridGen/trunk/
OMP_INCL=/usr/lib64/openmpi/1.4-gcc/include/
GSL_LIB=/home/e495299/codes/
# -----------------------------------------------------------
# Target
# -----------------------------------------------------------
SOLVER=StretchedStaticCylindricalGridSolver
MAIN=StretchedSolver
TARGET=stretchedSolve
HEXA=../../../gridGen/trunk/Hexa8p3
MOM_SRC=NewMomentumSource
DEFINITION=/home/e495299/codes/gridGen/trunk/stretchedGrid/
# -----------------------------------------------------------
# Compiler Options
# -----------------------------------------------------------
CXX=g++
OPT_FLAGS=-Wno-deprecated -O3 
DEBUG_FLAGS=-g -DBZ_DEBUG -DTVMET_DEBUG -Wno-deprecated
LINK=-L/usr/lib64/openmpi/1.4-gcc/lib/ -lm -lblitz -lgsl -lgslcblas -lmpi
# -----------------------------------------------------------
# Optimization Flags
# -----------------------------------------------------------
#FLAGS=$(DEBUG_FLAGS)
#FLAGS=$(DEBUG_FLAGS) -DUSE_MPI -fopenmp -pg
FLAGS=-ggdb -O3 -DUSE_MPI -pg -fopenmp -I$(DEFINITION) -I$(OMP_INCL)
# -----------------------------------------------------------
# Dependency and compilation 
# -----------------------------------------------------------
#$(TARGET): $(MAIN).o $(SOLVER).o $(DET_DT) $(TEST_EX) $(HEXA).cpp $(HEXA).h $(MOM_SRC).o
$(TARGET): $(MAIN).o $(SOLVER).o $(DET_DT) $(TEST_EX) $(HEXA).cpp $(HEXA).h 
	$(CXX) $(FLAGS) -L./lib/ -L$(GSL_LIB)/lib -L$(BZ_DIR)/lib -L$(TVMET_DIR)/lib -o $(TARGET) $(MOM_SRC).o $(SOLVER).o $(MAIN).o $(HEXA).o $(LINK)
$(MAIN).o: $(MAIN).cpp $(SOLVER).o $(MOM_SRC).o 
	$(CXX) -c $(FLAGS) -I$(HEXA8P3) -I$(GSL_LIB)/include -I$(BZ_DIR)/include -I$(TVMET_DIR)/include -c $(MAIN).cpp
$(SOLVER).o: $(SOLVER).cpp $(SOLVER).h
	$(CXX) -c $(FLAGS) -I$(HEXA8P3) -I$(GSL_LIB)/include -I$(BZ_DIR)/include -I$(TVMET_DIR)/include  -c $(SOLVER).cpp
$(MOM_SRC).o: $(MOM_SRC).cpp
	$(CXX) -c $(FLAGS) -I$(HEXA8P3) -I$(GSL_LIB)/include -I$(BZ_DIR)/include -I$(TVMET_DIR)/include  -c $(MOM_SRC).cpp
all: $(TARGET) 
clean:
	rm *.o

