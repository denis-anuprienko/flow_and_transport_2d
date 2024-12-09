#ifndef HEADER_H
#define HEADER_H

#include "inmost.h"
#include <iostream>

//#include "Process/process.h"
//#include "ProcessData/processdata.h"
//#include "Discretization/discretization.h"
//#include "Utils/utils.h"

enum{
    T_ASSEMBLE = 0, // assembly of residual
    T_PRECOND,      // preconditioner computation
    T_ITER,         // Bi-CGSTAB iterations
    T_SOLVE,        // solution time not including init (Problem::run())
    T_IO,           // reading mesh and saving result
    T_INIT,         // every preparatory move in Problem, including discretization building
    T_UPDATE,       // update of solution, may include flux and stress calculation
};

enum{
	BC_FLOW_DIR = 0,
	BC_FLOW_NEUM
};

enum{
	BC_MECH_DIR = 0,
	BC_MECH_NEUM,
	BC_MECH_ROLLER
};

#endif
