#ifndef PROCESS_H
#define PROCESS_H

#include "inmost.h"
#include "../ProcessData/processdata.h"
#include "../Discretization/discretization.h"

#define pprintf(ARG) if(rank == 0) printf((ARG))
#define Cout if(rank == 0) std::cout

using namespace INMOST;

// Base class for physical processes
// Main purpose is to fill residual, except for time derivatives
//
//class Process
//{
//protected:
//    // Pointer to the mesh
//    Mesh *m;
//    /// Current time
//    double t;
//    /// Time step size
//    double dt;

//public:
//    Process(Mesh *m_);
//    virtual ~Process();
//    virtual void prepareVars(Automatizator &aut) = 0;
//    virtual void fillResidual(Residual &, dynamic_variable &) = 0;
//    void SetTime(double t_) { t = t_; }
//    void SetTimeStep(double dt_) { dt = dt_; }
//};

// Confined flow in 2D
//
class Process_Flow2D : public AbstractSubModel
{
protected:
    /// Mesh name
    std::string meshName;
    /// Mesh
    Mesh *m;
    /// Water head entry
    SingleEntry entryH;
    /// Water head tag
    Tag tagH;
    /// Flux tag on cells
    Tag tagFlux;
    /// Normal flux tag on faces
    Tag tagFluxF;
    /// Water head dynamic_variable
    dynamic_variable varH;
    /// PD object that provides discr. with tensor, BCs, etc.
    ProcessData_Flow2D *pd;
    /// TPFA finite volume discretization
    FV_Diffusion2D *fv;
    /// Time step size
    double dt;
    /// Current time
    double t;
	/// Rank (for output)
	int rank;
public:
    Process_Flow2D(){}
    Process_Flow2D(Mesh *m_, std::string fv_name = "TPFA");
    virtual ~Process_Flow2D();
    bool PrepareEntries(Model & model);
    bool Initialize(Model & model);
    bool PrepareIterations();
    bool FillResidual(Residual & R) const;
    bool UpdateSolution(const Sparse::Vector & sol, double alpha);
    bool UpdateTimeStep();
    bool SetTimeStep(double dt_);
    bool SetTime(double t_);
    bool RestoreTimeStep();
    double UpdateMultiplier(const Sparse::Vector & sol) const;
    double AdjustTimeStep(double dt) const;
	void computeFlux();
};

class Process_Advection2D : public AbstractSubModel
{

protected:
    /// Mesh name
    std::string meshName;
    /// Mesh
    Mesh *m;
    /// Concentration entry
    SingleEntry entryC;
    /// Concentration tag
    Tag tagC;
    /// Concentration tag on previous step
    Tag tagCold;
    /// Concentration dynamic_variable
    dynamic_variable varC;
    /// PD object that provides discr. with tensor, BCs, etc.
    ProcessData_Advection2D *pd;
    /// TPFA finite volume discretization
    FV_Diffusion2D *fv;
    /// Time step size
    double dt;
    /// Current time
    double t;
    /// Rank (for output)
    int rank;
public:
    Process_Advection2D(){}
    Process_Advection2D(Mesh *m_);
    virtual ~Process_Advection2D();
    bool PrepareEntries(Model & model);
    bool Initialize(Model & model);
    bool PrepareIterations();
    bool FillResidual(Residual & R) const;
    bool UpdateSolution(const Sparse::Vector & sol, double alpha);
    bool UpdateTimeStep();
    bool SetTimeStep(double dt_);
    bool RestoreTimeStep();
    double UpdateMultiplier(const Sparse::Vector & sol) const;
    double AdjustTimeStep(double dt) const;
    double getOutflow();
};

#endif
