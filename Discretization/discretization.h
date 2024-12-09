#ifndef DISCRETIZATION_H
#define DISCRETIZATION_H

#include "inmost.h"
#include "../ProcessData/processdata.h"

#define pprintf(ARG) if(rank == 0) printf((ARG))
#define Cout if(rank == 0) std::cout

// Base class for discretization methods
//
class Discretization
{
protected:
    INMOST::Mesh *m;
	int rank;

public:
    Discretization(ProcessData *pd_);
    ~Discretization();
    virtual void build() = 0;
};

class FV_Diffusion2D : public Discretization
{
private:
public:
    FV_Diffusion2D(ProcessData *pd_);
    virtual ~FV_Diffusion2D();
    virtual INMOST::variable getDgrad(const INMOST::Face &f,
                                      INMOST::dynamic_variable var) = 0;
	virtual void rebuildBCs(void) = 0;
};

class FV_Diffusion2D_TPFA : public FV_Diffusion2D
{
private:
    INMOST::Tag tagT; // TPFA coefficient tag
	ProcessData_Diffusion2D *pd;
public:
	FV_Diffusion2D_TPFA(ProcessData_Diffusion2D *pd_);
    ~FV_Diffusion2D_TPFA();
	void build(void);
    INMOST::variable getDgrad(const INMOST::Face &f,
                              INMOST::dynamic_variable var);
	void rebuildBCs(void);
};


#endif
