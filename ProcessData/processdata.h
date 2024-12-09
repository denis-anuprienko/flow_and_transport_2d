#ifndef PROCESSDATA_H
#define PROCESSDATA_H

#include "inmost.h"
#include "../header.h"

#define pprintf(ARG) if(rank == 0) printf((ARG))

class ProcessData
{
protected:
    INMOST::Mesh *m;
	int rank;
    double t; // current time
public:
    ProcessData(INMOST::Mesh *m_);
    ~ProcessData();
    INMOST::Mesh *getMesh();
    void SetTime(const double t_) { t = t_; }
};

class ProcessData_Diffusion2D : public ProcessData
{
private:
public:
	ProcessData_Diffusion2D(INMOST::Mesh *m_) : ProcessData (m_) { (void) m; }
	~ProcessData_Diffusion2D(){}
	virtual int getDiffusionBCtype(const INMOST::Face &f) = 0;
	virtual void getDiffusionBC(const INMOST::Face &f, double *res) = 0;
	virtual void getDiffusionTensor(const INMOST::Cell &c, double *res) = 0;
	virtual void getSourceTerm(const INMOST::Cell &c, double *res) = 0;
};

class ProcessData_Flow2D : public ProcessData_Diffusion2D
{
private:
	INMOST::Tag tagK;
	INMOST::Tag tagBC;
	INMOST::Tag tagSource;
	bool haveK;
	bool haveBC;
	bool haveSource;
public:
    ProcessData_Flow2D(INMOST::Mesh *m_);
	~ProcessData_Flow2D();
	int getDiffusionBCtype(const INMOST::Face &f);
    void getDiffusionBC(const INMOST::Face &f, double *res);
    void getDiffusionTensor(const INMOST::Cell &c, double *res);
    void getSourceTerm(const INMOST::Cell &c, double *res);
};

class ProcessData_Advection2D : public ProcessData
{
private:
    INMOST::Tag tagSourceC;
    const std::string tagnameFlux = "FluxF";
    INMOST::Tag tagFlux;
	INMOST::Tag tagR;
    bool gotFlux;
	bool gotR;
public:
    ProcessData_Advection2D(INMOST::Mesh *m_);
    ~ProcessData_Advection2D();
	// Get flux from tag "Flux"
    double getFlux(const INMOST::Face &f);
	double getRetardation(const INMOST::Cell &c);
    void getSourceTerm(const INMOST::Cell &c, double *res);
};

#endif
