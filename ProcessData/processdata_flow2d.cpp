#include "processdata.h"

using namespace INMOST;
using namespace std;

ProcessData_Flow2D::ProcessData_Flow2D(Mesh *m_)
	: ProcessData_Diffusion2D (m_)
{
	haveK = false;
	haveBC = false;
	haveSource = false;
	if(m->HaveTag("K")){
		cout << "K is present on the mesh" << endl;
		haveK = true;
		tagK = m->GetTag("K");
	}
	if(m->HaveTag("BOUNDARY_CONDITION")){
		cout << "BC tag is present on the mesh" << endl;
		haveBC = true;
		tagBC = m->GetTag("BOUNDARY_CONDITION");
	}
	if(m->HaveTag("SOURCE")){
		cout << "Source tag is present on the mesh" << endl;
		haveSource = true;
		tagSource = m->GetTag("SOURCE");
	}
}

ProcessData_Flow2D::~ProcessData_Flow2D()
{

}

int ProcessData_Flow2D::getDiffusionBCtype(const INMOST::Face &f)
{
	if(haveBC){
		return static_cast<int>(f.RealArray(tagBC)[1]);
	}
    return BC_FLOW_NEUM;
}

void ProcessData_Flow2D::getDiffusionBC(const INMOST::Face &f, double *res)
{
    double x[2];
    f.Barycenter(x);
	if(haveBC)
		res[0] = f.RealArray(tagBC)[2];
	else
        res[0] = 0;
}

void ProcessData_Flow2D::getDiffusionTensor(const INMOST::Cell &c, double *res)
{

	if(haveK){
        for(unsigned i = 0; i < 4; i++)
			res[i] = c.RealArray(tagK)[i];
	}
	else{
		double k = 1e0;
		res[0] = k;
		res[1] = 0.0;
		res[2] = 0.0;
        res[3] = k;
        //res[4] = k;
        //res[5] = 0.0;
        //res[6] = 0.0;
        //res[7] = 0.0;
        //res[8] = k;
	}
}

void ProcessData_Flow2D::getSourceTerm(const INMOST::Cell &c, double *res)
{
    *res = 0;
    double x[3];
    c.Barycenter(x);
	if(haveSource){
		*res = c.Real(tagSource);
		return;
    }
    //printf("cell %d: src = %e\n", c.LocalID(), *res);
}
