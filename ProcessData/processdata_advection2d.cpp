#include "processdata.h"

using namespace INMOST;
using namespace std;

ProcessData_Advection2D::ProcessData_Advection2D(Mesh *m_)
    : ProcessData (m_)
{
    gotFlux = false;
}

ProcessData_Advection2D::~ProcessData_Advection2D()
{

}


void ProcessData_Advection2D::getSourceTerm(const INMOST::Cell &c, double *res)
{
    double x[3];
    c.Barycenter(x);
    *res = 0;
}

double ProcessData_Advection2D::getFlux(const INMOST::Face &f)
{
    if(!gotFlux){
        if(!m->HaveTag(tagnameFlux)){
            cout << "No flux tag!" << endl;
            exit(97);
        }
        tagFlux = m->GetTag(tagnameFlux);
        gotFlux = true;
    }
    return f.Real(tagFlux);
}

double ProcessData_Advection2D::getRetardation(const INMOST::Cell &c)
{
    if (!gotR) {
        if (!m->HaveTag("Retardation")) {
            cout << "No alpha tag!" << endl;
            exit(97);
        }
        tagR = m->GetTag("Retardation");
        gotR = true;
    }
    return c.Real(tagR);
}


