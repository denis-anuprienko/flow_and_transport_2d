#include "process.h"

using namespace INMOST;
using namespace std;

Process_Advection2D::Process_Advection2D(Mesh *m_)
{
    m = m_;
    pd = new ProcessData_Advection2D(m_);
    meshName = "mesh";
}

Process_Advection2D::~Process_Advection2D()
{
    if(pd != nullptr)
       ;//delete pd;
}

bool Process_Advection2D::PrepareEntries(Model &model)
{
    m = model.GetMesh(meshName);
    rank = m->GetProcessorRank();
    if(m->HaveTag("Conc"))
        tagC = m->GetTag("Conc");
    else
        tagC = m->CreateTag("Conc", DATA_REAL, CELL, NONE, 1);
    tagCold = m->CreateTag("Conc_Old", DATA_REAL, CELL, NONE, 1);
    entryC.SetTag(tagC);
    entryC.SetElementType(CELL);
    model.AddEntry("Conc", entryC);
    varC = dynamic_variable(&entryC);
    return true;
}

bool Process_Advection2D::Initialize(Model &model)
{
    Tag tagId = m->GetTag("CellEntityIds");
    if(!tagId.isValid()){
        Cout << "Something wrong with tag CellEntityIds" << endl;
        exit(98);
    }
    for(auto icell = m->BeginCell(); icell != m->EndCell(); icell++){
        Cell c = icell->getAsCell();
        c.Real(tagC) = 0.0;

        if(c.Integer(tagId) == 8)
            c.Real(tagC) = 1;

        c.Real(tagCold) = c.Real(tagC);
    }
    return true;
}

bool Process_Advection2D::PrepareIterations()
{
    return true;
}

bool Process_Advection2D::FillResidual(Residual &R) const
{
    for(auto icell = m->BeginCell(); icell != m->EndCell(); icell++){
        Cell cell = icell->getAsCell();
        double V = cell.Volume();
        double lambda = 1;
        R[varC.Index(cell)] -=
            (
                (varC(cell) - cell.Real(tagCold)) / dt
            ) * V;

        auto faces = cell.getFaces();
        for(auto iface = faces.begin(); iface != faces.end(); iface++){
            Face f = iface->getAsFace();
            variable flux = pd->getFlux(f);
            Cell cp = f.BackCell(), cm = f.FrontCell();
            if(cm.isValid()){
                if(flux.GetValue() > 0.)
                    flux *= varC(cp);
                else
                    flux *= varC(cm);
                if(cell == cm)
                    flux *= -1;
            }
            else
                flux *= varC(cp);
            R[varC.Index(cell)] -= flux;
        }
    }
    return true;
}

bool Process_Advection2D::UpdateSolution(const Sparse::Vector &sol, double alpha)
{
    for(auto icell = m->BeginCell(); icell != m->EndCell(); icell++){
        Cell c = icell->getAsCell();
        c.Real(tagC) -= sol[varC.Index(c)];
    }
    return true;
}

bool Process_Advection2D::UpdateTimeStep()
{
    for(auto icell = m->BeginCell(); icell != m->EndCell(); icell++){
        Cell c = icell->getAsCell();
        c.Real(tagCold) = c.Real(tagC);
    }
    return true;
}

bool Process_Advection2D::SetTimeStep(double dt_)
{
    dt = dt_;
    return true;
}

bool Process_Advection2D::RestoreTimeStep()
{
    return true;
}

double Process_Advection2D::UpdateMultiplier(const Sparse::Vector &sol) const
{
    return 1.0;
}

double Process_Advection2D::AdjustTimeStep(double dt) const
{
    return dt;
}

double Process_Advection2D::getOutflow()
{
    double outflow = 0.0;
    for(auto iface = m->BeginFace(); iface != m->EndFace(); iface++){
        Face f = iface->getAsFace();
        if(!f.Boundary())
            continue;
        double nf[2];
        f.UnitNormal(nf);
        if(fabs(nf[1] - 1.0) > 1e-10)
            continue;
        // Upper face
        outflow += pd->getFlux(f) * f.BackCell().Real(tagC);
    }
    return outflow;
}
