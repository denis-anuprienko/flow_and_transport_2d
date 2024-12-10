#include "process.h"

using namespace INMOST;
using namespace std;

Process_Flow2D::Process_Flow2D(Mesh *m_, string fv_name)

{
    if(m_ == nullptr){
        exit(2);
    }
    //m = m_;
    cout << "Creating PD" << endl;
    pd = new ProcessData_Flow2D(m_);
	Cout << "FV type: " << fv_name << endl;
	if(fv_name == "TPFA")
		fv = new FV_Diffusion2D_TPFA(pd);
    else {
        cout << "Unidentified FV name!" << endl;
        exit(3);
    }
	fv->build();
    meshName = "mesh";
}

Process_Flow2D::~Process_Flow2D()
{
	delete pd;
	delete fv;
}

bool Process_Flow2D::PrepareEntries(Model & model)
{
    printf("PREPARE\n");
    m = model.GetMesh(meshName);
	rank = m->GetProcessorRank();
	if(m->HaveTag("Water_Head"))
		tagH = m->GetTag("Water_Head");
	else
		tagH = m->CreateTag("Water_Head", DATA_REAL, CELL, NONE, 1);
    entryH.SetTag(tagH);
    entryH.SetElementType(CELL);
    model.AddEntry("Water_Head", entryH);
    auto names = model.GetEntriesNames();
    varH = dynamic_variable(&entryH);
    tagFlux = m->CreateTag("Flux", DATA_REAL, CELL, NONE, 3);
    tagFluxF = m->CreateTag("FluxF", DATA_REAL, FACE, NONE, 1);
    return true;
}

bool Process_Flow2D::Initialize(Model & model)
{
    (void) model;
    for(auto icell = m->BeginCell(); icell != m->EndCell(); icell++){
        Cell c = icell->getAsCell();
        double x[2], h;
        c.Barycenter(x);
    }
	Cout << "Initialized: Flow" << endl;
	return true;
}

bool Process_Flow2D::PrepareIterations()
{
    return true;
}

bool Process_Flow2D::FillResidual(Residual &R) const
{
    for(auto iface = m->BeginFace(); iface != m->EndFace(); iface++){
        Face f = iface->getAsFace();
        Cell cp = f.BackCell(), cm = f.FrontCell();
        variable q = -1. * fv->getDgrad(f, varH);
        R[varH.Index(cp)] += q ;
		if(cm.isValid())
            R[varH.Index(cm)] -= q ;
    }
    for(auto icell = m->BeginCell(); icell != m->EndCell(); icell++){
        Cell c = icell->getAsCell();
        double src;
        pd->getSourceTerm(c, &src);
        R[varH.Index(c)] -= src * c.Volume();
        //printf("Cell %d add: %e\n", c.LocalID(), src*c.Volume());
	}
    return true;
}

bool Process_Flow2D::UpdateSolution(const Sparse::Vector & sol, double alpha)
{
	(void) alpha;
    for(auto icell = m->BeginCell(); icell != m->EndCell(); icell++){
        Cell c = icell->getAsCell();
		entryH.Value(c) -= sol[entryH.Index(c)];
    }

    computeFlux();
	m->ExchangeData(tagH, CELL);
    return true;
}

bool Process_Flow2D::UpdateTimeStep()
{
    return true;
}

bool Process_Flow2D::SetTimeStep(double dt_)
{
    dt = dt_;
    return true;
}

bool Process_Flow2D::RestoreTimeStep()
{
    return true;
}

double Process_Flow2D::UpdateMultiplier(const Sparse::Vector & sol) const
{
	(void) sol;
    return 1.0;
}

double Process_Flow2D::AdjustTimeStep(double dt) const
{
    return dt;
}

bool Process_Flow2D::SetTime(double t_)
{
    t = t_;
    pd->SetTime(t_);
    return true;
}

void Process_Flow2D::computeFlux()
{
	for(auto icell = m->BeginCell(); icell != m->EndCell(); icell++){
		auto faces = icell->getFaces();
		unsigned m = static_cast<unsigned>(faces.size());
        double cflux[2] = {0., 0.};
		for(unsigned k = 0; k < m; k++){
			double flux = (-1 * fv->getDgrad(faces[k],varH)).GetValue() / faces[k].Area();
            faces[k].Real(tagFluxF) = flux;
            double nor[2];
			faces[k].UnitNormal(nor);
            for(int i = 0; i < 2; i++)
				cflux[i] += nor[i] * flux;
		}
        for(unsigned i = 0; i < 2; i++)
			icell->RealArray(tagFlux)[i] = cflux[i]/m;
        icell->RealArray(tagFlux)[2] = 0.0;
	}
	m->ExchangeData(tagFlux, CELL);
}
