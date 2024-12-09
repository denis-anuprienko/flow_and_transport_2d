#include "discretization.h"

using namespace INMOST;
using namespace std;

FV_Diffusion2D_TPFA::FV_Diffusion2D_TPFA(ProcessData_Diffusion2D *pd_)
    : FV_Diffusion2D (pd_)
{
    pd = pd_;
    m = pd->getMesh();
}

FV_Diffusion2D_TPFA::~FV_Diffusion2D_TPFA()
{

}

void FV_Diffusion2D_TPFA::build()
{
    tagT = m->CreateTag("TPFA_trans", DATA_REAL, FACE, NONE, 1);
    for(auto iface = m->BeginFace(); iface != m->EndFace(); iface++){
		Face f = iface->getAsFace();
		// Don't build for Neumann
		if(f.Boundary() && pd->getDiffusionBCtype(f) == 1)
			continue;
//		printf("Building TPFA for face %d\n", f.DataLocalID());
//		fflush(stdout);
        double xf[2];
		f.Barycenter(xf);

		if(f.Boundary()){
			if(pd->getDiffusionBCtype(f) == -1){
				Cout << "Undefined BC type for face!" << endl;
                exit(4);
			}
			// Here 'p' and 'm' refer to '+' and '-'
			Cell cp = f.BackCell();
            double xp[2];
			cp.Barycenter(xp);

            rMatrix Dp(2,2), ne(2,1), lp(2,1);
			// initialize diffusion tensors
			pd->getDiffusionTensor(cp, Dp.data());
			//Dp.Print();

			// Get unit normal for face
			f.UnitNormal(ne.data());

			// Compute l's
            for(unsigned k = 0; k < 2; k++)
				lp(k,0) = xf[k] - xp[k];
			double maglp = lp.FrobeniusNorm();
			lp /= (maglp * maglp);
			//lp /= (xf[0] - xp[0])*(xf[0] - xp[0]) + (xf[1] - xp[1])*(xf[1] - xp[1]) + (xf[2] - xp[2])*(xf[2] - xp[2]);

			double coef = (Dp*lp).DotProduct(ne);
			coef *= f.Area();
            f.Real(tagT) = 1;//shit coef;
			//printf("b.face %d (%lf %lf %lf): T = %e\n", f.LocalID(), xf[0],xf[1],xf[2], coef);
		}
		else{ // internal face
			// Here 'p' and 'm' refer to '+' and '-'
			Cell cp = f.BackCell(), cm = f.FrontCell();
            double xp[2], xm[2];
			cp.Barycenter(xp);
			cm.Barycenter(xm);

            rMatrix Dp(2,2), Dm(2,2), ne(2,1), lp(2,1), lm(2,1);

			// initialize diffusion tensors
			pd->getDiffusionTensor(cp, Dp.data());
			pd->getDiffusionTensor(cm, Dm.data());

            if(f.DataLocalID() == 10){
                Dp.Print();
                Dm.Print();
            }
			// Get unit normal for face
			f.UnitNormal(ne.data());

			// Compute l's
            for(unsigned k = 0; k < 2; k++){
				lp(k,0) = xf[k] - xp[k];
				lm(k,0) = xf[k] - xm[k];
			}
			double lpmag = lp.FrobeniusNorm();
			double lmmag = lm.FrobeniusNorm();
			lp /= (lpmag * lpmag);
			lm /= (lmmag * lmmag);

			double coef = (Dp*lp).DotProduct(ne) * (Dm*lm).DotProduct(ne);
			coef /= ((Dp*lp).DotProduct(ne) - (Dm*lm).DotProduct(ne));
			coef *= f.Area();
            f.Real(tagT) = -coef;
            //printf("face %d: T = %e\n", f.LocalID(), coef);
		}
    }
	m->ExchangeData(tagT, FACE);
	Cout << "TPFA is built" << endl;
}

void FV_Diffusion2D_TPFA::rebuildBCs()
{

}

variable FV_Diffusion2D_TPFA::getDgrad(const Face & f, dynamic_variable var)
{
    variable res;
	if(f.Boundary()){
		// Neumann case
		if(pd->getDiffusionBCtype(f) == 1){
			double val;
			// val = -K*grad(h), K*grad(h) = -val
			pd->getDiffusionBC(f, &val);
			return -1*val;
		}

		// Dirichlet
        Cell cp = f.BackCell();
        double BCval;
        pd->getDiffusionBC(f, &BCval);
        //printf("face %d, BC = %lf\n", f.DataLocalID(), BCval);
		res = f.Real(tagT) * (BCval - var(cp));

        double xf[2];
        f.Barycenter(xf);
//        if(xf[2] > 1-1e-10 || xf[2] < 1e-10)
//            res = 0;
//        printf("Face %d: flux = %e\n", f.LocalID(), res.GetValue());
    }
    else{
        Cell cp = f.BackCell(), cm = f.FrontCell();
        res = f.Real(tagT) * (var(cm) - var(cp));
    }
    return res;
}
