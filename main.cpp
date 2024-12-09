#include "header.h"
#include "Process/process.h"

// TODO:
// 1. Integrate in meshDiam()
// 2. Output extension
// 3. Compute rank
// 4. Dirichlet nodes detection in parallel

using namespace std;
using namespace INMOST;

#define pprintf(ARG) if(rank == 0) printf((ARG))
#define Cout std::cout

class Problem
{
private:

    // Mesh
    Mesh *mesh;

    // Automatizator
    Automatizator *aut;

    // Model
    Model *model;

    // Processes
    Process_Flow2D *Flow;
    Process_Advection2D *Tran;

    // Linear solvers
    Solver *solverFlow;
    Solver *solverTran;

    // Solution parameters
    double dT, T;
    int nt;
    string params_path;          // parameter file path
    string mesh_path;            // mesh path
    string save_dir;             // results folder
    double newton_flow_rtol;     // newton relative tolerance
    double newton_flow_atol;     // newton absolute tolerance
    int    newton_flow_maxit;
    double newton_tran_rtol;
    double newton_tran_atol;
    int    newton_tran_maxit;
    string solver_type_flow;     // linear solver type
    string solver_type_tran;
    double droptol_flow;
    double droptol_tran;
    double lin_atol_flow;
    double lin_atol_tran;
    double lin_rtol_flow;
    double lin_rtol_tran;
    string fv_type;              // can be "TPFA"
    int save_intensity;          // save VTK every ... step
    bool save_sol;               // if we save solution as (P)VTK
    string output_extension;     // vtk or pvtk
    unsigned total_iter_lin;

    // Auxiliary
    double times[7];
    double tglob;                // global timer
    int iter_linear;             // total number of linear iters
    Tag tagGhost;
    ofstream outC;               // concentration outflow

public:
    Problem(string param_path);
    ~Problem();
    void setDefaultParams();
    void readParams();
    void readMesh();
    void run();
    void timeStepTransport();
    void printTimes();
};

Problem::Problem(string param_path)
{
    tglob = Timer();
    params_path = param_path;
    iter_linear = 0;
    setDefaultParams();
    readParams();
    nt = static_cast<int>(T/dT);
    std::fill_n(times, 6, 0.0);

    // Init null pointers
    aut = nullptr;
    mesh = nullptr;
    model = nullptr;
    Flow = nullptr;
    Tran = nullptr;
    solverFlow = nullptr;
    solverTran = nullptr;

    mesh = new Mesh;
    output_extension = ".vtk";
    readMesh();

    double ttt = Timer();

    aut = new Automatizator();
    model = new Model(*aut);
    model->AddMesh("mesh", *mesh);
    cout << "Before create procs" << endl;
    Flow = new Process_Flow2D(mesh, fv_type);
    cout << "Created flow" << endl;
    Tran = new Process_Advection2D(mesh);
    cout << "Created procs" << endl;
    model->AddSubModel("Flow2D", *Flow);
    model->AddSubModel("Tran2D", *Tran);
    model->PrepareEntries();
    model->Initialize();
    model->SetTimeStep(dT);
    model->PrepareIterations();
    printf("IND: %d %d\n", aut->GetFirstIndex(), aut->GetLastIndex());
    fflush(stdout);

    solverFlow = new Solver(solver_type_flow);
    solverTran = new Solver(solver_type_tran);

    solverFlow->SetParameterReal("relative_tolerance", lin_rtol_flow);
    solverFlow->SetParameterReal("absolute_tolerance", lin_atol_flow);
    solverFlow->SetParameterReal("drop_tolerance", droptol_flow);

    solverTran->SetParameterReal("relative_tolerance", lin_rtol_tran);
    solverTran->SetParameterReal("absolute_tolerance", lin_atol_tran);
    solverTran->SetParameterReal("drop_tolerance", droptol_tran);

    total_iter_lin = 0;

    outC.open("outC.txt");

    times[T_INIT] += Timer() - ttt;
}

Problem::~Problem()
{
    printTimes();
    if(Flow != nullptr)
        delete Flow;
    if(Tran != nullptr)
        delete Tran;
    if(aut != nullptr)
        delete aut;
    if(model != nullptr)
        delete model;
    if(mesh != nullptr)
        delete mesh;
    if(solverFlow != nullptr)
        delete solverFlow;
    if(solverTran != nullptr)
        delete solverTran;
}

void Problem::setDefaultParams()
{
    mesh_path = "sol0.vtk";
    dT = 1.0;
    T  = dT * 4;
    save_dir = ".";
    newton_flow_rtol = 1e-5;     // for subproblems in splitting
    newton_flow_atol = 1e-5;
    newton_flow_maxit = 15;
    newton_tran_rtol = 1e-6;
    newton_tran_atol = 1e-6;
    newton_tran_maxit = 5;
    solver_type_flow = "inner_mptiluc";
    solver_type_tran = "inner_ilu2";
    droptol_flow = 1e-2;
    droptol_tran = 1e-2;
    lin_atol_flow = 1e-9;
    lin_atol_tran = 1e-9;
    lin_rtol_flow = 1e-6;
    lin_rtol_tran = 1e-6;
    fv_type = "TPFA";
    save_intensity = 1;
    save_sol = true;
}

void Problem::readParams()
{
    double t = Timer();
    std::ifstream file(params_path);
    if(file.fail()){
        cout << "Could not open param file " << params_path << endl;
        exit(1);
    }
    std::string line;
    while (getline(file, line))
    {
        std::istringstream iss(line);
        std::string firstword;
        iss >> firstword;
        if(firstword[0] == '#'){
            //std::cout << "Found a comment line" << std::endl;
            continue;
        }

        if(firstword == "dt")
            iss >> dT;
        if(firstword == "T")
            iss >> T;
        if(firstword == "mesh")
            iss >> mesh_path;
        if(firstword == "save_dir")
            iss >> save_dir;
        if(firstword == "newton_flow_rtol")
            iss >> newton_flow_rtol;
        if(firstword == "newton_flow_atol")
            iss >> newton_flow_atol;
        if(firstword == "newton_flow_maxit")
            iss >> newton_flow_maxit;
        if(firstword == "newton_tran_rtol")
            iss >> newton_tran_rtol;
        if(firstword == "newton_tran_atol")
            iss >> newton_tran_atol;
        if(firstword == "newton_tran_maxit")
            iss >> newton_tran_maxit;
        if(firstword == "solver_type_flow")
            iss >> solver_type_flow;
        if(firstword == "solver_type_tran")
            iss >>solver_type_tran;
        if(firstword == "droptol_flow")
            iss >> droptol_flow;
        if(firstword == "droptol_tran")
            iss >> droptol_tran;
        if(firstword == "lin_atol_flow")
            iss >> lin_atol_flow;
        if(firstword == "lin_atol_tran")
            iss >> lin_atol_tran;
        if(firstword == "lin_rtol_flow")
            iss >> lin_rtol_flow;
        if(firstword == "lin_rtol_tran")
            iss >> lin_rtol_tran;
        if(firstword == "fv_type")
            iss >> fv_type;
        if(firstword == "save_intensity")
            iss >> save_intensity;
        if(firstword == "save_sol")
            iss >> save_sol;
    }
    times[T_IO] += Timer() - t;
    Cout << "Finished reading params from " << params_path << endl;
}

void Problem::readMesh()
{
    Cout << "Loading mesh " << mesh_path << endl;
    double ttt = Timer(); // T_IO
    mesh->Load(mesh_path);
    Cout << "Number of cells: " << mesh->NumberOfCells() << endl;
    Cout << "Number of faces: " << mesh->NumberOfFaces() << endl;
    Cout << "Number of edges: " << mesh->NumberOfEdges() << endl;
    Cout << "Number of nodes: " << mesh->NumberOfNodes() << endl;

    times[T_IO] += Timer() - ttt;

    ttt = Timer(); // T_INIT
    Mesh::GeomParam table;
    table[MEASURE] = CELL | FACE;
    table[ORIENTATION] = FACE;
    table[NORMAL] = FACE;
    table[BARYCENTER] = CELL | FACE;
    mesh->PrepareGeometricData(table);
    for(auto iface = mesh->BeginFace(); iface != mesh->EndFace(); iface++)
        iface->FixNormalOrientation();
    mesh->ReorderEmpty(CELL|FACE|EDGE|NODE);
    mesh->AssignGlobalID(CELL|FACE|EDGE|NODE);

    table[MEASURE] = CELL | FACE;
    table[ORIENTATION] = FACE;
    table[NORMAL] = FACE;
    table[BARYCENTER] = CELL | FACE | EDGE | NODE;
    table[CENTROID] = CELL | FACE | EDGE | NODE;
    mesh->PrepareGeometricData(table);

    for(auto iface = mesh->BeginFace(); iface != mesh->EndFace(); iface++){
        iface->FixNormalOrientation();
    }
    
    // Name: GERA_K
    // Data type: double
    // Defined on: cells
    // Sparse on: nowhere
    // 4 numbers for each cell
    // [ Kx  Kxy]
    // [ Kxy Ky ]
    Tag tagK = mesh->CreateTag("K", DATA_REAL, CELL, NONE, 4);
    Tag tagR = mesh->CreateTag("Retardation", DATA_REAL, CELL, NONE, 1);
    for(auto icell = mesh->BeginCell(); icell != mesh->EndCell(); icell++){
        Cell c = icell->getAsCell();
        double K = 1.0;
        double R = 0.0;

        c.RealArray(tagK)[0] = K;
        c.RealArray(tagK)[1] = 0;
        c.RealArray(tagK)[2] = 0;
        c.RealArray(tagK)[3] = K;

        c.RealArray(tagR)[0] = R;
    }

    Tag tagBC = mesh->CreateTag("BOUNDARY_CONDITION", DATA_REAL, FACE, FACE, 3);
    for(auto iface = mesh->BeginFace(); iface != mesh->EndFace(); iface++){
        Face f = iface->getAsFace();
        if(!f.Boundary())
            continue;
        f.RealArray(tagBC)[0] = 0.0;
        f.RealArray(tagBC)[1] = 1.0;
        f.RealArray(tagBC)[2] = 0.0;

        double nf[2];
        f.UnitNormal(nf);
        // Top boundary: water head = 0
        if(fabs(nf[1]-1.0) < 1e-10){
            f.RealArray(tagBC)[0] = 1.0;
            f.RealArray(tagBC)[1] = 0.0;
            f.RealArray(tagBC)[2] = 0.0;
        }
        // Bottom boundary: water head = 100
        if(fabs(nf[1]+1.0) < 1e-10){
            f.RealArray(tagBC)[0] = 1.0;
            f.RealArray(tagBC)[1] = 0.0;
            f.RealArray(tagBC)[2] = 100.0;
        }
    }

    times[T_INIT] += Timer() - ttt;
}

void Problem::printTimes()
{
    printf("Total linear iterations: %u\n", total_iter_lin);
    printf("\n+=========================\n");
    printf("| T_solve    = %lf\n", times[T_SOLVE]);
    printf("| T_assemble = %lf\n", times[T_ASSEMBLE]);
    printf("| T_precond  = %lf\n", times[T_PRECOND]);
    printf("| T_iter     = %lf\n", times[T_ITER]);
    printf("| T_IO       = %lf\n", times[T_IO]);
    printf("| T_update   = %lf\n", times[T_UPDATE]);
    printf("| T_init     = %lf\n", times[T_INIT]);
    printf("+-------------------------\n");
    double sum = 0.0;
    for(int i = 0; i < 6; i++) sum += times[i];
    printf("| T_total    = %lf\n", Timer() - tglob);
    printf("+=========================\n");
}

void Problem::run()
{
    double trun = Timer();
    Cout << "Solving groundwater flow problem" << endl;
    Cout << "Solver: " << solver_type_tran << " (" << droptol_tran;
    Cout << ", " << lin_rtol_tran << ", " << lin_atol_tran << ")" << endl;
    // Activate Flow and deactivate Tran
    model->ActivateEntry("Water_Head");
    model->ActivateSubModel("Flow2D");
    model->DeactivateEntry("Conc");
    model->DeactivateSubModel("Tran2D");
    model->ToggleEntryState();
    // Residual contains system of eqautions
    Residual R("", aut->GetFirstIndex(), aut->GetLastIndex());
    //printf("IND: %d %d\n", aut->GetFirstIndex(), aut->GetLastIndex());
    //fflush(stdout);

    double ttt = Timer();
    model->FillResidual(R);
    times[T_ASSEMBLE] += Timer() - ttt;

    //Cout << "R.Norm = " << R.Norm() << endl;
    Cout << "R.Size = " << R.GetResidual().Size() << endl;

    Sparse::Vector v("", aut->GetFirstIndex(), aut->GetLastIndex());
    ttt = Timer();
    // Set Jacobi matrix
    solverFlow->SetMatrix(R.GetJacobian());
    times[T_PRECOND] += Timer() - ttt;

    ttt = Timer();
    bool solved = solverFlow->Solve(R.GetResidual(), v);
    times[T_ITER] += Timer() - ttt;
    total_iter_lin += solverFlow->Iterations();
    Cout << "Residual: " << solverFlow->Residual() << endl;
    if(!solved){
        Cout << "Linear solver failed: " << solverFlow->GetReason() << endl;
        exit(1);
    }
    Cout << "Lin.it: " << solverFlow->Iterations() << endl;
    Cout << "Success" << endl;
    ttt = Timer(); // T_UPDATE
    model->UpdateSolution(v, 1.0);
    times[T_UPDATE] += Timer() - ttt;

    if(save_sol)
        mesh->Save("Head" + output_extension);

    model->DeactivateSubModel("Flow2D");
    model->DeactivateEntry("Water_Head");

    model->ActivateEntry("Conc");
    model->ActivateSubModel("Tran2D");
    model->ToggleEntryState();

    double T = 0.0;
    mesh->Save("sol0" + output_extension);
    for(int itime = 1; itime < nt; itime++){
        if(itime%save_intensity == 0)
            Cout << "\n\n===== Time step " << itime << ", tbeg = " << T << " =====" << endl;
        timeStepTransport();
        if(itime%save_intensity == 0)
            mesh->Save("sol" + to_string(itime/save_intensity) + output_extension);
        T += dT;
        double outc = Tran->getOutflow();
        outC << T << " " << outc << endl;
    }

    times[T_SOLVE] += Timer() - trun;
}

void Problem::timeStepTransport()
{
    Residual Rtran("", aut->GetFirstIndex(), aut->GetLastIndex());
    double ttt = Timer();
    model->FillResidual(Rtran);
    times[T_ASSEMBLE] += Timer() - ttt;

    Sparse::Vector v("", aut->GetFirstIndex(), aut->GetLastIndex());
    ttt = Timer();
    //solverTran->SetParameter("verbosity", "3");
    solverTran->SetMatrix(Rtran.GetJacobian());
    times[T_PRECOND] += Timer() - ttt;

    ttt = Timer();
    bool solved = solverTran->Solve(Rtran.GetResidual(), v);
    times[T_ITER] += Timer() - ttt;
    total_iter_lin += solverTran->Iterations();
    //Cout << "Residual: " << solverTran->Residual() << endl;
    if(!solved){
        Cout << "Linear solver failed: " << solverTran->GetReason() << endl;
        exit(1);
    }
    //Cout << "Lin.it: " << solverTran->Iterations() << endl;
    //Cout << "Success" << endl;
    ttt = Timer(); // T_UPDATE
    model->UpdateSolution(v, 1.0);
    times[T_UPDATE] += Timer() - ttt;

    model->UpdateTimeStep();
}

int main(int argc, char *argv[])
{
    if(argc != 2){
            cout << "Usage: main <params_file>" << endl;
        return 1;
    }

    // Intialize INMOST
    Solver::Initialize(&argc, &argv);
    Mesh::Initialize(&argc, &argv);
    Partitioner::Initialize(&argc, &argv);

    {
        Problem P(argv[1]);
        //cout << "Ready to run" << endl;
        P.run();
    }

    Solver::Finalize();
    Mesh::Finalize();
    Partitioner::Finalize();
}
