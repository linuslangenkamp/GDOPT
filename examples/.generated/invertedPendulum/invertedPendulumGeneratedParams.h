//defines

#define INIT_VARS InitVars::CONST
#define RADAU_INTEGRATOR IntegratorSteps::Steps1
#define INTERVALS 1000
#define FINAL_TIME 12
#define LINEAR_SOLVER LinearSolver::MUMPS
#define MESH_ALGORITHM MeshAlgorithm::L2_BOUNDARY_NORM
#define MESH_ITERATIONS 0
#define TOLERANCE 1e-14
#define EXPORT_OPTIMUM_PATH "/tmp"

// values for runtime parameters
#define PARAMETER_MS_VALUE 1.5
#define PARAMETER_MP_VALUE 0.5
#define PARAMETER_R_VALUE 1
