//defines

#define INIT_VARS InitVars::SOLVE
#define RADAU_INTEGRATOR IntegratorSteps::Steps1
#define INTERVALS 5000
#define FINAL_TIME 8
#define LINEAR_SOLVER LinearSolver::MA57
#define MESH_ALGORITHM MeshAlgorithm::L2_BOUNDARY_NORM
#define MESH_ITERATIONS 0
#define TOLERANCE 1e-14
#define EXPORT_OPTIMUM_PATH "/tmp"
#define INITIAL_STATES_PATH "/tmp"

// values for runtime parameters
#define MS_VALUE 1.5
#define MP_VALUE 0.5
#define R_VALUE 1
