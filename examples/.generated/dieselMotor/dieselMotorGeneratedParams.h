//defines

#define INIT_VARS InitVars::SOLVE
#define RADAU_INTEGRATOR IntegratorSteps::Steps3
#define INTERVALS 100
#define FINAL_TIME 0.5
#define LINEAR_SOLVER LinearSolver::MUMPS
#define MESH_ALGORITHM MeshAlgorithm::L2_BOUNDARY_NORM
#define MESH_ITERATIONS 5
#define TOLERANCE 1e-14
#define EXPORT_OPTIMUM_PATH "/tmp"
#define INITIAL_STATES_PATH "/tmp"
