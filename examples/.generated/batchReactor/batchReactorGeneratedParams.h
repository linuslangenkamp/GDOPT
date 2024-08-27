//defines

#define INIT_VARS InitVars::CONST
#define RADAU_INTEGRATOR IntegratorSteps::Steps3
#define INTERVALS 25
#define FINAL_TIME 1
#define LINEAR_SOLVER LinearSolver::MA57
#define MESH_ALGORITHM MeshAlgorithm::L2_BOUNDARY_NORM
#define MESH_ITERATIONS 10
#define TOLERANCE 1e-14
#define EXPORT_OPTIMUM_PATH "/tmp"
#define EXPORT_JACOBIAN_PATH "/tmp"

// values for runtime parameters
#define REACT_SPEED_VALUE 1
#define DEPLETION_SPEED_VALUE 0.3
