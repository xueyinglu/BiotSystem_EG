#ifndef BIOT_SYSTEM_H_
#define BIOT_SYSTEM_H_
#include "DealiiHeader.h"
#include "ParameterReader.h"
#include "InitialPressure.h"
#include "RightHandSide.h"
#include "PermFunction.h"
#include "LambdaFunction.h"
#include "MuFunction.h"
namespace LA
{
using namespace dealii::LinearAlgebraTrilinos;
}
using namespace dealii;
using namespace std;
class BiotSystem
{
    /* Biot system with EG for flow and CG for mechanics*/
public:
    BiotSystem(ParameterHandler &prm);
    // virtual BiotSystem();
    void run_fixed_stress();
    void check_disp_solver_convergence();

    enum TestCase{
        none, 
        benchmark,
        mandel,
        terzaghi,
        heterogeneous
    };


private:
    MPI_Comm mpi_com;
    ParameterHandler &prm;
    double del_t = 0.01;
    
    double T0 =0.0;     // initial time
    double T = 1;       // final time
    double t = 0;
    int timestep = 0;
    int num_global_refinement = 5;
    double h = 1./pow(2,num_global_refinement);
    bool adaptivity = false;
    // parallel::distributed::Triangulation<dim> triangulation;
    Triangulation<dim> triangulation;
    /* EG pressure solution */
    FESystem<dim> fe_pressure;

    DoFHandler<dim> dof_handler_pressure;
    BlockSparsityPattern sparsity_pattern_pressure;
    ConstraintMatrix constraints_pressure;

    LA::MPI::BlockSparseMatrix system_matrix_pressure;
    LA::MPI::BlockVector solution_pressure, initial_pressure, prev_timestep_sol_pressure, prev_fs_sol_pressure;
    LA::MPI::BlockVector system_rhs_pressure;
    std::vector<IndexSet> partition_pressure;
    std::vector<IndexSet> partition_relevant_pressure;
    IndexSet relevant_set_pressure;
    LA::MPI::PreconditionAMG preconditioner_pressure_cg;
    LA::MPI::PreconditionAMG preconditioner_pressure_dg;

    // displacement solution
    FESystem<dim> fe_displacement;
    ConstraintMatrix constraints_displacement;
    DoFHandler<dim> dof_handler_displacement;

    SparsityPattern sparsity_pattern_displacement;
    SparseMatrix<double> system_matrix_displacement;

    Vector<double> solution_displacement, prev_timestep_sol_displacement, prev_fs_sol_displacement;
    Vector<double> initial_displacement;
    Vector<double> system_rhs_displacement;
    /*
    LA::MPI::BlockVector solution_displacement, prev_timestep_sol_displacement, prev_fs_sol_displacement;
    LA::MPI::BlockVector initial_displacement;
    LA::MPI::BlockVector system_rhs_displacement;
    std::vector<IndexSet> partition_displacement;
    std::vector<IndexSet> partition_relevant_displacement;
    IndexSet relevant_set_displacement;
    */
    ConvergenceTable convergence_table;
    vector<double> l2_error_p;
    vector<double> l2_error_u;
    vector<double> energy_error_u;
    vector<double> h_error_p_sq;
    vector<int> num_fs;

    // Data
    bool b_p_mult; //flag for permeability multiplier
    double mu_f = 1; // fluid viscosity
    RightHandSide right_hand_side; // mechanics equation body force
    ConstantFunction<dim> permeability;
    ConstantFunction<dim> lambda, mu;
    PermFunction perm_function;
    LambdaFunction lambda_function;
    MuFunction mu_function;
    double pressure_dirichlet_bc;
    double initial_pressure_value;
    Tensor<1,dim> traction_bc;
    TestCase test_case= benchmark;
    // coupling
    int criteria = 3; // 1: change in mean stress; 2: change in relative mean stress; 3: a posteriori
    double biot_alpha = 0.75;
    double K_b = 7./12; //K_b = lambda +2/3*mu
    double biot_inv_M = 3./28;
    double tol_fixed_stress = 1e-5;

    // EG for flow
    int degree = 1;
    double gamma_penal;
    double min_cell_diameter;
    bool bCG_WeaklyBD = true;
    double d_SForm = 0; //SIPG

    /* element-wise a posteriori error indicators*/
    DoFHandler<dim> dof_handler_output;
    FESystem<dim> fe_output;
    Vector<double> cell_eta_p;
    Vector<double> cell_eta_u;
    Vector<double> cell_stress_xx;
    Vector<double> cell_stress_xy;
    Vector<double> cell_stress_yy;
    Vector<double> cell_strain_xx;
    Vector<double> cell_strain_xy;
    Vector<double> cell_strain_yy;
    Vector<double> cell_vstrain;

    /* global a posteriori error estimators (recorded for each time step) */
    vector<double> eta_fs;
    vector<double> eta_alg;
    vector<double> eta_time;
    vector<double> eta_flow;
    vector<double> eta_jump;
    vector<double> eta_pen;
    vector<double> eta_partial_p_J;
    vector<double> eta_p_J;

    vector<double> eta_face_partial_sigma_n; // the errors on the tensor's time derivative for time step n
    vector<double> eta_N_partial_sigma_n; // the errors on the tensor's time derivative on the traction boundary for time step n
    vector<double> eta_face_partial_sigma; // the (cumulative in time ) errors on the tensor's time derivative_form
    vector<double> eta_partial_u_n; // the errors on the displacement's time derivative for time step n
    vector<double> eta_partial_u; // the (cumulative in time) errors on the displacement's time derivative
    vector<double> eta_face_sigma_n; //the errors on the tensor at time t_n;
    vector<double> eta_N_sigma_n; // the errors on the tensor at the traction bc at time step n
    vector<double> eta_face_sigma; //the errors on the tensor at final time
    vector<double> eta_u_n; // the errors on the displacement at time t_n; 
    vector<double> eta_u; // the errors on the displacement at final time
    vector<double> eta_sum; // the sum of all error indicators
    ConvergenceTable p_indicators_table;
    ConvergenceTable u_indicators_table;
    ConvergenceTable efficiency_table;

    string filename_base;
    int output_frequency;

    void set_control_parameters();
    void make_grid();
    // void setup_system();
    void setup_system_eg();
    void set_material_properties();
    void refine_mesh();
    // void assemble_system_pressure();
    void assemble_system_pressure_eg();
    void assemble_system_displacement();

    void solve_pressure_eg();
    void solve_displacement();

    void fixed_stress_iteration();

    vector<double> check_fs_convergence(); // check the convergence of fixed-stress iteration

    void output_displacement(int timestep, int fs_count) const;
    void output_pressure(int timestep, int fs_count) const;
    void output_error();
    void calc_error(); // compute the errors
    //void process_solution(int fs_count); // compute the errors
    void plot_error() const;

    void calc_strain_stress();    
    void calc_a_posteriori_indicators_p_eg();
    void calc_a_posteriori_indicators_u();

    double calc_u_energy_norm();
    void calc_p_h_norm();
    void calc_efficiency_indices();

};

#endif
