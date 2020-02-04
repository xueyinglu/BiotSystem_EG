#ifndef BIOT_SYSTEM_H_
#define BIOT_SYSTEM_H_
#include "DealiiHeader.h"
#include "InitialPressure.h"
#include "RightHandSide.h"
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
    BiotSystem();
    BiotSystem(int _num_global_refinement, double _del_t, double _T, double _fs_tol);
    BiotSystem(string testcase, int _num_global_refinement, double _del_t, double _T);
    // virtual BiotSystem();
    void run_fixed_stress();
    void check_disp_solver_convergence();

    enum TestCase
    {
        none,
        simple,
        mandel
    };

private:
    double del_t = 0.01;
    double T = 1;
    double t = 0;
    int timestep = 0;
    int num_global_refinement = 5;
    double h = 1. / pow(2, num_global_refinement);

    Triangulation<dim> triangulation;
    /* EG pressure solution */
    FESystem<dim> fe_pressure;

    DoFHandler<dim> dof_handler_pressure;
    BlockSparsityPattern sparsity_pattern_pressure;
    ConstraintMatrix constraints_pressure;

    LA::MPI::BlockSparseMatrix system_matrix_pressure;
    LA::MPI::BlockVector solution_pressure, prev_timestep_sol_pressure, prev_fs_sol_pressure;
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

    Vector<double> solution_displacement;
    Vector<double> system_rhs_displacement;
    ConvergenceTable convergence_table;

    vector<double> l2_error_p;
    vector<double> l2_error_u;
    vector<double> energy_error_u;

    // Data
    double mu_f = 1;               // fluid viscosity
    RightHandSide right_hand_side; // mechanics equation body force
    InitialPressure initial_pressure;
    ConstantFunction<dim> permeability;
    ConstantFunction<dim> lambda, mu;

    TestCase test_case;
    // coupling

    double biot_alpha = 0.75;
    double K_b = 7. / 12; //K_b = lambda +2/3*mu
    double biot_inv_M = 3. / 28;
    double tol_fixed_stress = 1e-5;
    Vector<double> prev_timestep_sol_displacement;
    Vector<double> prev_fs_sol_displacement;

    // EG for flow
    int degree = 1;
    double gamma_penal;
    double min_cell_diameter;
    bool bNeaumannBD = false;
    bool bCG_WeaklyBD = false;
    double d_SForm = 0; //SIPG

    /* element-wise a posteriori error indicators*/
    DoFHandler<dim> dof_handler_output;
    FESystem<dim> fe_output;
    Vector<double> cell_eta_time;
    Vector<double> cell_eta_E_p;
    Vector<double> cell_eta_E_partial_u;
    Vector<double> cell_eta_E_u;
    /* global a posteriori error estimators (recorded for each time step) */
    vector<double> eta_fs;
    vector<double> eta_alg;
    vector<double> eta_time;
    // vector<double> eta_p_residual; // error of flow residual at time t_n
    // vector<double> eta_flux_jump;  // error of flux jump at time t_n
    vector<double> eta_flow;
    vector<double> eta_jump;
    vector<double> eta_pen;
    vector<double> eta_partial_p_J;
    vector<double> eta_p_J;

    vector<double> eta_face_partial_sigma_n; // the errors on the tensor's time derivative for time step n
    vector<double> eta_face_partial_sigma;   // the (cumulative in time ) errors on the tensor's time derivative_form
    vector<double> eta_partial_u_n;          // the errors on the displacement's time derivative for time step n
    vector<double> eta_partial_u;            // the (cumulative in time) errors on the displacement's time derivative
    vector<double> eta_face_sigma_n;         //the errors on the tensor at time t_n;
    vector<double> eta_face_sigma;           //the errors on the tensor at final time
    vector<double> eta_u_n;                  // the errors on the displacement at time t_n;
    vector<double> eta_u;                    // the errors on the displacement at final time
    vector<double> eta_sum;                  // the sum of all error indicators
    ConvergenceTable p_indicators_table;
    ConvergenceTable u_indicators_table;
    ConvergenceTable efficiency_table;

    void make_grid();
    // void setup_system();
    void setup_system_eg();

    // void assemble_system_pressure();
    void assemble_system_pressure_eg();
    void assemble_system_displacement();

    void set_newton_bc_pressure();
    void solve_pressure_eg();
    void solve_displacement();

    void fixed_stress_iteration();

    double check_fs_convergence(); // check the convergence of fixed-stress iteration

    void output_displacement(int timestep, int fs_count) const;
    void output_pressure(int timestep, int fs_count) const;
    void output_error();
    void calc_error(); // compute the errors
    //void process_solution(int fs_count); // compute the errors
    //  void plot_error() const;

    void calc_a_posteriori_indicators_p_eg();
    void calc_a_posteriori_indicators_u();

    double calc_u_energy_norm();
    void calc_efficiency_indices();

    MPI_Comm mpi_com;
};

#endif