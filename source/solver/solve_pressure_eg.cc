#include "BiotSystem.h"
#include "AuxTools.h"
void BiotSystem::solve_pressure_eg(){
    LA::MPI::BlockVector sol(system_rhs_pressure);

    SolverControl solver_control(10000, 1e-10 * system_rhs_pressure.l2_norm());
    // SolverControl solver_control(10000, 1e-12);
    SolverGMRES<LA::MPI::BlockVector> solver(solver_control);

    {
        LA::MPI::PreconditionAMG::AdditionalData data;
        data.elliptic = true;
        data.higher_order_elements = false;
        data.smoother_sweeps = 2;
        data.aggregation_threshold = 0.01;
        preconditioner_pressure_cg.initialize(system_matrix_pressure.block(0, 0), data);
        preconditioner_pressure_dg.initialize(system_matrix_pressure.block(1, 1), data);
    }

    BlockDiagonalPreconditioner<LA::MPI::PreconditionAMG, LA::MPI::PreconditionAMG>
        preconditioner(system_matrix_pressure,
                       preconditioner_pressure_cg, preconditioner_pressure_dg);

    solver.solve(system_matrix_pressure, solution_pressure,
                 system_rhs_pressure, preconditioner);

    constraints_pressure.distribute(solution_pressure);
    // newton_update_pressure = sol;

    //solver_dg.solve(system_matrix_pressure.block(1,1), sol_inc.block(1), Residual_Vector.block(1));
    //solver_cg.solve(system_matrix_pressure.block(0,0), sol_inc.block(0), Residual_Vector.block(0));

}