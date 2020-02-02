#include "BiotSystem.h"
void BiotSystem::solve_displacement(){
    SolverControl solver_control(1000, 1e-12);
  SolverCG<> cg(solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix_displacement, 1.2);

  cg.solve(system_matrix_displacement, solution_displacement, system_rhs_displacement,
           preconditioner);

  constraints_displacement.distribute(solution_displacement);
}