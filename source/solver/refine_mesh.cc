#include "BiotSystem.h"

void BiotSystem::refine_mesh()
{

    /*
    GridRefinement::refine_and_coarsen_fixed_fraction(triangulation,
                                                      cell_eta_u,
                                                      0.6,
                                                      0.4);
                                                      */
    //GridRefinement::refine_and_coarsen_fixed_number(triangulation, cell_eta_p,0.1, 0.4);
    Vector<double> cell_eta_refine = cell_eta_p;
    cell_eta_refine /= cell_eta_p.linfty_norm();
    Vector<double> dum = cell_eta_u;
    dum /= cell_eta_u.linfty_norm();
    cell_eta_refine += dum;
    parallel::distributed::GridRefinement::
        refine_and_coarsen_fixed_number(triangulation,
                                        cell_eta_refine,
                                        0.1, 0.2);
    const unsigned int max_grid_level = 9;
    const unsigned int min_grid_level = 3;
    if (triangulation.n_levels() > max_grid_level)
    {
        for (const auto &cell :
             triangulation.active_cell_iterators_on_level(max_grid_level))
            cell->clear_refine_flag();
        for (const auto &cell :
             triangulation.active_cell_iterators_on_level(min_grid_level))
            cell->clear_coarsen_flag();
    }
    cout << "prepare to refine mesh " << endl;
    parallel::distributed::SolutionTransfer<dim, LA::MPI::BlockVector>
        solution_trans_p(dof_handler_pressure);
    SolutionTransfer<dim> solution_trans_u(dof_handler_displacement);
    triangulation.prepare_coarsening_and_refinement();
    LA::MPI::BlockVector prev_sol_p = solution_pressure;
    Vector<double> prev_sol_u = solution_displacement;
    cout << "line 39 " << endl;
    solution_trans_p.prepare_for_coarsening_and_refinement(prev_sol_p);
    cout << "line 40 " << endl;
    solution_trans_u.prepare_for_coarsening_and_refinement(prev_sol_u);
    triangulation.execute_coarsening_and_refinement();
    cout << "line 43 " << endl;
    setup_system_eg();
    cout << "line 42 " << endl;
    LA::MPI::BlockVector distributed_solution(partition_pressure);
    solution_trans_p.interpolate(distributed_solution);
    solution_pressure = distributed_solution;
    solution_trans_u.interpolate(prev_sol_u, solution_displacement);
    constraints_pressure.distribute(solution_pressure);
    constraints_displacement.distribute(solution_displacement);

    typename DoFHandler<dim>::active_cell_iterator cell =
                                                       dof_handler_pressure.begin_active(),
                                                   endc = dof_handler_pressure.end();

    for (; cell != endc; ++cell)
        if (cell->is_locally_owned())
        {
            min_cell_diameter = std::min(cell->diameter(), min_cell_diameter);
        }

    min_cell_diameter = -Utilities::MPI::max(-min_cell_diameter, mpi_com);
}