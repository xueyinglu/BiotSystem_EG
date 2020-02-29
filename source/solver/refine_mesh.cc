#include "BiotSystem.h"

void BiotSystem::refine_mesh()
{

    Vector<double> cell_eta_refine = cell_eta_p;
    cout << "cell_eta_p.L_infty =" << cell_eta_p.linfty_norm() << endl;
    cell_eta_refine /= cell_eta_p.linfty_norm();
    Vector<double> dum = cell_eta_u;
    dum /= cell_eta_u.linfty_norm();
    cout << "cell_eta_u.L_infty =" << cell_eta_u.linfty_norm() << endl;
    cell_eta_refine += dum;

    for (const auto &n : cell_eta_refine)
    {
        if (n < 0)
        {
            cout << n << endl;
        }
    }
    /*
    parallel::distributed::GridRefinement::
        refine_and_coarsen_fixed_number(triangulation,
                                        cell_eta_refine,
                                        0.1, 0.);
    */
    GridRefinement::refine_and_coarsen_fixed_number(triangulation, cell_eta_refine, 0.1, 0.2);
    const unsigned int max_grid_level = 9;
    const unsigned int min_grid_level = 3;
    if (triangulation.n_levels() > max_grid_level)
    {
        for (const auto &cell :
             triangulation.active_cell_iterators_on_level(max_grid_level))
            if (cell->is_locally_owned())
                cell->clear_refine_flag();
        for (const auto &cell :
             triangulation.active_cell_iterators_on_level(min_grid_level))
            if (cell->is_locally_owned())
                cell->clear_coarsen_flag();
    }
    SolutionTransfer<dim, BlockVector<double>> solution_trans_p(dof_handler_pressure);
    SolutionTransfer<dim> solution_trans_u(dof_handler_displacement);
    triangulation.prepare_coarsening_and_refinement();
    LA::MPI::BlockVector prev_sol_p = solution_pressure;
    // workaround for converting LA::MPI::BlockVector to BlockVector to use ordinary GridRefinement.
    BlockVector<double> test(prev_sol_p);
    Vector<double> prev_sol_u = solution_displacement;
    cout << "line 39 " << endl;
    solution_trans_p.prepare_for_coarsening_and_refinement(test);
    cout << "line 40 " << endl;
    solution_trans_u.prepare_for_coarsening_and_refinement(prev_sol_u);
    triangulation.execute_coarsening_and_refinement();
    cout << "line 43 " << endl;
    setup_system_eg();
    cout << "line 60 " << endl;
    /*
    LA::MPI::BlockVector distributed_solution(partition_pressure);
    solution_trans_p.interpolate(distributed_solution);
    solution_pressure = distributed_solution;
*/
    BlockVector<double> dum2(solution_pressure);
    solution_trans_p.interpolate(test, dum2 );
    cout << "line 65 " << endl;
    solution_pressure = dum2;
    cout << "line 71 " << endl;
    solution_trans_u.interpolate(prev_sol_u, solution_displacement);
    cout << "line 73 " << endl;
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
