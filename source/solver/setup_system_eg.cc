#include "BiotSystem.h"

void BiotSystem::setup_system_eg()
{
  /********************************** pressure *********************************/
  system_matrix_pressure.clear();
  dof_handler_pressure.distribute_dofs(fe_pressure);

  DoFRenumbering::component_wise(dof_handler_pressure); //, block_component);

  std::vector<types::global_dof_index> dofs_per_component(2);
  DoFTools::count_dofs_per_component(dof_handler_pressure,
                                     dofs_per_component);

  const unsigned int n_p_cg = dofs_per_component[0];
  const unsigned int n_p_dg = dofs_per_component[1];
  cout << "====== # of Pressure Dofs ==========================" << std::endl;
  cout << "DoFs(CG): " << n_p_cg << " (DG) + " << n_p_dg << " = "
       << n_p_cg + n_p_dg << std::endl;

  // number_of_dof_cycle_pres_CG.push_back(n_p_cg);
  // number_of_dof_cycle_pres_DG.push_back(n_p_dg);
  // number_of_dof_cycle_pres_EG.push_back(n_p_cg + n_p_dg);

  // MPI
  partition_pressure.clear();

  partition_pressure.push_back(dof_handler_pressure.locally_owned_dofs().get_view(0, n_p_cg));
  partition_pressure.push_back(dof_handler_pressure.locally_owned_dofs().get_view(n_p_cg, n_p_cg + n_p_dg));
  DoFTools::extract_locally_relevant_dofs(dof_handler_pressure, relevant_set_pressure);
  partition_relevant_pressure.clear();

  partition_relevant_pressure.push_back(relevant_set_pressure.get_view(0, n_p_cg));
  partition_relevant_pressure.push_back(relevant_set_pressure.get_view(n_p_cg, n_p_cg + n_p_dg));

  constraints_pressure.clear();
  constraints_pressure.reinit(relevant_set_pressure);

  if (bCG_WeaklyBD == false )
    set_newton_bc_pressure();

  DoFTools::make_hanging_node_constraints(dof_handler_pressure,
                                          constraints_pressure);

  constraints_pressure.close();

  TrilinosWrappers::BlockSparsityPattern csp(partition_pressure, mpi_com);

  DoFTools::make_flux_sparsity_pattern(dof_handler_pressure,
                                       csp,
                                       constraints_pressure,
                                       false,
                                       Utilities::MPI::this_mpi_process(mpi_com));
  csp.compress();

  system_matrix_pressure.reinit(csp);

  solution_pressure.reinit(partition_relevant_pressure, mpi_com);

  prev_timestep_sol_pressure.reinit(partition_relevant_pressure, mpi_com);
  prev_fs_sol_pressure.reinit(partition_relevant_pressure, mpi_com);

  system_rhs_pressure.reinit(partition_pressure, mpi_com);
  /*************************** displacement **********************************/
  dof_handler_displacement.distribute_dofs(fe_displacement);
  std::cout << "Number of degrees of freedom for displacement: " << dof_handler_displacement.n_dofs() << std::endl;
  constraints_displacement.clear();
  DoFTools::make_hanging_node_constraints(dof_handler_displacement,
                                          constraints_displacement);
  constraints_displacement.close();

  DynamicSparsityPattern dsp_displacement(dof_handler_displacement.n_dofs(), dof_handler_displacement.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler_displacement, dsp_displacement,
                                  constraints_displacement,
                                  /*keep_constrained_dofs = */ true);
  sparsity_pattern_displacement.copy_from(dsp_displacement);

  system_matrix_displacement.reinit(sparsity_pattern_displacement);

  solution_displacement.reinit(dof_handler_displacement.n_dofs());
  system_rhs_displacement.reinit(dof_handler_displacement.n_dofs());

  dof_handler_output.distribute_dofs(fe_output);
  cell_eta_time.reinit(dof_handler_output.n_dofs());
  cell_eta_E_p.reinit(dof_handler_output.n_dofs());
  cell_eta_E_partial_u.reinit(dof_handler_output.n_dofs());
  cell_eta_E_u.reinit(dof_handler_output.n_dofs());
}