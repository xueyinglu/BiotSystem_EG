#include "BiotSystem.h"

BiotSystem::BiotSystem(ParameterHandler &param) : mpi_com(MPI_COMM_WORLD),
                                                 prm(param),
                                              fe_pressure(FE_Q<dim>(1), 1,
                                                            FE_DGQ<dim>(0), 1),
                                              dof_handler_pressure(triangulation),
                                              fe_displacement(FE_Q<dim>(1), dim),
                                              dof_handler_displacement(triangulation),
                                              permeability(0.05),
                                              lambda(0.5),
                                              mu(0.125),
                                              dof_handler_output(triangulation),
                                              fe_output(FE_DGQ<dim>(0), 1)

{
}
