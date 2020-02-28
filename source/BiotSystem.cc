#include "BiotSystem.h"
BiotSystem::BiotSystem() : fe_pressure(FE_Q<dim>(1), 1),
                           dof_handler_pressure(triangulation),
                           fe_displacement(FE_Q<dim>(1), dim),
                           dof_handler_displacement(triangulation),
                           permeability(0.05),
                           lambda(0.5),
                           mu(0.125),
                           dof_handler_output(triangulation),
                           fe_output(FE_DGQ<dim>(0), 1)
{
    test_case = TestCase::terzaghi;
}
BiotSystem::BiotSystem(int _num_global_refinement, double _del_t, double _T, double _fs_tol) : fe_pressure(FE_Q<dim>(1), 1,
                                                                                                FE_DGQ<dim>(0),1),
                                                                                               dof_handler_pressure(triangulation),
                                                                                               fe_displacement(FE_Q<dim>(1), dim),
                                                                                               dof_handler_displacement(triangulation),
                                                                                               permeability(0.05),
                                                                                               lambda(0.5),
                                                                                               mu(0.125),
                                                                                               dof_handler_output(triangulation),
                                                                                               fe_output(FE_DGQ<dim>(0), 1),
                                                                                               mpi_com(MPI_COMM_WORLD)
{
    num_global_refinement = _num_global_refinement;
    del_t = _del_t;
    T = _T;
    h = 1. / std::pow(2, num_global_refinement);
    // test_case = TestCase::benchmark;
    test_case = TestCase::heterogeneous;
    tol_fixed_stress = _fs_tol;
    min_cell_diameter = h;
    gamma_penal = 1;
    criteria = 1;
    adaptivity = true;
}

BiotSystem::BiotSystem(int _num_global_refinement, double _del_t, double _T, double _fs_tol, int _criteria) : fe_pressure(FE_Q<dim>(1), 1,
                                                                                                FE_DGQ<dim>(0),1),
                                                                                                dof_handler_pressure(triangulation),
                                                                                                fe_displacement(FE_Q<dim>(1), dim),
                                                                                                dof_handler_displacement(triangulation),
                                                                                                permeability(0.05),
                                                                                                lambda(0.5),
                                                                                                mu(0.125),
                                                                                                dof_handler_output(triangulation),
                                                                                                fe_output(FE_DGQ<dim>(0), 1),
                                                                                                mpi_com(MPI_COMM_WORLD)
{
    num_global_refinement = _num_global_refinement;
    del_t = _del_t;
    T = _T;
    h = 1. / std::pow(2, num_global_refinement);
    // test_case = TestCase::benchmark;
    test_case = TestCase::terzaghi;
    tol_fixed_stress = _fs_tol;
    criteria = _criteria;
}
