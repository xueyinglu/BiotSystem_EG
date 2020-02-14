#include "BiotSystem.h"

void BiotSystem::set_newton_bc_pressure()
{
    if (test_case == TestCase::benchmark)
    {

        cout << "    :  Set Benchmark BC PRESSURE " << std::endl;

        std::vector<bool> component_mask(2, false);
        component_mask[0] = true;
        component_mask[1] = false;

        VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                                 0,
                                                 ZeroFunction<dim>(2),
                                                 constraints_pressure,
                                                 component_mask);
    }
    else if (test_case == TestCase::terzaghi || test_case == TestCase::heterogeneous)
    {
        cout << " Set Terzaghi/Heterogeneous PRESSURE BC" << std::endl;
        vector<double> eg_p_bc;
        eg_p_bc.push_back(pressure_dirichlet_bc);
        eg_p_bc.push_back(0.);

        VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                                 2,
                                                 ConstantFunction<dim>(eg_p_bc),
                                                 constraints_pressure);
        // Shall we set the DG component to zero?
        /*
        component_mask[0] = false;
        component_mask[1] = true;

        VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                                 2,
                                                 ZeroFunction<dim>(2),
                                                 constraints_pressure,
                                                 component_mask);
        */
        cout << " End of set Terzaghi PRESSURE BC" << std::endl;
    }
}