#include "BiotSystem.h"

void BiotSystem::set_newton_bc_pressure()
{

    cout << "    :  Set Newton BC PRESSURE " << std::endl;

    std::vector<bool> component_mask(2, false);
    component_mask[0] = true;
    component_mask[1] = false;

    VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                             0,
                                             ZeroFunction<dim>(2),
                                             constraints_pressure,
                                             component_mask);

    cout << "    :  Set Newton BC PRESSURE - END" << std::endl;
}