#include "BiotSystem.h"
void BiotSystem::output_pressure(int timestep, int fs_count) const {
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler_pressure);
    data_out.add_data_vector(solution_pressure, "solution");
    data_out.build_patches();
    //std::ofstream output("pressure-" + std::to_string(timestep) +"-"+std::to_string(fs_count)+ ".vtk");
    std::ofstream output("visual/pressure-" + std::to_string(timestep) + ".vtk");
    data_out.write_vtk(output);
}