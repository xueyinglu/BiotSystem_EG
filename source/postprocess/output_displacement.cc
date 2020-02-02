#include "BiotSystem.h"
#include <iostream>
#include <fstream>
using namespace std;

void BiotSystem::output_displacement(int time, int fs_count) const{
    // ofstream output("displacment-" + to_string(time) + "-" + to_string(fs_count) +".vtk");
    ofstream output("visual/displacment-" + to_string(time) +".vtk");
    vector<string> sol_names;
    sol_names.push_back("x_disp");
    sol_names.push_back("y_disp");
    
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler_displacement);
    data_out.add_data_vector(solution_displacement, sol_names);
    data_out.build_patches();
    data_out.write_vtk(output);
}