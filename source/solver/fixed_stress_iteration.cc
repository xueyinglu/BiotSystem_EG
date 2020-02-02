#include "BiotSystem.h"
using namespace std;
void BiotSystem::fixed_stress_iteration()
{
    bool iteration = true;
    int fs_count = 0;
    double l2_norm;
    while (iteration)
    {
        prev_fs_sol_pressure = solution_pressure;
        prev_fs_sol_displacement = solution_displacement;
        fs_count++;
        cout << "fixed stress no. " << fs_count << endl;
        cout << "solve for pressure" << endl;
        assemble_system_pressure_eg();
        solve_pressure_eg();
        cout << "solve for displacement" << endl;
        assemble_system_displacement();
        solve_displacement();
        l2_norm = check_fs_convergence();
        iteration = (l2_norm > tol_fixed_stress);

        // process_solution(fs_count);
        // plot_error(fs_count);
    }

    cout << "fixed stress iteration converged!" << endl;
    //TODO: not correct
    eta_fs.push_back(1 / sqrt(del_t) * l2_norm);
}