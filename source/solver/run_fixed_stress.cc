#include "BiotSystem.h"
#include "PressureSolution.h"
#include "PressureSolutionEG.h"
using namespace std;
void BiotSystem::run_fixed_stress()
{
    set_control_parameters();
    set_material_properties();
    make_grid();
    setup_system_eg();
    if (test_case == TestCase::benchmark)
    {
        /* Interpolate initial EG pressure */
        // solution_pressure = 0.;
        // FEValuesExtractors::Scalar pressure_cg(0);
        // ComponentMask pressure_cg_mask = fe_pressure.component_mask(pressure_cg);

        VectorTools::interpolate(dof_handler_pressure,
                                 PressureSolutionEG(0),
                                 solution_pressure);

        prev_timestep_sol_pressure = solution_pressure;
        // Initialize u_0
        cout << "Initializing u_0" << endl;
        assemble_system_displacement();
        solve_displacement();
        prev_timestep_sol_displacement = solution_displacement;
    }
    else if (test_case == TestCase::terzaghi || test_case == TestCase::mandel)
    { // p_0 = 0; u_0 = 0;
        cout << "Benchmark Terzaghi : p_0 =0; u_0 = 0" << endl;
        VectorTools::interpolate(dof_handler_pressure,
                                 ZeroFunction<dim>(2),
                                 solution_pressure);
        prev_timestep_sol_pressure = solution_pressure;
        VectorTools::interpolate(dof_handler_displacement,
                                 ZeroFunction<dim>(dim),
                                 solution_displacement);
        prev_timestep_sol_displacement = solution_displacement;
    }
    else if (test_case == TestCase::heterogeneous)
    {
        cout << "heterogeneous test" << endl;
        vector<double> eg_initial_p;
        eg_initial_p.push_back(initial_pressure_value);
        eg_initial_p.push_back(0.);
        VectorTools::interpolate(dof_handler_pressure,
                                 ConstantFunction<dim>(eg_initial_p),
                                 solution_pressure);
        prev_timestep_sol_pressure = solution_pressure;
        cout << "Initializing u_0" << endl;
        assemble_system_displacement();
        solve_displacement();
        prev_timestep_sol_displacement = solution_displacement;
    }

    for (timestep = 1; timestep < ((T + 1e-5) / del_t); timestep++)
    {
        cout << "timestep = " << timestep << endl;
        t += del_t;
        fixed_stress_iteration();

        if (timestep % output_frequency == 0)
        {
            plot_error();
            if (test_case == TestCase::benchmark || test_case == TestCase::terzaghi || test_case == TestCase::mandel)
            {
                calc_error();
            }

            if (criteria != 3)
            {
                calc_a_posteriori_indicators_p_eg();
                calc_a_posteriori_indicators_u();
            }
            if (test_case == TestCase::benchmark || test_case == TestCase::mandel)
            {
                calc_efficiency_indices();
            }
            calc_strain_stress();
        }
        if (adaptivity == true)
        {
            refine_mesh();
        }

        prev_timestep_sol_displacement = solution_displacement;
        prev_timestep_sol_pressure = solution_pressure;
    }
    output_error();
}
