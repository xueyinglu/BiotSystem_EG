#include "BiotSystem.h"
#include "PressureSolutionEG.h"
#include "TerzaghiPressureEG.h"
#include "MandelPressureEG.h"
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
    else if (test_case == TestCase::terzaghi)
    { // p_0 = 0; u_0 = 0;
        cout << "Benchmark Terzaghi : T_0 = " << T0 << endl;
        VectorTools::interpolate(dof_handler_pressure,
                                 TerzaghiPressureEG(T0),
                                 solution_pressure);
        prev_timestep_sol_pressure = solution_pressure;
        // Initialize u_0
        cout << "Terzaghi probelm: solving for u_0" << endl;
        assemble_system_displacement();
        solve_displacement();
        prev_timestep_sol_displacement = solution_displacement;
    }
    else if (test_case == TestCase::mandel)
    {

        cout << "Benchmark Mandel : T0 = " << T0 << endl;
        VectorTools::interpolate(dof_handler_pressure,
                                 MandelPressureEG(T0),
                                 solution_pressure);

        prev_timestep_sol_pressure = solution_pressure;
        // Initialize u_0
        cout << "Mandel probelm: solving for u_0" << endl;
        assemble_system_displacement();
        solve_displacement();
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

    for (timestep = 1; timestep < ((T - T0) / del_t + 1); timestep++)
    {
        t += del_t;
        cout << "--------------------- timestep = " << timestep << ", t = " << t << "------------------------" << endl;
        fixed_stress_iteration();

        if (criteria != 3)
        {
            calc_a_posteriori_indicators_p_eg();
            calc_a_posteriori_indicators_u();
        }
        if (timestep % output_frequency == 0)
        {
            plot_error();
            if (test_case == TestCase::benchmark || test_case == TestCase::terzaghi || test_case == TestCase::mandel)
            {
                calc_error();
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
