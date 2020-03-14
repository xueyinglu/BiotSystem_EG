#include "BiotSystem.h"

void BiotSystem::set_control_parameters()
{
    prm.enter_subsection("Global parameters");
    num_global_refinement = prm.get_integer("Global refinement numbers");
    h = 1. / pow(2, num_global_refinement);
    min_cell_diameter = h;
    cout <<"h=" << h << endl;
    T = prm.get_double("Final time");
    del_t = prm.get_double("Timestep size");
    criteria = prm.get_integer("Fixed-stress criteria");
    tol_fixed_stress = prm.get_double("Fixed-stress tolerance");
    filename_base = prm.get("Output filename");
    output_frequency = prm.get_integer("Output frequency");
    adaptivity = prm.get_bool("Adaptivity flag");
    b_p_mult = prm.get_bool("Permeability multiplier flag");
    if (prm.get("Test case") == "benchmark")
    {
        test_case = TestCase::benchmark;
    }
    else if (prm.get("Test case") == "terzaghi")
    {
        test_case = TestCase::terzaghi;
    }
    else if (prm.get("Test case") == "mandel")
    {
        test_case = TestCase::mandel;
    }
    else if (prm.get("Test case") == "heterogeneous")
    {
        test_case = TestCase::heterogeneous;
    }
    gamma_penal = prm.get_double("Penalization parameter");
    cout << "gamma_penal = " << gamma_penal << endl;
    prm.leave_subsection();
}