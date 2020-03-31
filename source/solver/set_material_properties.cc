#include "BiotSystem.h"

void BiotSystem::set_material_properties()
{
    cout << "set material properties" << endl;
    if (test_case == TestCase::terzaghi)
    {
        cout << " Terzaghi Phillips Phillips Test Case 1" << endl;
        double perm = 1e-3;
        double E = 1e5;  // Young's modulus
        double nu = 0.2; // Poisson's ratio

        double lame_lambda = E * nu / (1 + nu) / (1 - 2 * nu);
        double lame_mu = E / 2 / (1 + nu);

        traction_bc[0] = 0;
        traction_bc[1] = 1e3;
        pressure_dirichlet_bc = 0;
        initial_pressure_value = 0;
        biot_alpha = 1;
        mu_f = 1;
        K_b = lame_lambda + 2. / 3 * lame_mu;
        biot_inv_M = 0.1;
        permeability = ConstantFunction<dim>(perm);
        lambda = ConstantFunction<dim>(lame_lambda);
        mu = ConstantFunction<dim>(lame_mu);
    }

    if (test_case == TestCase::heterogeneous)
    {
        cout << " Heterogeneous test" << endl;
        double perm = 1e-16;
        mu_f = 1e-6;
        double E = 1e5;  // Young's modulus
        double nu = 0.2; // Poisson's ratio

        double lame_lambda = E * nu / (1 + nu) / (1 - 2 * nu);
        double lame_mu = E / 2 / (1 + nu);

        traction_bc[0] = 0;
        // traction_bc[1] = 1000;
        // pressure_dirichlet_bc = 0;
        // initial_pressure_value = 1000;
        traction_bc[1] = 2000;
        pressure_dirichlet_bc = 2000;
        initial_pressure_value = 2000;
        biot_alpha = 1;
        // biot_alpha = 0;
        K_b = lame_lambda + 2. / 3 * lame_mu;
        biot_inv_M = 0.1;
        permeability = ConstantFunction<dim>(perm);
        lambda = ConstantFunction<dim>(lame_lambda);
        mu = ConstantFunction<dim>(lame_mu);
    }

    if (test_case == TestCase:: mandel){
        /*
        // Bin Wang Data
        double perm = 1e-13;
        mu_f =1e-3;
        double F = 5.94e8;
        double E = 5.94e9;
        double nu =0.2;
        biot_alpha =1.0;
        biot_inv_M = 1.0 / 1.65e10;
        */
        // Phillips Phillips Data 
        double perm = 1.0e-2;
        mu_f =1.0;
        double F = 2.0e3;
        double E = 1.0e4;
        double nu =0.2;
        biot_alpha =1.0;
        biot_inv_M = 1e-4;
        
        permeability = ConstantFunction<dim>(perm);
        double lame_lambda = E * nu / (1 + nu) / (1 - 2 * nu);
        double lame_mu = E / 2 / (1 + nu);
        K_b = lame_lambda + 2. / 3 * lame_mu;
        // biot_inv_M = K_b/(1-biot_alpha)/(biot_alpha - phi0);
        lambda = ConstantFunction<dim>(lame_lambda);
        mu = ConstantFunction<dim>(lame_mu);
        traction_bc[0] = 0;
        traction_bc[1] = -2*F; 
        initial_pressure_value = 0;
        pressure_dirichlet_bc =0;
    }
}
