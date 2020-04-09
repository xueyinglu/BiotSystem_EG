
#include "BiotSystem.h"
#include "AuxTools.h"
#include "DisplacementSolution.h"
#include "MandelDisplacement.h"
using namespace std;

vector<double> BiotSystem::calc_u_energy_norm()
{

    QGauss<dim> quadrature_displacement(fe_displacement.degree + 2);
    FEValues<dim> fe_value_displacement(fe_displacement,
                                        quadrature_displacement, update_values | update_quadrature_points | update_gradients | update_JxW_values);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_displacement.begin_active(),
                                                   endc = dof_handler_displacement.end();
    const FEValuesExtractors::Vector displacements(0);
    const unsigned int n_q_points = quadrature_displacement.size();
    vector<Vector<double>> sol_u_values(n_q_points, Vector<double>(dim));
    vector<Tensor<1, dim>> ana_u_values(n_q_points);
    vector<vector<Tensor<1, dim>>> grad_u_values(n_q_points, vector<Tensor<1, dim>>(dim));
    DisplacementSolution disp(t);
    MandelDisplacement mandel_u(t);
    vector<double> lambda_values(n_q_points);
    vector<double> mu_values(n_q_points);

    double u_energy_norm = 0;
    double l2_strain = 0;
    double h_div_u = 0;
    for (; cell != endc; ++cell)
    {
        fe_value_displacement.reinit(cell);
        fe_value_displacement.get_function_values(solution_displacement, sol_u_values);
        fe_value_displacement.get_function_gradients(solution_displacement, grad_u_values);
        lambda.value_list(fe_value_displacement.get_quadrature_points(), lambda_values);
        mu.value_list(fe_value_displacement.get_quadrature_points(), mu_values);

        for (unsigned int q = 0; q < n_q_points; q++)
        {
            Tensor<2, dim> true_grad;
            disp.gradient_value(fe_value_displacement.quadrature_point(q), true_grad);
            if (test_case == TestCase::mandel)
            {
                mandel_u.gradient_value(fe_value_displacement.quadrature_point(q), true_grad);
            }
            Tensor<2, dim> grad_u = Tensors::get_grad_u<dim>(q, grad_u_values);
            Tensor<2, dim> e_strain = 0.5 * ((grad_u - true_grad) + transpose(grad_u - true_grad));
            double e_div = Tensors::get_divergence_u(grad_u - true_grad);
            l2_strain += 4 * mu_values[q] * mu_values[q] * e_strain.norm_square() * fe_value_displacement.JxW(q);
            h_div_u += lambda_values[q] * lambda_values[q] * e_div * e_div * fe_value_displacement.JxW(q);
            u_energy_norm += (2 * mu_values[q] * e_strain.norm_square() + lambda_values[q] * e_div * e_div) * fe_value_displacement.JxW(q);
        }
    }
    u_energy_norm = sqrt(u_energy_norm);
    double energy_u_1 = sqrt(l2_strain) + sqrt(h_div_u);
    vector<double> energy_u;
    energy_u.push_back(u_energy_norm);
    energy_u.push_back(energy_u_1);
    return energy_u;
}