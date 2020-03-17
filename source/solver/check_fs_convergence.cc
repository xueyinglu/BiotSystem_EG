#include "BiotSystem.h"
#include "AuxTools.h"
using namespace std;
// check the convergence of fixed stress, return the residual
vector<double> BiotSystem::check_fs_convergence()
{
    QGauss<dim> quadrature_pressure(degree + 1);
    QGauss<dim> quadrature_displacement(fe_displacement.degree + 1);
    FEValues<dim> fe_value_pressure(fe_pressure,
                                    quadrature_pressure, update_values | update_quadrature_points | update_gradients | update_JxW_values);
    FEValues<dim> fe_value_displacement(fe_displacement,
                                        quadrature_displacement, update_values | update_quadrature_points | update_gradients | update_JxW_values);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_pressure.begin_active(),
                                                   endc = dof_handler_pressure.end();

    typename DoFHandler<dim>::active_cell_iterator
        cell_displacement = dof_handler_displacement.begin_active();
    const unsigned int n_q_points = quadrature_pressure.size();

    vector<vector<Tensor<1, dim>>> prev_fs_sol_grad_u_values(n_q_points, vector<Tensor<1, dim>>(dim));
    vector<vector<Tensor<1, dim>>> grad_u_values(n_q_points, vector<Tensor<1, dim>>(dim));
    vector<Vector<double>> prev_fs_sol_pressure_values(n_q_points, Vector<double>(2));
    vector<Vector<double>> pressure_values(n_q_points, Vector<double>(2));

    Vector<double> difference_in_u = solution_displacement;
    difference_in_u -= prev_fs_sol_displacement;
    cout << "difference in u = " << difference_in_u.l2_norm() << endl;

    double mean_stress;
    double prev_fs_mean_stress;
    double residual = 0.0;
    double l2square_mean_stress = 0.0;

    for (; cell != endc; ++cell, ++cell_displacement)
    {
        fe_value_pressure.reinit(cell);
        fe_value_displacement.reinit(cell_displacement);
        fe_value_pressure.get_function_values(solution_pressure, pressure_values);
        fe_value_pressure.get_function_values(prev_fs_sol_pressure, prev_fs_sol_pressure_values);
        fe_value_displacement.get_function_gradients(solution_displacement, grad_u_values);
        fe_value_displacement.get_function_gradients(prev_fs_sol_displacement, prev_fs_sol_grad_u_values);

        for (unsigned int q = 0; q < n_q_points; q++)
        {
            const Tensor<2, dim> grad_u = Tensors ::get_grad_u<dim>(q, grad_u_values);
            const double div_u = Tensors ::get_divergence_u<dim>(grad_u);
            const Tensor<2, dim> prev_fs_grad_u = Tensors ::get_grad_u<dim>(q, prev_fs_sol_grad_u_values);
            const double prev_fs_div_u = Tensors ::get_divergence_u<dim>(prev_fs_grad_u);
            mean_stress = K_b * div_u - biot_alpha * (pressure_values[q][0] + pressure_values[q][1]);
            prev_fs_mean_stress = K_b * prev_fs_div_u - biot_alpha * (prev_fs_sol_pressure_values[q][0] + prev_fs_sol_pressure_values[q][1]);
            residual += (mean_stress - prev_fs_mean_stress) * (mean_stress - prev_fs_mean_stress) * fe_value_pressure.JxW(q);
            l2square_mean_stress += mean_stress * mean_stress * fe_value_pressure.JxW(q);
        }
    }
    //cout << "fixed stress iteration convergence criteria = " << sqrt(residual/l2square_mean_stress) << endl;

    double change_ms = sqrt(residual);;
    double rel_change_ms = sqrt(residual / l2square_mean_stress);
    if (criteria != 2)
    {
        
        cout << "fixed stress iteration convergence criteria 1 = " << change_ms << endl;
    }
    else if (criteria == 2)
    {
        
        cout << "fixed stress iteration convergence criteria 2 = " << change_ms << endl;
    }
    vector<double> results;
    results.push_back(change_ms);
    results.push_back(rel_change_ms);
    return results;
}
