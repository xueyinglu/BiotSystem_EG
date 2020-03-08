#include "BiotSystem.h"
#include "PressureSolution.h"
void BiotSystem::calc_p_h_norm()
{

    QGauss<dim> quadrature_pressure(fe_pressure.degree + 1);
    FEValues<dim> fe_value_pressure(fe_pressure, quadrature_pressure,
                                    update_values | update_quadrature_points | update_gradients | update_hessians | update_JxW_values);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_pressure.begin_active(),
                                                   endc = dof_handler_pressure.end();
    double h_norm_p_sq = 0;
    const unsigned int n_q_points = quadrature_pressure.size();
    vector<Tensor<1, dim>> grad_p_values(n_q_points);
    vector<double> permeability_values(n_q_points);
    for (; cell != endc; ++cell)
    {
        fe_value_pressure.reinit(cell);
        permeability.value_list(fe_value_pressure.get_quadrature_points(), permeability_values);
        if (test_case == TestCase::heterogeneous)
        {
            perm_function.value_list(fe_value_pressure.get_quadrature_points(), permeability_values);
        }
        for (unsigned int q = 0; q < n_q_points; q++)
        {
            Tensor<1, dim> true_grad_p;
            PressureSolution(t).gradient_value(fe_value_pressure.quadrature_point(q), true_grad_p);
            h_norm_p_sq += permeability_values[q] * (grad_p_values[q] - true_grad_p).norm_square() *fe_value_pressure.JxW(q);
        }
    }
    h_error_p_sq.push_back(h_norm_p_sq);

}