#include "BiotSystem.h"
#include "AuxTools.h"

void BiotSystem::calc_strain_stress()
{

    QGauss<dim> quadrature(degree + 1);
    FEValues<dim> fe_value_displacement(fe_displacement, quadrature,
                                        update_values | update_quadrature_points | update_gradients | update_hessians | update_JxW_values);
    typename DoFHandler<dim>::active_cell_iterator
        cell_u = dof_handler_displacement.begin_active(),
        endc = dof_handler_displacement.end();
    typename DoFHandler<dim>::active_cell_iterator cell_output = dof_handler_output.begin_active();
    const unsigned int n_q_points = quadrature.size();
    vector<double> lambda_values(n_q_points);
    vector<double> mu_values(n_q_points);
    vector<vector<Tensor<1, dim>>> grad_u_values(n_q_points, vector<Tensor<1, dim>>(dim));
    Tensor<2, dim> identity = Tensors::get_Identity<dim>();
    vector<types::global_dof_index> output_dofs(dof_handler_output.get_fe().dofs_per_cell);

    cell_stress_xx = 0.;
    cell_stress_xy = 0.;
    cell_stress_yy = 0.;
    cell_strain_xx = 0.;
    cell_strain_xy = 0.;
    cell_strain_yy = 0.;
    cell_vstrain = 0.;

    for (; cell_u != endc; ++cell_u, ++cell_output)
    {
        cell_output->get_dof_indices(output_dofs);
        fe_value_displacement.reinit(cell_u);
        lambda.value_list(fe_value_displacement.get_quadrature_points(), lambda_values);
        mu.value_list(fe_value_displacement.get_quadrature_points(), mu_values);
        if (test_case == TestCase::heterogeneous)
        {
            lambda_function.value_list(fe_value_displacement.get_quadrature_points(), lambda_values);
            mu_function.value_list(fe_value_displacement.get_quadrature_points(), mu_values);
        }
        fe_value_displacement.get_function_gradients(solution_displacement, grad_u_values);
        for (unsigned int q = 0; q < n_q_points; q++)
        {
            Tensor<2, dim> grad_u = Tensors::get_grad_u<dim>(q, grad_u_values);
            Tensor<2, dim> E = 0.5 * (grad_u + transpose(grad_u));
            Tensor<2, dim> sigma = 2 * mu_values[q] * E + lambda_values[q] * trace(E) * identity;
            double vstrain = trace(E);

            cell_strain_xx[output_dofs[0]] += E[0][0];
            cell_strain_xy[output_dofs[0]] += E[0][1];
            cell_strain_yy[output_dofs[0]] += E[1][1];
            cell_stress_xx[output_dofs[0]] += sigma[0][0];
            cell_stress_xy[output_dofs[0]] += sigma[0][1];
            cell_stress_yy[output_dofs[0]] += sigma[1][1];
            cell_vstrain[output_dofs[0]] += vstrain;
        }
    }
    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler_output);
    data_out.add_data_vector(cell_stress_xx, "stress_xx", DataOut<dim>::type_dof_data);
    data_out.add_data_vector(cell_stress_xy, "stress_xy", DataOut<dim>::type_dof_data);
    data_out.add_data_vector(cell_stress_yy, "stress_yy", DataOut<dim>::type_dof_data);
    data_out.add_data_vector(cell_strain_xx, "strain_xx", DataOut<dim>::type_dof_data);
    data_out.add_data_vector(cell_strain_xy, "strain_xy", DataOut<dim>::type_dof_data);
    data_out.add_data_vector(cell_strain_yy, "strain_yy", DataOut<dim>::type_dof_data);
    data_out.add_data_vector(cell_vstrain, "vstrain", DataOut<dim>::type_dof_data);
    data_out.build_patches();
    ofstream output("visual/stress" + to_string(timestep) + ".vtk");
    data_out.write_vtk(output);
}