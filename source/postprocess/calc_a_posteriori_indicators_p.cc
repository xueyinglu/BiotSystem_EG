#include "BiotSystem.h"
#include <fstream>
#include <iostream>
using namespace std;

void BiotSystem::calc_a_posteriori_indicators_p()
{
    /* calculate eta_alg */
    double dum1 = 0;
    double dum2 = 0;
    for (auto &n : eta_fs)
    {
        dum1 += n;
        dum2 += n * n;
    }
    //eta_alg.push_back(del_t * dum1 * dum1 + h * h * dum2);
    eta_alg.push_back(del_t * dum1 * dum1);

    /* calculate eta_time */
    double eta_t_p_n = 0; // = Del_t /3 * \| k^{1/2} \nabla (p_h^{n,l} - p_h^{n-1})\|_{L2(\Omega)}
    QGauss<dim> quadrature_pressure(fe_pressure.degree + 1);
    FEValues<dim> fe_value_pressure(fe_pressure, quadrature_pressure,
                                    update_values | update_quadrature_points | update_gradients | update_hessians | update_JxW_values);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_pressure.begin_active(),
                                                   endc = dof_handler_pressure.end();
    typename DoFHandler<dim>::active_cell_iterator cell_output = dof_handler_output.begin_active();
    const unsigned int n_q_points = quadrature_pressure.size();
    vector<Tensor<1, dim>> grad_p_values(n_q_points);
    vector<Tensor<1, dim>> prev_timestep_grad_p_values(n_q_points);
    vector<double> permeability_values(n_q_points);

    vector<types::global_dof_index> output_dofs(dof_handler_output.get_fe().dofs_per_cell);
    cell_eta_time = 0;
    cell_eta_E_p = 0;
    for (; cell != endc; ++cell, ++cell_output)
    {
        fe_value_pressure.reinit(cell);
        fe_value_pressure.get_function_gradients(solution_pressure, grad_p_values);
        cell_output->get_dof_indices(output_dofs);
        fe_value_pressure.get_function_gradients(prev_timestep_sol_pressure, prev_timestep_grad_p_values);
        permeability.value_list(fe_value_pressure.get_quadrature_points(), permeability_values);
        for (unsigned int q = 0; q < n_q_points; q++)
        {
            Tensor<1, dim> cell_difference = grad_p_values[q] - prev_timestep_grad_p_values[q];
            eta_t_p_n += permeability_values[q] *
                         cell_difference.norm_square() *
                         fe_value_pressure.JxW(q);
            cell_eta_time[output_dofs[0]] += permeability_values[q] *
                                             cell_difference.norm_square() *
                                             fe_value_pressure.JxW(q);
        }
    }
    eta_t_p_n = eta_t_p_n * del_t / 3;
    if (timestep == 1)
    {
        eta_time.push_back(eta_t_p_n);
    }
    else
    {
        eta_time.push_back(eta_time.back() + eta_t_p_n);
    }

    /* calculate eta_p_residual, eta_flux_jump, eta_flow */

    double eta_flow_n = 0;
    double eta_E_p = 0;
    // residual at cells
    QGauss<dim> quadrature_displacement(fe_displacement.degree + 1);
    FEValues<dim> fe_value_displacement(fe_displacement, quadrature_displacement,
                                        update_values | update_quadrature_points | update_gradients | update_JxW_values);
    vector<vector<Tensor<1, dim>>> prev_timestep_grad_u_values(n_q_points, vector<Tensor<1, dim>>(dim));
    vector<vector<Tensor<1, dim>>> grad_u_values(n_q_points, vector<Tensor<1, dim>>(dim));

    vector<double> prev_timestep_p_values(n_q_points);
    vector<double> p_values(n_q_points);
    vector<double> laplacian_p_values(n_q_points);
    cell = dof_handler_pressure.begin_active();
    cell_output = dof_handler_output.begin_active();
    typename DoFHandler<dim>::active_cell_iterator
        cell_displacement = dof_handler_displacement.begin_active();
    for (; cell != endc; ++cell, ++cell_displacement, ++cell_output)
    {
        fe_value_pressure.reinit(cell);
        fe_value_displacement.reinit(cell_displacement);
        fe_value_pressure.get_function_values(solution_pressure, p_values);
        fe_value_pressure.get_function_values(prev_timestep_sol_pressure, prev_timestep_p_values);
        fe_value_pressure.get_function_laplacians(solution_pressure, laplacian_p_values);
        fe_value_displacement.get_function_gradients(solution_displacement, grad_u_values);
        fe_value_displacement.get_function_gradients(prev_timestep_sol_displacement, prev_timestep_grad_u_values);
        permeability.value_list(fe_value_pressure.get_quadrature_points(), permeability_values);

        cell_output->get_dof_indices(output_dofs);
        for (unsigned int q = 0; q < n_q_points; q++)
        {
            // TODO: extension to non constant permeability
            double cell_residual = 1 / mu_f * permeability_values[q] * laplacian_p_values[q] - biot_inv_M / del_t * (p_values[q] - prev_timestep_p_values[q]) - biot_alpha / del_t * (grad_u_values[q][0][0] + grad_u_values[q][1][1] - prev_timestep_grad_u_values[q][0][0] - prev_timestep_grad_u_values[q][1][1]);
            // cout << "cell residual = " << cell_residual << endl;
            eta_E_p += cell_residual * cell_residual * fe_value_pressure.JxW(q);
            cell_eta_E_p[output_dofs[0]] += cell_residual * cell_residual * fe_value_pressure.JxW(q);
        }
    }

    eta_E_p = eta_E_p * h * h * del_t;
    eta_p_residual.push_back(eta_E_p);

    // jump of flux at cell boundaries.
    double eta_flux_e = 0;
    QGauss<dim - 1> face_quadrature(fe_pressure.degree + 1);
    // const unsigned int n_face_q_points = face_quadrature.size();
    FEFaceValues<dim> fe_face_values(fe_pressure, face_quadrature,
                                     update_values | update_normal_vectors | update_gradients | update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_values_neighbor(fe_pressure, face_quadrature,
                                              update_values | update_gradients | update_quadrature_points | update_JxW_values);

    // vector<Tensor<1, dim>> face_grad_p_values(n_face_q_points);
    cell = dof_handler_pressure.begin_active();
    for (; cell != endc; ++cell)
    {
        fe_value_pressure.reinit(cell);
        for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
        {
            const auto face = cell->face(face_no);
            if (!face->at_boundary())
            {
                Assert(cell->neighbor(face_no).state() == IteratorState::valid, ExcInternalError());
                const auto neighbor = cell->neighbor(face_no);
                vector<Tensor<1, dim>> face_grad_p_values(fe_face_values.n_quadrature_points);
                vector<Tensor<1, dim>> neighbor_grad_p_values(fe_face_values.n_quadrature_points);
                vector<double> face_perm_values(fe_face_values.n_quadrature_points);
                const unsigned int neighbor_face = cell->neighbor_of_neighbor(face_no);
                fe_face_values.reinit(cell, face_no);
                fe_face_values_neighbor.reinit(neighbor, neighbor_face);
                fe_face_values.get_function_gradients(solution_pressure, face_grad_p_values);
                fe_face_values_neighbor.get_function_gradients(solution_pressure, neighbor_grad_p_values);
                permeability.value_list(fe_face_values.get_quadrature_points(), face_perm_values);
                // vector<Point<dim>> v_normal1 =fe_face_values.get_normal_vectors();
                // vector<Point<dim>> v_normal2 =fe_face_values_neighbor.get_normal_vectors();
                for (unsigned int q = 0; q < fe_face_values.n_quadrature_points; q++)
                {
                    const Tensor<1, dim> &n = fe_face_values.normal_vector(q);
                    // TODO: extension to when permeability changes over the face
                    double jump = face_perm_values[q] * (face_grad_p_values[q] * n - neighbor_grad_p_values[q] * n);
                    eta_flux_e += jump * jump * fe_face_values.JxW(q);
                }
            }
        }
    }
    eta_flux_e = eta_flux_e * h * del_t;
    eta_flux_jump.push_back(eta_flux_e);
    eta_flow_n = eta_E_p + eta_flux_e;
    if (timestep == 1)
    {
        eta_flow.push_back(eta_flow_n);
    }
    else
    {
        eta_flow.push_back(eta_flow.back() + eta_flow_n);
    }

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler_output);
    data_out.add_data_vector(cell_eta_time, "eta_time", DataOut<dim>::type_dof_data);
    data_out.add_data_vector(cell_eta_E_p, "eta_E_p", DataOut<dim>::type_dof_data);
    data_out.build_patches();
    ofstream output("visual/indicators-p" + to_string(timestep) + ".vtk");
    data_out.write_vtk(output);

    p_indicators_table.add_value("time", t);
    // p_indicators_table.add_value("eta_fs", eta_fs.back());
    p_indicators_table.add_value("eta_alg", eta_alg.back());
    p_indicators_table.add_value("eta_time", eta_time.back());
    p_indicators_table.add_value("eta_flow", eta_flow.back());
    // p_indicators_table.add_value("eta_p_residual", eta_p_residual.back());
    // p_indicators_table.add_value("eta_flux_jump", eta_flux_jump.back());
    p_indicators_table.add_value("sum", eta_alg.back() + eta_time.back() +eta_flow.back());
    p_indicators_table.add_value("error", l2_error_p.back());
    p_indicators_table.add_value("eff", (eta_alg.back() + eta_time.back() +eta_flow.back())/ l2_error_p.back());
}