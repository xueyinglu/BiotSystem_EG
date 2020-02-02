#include "BiotSystem.h"
#include "AuxTools.h"
using namespace std;
#include <numeric>
void BiotSystem::calc_a_posteriori_indicators_u()
{

    /* calculate eta_face_partial_sigma_n, eta_face_partial_sigma, eta_face_sigma_n, eta_face_sigma */
    double eta_e_partial_sigma = 0;
    double eta_e_sigma = 0;
    QGauss<dim - 1> face_quadrature(fe_pressure.degree + 1);
    // const unsigned int n_face_q_points = face_quadrature.size();
    FEFaceValues<dim> fe_face_p(fe_pressure, face_quadrature,
                                update_values | update_normal_vectors | update_gradients | update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_neighbor_p(fe_pressure, face_quadrature,
                                         update_values | update_gradients | update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_u(fe_displacement, face_quadrature,
                                update_values | update_normal_vectors | update_gradients | update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_neighbor_u(fe_displacement, face_quadrature,
                                         update_values | update_gradients | update_quadrature_points | update_JxW_values);
    const unsigned int dofs_per_cell = fe_displacement.dofs_per_cell;
    Tensor<2, dim> identity = Tensors::get_Identity<dim>();
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_pressure.begin_active(),
                                                   endc = dof_handler_pressure.end();
    typename DoFHandler<dim>::active_cell_iterator
        cell_u = dof_handler_displacement.begin_active();

    for (; cell != endc; ++cell, ++cell_u)
    {
        for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
        {
            typename DoFHandler<dim>::face_iterator face_p = cell->face(face_no);
            typename DoFHandler<dim>::face_iterator face_u = cell_u->face(face_no);
            if (!face_p->at_boundary())
            //if (! cell -> at_boundary(face_no))
            {
                Assert(cell->neighbor(face_no).state() == IteratorState::valid, ExcInternalError());
                const typename DoFHandler<dim>::cell_iterator neighbor_p = cell->neighbor(face_no);
                const typename DoFHandler<dim>::cell_iterator neighbor_u = cell_u->neighbor(face_no);
                vector<vector<Tensor<1, dim>>> face_grad_u_values(fe_face_u.n_quadrature_points, vector<Tensor<1, dim>>(dim));
                vector<vector<Tensor<1, dim>>> neighbor_grad_u_values(fe_face_u.n_quadrature_points, vector<Tensor<1, dim>>(dim));
                vector<vector<Tensor<1, dim>>> prev_timestep_face_grad_u_values(fe_face_u.n_quadrature_points, vector<Tensor<1, dim>>(dim));
                vector<vector<Tensor<1, dim>>> prev_timestep_neighbor_grad_u_values(fe_face_u.n_quadrature_points, vector<Tensor<1, dim>>(dim));
                vector<double> face_p_values(fe_face_p.n_quadrature_points);
                vector<double> neighbor_p_values(fe_face_p.n_quadrature_points);
                vector<double> prev_timestep_face_p_values(fe_face_p.n_quadrature_points);
                vector<double> prev_timestep_neighbor_p_values(fe_face_p.n_quadrature_points);
                vector<double> lambda_values(fe_face_p.n_quadrature_points);
                vector<double> mu_values(fe_face_p.n_quadrature_points);

                const unsigned int neighbor_face_p = cell->neighbor_of_neighbor(face_no);
                const unsigned int neighbor_face_u = cell_u->neighbor_of_neighbor(face_no);
                fe_face_p.reinit(cell, face_no);
                fe_face_neighbor_p.reinit(neighbor_p, neighbor_face_p);
                fe_face_u.reinit(cell_u, face_no);
                fe_face_neighbor_u.reinit(neighbor_u, neighbor_face_u);
                fe_face_p.get_function_values(solution_pressure, face_p_values);
                fe_face_p.get_function_values(prev_timestep_sol_pressure, prev_timestep_face_p_values);
                fe_face_neighbor_p.get_function_values(solution_pressure, neighbor_p_values);
                fe_face_neighbor_p.get_function_values(prev_timestep_sol_pressure, prev_timestep_neighbor_p_values);

                fe_face_u.get_function_gradients(solution_displacement, face_grad_u_values);
                fe_face_u.get_function_gradients(prev_timestep_sol_displacement, prev_timestep_face_grad_u_values);
                fe_face_neighbor_u.get_function_gradients(solution_displacement, neighbor_grad_u_values);
                fe_face_neighbor_u.get_function_gradients(prev_timestep_sol_displacement, prev_timestep_neighbor_grad_u_values);

                lambda.value_list(fe_face_p.get_quadrature_points(), lambda_values);
                mu.value_list(fe_face_p.get_quadrature_points(), mu_values);
                //vector<Point<dim>> v_normal1 = fe_face_p.get_normal_vectors();
                // vector<Point<dim>> v_normal2 = fe_face_neighbor_p.get_normal_vectors();
                for (unsigned int q = 0; q < fe_face_p.n_quadrature_points; q++)
                {
                    Tensor<2, dim> face_grad_u = Tensors::get_grad_u<dim>(q, face_grad_u_values) - Tensors::get_grad_u<dim>(q, prev_timestep_face_grad_u_values);
                    Tensor<2, dim> face_E = 0.5 * (face_grad_u + transpose(face_grad_u));
                    Tensor<2, dim> face_sigma = 2 * mu_values[q] * face_E + lambda_values[q] * trace(face_E) * identity;

                    Tensor<2, dim> neighbor_grad_u = Tensors::get_grad_u<dim>(q, neighbor_grad_u_values) - Tensors::get_grad_u<dim>(q, prev_timestep_neighbor_grad_u_values);
                    Tensor<2, dim> neighbor_E = 0.5 * (neighbor_grad_u + transpose(neighbor_grad_u));
                    Tensor<2, dim> neighbor_sigma = 2 * mu_values[q] * neighbor_E + lambda_values[q] * trace(neighbor_E) * identity;

                    const Tensor<1, dim> &n = fe_face_p.normal_vector(q);

                    Tensor<1, dim> dum = (face_sigma + biot_alpha * (face_p_values[q] - prev_timestep_face_p_values[q]) * identity) * n - (neighbor_sigma + biot_alpha * (neighbor_p_values[q] - prev_timestep_neighbor_p_values[q]) * identity) * n;
                    eta_e_partial_sigma += dum.norm_square() * fe_face_p.JxW(q);

                    face_grad_u = Tensors::get_grad_u<dim>(q, face_grad_u_values);
                    face_E = 0.5 * (face_grad_u + transpose(face_grad_u));
                    face_sigma = 2 * mu_values[q] * face_E + lambda_values[q] * trace(face_E) * identity;

                    neighbor_grad_u = Tensors::get_grad_u<dim>(q, neighbor_grad_u_values);
                    neighbor_E = 0.5 * (neighbor_grad_u + transpose(neighbor_grad_u));
                    neighbor_sigma = 2 * mu_values[q] * neighbor_E + lambda_values[q] * trace(neighbor_E) * identity;
                    Tensor<1, dim> dum6 = (face_sigma + biot_alpha * face_p_values[q] * identity) * n - (neighbor_sigma + biot_alpha * neighbor_p_values[q] * identity) * n;
                    // Tensor<1, dim> dum6 = (face_sigma - biot_alpha * face_p_values[q] * identity) * n - (neighbor_sigma - biot_alpha * neighbor_p_values[q] * identity) * n;
                    /*debug

                    Tensor<1,dim> debug1 =(face_sigma + biot_alpha * face_p_values[q] * identity) * n;
                    Tensor<1,dim> debug2 =(face_sigma - biot_alpha * face_p_values[q] * identity) * n;
                    Tensor<1,dim> debug3 = biot_alpha*(face_p_values[q] - neighbor_p_values[q])*identity*n;
                    cout << "face_sigma = " << face_sigma << "neighbor_sigma = " << neighbor_sigma <<endl;
                    cout<< "dum6=" <<dum6<< " dum7="<<dum7<< endl;
                    cout<< "debug1=" <<debug1<< " debug2="<<debug2<<"debug3=" <<debug3<< endl;
                    */
                    eta_e_sigma += dum6.norm_square() * fe_face_p.JxW(q);
                    
                }
            }
        }
    }

    eta_e_partial_sigma = sqrt(h * eta_e_partial_sigma);
    eta_face_partial_sigma_n.push_back(eta_e_partial_sigma);
    cout << "eta_e_partial_sigma = " << eta_e_partial_sigma << endl;
    double dum2 = 0;
    for (auto &n : eta_face_partial_sigma_n)
    {
        dum2 += n;
    }
    eta_face_partial_sigma.push_back(dum2 * dum2);

    eta_e_sigma *= h;
    cout << "eta_e_sigma = " << eta_e_sigma << endl;
    if (timestep == 1)
    {
        eta_face_sigma.push_back(eta_e_sigma);
    }
    else
    {
        eta_face_sigma.push_back(eta_face_sigma_n.back() + eta_e_sigma);
    }
    eta_face_sigma_n.push_back(eta_e_sigma);

    /* calculate eta_partial_u_n, eta_partial_u, eta_u */

    QGauss<dim> quadrature(fe_pressure.degree + 1);
    FEValues<dim> fe_value_pressure(fe_pressure, quadrature,
                                    update_values | update_quadrature_points | update_gradients | update_hessians | update_JxW_values);
    FEValues<dim> fe_value_displacement(fe_displacement, quadrature,
                                        update_values | update_quadrature_points | update_gradients | update_hessians | update_JxW_values);
    cell = dof_handler_pressure.begin_active();

    cell_u = dof_handler_displacement.begin_active();

    typename DoFHandler<dim>::active_cell_iterator cell_output = dof_handler_output.begin_active();
    const unsigned int n_q_points = quadrature.size();
    vector<Tensor<1, dim>> grad_p_values(n_q_points);
    vector<Tensor<1, dim>> prev_timestep_grad_p_values(n_q_points);
    vector<vector<Tensor<2, dim>>> hessian_u_values(n_q_points, vector<Tensor<2, dim>>(dim));
    vector<vector<Tensor<2, dim>>> prev_timestep_hessian_u_values(n_q_points, vector<Tensor<2, dim>>(dim));

    vector<double> lambda_values(n_q_points);
    vector<double> mu_values(n_q_points);

    vector<types::global_dof_index> output_dofs(dof_handler_output.get_fe().dofs_per_cell);
    cell_eta_E_partial_u = 0;
    cell_eta_E_u = 0;
    double eta_E_partial_u = 0;
    double eta_E_u = 0;
    for (; cell != endc; ++cell, ++cell_u, ++cell_output)
    {
        fe_value_pressure.reinit(cell);
        fe_value_displacement.reinit(cell_u);

        fe_value_pressure.get_function_gradients(solution_pressure, grad_p_values);
        fe_value_pressure.get_function_gradients(prev_timestep_sol_pressure, prev_timestep_grad_p_values);

        fe_value_displacement.get_function_hessians(solution_displacement, hessian_u_values);
        fe_value_displacement.get_function_hessians(prev_timestep_sol_displacement, prev_timestep_hessian_u_values);

        lambda.value_list(fe_value_displacement.get_quadrature_points(), lambda_values);
        mu.value_list(fe_value_displacement.get_quadrature_points(), mu_values);

        cell_output->get_dof_indices(output_dofs);
        for (unsigned int q = 0; q < n_q_points; q++)
        {
            Tensor<1, dim> dum3;
            dum3[0] = (2 * mu_values[q] + lambda_values[q]) * (hessian_u_values[q][0][0][0] - prev_timestep_hessian_u_values[q][0][0][0]) 
                    + (lambda_values[q] + mu_values[q]) * (hessian_u_values[q][1][0][1] - prev_timestep_hessian_u_values[q][1][0][1]) 
                    + mu_values[q] * (hessian_u_values[q][0][1][1] - prev_timestep_hessian_u_values[q][0][1][1]) 
                    - biot_alpha * (grad_p_values[q][0] - prev_timestep_grad_p_values[q][0]);
            dum3[1] = (2 * mu_values[q] + lambda_values[q]) * (hessian_u_values[q][1][1][1] - prev_timestep_hessian_u_values[q][1][1][1]) 
                    + (lambda_values[q] + mu_values[q]) * (hessian_u_values[q][0][0][1] - prev_timestep_hessian_u_values[q][0][0][1]) 
                    + mu_values[q] * (hessian_u_values[q][1][0][0] - prev_timestep_hessian_u_values[q][1][0][0]) 
                    - biot_alpha * (grad_p_values[q][1] - prev_timestep_grad_p_values[q][1]);

            eta_E_partial_u += dum3.norm_square() * fe_value_displacement.JxW(q);
            cell_eta_E_partial_u[output_dofs[0]] += dum3.norm_square() * fe_value_displacement.JxW(q);
            Tensor<1, dim> dum5;
            dum5[0] = (2 * mu_values[q] + lambda_values[q]) * hessian_u_values[q][0][0][0] + (lambda_values[q] + mu_values[q]) * hessian_u_values[q][1][0][1] + mu_values[q] * hessian_u_values[q][0][1][1] - biot_alpha * grad_p_values[q][0];
            dum5[1] = (2 * mu_values[q] + lambda_values[q]) * hessian_u_values[q][1][1][1] + (lambda_values[q] + mu_values[q]) * hessian_u_values[q][0][0][1] + mu_values[q] * hessian_u_values[q][1][0][0] - biot_alpha * grad_p_values[q][1];
            eta_E_u += dum5.norm_square() * fe_value_displacement.JxW(q);
            cell_eta_E_u[output_dofs[0]] += dum5.norm_square() * fe_value_displacement.JxW(q);
        }
    }

    eta_E_partial_u = sqrt(eta_E_partial_u) * h;
    eta_partial_u_n.push_back(eta_E_partial_u);
    double dum4 = 0;
    for (auto &n : eta_partial_u_n)
    {
        dum4 += n;
    }
    eta_partial_u.push_back(dum4 * dum4);

    eta_E_u = eta_E_u * h * h;
    if (timestep == 1)
    {
        eta_u.push_back(eta_E_u);
    }
    else
    {
        eta_u.push_back(eta_u_n.back() + eta_E_u);
    }
    eta_u_n.push_back(eta_E_u);

    u_indicators_table.add_value("time", t);
    u_indicators_table.add_value("eta_face_partial_sigma", eta_face_partial_sigma.back());
    u_indicators_table.add_value("eta_face_sigma", eta_face_sigma.back());
    u_indicators_table.add_value("eta_partial_u", eta_partial_u.back());
    u_indicators_table.add_value("eta_u", eta_u.back());
    u_indicators_table.add_value("sum", eta_face_partial_sigma.back() + eta_face_sigma.back() + eta_partial_u.back()+eta_u.back());
    u_indicators_table.add_value("error", l2_error_u.back());
    u_indicators_table.add_value("eff", (eta_face_partial_sigma.back() + eta_face_sigma.back() + eta_partial_u.back()+eta_u.back())/ l2_error_u.back() );

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler_output);
    data_out.add_data_vector(cell_eta_E_partial_u, "eta_E_partial_u", DataOut<dim>::type_dof_data);
    data_out.add_data_vector(cell_eta_E_u, "eta_E_u", DataOut<dim>::type_dof_data);
    data_out.build_patches();
    ofstream output("output/indicators-u" + to_string(timestep) + ".vtk");
    data_out.write_vtk(output);
}