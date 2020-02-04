#include "BiotSystem.h"
#include <fstream>
#include <iostream>
using namespace std;

void BiotSystem::calc_a_posteriori_indicators_p_eg()
{
    /************************** calculate eta_alg ****************************/
    double dum1 = 0;
    double dum2 = 0;
    for (auto &n : eta_fs)
    {
        dum1 += n;
        dum2 += n * n;
    }
    // for EG flow
    eta_alg.push_back(del_t * dum1 * dum1 + h * h * dum2);
    // for CG flow
    // eta_alg.push_back(del_t * dum1 * dum1);

    /*************************** calculate integrals on the elements ***************************/
    // see notes for details of these local notations
    double eta_t_p_n = 0; // = Del_t /3 * \| k^{1/2} \nabla (p_h^{n,l} - p_h^{n-1})\|_{L2(\Omega)}
    double eta_E_p_n = 0;
    QGauss<dim> quadrature_pressure(degree + 2);
    FEValues<dim> fe_value_pressure(fe_pressure, quadrature_pressure,
                                    update_values | update_quadrature_points | update_gradients | update_hessians | update_JxW_values);
    
    QGauss<dim> quadrature_displacement(degree + 2);
    FEValues<dim> fe_value_displacement(fe_displacement, quadrature_displacement,
                                        update_values | update_quadrature_points | update_gradients | update_JxW_values);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_pressure.begin_active(),
                                                   endc = dof_handler_pressure.end();
    typename DoFHandler<dim>::active_cell_iterator
        cell_displacement = dof_handler_displacement.begin_active();
    typename DoFHandler<dim>::active_cell_iterator cell_output = dof_handler_output.begin_active();
    const unsigned int n_q_points = quadrature_pressure.size();

    vector<vector<Tensor<1, dim>>> grad_p_values(n_q_points, vector<Tensor<1, dim>>(2));
    vector<vector<Tensor<1, dim>>> prev_timestep_grad_p_values(n_q_points, vector<Tensor<1, dim>>(2));
    vector<double> permeability_values(n_q_points);
    
    vector<vector<Tensor<1, dim>>> prev_timestep_grad_u_values(n_q_points, vector<Tensor<1, dim>>(dim));
    vector<vector<Tensor<1, dim>>> grad_u_values(n_q_points, vector<Tensor<1, dim>>(dim));

    vector<Vector<double>> p_values(n_q_points, Vector<double> (2));
    vector<Vector<double>> prev_timestep_p_values(n_q_points, Vector<double> (2));
    vector<Vector<double>> laplacian_p_values(n_q_points, Vector<double> (2));

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
            Tensor<1, dim> cell_difference = grad_p_values[q][0]+ grad_p_values[q][1] - prev_timestep_grad_p_values[q][0]-prev_timestep_grad_p_values[q][1] ;
            eta_t_p_n += permeability_values[q] *
                         cell_difference.norm_square() *
                         fe_value_pressure.JxW(q);
            cell_eta_time[output_dofs[0]] += permeability_values[q] *
                                             cell_difference.norm_square() * fe_value_pressure.JxW(q);

            double cell_residual = (1 / mu_f * permeability_values[q] * (laplacian_p_values[q][0] + laplacian_p_values[q][1] )
                                    -  biot_inv_M / del_t * (p_values[q][0]+p_values[q][1] - prev_timestep_p_values[q][0] - prev_timestep_p_values[q][1]) 
                                    - biot_alpha / del_t * (grad_u_values[q][0][0] + grad_u_values[q][1][1] - prev_timestep_grad_u_values[q][0][0] - prev_timestep_grad_u_values[q][1][1]));
            // cout << "cell residual = " << cell_residual << endl;
            eta_E_p_n += cell_residual * cell_residual * fe_value_pressure.JxW(q);
            cell_eta_E_p[output_dofs[0]] += cell_residual * cell_residual * fe_value_pressure.JxW(q);
        }
    }
    eta_t_p_n = eta_t_p_n * del_t / 3;
    eta_E_p_n = eta_E_p_n * h * h * del_t;
   
   /*************************** calculate integrals on the edges ***************************/ 
    // see notes for the formula
    double eta_pen_n = 0;
    double eta_t_J_n = 0; 
    double eta_partial_p_J_n = 0;
    double eta_p_J_n =0;
    double eta_flux_e_n = 0;
    QGauss<dim - 1> face_quadrature(fe_pressure.degree + 1);
    FEFaceValues<dim> fe_face_p(fe_pressure, face_quadrature,
                                update_values | update_normal_vectors | update_gradients | update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_neighbor_p(fe_pressure, face_quadrature,
                                         update_values | update_gradients | update_quadrature_points | update_JxW_values);
    cell = dof_handler_pressure.begin_active();
    for (; cell!= endc; ++cell)
    {
        for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
        {
            typename DoFHandler<dim>::face_iterator face_p = cell->face(face_no);
            if (!face_p->at_boundary())
            {
                Assert(cell->neighbor(face_no).state() == IteratorState::valid, ExcInternalError());
                const typename DoFHandler<dim>::cell_iterator neighbor_p = cell->neighbor(face_no);
                vector<Vector<double>> face_p_values(fe_face_p.n_quadrature_points, Vector<double>(2));
                vector<Vector<double>> neighbor_p_values(fe_face_neighbor_p.n_quadrature_points, Vector<double>(2));
                vector<Vector<double>> prev_timestep_face_p_values(fe_face_p.n_quadrature_points, Vector<double>(2));
                vector<Vector<double>> prev_timestep_neighbor_p_values(fe_face_neighbor_p.n_quadrature_points, Vector<double>(2));
                vector<vector<Tensor<1, dim>>> face_grad_p_values(fe_face_p.n_quadrature_points, vector<Tensor<1, dim>>(2));
                vector<vector<Tensor<1, dim>>> neighbor_grad_p_values(fe_face_neighbor_p.n_quadrature_points, vector<Tensor<1, dim>>(2));
                const unsigned int neighbor_face_p = cell->neighbor_of_neighbor(face_no);
                fe_face_p.reinit(cell, face_no);
                fe_face_neighbor_p.reinit(neighbor_p, neighbor_face_p);
                fe_face_p.get_function_values(solution_pressure, face_p_values);
                fe_face_p.get_function_values(prev_timestep_sol_pressure, prev_timestep_face_p_values);
                fe_face_neighbor_p.get_function_values(solution_pressure, neighbor_p_values);
                fe_face_neighbor_p.get_function_values(prev_timestep_sol_pressure, prev_timestep_neighbor_p_values);
                fe_face_p.get_function_gradients(solution_pressure, face_grad_p_values);
                fe_face_neighbor_p.get_function_gradients(solution_pressure, neighbor_grad_p_values);

                for (unsigned int q = 0; q < fe_face_p.n_quadrature_points; q++)
                {
                    const Tensor<1, dim> &n = fe_face_p.normal_vector(q);
                    double jump_t = (face_p_values[q][1]-prev_timestep_face_p_values[q][1])
                                -(neighbor_p_values[q][1]-prev_timestep_neighbor_p_values[q][1]);
                    double jump =face_p_values[q][1] - neighbor_p_values[q][1] ;
                    eta_t_J_n += jump_t*jump_t*fe_face_p.JxW(q);
                    eta_pen_n += jump*jump*fe_face_p.JxW(q);
                    eta_partial_p_J_n += jump_t*jump_t*fe_face_p.JxW(q);
                    eta_p_J_n+= jump*jump*fe_face_p.JxW(q); 
                    double flux_jump = permeability.value(fe_face_p.quadrature_point(q),0) * ((face_grad_p_values[q][0] + face_grad_p_values[q][1])  * n 
                                    - (neighbor_grad_p_values[q][0] + neighbor_grad_p_values[q][1] )  * n);
                    eta_flux_e_n += flux_jump * flux_jump * fe_face_p.JxW(q);
                    //TODO add the face integrals to the visualization
                }

            }
        }

    }
    eta_t_J_n = del_t/3*gamma_penal/min_cell_diameter*eta_t_J_n;
    eta_pen_n = del_t * gamma_penal/min_cell_diameter * eta_pen_n;
    eta_partial_p_J_n = sqrt(gamma_penal * min_cell_diameter * eta_partial_p_J_n);
    eta_p_J_n = gamma_penal * min_cell_diameter * eta_p_J_n;
    eta_flux_e_n = eta_flux_e_n * h * del_t;


    if (timestep == 1)
    {
        eta_time.push_back(eta_t_p_n+eta_t_J_n);
    }
    else
    {
        eta_time.push_back(eta_time.back() + eta_t_p_n +eta_t_J_n);
    }

    eta_pen.push_back(eta_pen_n);
    eta_partial_p_J.push_back(eta_partial_p_J_n);
    eta_p_J.push_back(eta_p_J_n);


    if (timestep == 1)
    {
        eta_flow.push_back(eta_E_p_n + eta_flux_e_n);
    }
    else
    {
        eta_flow.push_back(eta_flow.back() + eta_E_p_n + eta_flux_e_n);
    }

    double dum3 = 0;
    double dum4 = 0;
    for (auto &n : eta_pen)
    {
        dum3 += n;
    }

    for (auto &n : eta_partial_p_J){
        dum4 += n;
    }
    if(timestep == 1){
        eta_jump.push_back(dum3 + dum4*dum4 + eta_p_J.back());
    }
    else{
        eta_jump.push_back(dum3 + dum4*dum4 + eta_p_J.end()[-1]+ eta_p_J.end()[-2] );
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
    p_indicators_table.add_value("eta_jump", eta_jump.back());
    // p_indicators_table.add_value("eta_p_residual", eta_p_residual.back());
    // p_indicators_table.add_value("eta_flux_jump", eta_flux_jump.back());
    p_indicators_table.add_value("sum", eta_alg.back() + eta_time.back() +eta_flow.back() + eta_jump.back());
    p_indicators_table.add_value("error", l2_error_p.back());
    // p_indicators_table.add_value("eff", (eta_alg.back() + eta_time.back() +eta_flow.back())/ l2_error_p.back());
}