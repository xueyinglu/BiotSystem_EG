#include "BiotSystem.h"
#include "DisplacementSolution.h"
#include "InitialPressure.h"
#include "AuxTools.h"
using namespace std;
void BiotSystem::assemble_system_displacement()
{
    system_matrix_displacement.reinit(sparsity_pattern_displacement);
    system_rhs_displacement.reinit(dof_handler_displacement.n_dofs());
    QGauss<dim> quadrature_formula(fe_displacement.degree + 2);
    QGauss<dim - 1> face_quadrature_formula(fe_displacement.degree + 2);

    FEValues<dim> fe_values(fe_displacement, quadrature_formula,
                            update_values | update_gradients |
                                update_quadrature_points | update_JxW_values);

    FEFaceValues<dim> fe_face_values(fe_displacement, face_quadrature_formula,
                                 update_values | update_gradients |
                                     update_quadrature_points | update_JxW_values);

    FEValues<dim> fe_values_pressure(fe_pressure, quadrature_formula,
                                     update_values | update_quadrature_points |
                                     update_JxW_values | update_gradients);

    const unsigned int dofs_per_cell = fe_displacement.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    std::vector<double> lambda_values(n_q_points);
    std::vector<double> mu_values(n_q_points);
    std::vector<Vector<double>> rhs_values(n_q_points,
                                           Vector<double>(dim));
    //std::vector<Vector<double>> grad_p_values(n_q_points,
    //                                       Vector<double>(dim));
    std::vector<Vector<double>> pore_pressure_values_eg(n_q_points, Vector<double>(2));
    // std::vector<Tensor<1, dim>> grad_p_values(n_q_points);
    Tensor<2, dim> identity = Tensors::get_Identity<dim>();

    InitialPressure initial_pressure;
    const FEValuesExtractors::Vector displacements(0);
    // Now we can begin with the loop over all cells:
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_displacement.begin_active(),
                                                   endc = dof_handler_displacement.end();

    typename DoFHandler<dim>::active_cell_iterator
        cell_pressure = dof_handler_pressure.begin_active();

    for (; cell != endc; ++cell, ++cell_pressure)
    {
        cell_matrix = 0;
        cell_rhs = 0;

        fe_values.reinit(cell);
        fe_values_pressure.reinit(cell_pressure);

        lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
        mu.value_list(fe_values.get_quadrature_points(), mu_values);
        if (test_case == heterogeneous)
        {
            lambda_function.value_list(fe_values.get_quadrature_points(), lambda_values);
            mu_function.value_list(fe_values.get_quadrature_points(), mu_values);
        }
        right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                                          rhs_values);

        fe_values_pressure.get_function_values(solution_pressure, pore_pressure_values_eg);

        // fe_values_pressure.get_function_gradients(solution_pressure, grad_p_values);
        // initial_pressure.value_list(fe_values_pressure.get_quadrature_points(), pore_pressure_values);

        // Assemble the cell matrix as in elasticity_cg
        for (unsigned int q = 0; q < n_q_points; ++q)
        {
            std::vector<Tensor<1, dim>> phi_i_u(dofs_per_cell);
            std::vector<Tensor<2, dim>> phi_i_grads_u(dofs_per_cell);
            std::vector<Tensor<2, dim>> E_phi(dofs_per_cell);
            std::vector<Tensor<2, dim>> sigma_phi(dofs_per_cell);
            // cout <<"fe_values at quadrature_point q = " << fe_values_pressure.quadrature_point(q) <<endl;
            // Compute and store desired quantities
            for (unsigned int k = 0; k < dofs_per_cell; ++k)
            {
                phi_i_u[k] = fe_values[displacements].value(k, q);
                phi_i_grads_u[k] = fe_values[displacements].gradient(k, q);
                E_phi[k] = 0.5 * (phi_i_grads_u[k] + transpose(phi_i_grads_u[k]));
                sigma_phi[k] = 2.0 * mu_values[q] * E_phi[k] + lambda_values[q] * trace(E_phi[k]) * identity;
            }

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                { // assemble cell level matrix
                    cell_matrix(j, i) += fe_values.JxW(q) * scalar_product(sigma_phi[i], E_phi[j]);
                }
                // assemble cell level rhs as in elasticity_cg
                cell_rhs(i) += biot_alpha * (pore_pressure_values_eg[q][0] + pore_pressure_values_eg[q][1]) * trace(phi_i_grads_u[i]) * fe_values.JxW(q);
                // cell_rhs(i) += 0.75 * pore_pressure_values[q] * trace(phi_i_grads_u[i]) *fe_values.JxW(q);
                // cell_rhs(i) -= biot_alpha * (grad_p_values[q]*phi_i_u[i])*fe_values.JxW(q);
            }

        } // end q_point

        /**************************** Neumann BC -- Traction BC *******************************/
        // Apply traction BC on the bottom (y=0)
        if (test_case == TestCase::terzaghi || test_case == TestCase::heterogeneous)
        {   
            for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
            {   
                if (cell->face(face)->at_boundary() &&
                    (cell->face(face)->boundary_id() == 2))
                {

                    fe_face_values.reinit(cell, face);
                    for (unsigned int q = 0; q < n_face_q_points; ++q)
                    {
                        for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                            cell_rhs(i) += traction_bc * fe_face_values[displacements].value(i, q) * fe_face_values.JxW(q);
                        }
                    }
                }
            }
        }

        cell->get_dof_indices(local_dof_indices);
        constraints_displacement.distribute_local_to_global(cell_matrix, local_dof_indices, system_matrix_displacement);
        constraints_displacement.distribute_local_to_global(cell_rhs, local_dof_indices, system_rhs_displacement);
    }

    /**************************** Dirichlet BC  *******************************/
    if (test_case == TestCase::benchmark)
    {
        /*
    vector<bool> component_mask;
    component_mask.push_back(false);
    component_mask.push_back(true);
    std::map<types::global_dof_index, double> boundary_values; 
    
    VectorTools::interpolate_boundary_values(dof_handler_displacement,
                                             1,
                                             ZeroFunction<dim>(dim),
                                             boundary_values,
                                             ComponentMask(component_mask));

    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix_displacement,
                                       solution_displacement,
                                       system_rhs_displacement);
    component_mask[0] = true;
    component_mask[1] = false;
    VectorTools::interpolate_boundary_values(dof_handler_displacement,
                                             0,
                                             ZeroFunction<dim>(dim),
                                             boundary_values,
                                             ComponentMask(component_mask));
    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix_displacement,
                                       solution_displacement,
                                       system_rhs_displacement);
    */
        // Apply Dirichlet BC of the true solution
        std::map<types::global_dof_index, double> boundary_values;
        VectorTools::interpolate_boundary_values(dof_handler_displacement,
                                                 0,
                                                 DisplacementSolution(t),
                                                 //ZeroFunction<dim>(dim),
                                                 boundary_values);
        MatrixTools::apply_boundary_values(boundary_values,
                                           system_matrix_displacement,
                                           solution_displacement,
                                           system_rhs_displacement);
    }

    else if (test_case == TestCase::terzaghi || test_case == TestCase::heterogeneous)
    {
        // apply roller BC on left and right boundary
        vector<bool> component_mask;
        component_mask.push_back(true);
        component_mask.push_back(false);
        std::map<types::global_dof_index, double> boundary_values;

        VectorTools::interpolate_boundary_values(dof_handler_displacement,
                                                 1,
                                                 ZeroFunction<dim>(dim),
                                                 boundary_values,
                                                 ComponentMask(component_mask));

        MatrixTools::apply_boundary_values(boundary_values,
                                           system_matrix_displacement,
                                           solution_displacement,
                                           system_rhs_displacement);
        // apply roller BC on the top (y=H)
        component_mask[0] = false;
        component_mask[1] = true;
        VectorTools::interpolate_boundary_values(dof_handler_displacement,
                                                 0,
                                                 ZeroFunction<dim>(dim),
                                                 boundary_values,
                                                 ComponentMask(component_mask));
        MatrixTools::apply_boundary_values(boundary_values,
                                           system_matrix_displacement,
                                           solution_displacement,
                                           system_rhs_displacement);
    }
}
