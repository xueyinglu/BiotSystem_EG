#include "BiotSystem.h"

#include "PressureSolution.h"
#include "DisplacementSolution.h"
using namespace std;
void BiotSystem::calc_error()
{
    QGauss<dim> quadrature_pressure(degree +2);
    FEValues<dim> fe_values_pressure(fe_pressure, quadrature_pressure,
                                    update_values | update_quadrature_points | update_JxW_values | update_gradients);
    const unsigned int n_q_points = quadrature_pressure.size();
    vector<Vector<double>> sol_p_values(n_q_points, Vector<double>(2));
    vector<double> exact_p_values(n_q_points);
    PressureSolution exact_p(t);
    double L2_p_EG = 0;
     typename DoFHandler<dim>::active_cell_iterator
        cell = dof_handler_pressure.begin_active(),
        endc = dof_handler_pressure.end();

    for (; cell != endc; ++cell)
        if (cell->is_locally_owned())
        {

            fe_values_pressure.reinit(cell);

            //exact_pressure_gradient.vector_value_list(fe_values_pressure.get_quadrature_points(), exact_pressure_grads);
            exact_p.value_list(fe_values_pressure.get_quadrature_points(), exact_p_values);

            fe_values_pressure.get_function_values(solution_pressure, sol_p_values);
            //fe_values_pressure.get_function_gradients(solution_pressure, solution_grads_pressure);

            for (unsigned int q = 0; q < n_q_points; ++q)
            {
                //local_error_CG = exact_pressure_values[q] - solution_values_pressure[q][0];
                //L2_norm_CG += local_error_CG * local_error_CG * fe_values_pressure.JxW(q);

                //local_error_DG = exact_pressure_values[q] - solution_values_pressure[q][1];
                //L2_norm_DG += local_error_DG * local_error_DG * fe_values_pressure.JxW(q);

                L2_p_EG += (exact_p_values[q] - sol_p_values[q][1] - sol_p_values[q][0])
                *(exact_p_values[q] - sol_p_values[q][1] - sol_p_values[q][0])
                * fe_values_pressure.JxW(q);


            }
        }
    L2_p_EG = sqrt(L2_p_EG);
    /*
    Vector<float> difference_per_cell_pressure(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler_pressure,
                                      solution_pressure,
                                      PressureSolution(t),
                                      difference_per_cell_pressure,
                                      QGauss<dim>(fe_pressure.degree + 2),
                                      VectorTools::L2_norm);
    double L2_norm_pressure = difference_per_cell_pressure.l2_norm();
        //VectorTools::compute_global_error(triangulation,
        //                                  difference_per_cell_pressure,
        //                                  VectorTools::L2_norm);
    */
    /* Calculate the L2 norm of displacement solution */

    Vector<float> difference_per_cell_displacement(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler_displacement,
                                      solution_displacement,
                                      DisplacementSolution(t),
                                      difference_per_cell_displacement,
                                      QGauss<dim>(fe_displacement.degree + 2),
                                      VectorTools::L2_norm);
    double L2_norm_displacement = difference_per_cell_displacement.l2_norm();

    /*
    QGauss<dim> quadrature_displacement(fe_displacement.degree + 2);
    FEValues<dim> fe_value_displacement(fe_displacement,
                                        quadrature_displacement, update_values | update_quadrature_points | update_gradients | update_JxW_values);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_displacement.begin_active(),
                                                   endc = dof_handler_displacement.end();
    const FEValuesExtractors::Vector displacements (0);
    const unsigned int n_q_points = quadrature_displacement.size();
    vector<Vector<double>> sol_u_values(n_q_points,Vector<double> (dim));
    vector<Tensor<1,dim>> ana_u_values(n_q_points);
    DisplacementSolution disp(t);

    double L2_norm_displacement = 0;
    for(; cell!= endc; ++cell){
        fe_value_displacement.reinit(cell);
        fe_value_displacement.get_function_values(solution_displacement, sol_u_values);
        disp.value_list(fe_value_displacement.get_quadrature_points(), ana_u_values);

        for (unsigned int q =0; q<n_q_points; q++){
            Tensor<1, dim> sol;
            sol[0] =sol_u_values[q][0];
            sol[1] = sol_u_values[q][1];
            Tensor<1,dim> cell_difference= sol - ana_u_values[q];
            L2_norm_displacement += cell_difference.norm_square() *fe_value_displacement.JxW(q);
        }
    }
    */

    convergence_table.add_value("time", t);
    convergence_table.add_value("1/h", 1./h);
    convergence_table.add_value("L2_p", L2_p_EG);
    convergence_table.add_value("L2_u", L2_norm_displacement);
    double energy_norm = calc_u_energy_norm();
    energy_norm = sqrt(energy_norm);
    convergence_table.add_value("energy_u", energy_norm);
    l2_error_p.push_back(L2_p_EG);
    l2_error_u.push_back(L2_norm_displacement);
    energy_error_u.push_back(energy_norm);
}