
#include "BiotSystem.h"
#include "AuxTools.h"
#include "DisplacementSolution.h"
using namespace std;

double BiotSystem::calc_u_energy_norm(){

    QGauss<dim> quadrature_displacement(fe_displacement.degree + 2);
    FEValues<dim> fe_value_displacement(fe_displacement,
                                        quadrature_displacement, update_values | update_quadrature_points | update_gradients | update_JxW_values);
    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_displacement.begin_active(),
                                                   endc = dof_handler_displacement.end();
    const FEValuesExtractors::Vector displacements (0);
    const unsigned int n_q_points = quadrature_displacement.size();
    vector<Vector<double>> sol_u_values(n_q_points,Vector<double> (dim));
    vector<Tensor<1,dim>> ana_u_values(n_q_points);
    vector<vector<Tensor<1,dim>>> grad_u_values(n_q_points, vector<Tensor<1, dim>>(dim));
    DisplacementSolution disp(t);
    vector<double> lambda_values(n_q_points);
    vector<double> mu_values(n_q_points);

    double u_energy_norm_square = 0;
    for(; cell!= endc; ++cell){
        fe_value_displacement.reinit(cell);
        fe_value_displacement.get_function_values(solution_displacement, sol_u_values);
        fe_value_displacement.get_function_gradients(solution_displacement, grad_u_values);
        lambda.value_list(fe_value_displacement.get_quadrature_points(), lambda_values);
        mu.value_list(fe_value_displacement.get_quadrature_points(), mu_values);

        for (unsigned int q =0; q<n_q_points; q++){
            Tensor<2, dim> true_grad;
            disp.gradient_value(fe_value_displacement.quadrature_point(q), true_grad);
            Tensor<2,dim> grad_u = Tensors::get_grad_u<dim>(q, grad_u_values);
            Tensor<2,dim> e_strain = 0.5*((grad_u - true_grad) + transpose(grad_u - true_grad));
            double e_div = Tensors::get_divergence_u(grad_u - true_grad);
            u_energy_norm_square += (mu_values[q]/2*e_strain.norm_square() + lambda_values[q]/4* e_div*e_div) *fe_value_displacement.JxW(q);
        }
    }
    return u_energy_norm_square;
}