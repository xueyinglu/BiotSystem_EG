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
    vector<vector<Tensor<1, dim>>> grad_p_values(n_q_points, vector<Tensor<1, dim>>(2));
    vector<double> permeability_values(n_q_points);
    for (; cell != endc; ++cell)
    {
        fe_value_pressure.reinit(cell);
        fe_value_pressure.get_function_gradients(solution_pressure, grad_p_values);
        permeability.value_list(fe_value_pressure.get_quadrature_points(), permeability_values);
        if (test_case == TestCase::heterogeneous)
        {
            perm_function.value_list(fe_value_pressure.get_quadrature_points(), permeability_values);
        }
        for (unsigned int q = 0; q < n_q_points; q++)
        {
            Tensor<1, dim> true_grad_p;
            PressureSolution(t).gradient_value(fe_value_pressure.quadrature_point(q), true_grad_p);
            h_norm_p_sq += permeability_values[q] * (grad_p_values[q][0] + grad_p_values[q][1] - true_grad_p).norm_square() * fe_value_pressure.JxW(q);
        }
    }

    QGauss<dim - 1> face_quadrature(degree + 2);
    FEFaceValues<dim> fe_face_p(fe_pressure, face_quadrature,
                                update_values | update_normal_vectors | update_gradients | update_quadrature_points | update_JxW_values);
    FESubfaceValues<dim> fe_subface_p(fe_pressure, face_quadrature,
                                      update_values | update_normal_vectors | update_gradients | update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_neighbor_p(fe_pressure, face_quadrature,
                                         update_values | update_gradients | update_quadrature_points | update_JxW_values);

    cell = dof_handler_pressure.begin_active();

    for (; cell != endc; ++cell)
    {
        for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
        {
            typename DoFHandler<dim>::face_iterator face_p = cell->face(face_no);
            if (!face_p->at_boundary())
            {
                Assert(cell->neighbor(face_no).state() == IteratorState::valid, ExcInternalError());
                if (cell->face(face_no)->has_children())
                {
                    const unsigned int neighbor2 = cell->neighbor_face_no(face_no);

                    for (unsigned int subface_no = 0;
                         subface_no < cell->face(face_no)->number_of_children();
                         ++subface_no)
                    {
                        typename DoFHandler<dim>::cell_iterator neighbor_child = cell->neighbor_child_on_subface(face_no, subface_no);

                        Assert(!neighbor_child->has_children(), ExcInternalError());

                        fe_subface_p.reinit(cell, face_no, subface_no);
                        fe_face_neighbor_p.reinit(neighbor_child, neighbor2);
                        vector<Vector<double>> face_p_values(fe_subface_p.n_quadrature_points, Vector<double>(2));
                        vector<Vector<double>> neighbor_p_values(fe_face_neighbor_p.n_quadrature_points, Vector<double>(2));
                        fe_subface_p.get_function_values(solution_pressure, face_p_values);
                        fe_face_neighbor_p.get_function_values(solution_pressure, neighbor_p_values);
                        for (unsigned int q = 0; q < fe_subface_p.n_quadrature_points; q++)
                        {
                            double jump = face_p_values[q][1] - neighbor_p_values[q][1];
                            h_norm_p_sq+= jump *jump *fe_subface_p.JxW(q);
                        }
                    }
                }

                else if (!cell->neighbor_is_coarser(face_no))
                {
                    const typename DoFHandler<dim>::cell_iterator neighbor_p = cell->neighbor(face_no);
                    vector<Vector<double>> face_p_values(fe_face_p.n_quadrature_points, Vector<double>(2));
                    vector<Vector<double>> neighbor_p_values(fe_face_neighbor_p.n_quadrature_points, Vector<double>(2));
                    const unsigned int neighbor_face_p = cell->neighbor_of_neighbor(face_no);
                    fe_face_p.reinit(cell, face_no);
                    fe_face_neighbor_p.reinit(neighbor_p, neighbor_face_p);
                    fe_face_p.get_function_values(solution_pressure, face_p_values);
                    fe_face_neighbor_p.get_function_values(solution_pressure, neighbor_p_values);
                    for (unsigned int q = 0; q < fe_face_p.n_quadrature_points; q++)
                    {
                        double jump = face_p_values[q][1] - neighbor_p_values[q][1];
                        h_norm_p_sq+= jump *jump *fe_face_p.JxW(q);
                    }
                }
                else
                {
                    const auto neighbor = cell->neighbor(face_no);
                    std::pair<unsigned int, unsigned int> neighbor_face_subface = cell->neighbor_of_coarser_neighbor(face_no);

                    Assert(neighbor_face_subface.first < GeometryInfo<dim>::faces_per_cell, ExcInternalError());
                    Assert(neighbor_face_subface.second < neighbor->face(neighbor_face_subface.first)->number_of_children(), ExcInternalError());
                    Assert(neighbor->neighbor_child_on_subface(neighbor_face_subface.first, neighbor_face_subface.second) == cell, ExcInternalError());

                    fe_face_p.reinit(cell, face_no);

                    fe_subface_p.reinit(neighbor, neighbor_face_subface.first,
                                        neighbor_face_subface.second);
                    vector<Vector<double>> face_p_values(fe_face_p.n_quadrature_points, Vector<double>(2));
                    vector<Vector<double>> neighbor_p_values(fe_subface_p.n_quadrature_points, Vector<double>(2));
                    fe_face_p.get_function_values(solution_pressure, face_p_values);
                    fe_subface_p.get_function_values(solution_pressure, neighbor_p_values);
                    for (unsigned int q = 0; q < fe_face_p.n_quadrature_points; q++)
                    {
                        double jump = face_p_values[q][1] - neighbor_p_values[q][1];
                        h_norm_p_sq+= jump *jump *fe_face_p.JxW(q);
                    }
                }
            }
        }
    }

    h_error_p_sq.push_back(h_norm_p_sq);
}