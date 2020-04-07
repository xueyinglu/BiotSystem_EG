#include "BiotSystem.h"
#include "AuxTools.h"
using namespace std;
void BiotSystem::assemble_system_pressure_eg()
{
    system_matrix_pressure = 0;
    system_rhs_pressure = 0;
    QGauss<dim> quadrature(degree + 2);
    FEValues<dim> fe_value(fe_pressure,
                           quadrature, update_values | update_quadrature_points | update_gradients | update_JxW_values);
    FEValues<dim> fe_value_displacement(fe_displacement,
                                        quadrature, update_values | update_quadrature_points | update_gradients | update_JxW_values);
    // For face integrals
    QGauss<dim - 1> face_quadrature(degree + 2);
    FEFaceValues<dim> fe_face_values(fe_pressure, face_quadrature,
                                     update_values | update_normal_vectors |
                                         update_gradients |
                                         update_quadrature_points | update_JxW_values);
    FESubfaceValues<dim> fe_subface_values(fe_pressure, face_quadrature,
                                           update_values | update_normal_vectors |
                                               update_gradients |
                                               update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_face_values_neighbor(fe_pressure, face_quadrature,
                                              update_values | update_normal_vectors |
                                                  update_gradients |
                                                  update_quadrature_points | update_JxW_values);
    const unsigned int dofs_per_cell = fe_pressure.dofs_per_cell;
    const unsigned int n_q_points = quadrature.size();
    const unsigned int n_face_q_points = face_quadrature.size();
    // cout << "dofs_per_cell = " << dofs_per_cell << endl;
    // cout << "n_q_points = " << n_q_points << endl;
    // cout << "n_face_q_points = " << n_face_q_points << endl;
    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> cell_matrix_face(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);
    Vector<double> cell_rhs_face(dofs_per_cell);
    vector<unsigned int> local_dof_indices(dofs_per_cell);
    vector<unsigned int> local_dof_indices_face(dofs_per_cell);

    const FEValuesExtractors::Scalar pressure_cg(0);
    const FEValuesExtractors::Scalar pressure_dg(1);
    // Declaring test functions:
    std::vector<double> phi_i_p_cg(dofs_per_cell);
    std::vector<Tensor<1, dim>> phi_i_grads_p_cg(dofs_per_cell);

    std::vector<double> phi_i_p_dg(dofs_per_cell);
    std::vector<Tensor<1, dim>> phi_i_grads_p_dg(dofs_per_cell);

    std::vector<double> phi_i_p_face(dofs_per_cell);
    std::vector<double> phi_i_p_face_neighbor(dofs_per_cell);

    std::vector<Tensor<1, dim>> phi_i_grads_p_face(dofs_per_cell);
    std::vector<Tensor<1, dim>> phi_i_grads_p_face_neighbor(dofs_per_cell);

    std::vector<double> phi_i_p_face_cg(dofs_per_cell);
    std::vector<Tensor<1, dim>> phi_i_grads_p_face_cg(dofs_per_cell);
    std::vector<Tensor<1, dim>> phi_i_grads_p_face_neighbor_cg(dofs_per_cell);

    double penalty_term = gamma_penal / min_cell_diameter;
    double d_Big_K = 1;
    double eps = 1e-5;

    vector<double> permeability_values(n_q_points);
    vector<Vector<double>> prev_timestep_sol_pressure_values(n_q_points, Vector<double>(2));
    vector<Vector<double>> prev_fs_sol_pressure_values(n_q_points, Vector<double>(2));
    vector<vector<Tensor<1, dim>>> prev_timestep_sol_grad_u_values(n_q_points, vector<Tensor<1, dim>>(dim));
    vector<vector<Tensor<1, dim>>> prev_fs_sol_grad_u_values(n_q_points, vector<Tensor<1, dim>>(dim));
    double prev_timestep_mean_stress;
    double prev_fs_mean_stress;

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler_pressure.begin_active(),
                                                   endc = dof_handler_pressure.end();

    typename DoFHandler<dim>::active_cell_iterator
        cell_displacement = dof_handler_displacement.begin_active();
    for (; cell != endc; ++cell, ++cell_displacement)
    {
        cell->get_dof_indices(local_dof_indices);
        fe_value.reinit(cell);
        fe_value_displacement.reinit(cell_displacement);
        cell_matrix = 0;
        cell_rhs = 0;
        /* get the function values at current element */
        permeability.value_list(fe_value.get_quadrature_points(), permeability_values);
        if (test_case == TestCase::heterogeneous)
        {
            perm_function.value_list(fe_value.get_quadrature_points(), permeability_values);
        }
        fe_value.get_function_values(prev_timestep_sol_pressure, prev_timestep_sol_pressure_values);
        fe_value.get_function_values(prev_fs_sol_pressure, prev_fs_sol_pressure_values);
        fe_value_displacement.get_function_gradients(prev_timestep_sol_displacement, prev_timestep_sol_grad_u_values);
        fe_value_displacement.get_function_gradients(prev_fs_sol_displacement, prev_fs_sol_grad_u_values);
        /* assemble cell level matrix and rhs */
        for (unsigned int q = 0; q < n_q_points; q++)
        {
            for (unsigned int k = 0; k < dofs_per_cell; ++k)
            {
                phi_i_p_cg[k] = fe_value[pressure_cg].value(k, q);
                phi_i_grads_p_cg[k] = fe_value[pressure_cg].gradient(k, q);

                phi_i_p_dg[k] = fe_value[pressure_dg].value(k, q);
                phi_i_grads_p_dg[k] = fe_value[pressure_dg].gradient(k, q);
            }

            // calculate the mean stress values at the quadrature point
            const Tensor<2, dim> prev_time_grad_u = Tensors ::get_grad_u<dim>(q, prev_timestep_sol_grad_u_values);
            const double prev_time_div_u = Tensors ::get_divergence_u<dim>(prev_time_grad_u);
            const Tensor<2, dim> prev_fs_grad_u = Tensors ::get_grad_u<dim>(q, prev_fs_sol_grad_u_values);
            const double prev_fs_div_u = Tensors ::get_divergence_u<dim>(prev_fs_grad_u);
            prev_timestep_mean_stress = K_b * prev_time_div_u - biot_alpha * (prev_timestep_sol_pressure_values[q][0] + prev_timestep_sol_pressure_values[q][1]);
            prev_fs_mean_stress = K_b * prev_fs_div_u - biot_alpha * (prev_fs_sol_pressure_values[q][0] + prev_fs_sol_pressure_values[q][1]);
            // xueying: assign perm values
            d_Big_K = permeability.value(fe_value.quadrature_point(q), 0) / mu_f;
            if (test_case == TestCase::heterogeneous)
            {
                d_Big_K = perm_function.value(fe_value.quadrature_point(q), 0) / mu_f;
            }
            for (unsigned int i = 0; i < dofs_per_cell; i++)
            {
                for (unsigned int j = 0; j < dofs_per_cell; j++)
                {
                    // parabolic part
                    // CG part
                    cell_matrix(j, i) +=
                        ((biot_inv_M + biot_alpha * biot_alpha / K_b) / del_t * // (1/M + alpha^2/K_b)/del_t
                         phi_i_p_cg[i] * phi_i_p_cg[j] * fe_value.JxW(q));      // phi(x_q)*phi(x_q) dx
                                                                                // EG part
                    cell_matrix(j, i) +=
                        ((biot_inv_M + biot_alpha * biot_alpha / K_b) / del_t *
                         phi_i_p_dg[i] * phi_i_p_dg[j] * fe_value.JxW(q));
                    cell_matrix(j, i) +=
                        ((biot_inv_M + biot_alpha * biot_alpha / K_b) / del_t *
                         phi_i_p_cg[i] * phi_i_p_dg[j] * fe_value.JxW(q));
                    cell_matrix(j, i) +=
                        ((biot_inv_M + biot_alpha * biot_alpha / K_b) / del_t *
                         phi_i_p_dg[i] * phi_i_p_cg[j] * fe_value.JxW(q));
                    // elliptic part
                    cell_matrix(i, j) +=
                        (1.0 / mu_f * permeability_values[q] * // 1/mu_f * k
                         phi_i_grads_p_cg[i] * phi_i_grads_p_cg[j] *
                         fe_value.JxW(q)); // dx
                    cell_matrix(i, j) +=
                        (1.0 / mu_f * permeability_values[q] * // 1/mu_f * k
                         phi_i_grads_p_dg[i] * phi_i_grads_p_dg[j] *
                         fe_value.JxW(q)); // dx
                }
                // source term
                //cell_rhs(i) +=
                //    (fe_value.shape_value(i, q) * // phi_i(x_q)
                //     1 *                          // f(x_q)
                //     fe_value.JxW(q));            // dx

                if (test_case == TestCase::heterogeneous)
                {

                    // add a well term
                    /* ARMA paper */ /*
                    Point<2> well1 = Point<2>(0.5, 0.25 - 1. / 128);
                    Point<2> well2 = Point<2>(0.5, 0.5 - 1. / 128);
                    Point<2> well3 = Point<2>(0.5, 0.75 - 1. / 128);
                    */
                    Point<2> well1 = Point<2>(0.5, 21. / 64 + 1. / 128);
                    Point<2> well2 = Point<2>(0.5, 42. / 64 + 1. / 128);
                    if (fe_value.quadrature_point(q).distance(well1) < 1. / 128 ||
                        fe_value.quadrature_point(q).distance(well2) < 1. / 128) //||
                    //fe_value.quadrature_point(q).distance(well3) < 1. / 128)
                    {
                        cell_rhs(i) +=
                            (fe_value.shape_value(i, q) * // phi_i(x_q)
                             -20 *                        // f(x_q)
                             fe_value.JxW(q));            // dx
                    }
                }

                // prev time step
                cell_rhs(i) +=
                    ((biot_inv_M + biot_alpha * biot_alpha / K_b) / del_t * // (1/M + alpha^2/K_b)/del_t
                     (prev_timestep_sol_pressure_values[q][0] + prev_timestep_sol_pressure_values[q][1]) *
                     fe_value.shape_value(i, q) * fe_value.JxW(q));

                // change in mean stress
                cell_rhs(i) -=
                    (biot_alpha / K_b / del_t *
                     (prev_fs_mean_stress - prev_timestep_mean_stress) *
                     fe_value.shape_value(i, q) * fe_value.JxW(q));
            }
            // Integrals on the faces
            for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell; ++face_no)
            {

                double d_Big_K_neighbor = 0.;
                cell_matrix_face = 0;

                // On the Boundary
                if (cell->at_boundary(face_no))
                {
                    // Weakly impose Dirichlet BC where the boundary indicator is 0
                    if (
                        ((test_case == TestCase::terzaghi || test_case == TestCase::heterogeneous) && cell->face(face_no)->boundary_id() == 2) || (test_case == TestCase::mandel && cell->face(face_no)->boundary_id() == 1))
                    {
                        // cout << "weakly impose pressure dirichlet bc" << endl;
                        fe_face_values.reinit(cell, face_no);

                        double h_e = cell->face(face_no)->diameter();
                        penalty_term = gamma_penal / h_e;
                        for (unsigned int q = 0; q < n_face_q_points; ++q)
                        {
                            for (unsigned int k = 0; k < dofs_per_cell; ++k)
                            {
                                phi_i_p_face[k] = fe_face_values[pressure_dg].value(k, q);
                                phi_i_grads_p_face[k] = fe_face_values[pressure_dg].gradient(k, q);

                                phi_i_p_face_cg[k] = fe_face_values[pressure_cg].value(k, q);
                                phi_i_grads_p_face_cg[k] = fe_face_values[pressure_cg].gradient(k, q);
                            }

                            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                            {
                                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                                {

                                    if (bCG_WeaklyBD)
                                    {
                                        //(0,0) - CG Weakly imposed BD
                                        // I-1
                                        cell_matrix(i, j) += -d_Big_K * (phi_i_grads_p_face_cg[j]) * fe_face_values.normal_vector(q) * fe_face_values[pressure_cg].value(i, q) * fe_face_values.JxW(q);

                                        //(0,0) - CG Weakly imposed BD
                                        // I-2
                                        cell_matrix(i, j) += -d_SForm * d_Big_K * (phi_i_grads_p_face_cg[i]) * fe_face_values.normal_vector(q) * fe_face_values[pressure_cg].value(j, q) * fe_face_values.JxW(q);

                                        //(0,0) - CG Weakly imposed BD
                                        // I-4
                                        cell_matrix(i, j) += penalty_term * d_Big_K * (phi_i_p_face_cg[j]) * fe_face_values[pressure_cg].value(i, q) * fe_face_values.JxW(q);
                                    }
                                    //DG - Weakly imposed BD
                                    //(1,1)
                                    //These values are ZERO for EG
                                    cell_matrix(i, j) += -d_Big_K * fe_face_values[pressure_dg].gradient(j, q) * fe_face_values.normal_vector(q) * fe_face_values[pressure_dg].value(i, q) * fe_face_values.JxW(q);

                                    //These values are ZERO for EG
                                    cell_matrix(i, j) += -d_SForm * d_Big_K * fe_face_values[pressure_dg].gradient(i, q) * fe_face_values.normal_vector(q) * fe_face_values[pressure_dg].value(j, q) * fe_face_values.JxW(q);

                                    //DG - Weakly imposed BD
                                    //(1,1)
                                    // I-D-1
                                    cell_matrix(i, j) += penalty_term * d_Big_K * fe_face_values[pressure_dg].value(j, q) * fe_face_values[pressure_dg].value(i, q) * fe_face_values.JxW(q);

                                    // EG - Weakly Imposed BD.
                                    //(0,1)
                                    //I-E-1
                                    cell_matrix(i, j) += -d_Big_K * (phi_i_grads_p_face_cg[j]) * fe_face_values.normal_vector(q) * fe_face_values[pressure_dg].value(i, q) * fe_face_values.JxW(q);

                                    //I-E-2
                                    cell_matrix(i, j) += -d_SForm * d_Big_K * phi_i_grads_p_face_cg[i] * fe_face_values.normal_vector(q) * fe_face_values[pressure_dg].value(j, q) * fe_face_values.JxW(q);

                                    //I-E-4
                                    cell_matrix(i, j) += penalty_term * d_Big_K * fe_face_values[pressure_cg].value(j, q) * fe_face_values[pressure_dg].value(i, q) * fe_face_values.JxW(q);

                                    //I-E-5
                                    cell_matrix(i, j) += penalty_term * d_Big_K * fe_face_values[pressure_dg].value(j, q) * fe_face_values[pressure_cg].value(i, q) * fe_face_values.JxW(q);

                                } //int i.j

                                // I-E-2  &  I-E-3
                                cell_rhs(i) += -d_SForm * d_Big_K * pressure_dirichlet_bc // CG
                                               * fe_face_values.normal_vector(q) * fe_face_values[pressure_cg].gradient(i, q) * fe_face_values.JxW(q);

                                //(0,0) - EG Weakly imposed BD.
                                // I-E-4
                                cell_rhs(i) += penalty_term * d_Big_K * pressure_dirichlet_bc // - neighbor_pressure_boundary_term[q])
                                               * fe_face_values[pressure_dg].value(i, q) * fe_face_values.JxW(q);
                                cell_rhs(i) += penalty_term * d_Big_K * pressure_dirichlet_bc // - neighbor_pressure_boundary_term[q])
                                               * fe_face_values[pressure_cg].value(i, q) * fe_face_values.JxW(q);
                            }
                        } //for q
                    }
                } // cell - at_boundary
                // On the Interface
                else
                {

                    Assert(cell->neighbor(face_no).state() == IteratorState::valid, ExcInternalError());
                    const typename DoFHandler<dim>::cell_iterator neighbor = cell->neighbor(face_no);
                    // CASE A) HAS CHILDREN.
                    if (cell->face(face_no)->has_children())
                    {

                        // GET NEIGHBOR 2..
                        // NEIGHBOR FACE #NO
                        const unsigned int neighbor2 = cell->neighbor_face_no(face_no);

                        for (unsigned int subface_no = 0;
                             subface_no < cell->face(face_no)->number_of_children();
                             ++subface_no)
                        {
                            //Get Neighbor Child
                            typename DoFHandler<dim>::cell_iterator neighbor_child = cell->neighbor_child_on_subface(face_no, subface_no);

                            double h_e = cell->neighbor_child_on_subface(face_no, subface_no)->diameter();
                            Assert(!neighbor_child->has_children(), ExcInternalError());

                            fe_subface_values.reinit(cell, face_no, subface_no);
                            fe_face_values_neighbor.reinit(neighbor_child, neighbor2);

                            d_Big_K_neighbor = 1.;

                            //FOR EACH SUBSURFACE
                            cell_matrix_face = 0;
                            // IMPORTANT.
                            neighbor_child->get_dof_indices(local_dof_indices_face);

                            for (unsigned int q = 0; q < fe_subface_values.n_quadrature_points; ++q)
                            {
                                //for (unsigned int q=0; q<n_face_q_points; ++q){
                                // xueying : assign perm values
                                d_Big_K = permeability.value(fe_subface_values.quadrature_point(q), 0) / mu_f;
                                d_Big_K_neighbor = permeability.value(fe_face_values_neighbor.quadrature_point(q), 0) / mu_f;
                                if (test_case == TestCase::heterogeneous)
                                {
                                    d_Big_K = perm_function.value(fe_subface_values.quadrature_point(q) - eps * fe_subface_values.normal_vector(q), 0) / mu_f;
                                    d_Big_K_neighbor = perm_function.value(fe_face_values_neighbor.quadrature_point(q) - eps * fe_face_values_neighbor.normal_vector(q), 0) / mu_f;
                                }
                                // ADDED harmonic averaging of perm in the penalty term
                                double K_e = 2.0 * d_Big_K * d_Big_K_neighbor / (d_Big_K + d_Big_K_neighbor);
                                double beta_e = d_Big_K_neighbor / (d_Big_K + d_Big_K_neighbor);
                                if (beta_e != 0.5)
                                {
                                    gamma_penal = 1.0;
                                }
                                penalty_term = gamma_penal / h_e;
                                for (unsigned int k = 0; k < dofs_per_cell; ++k)
                                {
                                    phi_i_p_face_neighbor[k] = fe_face_values_neighbor[pressure_dg].value(k, q);
                                    phi_i_grads_p_face_neighbor[k] = fe_face_values_neighbor[pressure_dg].gradient(k, q);
                                    phi_i_grads_p_face_neighbor_cg[k] = fe_face_values_neighbor[pressure_cg].gradient(k, q);
                                    //phi_i_p_face[k]       = fe_face_values[pressure_dg].value (k, q);
                                    phi_i_grads_p_face[k] = fe_subface_values[pressure_dg].gradient(k, q);

                                    //phi_i_p_face_cg[k]       = fe_face_values[pressure_cg].value (k, q);
                                    phi_i_grads_p_face_cg[k] = fe_subface_values[pressure_cg].gradient(k, q);
                                }

                                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                                    {

                                        cell_matrix(i, j) += -(1 - beta_e) * d_Big_K * phi_i_grads_p_face_cg[j] * fe_subface_values.normal_vector(q) * fe_subface_values[pressure_dg].value(i, q) * fe_subface_values.JxW(q);

                                        cell_matrix(i, j) += -d_SForm * (1 - beta_e) * d_Big_K * phi_i_grads_p_face_cg[i] * fe_subface_values.normal_vector(q) * fe_subface_values[pressure_dg].value(j, q) * fe_subface_values.JxW(q);

                                        cell_matrix_face(i, j) += -beta_e * d_Big_K_neighbor * phi_i_grads_p_face_neighbor_cg[j] * fe_subface_values.normal_vector(q) * fe_subface_values[pressure_dg].value(i, q) * fe_subface_values.JxW(q);

                                        cell_matrix_face(i, j) += d_SForm * (1 - beta_e) * d_Big_K * phi_i_grads_p_face_cg[i] * fe_subface_values.normal_vector(q) * fe_face_values_neighbor[pressure_dg].value(j, q) * fe_subface_values.JxW(q);

                                        cell_matrix(i, j) += -(1 - beta_e) * d_Big_K * fe_subface_values[pressure_dg].gradient(j, q) * fe_subface_values.normal_vector(q) * fe_subface_values[pressure_dg].value(i, q) * fe_subface_values.JxW(q);

                                        cell_matrix(i, j) += -d_SForm * (1 - beta_e) * d_Big_K * fe_subface_values[pressure_dg].gradient(i, q) * fe_subface_values.normal_vector(q) * fe_subface_values[pressure_dg].value(j, q) * fe_subface_values.JxW(q);

                                        cell_matrix(i, j) += penalty_term * K_e * fe_subface_values[pressure_dg].value(j, q) * fe_subface_values[pressure_dg].value(i, q) * fe_subface_values.JxW(q);

                                        cell_matrix_face(i, j) += -penalty_term * K_e * fe_face_values_neighbor[pressure_dg].value(j, q) * fe_subface_values[pressure_dg].value(i, q) * fe_subface_values.JxW(q);

                                        cell_matrix_face(i, j) += -beta_e * d_Big_K_neighbor * phi_i_grads_p_face_neighbor[j] * fe_subface_values.normal_vector(q) * fe_subface_values[pressure_dg].value(i, q) * fe_subface_values.JxW(q);

                                        cell_matrix_face(i, j) += d_SForm * (1 - beta_e) * d_Big_K * fe_subface_values[pressure_dg].gradient(i, q) * fe_subface_values.normal_vector(q) * fe_face_values_neighbor[pressure_dg].value(j, q) * fe_subface_values.JxW(q);

                                    } //int J

                            } //For Q loop
                            constraints_pressure.distribute_local_to_global(cell_matrix_face,
                                                                            //local_dof_indices_face,
                                                                            local_dof_indices, //column
                                                                            local_dof_indices_face,
                                                                            system_matrix_pressure);

                        } //subface
                    }     // has children.
                    // CASE B) neighbor is not coarser
                    else if (!cell->neighbor_is_coarser(face_no))
                    {
                        const unsigned int neighbor_face = cell->neighbor_of_neighbor(face_no);

                        double h_e = cell->face(face_no)->diameter();
                        //DEBUG - necessary ?
                        fe_face_values.reinit(cell, face_no);
                        fe_face_values_neighbor.reinit(neighbor, neighbor_face);

                        //FOR EACH SUBSURFACE
                        cell_matrix_face = 0;
                        neighbor->get_dof_indices(local_dof_indices_face);

                        d_Big_K_neighbor = 1.;

                        for (unsigned int q = 0; q < n_face_q_points; ++q)
                        {
                            // xueying : assign perm values
                            d_Big_K = permeability.value(fe_face_values.quadrature_point(q), 0) / mu_f;
                            d_Big_K_neighbor = permeability.value(fe_face_values_neighbor.quadrature_point(q), 0) / mu_f;
                            if (test_case == TestCase::heterogeneous)
                            {
                                d_Big_K = perm_function.value(fe_face_values.quadrature_point(q) - eps * fe_face_values.normal_vector(q), 0) / mu_f;
                                d_Big_K_neighbor = perm_function.value(fe_face_values_neighbor.quadrature_point(q) - eps * fe_face_values_neighbor.normal_vector(q), 0) / mu_f;
                            }
                            double K_e = 2.0 * d_Big_K * d_Big_K_neighbor / (d_Big_K + d_Big_K_neighbor);
                            double beta_e = d_Big_K_neighbor / (d_Big_K + d_Big_K_neighbor);
                            if (beta_e != 0.5)
                            {
                                // cout << "d_Big_K = " << d_Big_K << endl;
                                // cout << "d_Big_K_neighbor = " << d_Big_K_neighbor << endl;
                                // cout << "beta_e = " << beta_e << endl;
                                gamma_penal = 1.0;
                            }
                            penalty_term = gamma_penal / h_e;
                            for (unsigned int k = 0; k < dofs_per_cell; ++k)

                            {
                                phi_i_p_face_neighbor[k] = fe_face_values_neighbor[pressure_dg].value(k, q);
                                phi_i_grads_p_face_neighbor[k] = fe_face_values_neighbor[pressure_dg].gradient(k, q);
                                phi_i_grads_p_face_neighbor_cg[k] = fe_face_values_neighbor[pressure_cg].gradient(k, q);
                                //phi_i_p_face[k]       = fe_face_values[pressure_dg].value (k, q);
                                phi_i_grads_p_face[k] = fe_face_values[pressure_dg].gradient(k, q);

                                //phi_i_p_face_cg[k]       = fe_face_values[pressure_cg].value (k, q);
                                phi_i_grads_p_face_cg[k] = fe_face_values[pressure_cg].gradient(k, q);
                            }

                            //DG DEBUG
                            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                                {

                                    cell_matrix(i, j) += -(1 - beta_e) * d_Big_K * phi_i_grads_p_face_cg[j] // + phi_i_grads_p_face_neighbor_cg[j])
                                                         * fe_face_values.normal_vector(q) * fe_face_values[pressure_dg].value(i, q) * fe_face_values.JxW(q);

                                    cell_matrix(i, j) += -d_SForm * (1 - beta_e) * d_Big_K * phi_i_grads_p_face_cg[i] // + phi_i_grads_p_face_neighbor_cg[j])
                                                         * fe_face_values.normal_vector(q) * fe_face_values[pressure_dg].value(j, q) * fe_face_values.JxW(q);

                                    cell_matrix_face(i, j) += -beta_e * d_Big_K_neighbor * phi_i_grads_p_face_neighbor_cg[j] * fe_face_values.normal_vector(q) * fe_face_values[pressure_dg].value(i, q) * fe_face_values.JxW(q);

                                    cell_matrix_face(i, j) += d_SForm * d_Big_K_neighbor * phi_i_grads_p_face_neighbor_cg[i] * fe_face_values.normal_vector(q) * fe_face_values[pressure_dg].value(j, q) * fe_face_values.JxW(q);

                                    cell_matrix(i, j) += -(1 - beta_e) * d_Big_K * fe_face_values[pressure_dg].gradient(j, q) * fe_face_values.normal_vector(q) * fe_face_values[pressure_dg].value(i, q) * fe_face_values.JxW(q);

                                    cell_matrix(i, j) += -d_SForm * (1 - beta_e) * d_Big_K * fe_face_values[pressure_dg].gradient(i, q) * fe_face_values.normal_vector(q) * fe_face_values[pressure_dg].value(j, q) * fe_face_values.JxW(q);

                                    cell_matrix(i, j) += penalty_term * K_e * fe_face_values[pressure_dg].value(j, q) * fe_face_values[pressure_dg].value(i, q) * fe_face_values.JxW(q);

                                    cell_matrix_face(i, j) += -penalty_term * K_e * fe_face_values_neighbor[pressure_dg].value(j, q) * fe_face_values[pressure_dg].value(i, q) * fe_face_values.JxW(q);

                                    cell_matrix_face(i, j) += -beta_e * d_Big_K_neighbor * phi_i_grads_p_face_neighbor[j] * fe_face_values.normal_vector(q) * fe_face_values[pressure_dg].value(i, q) * fe_face_values.JxW(q);

                                    cell_matrix_face(i, j) += d_SForm * (1 - beta_e) * d_Big_K * fe_face_values[pressure_dg].gradient(i, q) * fe_face_values.normal_vector(q) * fe_face_values_neighbor[pressure_dg].value(j, q) * fe_face_values.JxW(q);

                                } // For J loop

                        } //For Q loop
                        constraints_pressure.distribute_local_to_global(cell_matrix_face,
                                                                        //local_dof_indices_face,
                                                                        local_dof_indices, //column
                                                                        local_dof_indices_face,
                                                                        system_matrix_pressure);

                    } // else - no children no coarser

                    // CASE C) neighbor is coarser than cell
                    else
                    {

                        //i.e. neighbor is coarser than cell
                        std::pair<unsigned int, unsigned int> neighbor_face_subface = cell->neighbor_of_coarser_neighbor(face_no);

                        double h_e = cell->face(face_no)->diameter();
                        penalty_term = gamma_penal / h_e;
                        Assert(neighbor_face_subface.first < GeometryInfo<dim>::faces_per_cell, ExcInternalError());
                        Assert(neighbor_face_subface.second < neighbor->face(neighbor_face_subface.first)->number_of_children(), ExcInternalError());
                        Assert(neighbor->neighbor_child_on_subface(neighbor_face_subface.first, neighbor_face_subface.second) == cell, ExcInternalError());

                        //DEBUG - necessary ?
                        fe_face_values.reinit(cell, face_no);

                        fe_subface_values.reinit(neighbor, neighbor_face_subface.first,
                                                 neighbor_face_subface.second);

                        d_Big_K_neighbor = 1.;

                        //FOR EACH SUBSURFACE
                        cell_matrix_face = 0;
                        // IS this Correct ??
                        neighbor->get_dof_indices(local_dof_indices_face);

                        for (unsigned int q = 0; q < n_face_q_points; ++q)
                        {
                            //xueying : assign perm values
                            d_Big_K = permeability.value(fe_face_values.quadrature_point(q), 0) / mu_f;
                            d_Big_K_neighbor = permeability.value(fe_subface_values.quadrature_point(q), 0) / mu_f;
                            if (test_case == TestCase::heterogeneous)
                            {
                                d_Big_K = perm_function.value(fe_face_values.quadrature_point(q) - eps * fe_face_values.normal_vector(q), 0) / mu_f;
                                d_Big_K_neighbor = perm_function.value(fe_subface_values.quadrature_point(q) - eps * fe_subface_values.normal_vector(q), 0) / mu_f;
                            }
                            double K_e = 2.0 * d_Big_K * d_Big_K_neighbor / (d_Big_K + d_Big_K_neighbor);
                            double beta_e = d_Big_K_neighbor / (d_Big_K + d_Big_K_neighbor);
                            if (beta_e != 0.5)
                            {
                                gamma_penal = 1.0;
                            }
                            penalty_term = gamma_penal / h_e;
                            for (unsigned int k = 0; k < dofs_per_cell; ++k)
                            {
                                phi_i_p_face_neighbor[k] = fe_subface_values[pressure_dg].value(k, q);
                                phi_i_grads_p_face_neighbor[k] = fe_subface_values[pressure_dg].gradient(k, q);

                                phi_i_grads_p_face_neighbor_cg[k] = fe_subface_values[pressure_cg].gradient(k, q);

                                phi_i_grads_p_face[k] = fe_face_values[pressure_dg].gradient(k, q);
                                phi_i_grads_p_face_cg[k] = fe_face_values[pressure_cg].gradient(k, q);
                            }

                            //DG DEBUG
                            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                                for (unsigned int j = 0; j < dofs_per_cell; ++j)
                                {

                                    //(1,0)
                                    //H-1
                                    cell_matrix(i, j) += -(1 - beta_e) * d_Big_K * (phi_i_grads_p_face_cg[j]) // + phi_i_grads_p_face_neighbor_cg[j])
                                                         * fe_face_values.normal_vector(q) * fe_face_values[pressure_dg].value(i, q) * fe_face_values.JxW(q);

                                    //(0,1) - for SIPG
                                    //H-2
                                    //DEBUG May1. -- need to fix
                                    //cell_matrix(i,j) += - d_SForm
                                    //* 0.5 * d_Big_K
                                    //* fe_face_values[pressure_dg].gradient(i,q)
                                    //* fe_face_values.normal_vector(q)
                                    //* phi_i_p_face_cg[j]
                                    //* fe_face_values.JxW(q);
                                    cell_matrix(i, j) += -d_SForm * (1 - beta_e) * d_Big_K * fe_face_values[pressure_cg].gradient(i, q) * fe_face_values.normal_vector(q) * fe_face_values[pressure_dg].value(j, q) * fe_face_values.JxW(q);

                                    //Newly Debug
                                    //Different Matrix, Different Global Case
                                    // I = K, J=L

                                    //H-N-1
                                    cell_matrix_face(i, j) += -beta_e * d_Big_K_neighbor * phi_i_grads_p_face_neighbor_cg[j] * fe_face_values.normal_vector(q) * fe_face_values[pressure_dg].value(i, q) * fe_face_values.JxW(q);

                                    //H-N-2
                                    //DEBUG May1. -- need to fix
                                    //cell_matrix_face(i,j) +=  0.5 * d_SForm
                                    //*  d_Big_K
                                    //*  fe_face_values[pressure_dg].gradient(i,q)
                                    //*  fe_face_values.normal_vector(q)
                                    //*  fe_face_values_neighbor[pressure_cg].value(j,q) * fe_face_values.JxW(q);
                                    cell_matrix_face(i, j) += d_SForm * (1 - beta_e) * d_Big_K * fe_face_values[pressure_cg].gradient(i, q) * fe_face_values.normal_vector(q) * fe_subface_values[pressure_dg].value(j, q) * fe_face_values.JxW(q);

                                    //(1,1) - DG
                                    //These terms are ZERO in EG
                                    cell_matrix(i, j) += -(1 - beta_e) * d_Big_K * fe_face_values[pressure_dg].gradient(j, q) * fe_face_values.normal_vector(q) * fe_face_values[pressure_dg].value(i, q) * fe_face_values.JxW(q);

                                    //These terms are ZERO in EG
                                    cell_matrix(i, j) += -d_SForm * (1 - beta_e) * d_Big_K * fe_face_values[pressure_dg].gradient(i, q) * fe_face_values.normal_vector(q) * fe_face_values[pressure_dg].value(j, q) * fe_face_values.JxW(q);

                                    //(1,1)
                                    //H-3
                                    cell_matrix(i, j) += penalty_term * K_e * fe_face_values[pressure_dg].value(j, q) * fe_face_values[pressure_dg].value(i, q) * fe_face_values.JxW(q);

                                    // Newly DEBUG
                                    // Different Matrix, Diffeernt Global Case,
                                    // I = K, J = L
                                    //H-N-3
                                    cell_matrix_face(i, j) += -penalty_term * K_e * fe_subface_values[pressure_dg].value(j, q) * fe_face_values[pressure_dg].value(i, q) * fe_face_values.JxW(q);

                                    //These terms are ZERO in EG
                                    cell_matrix_face(i, j) += -beta_e * d_Big_K_neighbor * phi_i_grads_p_face_neighbor[j] * fe_face_values.normal_vector(q) * fe_face_values[pressure_dg].value(i, q) * fe_face_values.JxW(q);

                                    //These terms are ZERO in EG
                                    cell_matrix_face(i, j) += d_SForm * (1 - beta_e) * d_Big_K * fe_face_values[pressure_dg].gradient(i, q) * fe_face_values.normal_vector(q) * fe_subface_values[pressure_dg].value(j, q) * fe_face_values.JxW(q);

                                } // For J loop

                        } //For Q loop
                        constraints_pressure.distribute_local_to_global(cell_matrix_face,
                                                                        //local_dof_indices_face,
                                                                        local_dof_indices, //column
                                                                        local_dof_indices_face,
                                                                        system_matrix_pressure);

                    } //coarser
                }     // else- no boundary

            } // for face loop
            constraints_pressure.distribute_local_to_global(cell_matrix, local_dof_indices,
                                                            system_matrix_pressure);
            constraints_pressure.distribute_local_to_global(cell_rhs, local_dof_indices,
                                                            system_rhs_pressure);
        }
    }
    /*
    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler_pressure,
                                             0,
                                             ZeroFunction<2>(),
                                             boundary_values);
    MatrixTools::apply_boundary_values(boundary_values, system_matrix_pressure, solution_pressure, system_rhs_pressure);
    */
    system_matrix_pressure.compress(VectorOperation::add);
    system_rhs_pressure.compress(VectorOperation::add);
}
