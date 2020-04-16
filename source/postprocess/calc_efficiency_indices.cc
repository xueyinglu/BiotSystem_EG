#include "BiotSystem.h"
using namespace std;
void BiotSystem::calc_efficiency_indices()
{
    eta_sum.push_back(eta_alg.back() + eta_time.back() + eta_flow.back() + eta_jump.back() +
                      +eta_face_partial_sigma.back() + eta_face_sigma.back() + eta_partial_u.back() + eta_u.back());

    double error_sq = l2_error_p.back() * l2_error_p.back() * biot_inv_M / 4 + energy_error_u.back() * energy_error_u.back()/ 4;

    calc_p_h_norm();
    for (auto &n : h_error_p_sq)
    {
        error_sq += 1.0 / 2 / mu_f * del_t * n;
    }

    double eta_f = sqrt(eta_alg.back() + eta_time.back() + eta_flow.back() + eta_jump.back() );

    double I_eff1 = eta_f/ sqrt(error_sq);

    double eta_mechanics = sqrt(eta_face_partial_sigma.back() + eta_face_sigma.back() + eta_partial_u.back() + eta_u.back());

    double I_eff2 = eta_mechanics/(energy_error_u_1.back() + biot_alpha * l2_error_p.back());
    double lower_bound = eta_u_equilibrium.back() / (energy_error_u_1.back() + biot_alpha * l2_error_p.back());

    efficiency_table.add_value("time", t);
    efficiency_table.add_value("num_fs", num_fs.back());
    efficiency_table.add_value("eta_flow", eta_f );
    efficiency_table.add_value("error1", sqrt(error_sq)) ;
    efficiency_table.add_value("index1", I_eff1);
    efficiency_table.add_value("eta_mechanics", eta_mechanics );
    efficiency_table.add_value("error2", energy_error_u_1.back() + biot_alpha * l2_error_p.back());
    efficiency_table.add_value("index2", I_eff2);
    efficiency_table.add_value("eta_u_eq", eta_u_equilibrium.back());
    efficiency_table.add_value("lower_bound", lower_bound);

    efficiency_table.set_precision("eta_flow", 4);
    efficiency_table.set_precision("error1", 4);
    efficiency_table.set_precision("index1", 4);
    efficiency_table.set_precision("eta_mechanics", 4);
    efficiency_table.set_precision("error2", 4);
    efficiency_table.set_precision("index2", 4);
    efficiency_table.set_precision("eta_u_eq", 4);
    efficiency_table.set_precision("lowe_bound", 4);

    efficiency_table.set_scientific("time", true);
    efficiency_table.set_scientific("eta_flow", true);
    efficiency_table.set_scientific("error1", true);
    efficiency_table.set_scientific("index1", true);
    efficiency_table.set_scientific("eta_mechanics", true);
    efficiency_table.set_scientific("error2", true);
    efficiency_table.set_scientific("index2", true);
    efficiency_table.set_scientific("eta_u_eq", true);
    efficiency_table.set_scientific("lower_bound", true);

    efficiency_table.set_tex_table_caption("Efficiency Index: $h = 1/(2^" + to_string(num_global_refinement) + "), \\Delta t = " + to_string(del_t) + "$");
    ofstream efficiency_table_file(filename_base + "-efficiency-" + to_string(num_global_refinement) + "-" + to_string(del_t) + ".tex");
    efficiency_table.write_tex(efficiency_table_file, false);
}
