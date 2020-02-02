#include "BiotSystem.h"
using namespace std;
void BiotSystem::calc_efficiency_indices(){
    eta_sum.push_back(eta_alg.back()+eta_time.back() +eta_flow.back()
    +eta_face_partial_sigma.back() + eta_face_sigma.back() + eta_partial_u.back()+eta_u.back());

    double error = l2_error_p.back()*l2_error_p.back()*biot_inv_M/4 + energy_error_u.back()*energy_error_u.back();
    double I_eff = sqrt(eta_sum.back()/error);

    efficiency_table.add_value("time", t);
    efficiency_table.add_value("error", error);
    efficiency_table.add_value("eta_sum", eta_sum.back());
    efficiency_table.add_value("index", I_eff);
}