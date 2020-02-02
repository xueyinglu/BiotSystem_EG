#include "BiotSystem.h"
using namespace std;
#include <iostream>
#include <fstream>
void BiotSystem::output_error() {
    convergence_table.set_precision("L2_p", 4);
    convergence_table.set_precision("L2_u", 4);
    convergence_table.set_precision("energy_u", 4);
    convergence_table.set_scientific("L2_p", true);
    convergence_table.set_scientific("L2_u", true);
    convergence_table.set_scientific("energy_u", true);
    convergence_table.set_tex_caption("L2_p", "L^2-error p");
    convergence_table.set_tex_caption("L2_u", "L^2-error $\\mathbf{u}$");

    ofstream error_table_file("error-"+to_string(num_global_refinement) + "-" +to_string(del_t) +".tex");
    convergence_table.write_tex(error_table_file);

    /*
    // p_indicators_table.set_precision("eta_fs", 6);
    p_indicators_table.set_precision("eta_alg", 6);
    p_indicators_table.set_precision("eta_time", 6);
    p_indicators_table.set_precision("eta_flow", 6);
    // p_indicators_table.set_precision("eta_p_residual", 6);
    //p_indicators_table.set_precision("eta_flux_jump", 6);
    p_indicators_table.set_precision("sum", 6);
    p_indicators_table.set_precision("error", 6);
    p_indicators_table.set_precision("eff", 6);
    //p_indicators_table.set_scientific("eta_fs",true);
    p_indicators_table.set_scientific("eta_alg",true);
    p_indicators_table.set_scientific("eta_time", true);
    p_indicators_table.set_scientific("eta_flow", true);
    // p_indicators_table.set_scientific("eta_p_residual", true);
    // p_indicators_table.set_scientific("eta_flux_jump", true);
    p_indicators_table.set_scientific("sum", true);
    p_indicators_table.set_scientific("error", true);
    p_indicators_table.set_scientific("eff", true);
    //p_indicators_table.set_tex_caption("eta_fs", "$\\eta_{\\text{fs}}$");
    p_indicators_table.set_tex_caption("eta_alg", "$\\eta_{\\text{alg}}$");
    p_indicators_table.set_tex_caption("eta_time", "$\\eta_{\\text{time}}$");
    p_indicators_table.set_tex_caption("eta_flow", "$\\eta_{\\text{flow}}$");
    //p_indicators_table.set_tex_caption("eta_p_residual", "$\\eta_{\\text{p\\_residual}}$");
    //p_indicators_table.set_tex_caption("eta_flux_jump", "$\\eta_{\\text{flux\\_jump}}$");
    p_indicators_table.set_tex_table_caption("Flow indicators: $h = 1/(2^" + to_string(num_global_refinement)+ "), \\Delta t = " + to_string(del_t) +"$"); 
    
    u_indicators_table.set_precision("eta_face_partial_sigma", 6);
    u_indicators_table.set_precision("eta_face_sigma", 6);
    u_indicators_table.set_precision("eta_partial_u", 6);
    u_indicators_table.set_precision("eta_u", 6);
    u_indicators_table.set_precision("sum", 6);
    u_indicators_table.set_precision("error", 6);
    u_indicators_table.set_precision("eff", 6);
    u_indicators_table.set_scientific("eta_face_partial_sigma", true);
    u_indicators_table.set_scientific("eta_face_sigma", true);
    u_indicators_table.set_scientific("eta_partial_u", true);
    u_indicators_table.set_scientific("eta_u", true);
    u_indicators_table.set_scientific("sum", true);
    u_indicators_table.set_scientific("error", true);
    u_indicators_table.set_scientific("eff", true);
    u_indicators_table.set_tex_caption("eta_face_partial_sigma", "$\\eta_{\\mathcal{E}_{\\partial \\sigma}}$");
    u_indicators_table.set_tex_caption("eta_face_sigma", "$\\eta_{\\mathcal{E}_{\\sigma}}$");
    u_indicators_table.set_tex_caption("eta_partial_u", "$\\eta_{\\mathcal{T}_{\\partial u}}$");
    u_indicators_table.set_tex_caption("eta_u", "$\\eta_{\\mathcal{T}_{u}}$");
    u_indicators_table.set_tex_table_caption("Mechanics indicators: $h = 1/(2^" + to_string(num_global_refinement)+ "), \\Delta t = " + to_string(del_t) +"$"); 
    ofstream aposterior_table_file("aposteriori-"+to_string(num_global_refinement) + "-" +to_string(del_t)+".tex");
    p_indicators_table.write_tex(aposterior_table_file, false);
    u_indicators_table.write_tex(aposterior_table_file, false);

    efficiency_table.set_precision("error", 6);
    efficiency_table.set_precision("index", 6);
    efficiency_table.set_precision("eta_sum", 6);
    efficiency_table.set_scientific("error",true);
    efficiency_table.set_scientific("index",true);
    efficiency_table.set_scientific("eta_sum",true);
    efficiency_table.set_tex_table_caption("Efficiency Index: $h = 1/(2^" + to_string(num_global_refinement)+ "), \\Delta t = " + to_string(del_t) +"$");
    ofstream efficiency_table_file("efficiency-"+to_string(num_global_refinement) + "-" +to_string(del_t) +".tex");
    efficiency_table.write_tex(efficiency_table_file, false);
    */
}