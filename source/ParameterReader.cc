#include "ParameterReader.h"
ParameterReader::ParameterReader(
    ParameterHandler &paramhandler)
    : prm(paramhandler)
{
}

void ParameterReader::read_parameters(
    const std::string parameter_file)
{
    declare_parameters();

    if (!prm.read_input(parameter_file, true))
        AssertThrow(false, ExcMessage("could not read .prm!"));
}

void ParameterReader::declare_parameters()
{
    prm.enter_subsection("Global parameters");
    {
        prm.declare_entry("Global refinement numbers", "3",
                          Patterns::Integer(0));
        
        prm.declare_entry("Final time", "1.0", Patterns::Double(0));

        prm.declare_entry("Timestep size", "1.0", Patterns::Double(0));

        prm.declare_entry("Fixed-stress criteria", "1",
                          Patterns::Integer(0));

        prm.declare_entry("Fixed-stress tolerance", "1e-7", Patterns::Double(0));
        
        prm.declare_entry("Output filename", "solution_", Patterns::Anything());
        
        prm.declare_entry("Test case", "benchmark", Patterns::Anything());
        
        prm.declare_entry("Adaptivity flag", "false", Patterns::Bool());
        
        prm.declare_entry("Permeability multiplier flag", "false", Patterns::Bool());
        
        prm.declare_entry("Penalization parameter", "1000", Patterns::Double(0));
    }
    prm.leave_subsection();

}