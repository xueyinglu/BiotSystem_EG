#ifndef PARAMETERREADER_H_
#define PARAMETERREADER_H_
#include "DealiiHeader.h"

using namespace std;
using namespace dealii;
class ParameterReader : public Subscriptor
{
public:
    ParameterReader(
        ParameterHandler &);
    void
    read_parameters(
        const std::string);

private:
    void
    declare_parameters();
    ParameterHandler &prm;
};

#endif