#ifndef PRESSURE_SOLUTIONEG_H_
#define PRESSURE_SOLUTIONEG_H_
#include "DealiiHeader.h"
using namespace std;
using namespace dealii;
class PressureSolutionEG: public Function<dim>{
    public:
    PressureSolutionEG( double _t):
    Function<dim> (2) {t = _t;};

    virtual void vector_value(const Point<dim> &p, Vector<double> &values) const;

    private:
    double t;
};
#endif