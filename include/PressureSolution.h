#ifndef PRESSURE_SOLUTION_H_
#define PRESSURE_SOLUTION_H_
#include "DealiiHeader.h"
using namespace std;
using namespace dealii;
class PressureSolution: public Function<dim>{
    public:
    PressureSolution( double _t):
    Function<dim> (1) {t = _t;};

    virtual double value(const Point<dim> &p, const unsigned int component) const;

    private:
    double t;
};
#endif