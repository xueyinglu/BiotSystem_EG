#ifndef MANDEL_PRESSURE_H_
#define MANDEL_PRESSURE_H_
#include "DealiiHeader.h"
using namespace std;
using namespace dealii;

class MandelPressure : public Function<dim>
{
public:
    MandelPressure(double _t) : Function<dim>(1) { t = _t; };
    virtual double value(const Point<dim> &p, const unsigned int component) const;

private:
    double t;
};

#endif