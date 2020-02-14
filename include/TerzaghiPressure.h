#ifndef TERZAGHI_PRESSURE_H_
#define TERZAGHI_PRESSURE_H_
#include "DealiiHeader.h"
using namespace std;
using namespace dealii;
class TerzaghiPressure: public Function<dim>{
    public:
    TerzaghiPressure( double _t):
    Function<dim> (1) {t = _t;};

    virtual double value(const Point<dim> &p, const unsigned int component) const;

    private:
    double t;
};
#endif