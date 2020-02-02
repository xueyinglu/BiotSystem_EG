#ifndef INITIALPRESSURE_H_
#define INITIALPRESSURE_H_
// class for initial pressure for parabolic problems
#include "DealiiHeader.h"
using namespace std;
using namespace dealii;
class InitialPressure: public Function<dim>{
    public:
        InitialPressure():
        Function<dim>(2)
        {}

        virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const;
        virtual void grad_list(const vector<Point<dim>> &points, vector<Tensor<1, dim>> & values) const;

};
#endif