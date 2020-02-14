#ifndef PERM_FUNCTION_H
#define PERM_FUNCTION_H

#include "DealiiHeader.h"
using namespace std;
using namespace dealii;
class PermFunction: public Function<dim>{
    public:
        PermFunction():
        Function<dim>(1)
        {}

        virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const;
};
#endif