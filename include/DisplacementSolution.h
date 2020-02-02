#ifndef DISPLACEMENT_SOLUTION_H_
#define DISPLACEMENT_SOLUTION_H_
#include "DealiiHeader.h"
using namespace std;
using namespace dealii;

class DisplacementSolution: public Function<dim> {
    public :
    DisplacementSolution( double _t):
    Function<dim> (dim) {t = _t;};
    
    virtual void vector_value(const Point<dim> & point, Vector<double> & values) const;
    virtual void gradient_value(const Point<dim> & point, Tensor<2,dim> & gradient) const;
    virtual void vector_value_list(const vector<Point<dim> > & points, vector<Vector<double> > & value_list) const ;
    private :
    double t;
};

#endif