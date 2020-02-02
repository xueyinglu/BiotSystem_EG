#ifndef RIGHTHANDSIDE_H_
#define RIGHTHANDSIDE_H_
#include "DealiiHeader.h"
using namespace dealii;
  class RightHandSide :  public Function<dim>
  {
  public:
    RightHandSide ();

    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &values) const;

    virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                    std::vector<Vector<double> >   &value_list) const;
  };

#endif