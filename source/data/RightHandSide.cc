#include "../include/RightHandSide.h"
  RightHandSide::RightHandSide ()
    :
    Function<dim> (dim)
  {}


  inline
  void RightHandSide::vector_value (const Point<dim> &p,
                                         Vector<double>   &values) const
  {
    Assert (values.size() == dim,
            ExcDimensionMismatch (values.size(), dim));
    Assert (dim >= 2, ExcNotImplemented());

    Point<dim> point_1, point_2;
    point_1(0) = 0.5;
    point_2(0) = -0.5;

    if (((p-point_1).norm_square() < 0.2*0.2) ||
        ((p-point_2).norm_square() < 0.2*0.2))
      values(0) = 1;
    else
      values(0) = 0;

    if (p.norm_square() < 0.2*0.2)
      values(1) = 1;
    else
      values(1) = 0;
  }


  void RightHandSide::vector_value_list (const std::vector<Point<dim> > &points,
                                              std::vector<Vector<double> >   &value_list) const
  {
    Assert (value_list.size() == points.size(),
            ExcDimensionMismatch (value_list.size(), points.size()));

    const unsigned int n_points = points.size();

    for (unsigned int p=0; p<n_points; ++p)
      RightHandSide::vector_value (points[p],
                                        value_list[p]);
  }
