#include "DisplacementSolution.h"
using namespace std;
void DisplacementSolution::vector_value(const Point<dim> &p, Vector<double> & values) const
{

  double k = 0.05;
  double alpha = 0.75;
  double inv_M = 3. / 28;
  double PI = atan(1) * 4;
  double tx = PI *p(0);
  double ty = PI *p(1);
  double x = p(0);
  double y = p(1);
  double A = 2 * PI * PI * k / (alpha + inv_M);
  //double A = 0;
   values(0) = -exp(-A * t) / (2 * PI) * cos(tx) * sin(ty);
   values(1) = -exp(-A * t) / (2 * PI) * sin(tx) * cos(ty);

}

void DisplacementSolution::vector_value_list(const vector<Point<dim>> &points,
                                      vector<Vector<double>> &value_list) const
{
  Assert(value_list.size() == points.size(),
         ExcDimensionMismatch(value_list.size(), points.size()));

  const unsigned int n_points = points.size();

  for (unsigned int p = 0; p < n_points; ++p)
    vector_value(points[p], value_list[p]);
}

void DisplacementSolution:: gradient_value(const Point<dim> &p, Tensor<2,dim> & gradient) const {

  double k = 0.05;
  double alpha = 0.75;
  double inv_M = 3. / 28;
  double PI = atan(1) * 4;
  double tx = PI *p(0);
  double ty = PI *p(1);
  double x = p(0);
  double y = p(1);
  double A = 2 * PI * PI * k / (alpha + inv_M);
  gradient[0][0] = exp(-A * t) / 2  *sin(PI*x) *sin(PI*y);
  gradient[0][1] = - exp(-A * t) / 2  * cos(PI*x) * cos(PI*y);
  gradient[1][0] = - exp(-A * t) / 2 * cos(PI*x) * cos(PI*y);
  gradient[1][1] = exp(-A * t) / 2  *sin(PI*x) *sin(PI*y); 
}