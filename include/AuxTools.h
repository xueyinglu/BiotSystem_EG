#ifndef AUX_TOOLS_H_
#define AUX_TOOLS_H_
#include <fstream>
#include <sstream>
using namespace std;
using namespace dealii;

template <class PreconditionerA, class PreconditionerC>
class BlockDiagonalPreconditioner
{
public:
    BlockDiagonalPreconditioner(const LA::MPI::BlockSparseMatrix &M,
                                const PreconditionerA &pre_A, const PreconditionerC &pre_C)
        : matrix(M),
          prec_A(pre_A),
          prec_C(pre_C)
    {
    }

    void vmult(LA::MPI::BlockVector &dst,
               const LA::MPI::BlockVector &src) const
    {
        prec_A.vmult(dst.block(0), src.block(0));
        prec_C.vmult(dst.block(1), src.block(1));
    }

    const LA::MPI::BlockSparseMatrix &matrix;
    const PreconditionerA &prec_A;
    const PreconditionerC &prec_C;
};

// Define some tensors for cleaner notation later.
namespace Tensors
{
template <int dim>
inline Tensor<1, dim>
get_grad_p(unsigned int q,
           std::vector<std::vector<Tensor<1, dim>>> old_solution_grads)
{
  Tensor<1, dim> grad_p;
  grad_p[0] = old_solution_grads[q][0][0];
  grad_p[1] = old_solution_grads[q][0][1];
  if (dim == 3)
    grad_p[2] = old_solution_grads[q][0][2];

  return grad_p;
}

template <int dim>
inline double
get_deviator_norm(const Tensor<2, dim> deviator)
{

  return std::sqrt(deviator[0][0] * deviator[0][0] + deviator[0][1] * deviator[0][1] + deviator[1][0] * deviator[1][0] + deviator[1][1] * deviator[1][1]);
  if (dim == 3)
  {
    cout << " Wrong get deviator norm: to be fixed" << endl;
  }
}

template <int dim>
inline Tensor<1, dim>
get_grad_pf(
    unsigned int q,
    const std::vector<std::vector<Tensor<1, dim>>> &old_solution_grads)
{
  Tensor<1, dim> grad_pf;
  grad_pf[0] = old_solution_grads[q][dim][0];
  grad_pf[1] = old_solution_grads[q][dim][1];
  if (dim == 3)
    grad_pf[2] = old_solution_grads[q][dim][2];

  return grad_pf;
}

template <int dim>
inline Tensor<2, dim>
get_grad_u(unsigned int q,
           std::vector<std::vector<Tensor<1, dim>>> old_solution_grads)
{
  Tensor<2, dim> grad_u;
  grad_u[0][0] = old_solution_grads[q][0][0];
  grad_u[0][1] = old_solution_grads[q][0][1];

  grad_u[1][0] = old_solution_grads[q][1][0];
  grad_u[1][1] = old_solution_grads[q][1][1];
  if (dim == 3)
  {
    grad_u[0][2] = old_solution_grads[q][0][2];

    grad_u[1][2] = old_solution_grads[q][1][2];

    grad_u[2][0] = old_solution_grads[q][2][0];
    grad_u[2][1] = old_solution_grads[q][2][1];
    grad_u[2][2] = old_solution_grads[q][2][2];
  }

  return grad_u;
}

template <int dim>
inline Tensor<2, dim>
get_Identity()
{
  Tensor<2, dim> identity;
  identity[0][0] = 1.0;
  identity[1][1] = 1.0;
  if (dim == 3)
    identity[2][2] = 1.0;

  return identity;
}

template <int dim>
inline Tensor<1, dim>
get_u(unsigned int q,
      std::vector<Vector<double>> old_solution_values)
{
  Tensor<1, dim> u;
  u[0] = old_solution_values[q](0);
  u[1] = old_solution_values[q](1);
  if (dim == 3)
    u[2] = old_solution_values[q](2);

  return u;
}

template <int dim>
inline Tensor<1, dim>
get_u_LinU(const Tensor<1, dim> phi_i_u)
{
  Tensor<1, dim> tmp;
  tmp[0] = phi_i_u[0];
  tmp[1] = phi_i_u[1];
  if (dim == 3)
    tmp[2] = phi_i_u[2];
  return tmp;
}

template <int dim>
inline double
get_divergence_u(const Tensor<2, dim> grad_u)
{
  double tmp;
  if (dim == 2)
  {
    tmp = grad_u[0][0] + grad_u[1][1];
  }
  else if (dim == 3)
  {
    tmp = grad_u[0][0] + grad_u[1][1] + grad_u[2][2];
  }

  return tmp;
}

} // namespace Tensors

#endif