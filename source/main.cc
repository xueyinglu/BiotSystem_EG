#include "BiotSystem.h"
using namespace std;
int main(int argc, char* argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  try
  {
    int _numGlobalRefinement = atoi(argv[1]);
    double _delT = stod(argv[2]);
    double _T = stod(argv[3]);
    double _tol = stod(argv[4]);
    BiotSystem convergence_test(_numGlobalRefinement, _delT, _T, _tol);
    cout << "before run fixed stress" << endl;
    convergence_test.run_fixed_stress();
  //   BiotSystem u_convergence;
  //   u_convergence.check_disp_solver_convergence();
  }
  catch (std::exception &exc)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

    return 1;
  }
  catch (...)
  {
    std::cerr << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }

  return 0;
}