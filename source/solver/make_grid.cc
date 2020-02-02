#include "BiotSystem.h"

void BiotSystem::make_grid()
{
    switch (test_case)
    {
    case simple:
        GridGenerator::hyper_cube(triangulation, 0, 1);
        //itriangulation.begin_active()->face(0) -> set_boundary_id(1);
        //triangulation.begin_active()->face(1) -> set_boundary_id(1);
        triangulation.refine_global(num_global_refinement);
        std::cout << "Number of active cells:" << triangulation.n_active_cells() << std::endl;
    //case mandel:
    //    GridGenerator::hyper_rectangle(triangulation, Point<dim>(0, 0), Point<dim>(100, 10));
    //    triangulation.refine_global(num_global_refinement);
    //    std::cout << "Number of active cells:" << triangulation.n_active_cells() << std::endl;
    }
}