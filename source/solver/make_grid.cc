#include "BiotSystem.h"

void BiotSystem::make_grid()
{
    if (test_case == TestCase::benchmark)
    {
        GridGenerator::hyper_cube(triangulation, 0, 1);
        // triangulation.begin_active()->face(0)->set_boundary_id(1);
        // triangulation.begin_active()->face(1)->set_boundary_id(1);
        triangulation.refine_global(num_global_refinement);
        std::cout << "Number of active cells:" << triangulation.n_active_cells() << std::endl;
    }
    else if (test_case == TestCase::terzaghi || test_case == TestCase::heterogeneous)
    {
        GridGenerator::hyper_cube(triangulation, 0, 1);
        // top: id = 0;
        // left and right : id =1;
        // bottom: id = 2;
        triangulation.begin_active()->face(0)->set_boundary_id(1);
        triangulation.begin_active()->face(1)->set_boundary_id(1);
        triangulation.begin_active()->face(2)->set_boundary_id(2);
        triangulation.refine_global(num_global_refinement);
        std::cout << "Number of active cells:" << triangulation.n_active_cells() << std::endl;
    }
    else if (test_case = TestCase::mandel)
    {   
        cout << "---------------Mandel problem: Making grid----------------"<<endl;
        // Bin Wang Data
        // GridGenerator::hyper_rectangle(triangulation, Point<dim>(0, 0), Point<dim>(100, 10));
        // Phillips Phillips Data 
        GridGenerator::hyper_cube(triangulation, 0, 1);
        triangulation.begin_active()->face(0)->set_boundary_id(0);
        triangulation.begin_active()->face(1)->set_boundary_id(1);
        triangulation.begin_active()->face(2)->set_boundary_id(2);
        triangulation.begin_active()->face(3)->set_boundary_id(3);
        triangulation.refine_global(num_global_refinement);
        std::cout << "Number of active cells:" << triangulation.n_active_cells() << std::endl;
    }

    else if (test_case == TestCase::L_shape)
    {
        cout << "---------------L-shape: Making grid----------------"<<endl;
        GridGenerator::hyper_L(triangulation);
        triangulation.refine_global(num_global_refinement);
        std::cout << "Number of active cells:" << triangulation.n_active_cells() << std::endl;
    }
}
