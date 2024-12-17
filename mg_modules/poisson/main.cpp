
#include "poisson.hpp"
#include "../include/multigrid_base.hpp"
#include <iostream>
#include <vector>
#include <memory>
#include "../include/ndArray.h"

int main(int argc, char** argv)
{
        // Declare some option variables
    
        
        // Set up command-line options
       
        
        // Run through given aspect ratios
        std::cout << "Running..." << std::endl;

        // Solve problems   
 
    //    struct multigrid::Settings _set;

        std::array<double, 3> starts = {-0.5, -0.5, -0.5} ;
        std::array<double, 3> ends  = {0.5, 0.5, 0.5};
        multigrid::TransfertOperator rest_Op = multigrid::PointAverage;
        multigrid::TransfertOperator int_Op = multigrid::Linear;
           std::cout << "11111111111 Finished!" << std::endl;
        // std::unique_ptr<Poisson<double, 3>> problem(new Poisson<double,3>(settings,starts,ends,rest_Op,int_Op ));

        // multigrid::Poisson<double, 3> problem(settings, starts,ends, rest_Op, int_Op);
        multigrid::my_settings setss;
        // std::cout << "Finished!" << std::endl;
        multigrid::Poisson<double,3> problem(setss,starts,ends,rest_Op,int_Op);

        std::cout << setss.max_iter << "\n";
        // problem.solve();
        std::cout << "Finished!" << std::endl;
        return 0;

        // const nd::index_t size[1] = {5};
        // nd::ndArray<double, 1> A(new double[5], size, true);
        // A({1}) = 0;
        // A.info();
        // std::cout << A.size(0);

        // return 0;
    // std::cout << " poisson solver. " << "\n\n";
    // return 0;
    // try{
    //     // Declare some option variables
    
        
    //     // Set up command-line options
       
        
    //     // Run through given aspect ratios
    //     std::cout << "Running..." << std::endl;

    //     // Solve problems   
 
    // //    struct multigrid::Settings _set;

    //     std::array<double, 3> starts = {-0.5, -0.5, -0.5} ;
    //     std::array<double, 3> ends  = {0.5, 0.5, 0.5};
    //     multigrid::TransfertOperator rest_Op = multigrid::PointAverage;
    //     multigrid::TransfertOperator int_Op = multigrid::Linear;
    //     // std::unique_ptr<Poisson<double, 3>> problem(new Poisson<double,3>(settings,starts,ends,rest_Op,int_Op ));

    //     // multigrid::Poisson<double, 3> problem(settings, starts,ends, rest_Op, int_Op);
    //     multigrid::my_settings setss;
    //      multigrid::Poisson<double,3> problem(setss);

    //     std::cout << setss.max_iter << "\n";
    //     // problem.solve();
    //     std::cout << "Finished!" << std::endl;
    //     return 0;

    //     const nd::index_t size[1] = {5};
    //     nd::ndArray<double, 1> A(new double[5], size, true);
    //     A({1}) = 0;
    //     A.info();
    //     std::cout << A.size(0);
    // } catch (std::exception& e) {
	// 	std::cout << e.what() << std::endl;
    //     return 1;
    // }    
}