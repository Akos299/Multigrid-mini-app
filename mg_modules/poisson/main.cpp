
#include "poisson.hpp"
#include <iostream>
#include <vector>
#include <memory>

int main(int argc, char** argv)
{
    std::cout << " poisson solver. " << "\n\n";
    return 0;
    try{
        // Declare some option variables
    
        
        // Set up command-line options
       
        
        // Run through given aspect ratios
        std::cout << "Running..." << std::endl;

        // Solve problems   
 
        multigrid::Settings<double> settings;
        std::array<double, 3> starts = {-0.5, -0.5, -0.5} ;
        std::array<double, 3> ends  = {0.5, 0.5, 0.5};
        multigrid::TransfertOperator rest_Op = multigrid::PointAverage;
        multigrid::TransfertOperator int_Op = multigrid::Linear;
        std::unique_ptr<Poisson<double, 3>> problem(new Poisson<double,3>(settings,starts,ends,rest_Op,int_Op ));
        problem->solve();
        std::cout << "Finished!" << std::endl;
        return 0;
        
    } catch (std::exception& e) {
		std::cout << e.what() << std::endl;
        return 1;
    }    
}