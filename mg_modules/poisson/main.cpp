
#include "poisson.hpp"
#include "multigrid_base.hpp"
#include <iostream>
#include <vector>
#include <memory>
#include "ndArray.h"

int main(int argc, char** argv)
{
        // Declare some option variables
    
        
        // Set up command-line options
       
        
        // Run through given aspect ratios
        std::cout << "Running..." << std::endl;

//================== @@@@@@@ Solve problems  @@@@@@@ =================== 

        multigrid::my_settings init_setting;
        /* set the level numbers */
        init_setting.nb_levels = 3; // set the level numbers
        /* set finest grid geometry*/
        init_setting.nx_fine   = 16; init_setting.ny_fine   = 16; init_setting.nz_fine   = 16;
        /* set the maximum number of iteration*/
        init_setting.max_iter =  5;
        /* set the mgi cycle type */
        init_setting.cycle_type = multigrid::vCycle;
        /* set the number of pre and post smoothings*/
        init_setting.npre = 2; init_setting.npost = 2;
        /* set the coarser levels boundary conditions
        * for the mgi we used zero-fixed boundary condition for all coaser levels.
        * But for FAS or non-linear problems, we should used the physical boundary
        * conditions.
        */
        init_setting.bc_condition_type=multigrid::ZeroFixed;

        /* set the initial grid sizes */
        std::array<double, 3> starts = {-0.5, -0.5, -0.5} ;
        std::array<double, 3> ends  = {0.5, 0.5, 0.5};

        /* set transfert operators : restriction and prolongation */
        multigrid::TransfertOperator rest_Op = multigrid::PointAverage; // restriction operator
        multigrid::TransfertOperator int_Op = multigrid::Linear;  // interpolation/prolongation operator
        std::cout << "11111111111 Finished!" << std::endl;


        multigrid::Poisson<double,3> problem(init_setting,starts,ends,rest_Op,int_Op);

        std::cout << init_setting.max_iter << "\n";
        problem.solve();
        std::cout << "Finished!" << std::endl;
        return 0;
    
}