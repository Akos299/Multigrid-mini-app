#ifndef __STACK_HPP__
#define __STACK_HPP__

#include "../include/transfert_operators.hpp"
#include "../include/boundary_condition.hpp"

namespace multigrid {

/****************** Stack class *************************
*  This is a stack of all the grid hierarchy from the coarsest level to the finest one.
*/
    /** @brief a stack to track grids at all level */
    template <typename T, int Ndim>
    class Stack : public std::vector<nd::ndArray<T, Ndim>>
    {
        public:
            Stack(const Settings<T> &s)
            {
                nb_levels = s.nb_levels;
                nx_fine = s.nx_fine;
                ny_fine = s.ny_fine;
                nz_fine = s.nz_fine;
                finest_level = (nb_levels - 1);
                this->resize(nb_levels);
                // set finest grid size
                nd::index_t nx{nx_fine + 2*s.Nghost}, ny{ny_fine + 2*s.Nghost}, nz{nz_fine + 2*s.Nghost};
                if (Ndim == 2)
                    nz = 1;
                nd::index_t finest_size[3] = {nx,ny,nz};
                nd::index_t finest_numel   = nx*ny*nz;
                if(Ndim == 2)
                    (*this)[finest_level] = nd::ndArray<T,Ndim>(new T[finest_numel], {nx,ny}, true);
                else if(Ndim == 3)
                        (*this)[finest_level] = nd::ndArray<T,Ndim>(new T [finest_numel], {nx,ny,nz}, true);
                else
                {
                    std::cout << "###### Warning in Stack<T,Ndim>::Stack(const Settings &s) : " << "\n";
                    std::cout << "Only dimension = 2/3 can be supported." << "\n";
                }
                // Data set for all levels
                for (auto level = finest_level-1; level >= coarsest_level; level--)
                {
                    /* THIS IS GOOD FOR vertex-based multigrid-grid*/
                    // nx = 2 * (nx - 1) + 1;
                    // ny = 2 * (ny - 1) + 1;
                    // nz = (Ndim == 3) ? (2 * (nz - 1) + 1) : 0;
                    nx = int(nx/2) +1;
                    ny = int(ny/2) +1;
                    nz = (Ndim == 3) ? (int(nz/2)+1) : 0;
                    auto lev_numel = (Ndim == 2) ? (nx * ny) : (nx * ny * nz);
                    if(Ndim == 2)
                        (*this)[level] = nd::ndArray<T,Ndim>(new T [lev_numel], {nx,ny}, true);
                    else if(Ndim == 3)
                        (*this)[level] = nd::ndArray<T,Ndim>(new T [lev_numel], {nx,ny,nz}, true);
                    else
                    {
                        std::cout << "###### Warning in Stack<T,Ndim>::Stack(const Settings &s) : " << "\n";
                        std::cout << "Only dimension = 2/3 can be supported." << "\n";
                    }
                }

                // TODO Boundary conditions for all levels
                // Initialize Boundaries conditions

                for(auto level = coarsest_level; level < nb_levels; level++)
                    (*this)[level].set_zero();
            };
            virtual ~Stack() {};

            size_t finest_level;                    // finest level
            static const size_t coarsest_level = 0; // coarsest level. for future extension coarsestLevel can be different to zero.

            /** @brief function to coarse a grid from level l to level l-1 */
            template<TransfertOperator Op>
            inline void coarsen(size_t level);

            /** @brief function to coarse a grid from level l to level l-1 and store result in coarseLevelData */
            template<TransfertOperator Op>
            inline void coarsen(size_t level, nd::ndArray<T, Ndim> &coarse_level_data);

            /** @brief function to refine a grid from level l to level l+1 */
            template<TransfertOperator Op>
            inline void refine(size_t level);
    
            /** @brief function to refine a grid from level l to level l+1 and store result in fineLevelData */
            template<TransfertOperator Op>
            inline void refine(size_t level, nd::ndArray<T, Ndim> &fine_level_data);
            // BoundaryConditions bc; // boundary condition

        private:
            const size_t nb_levels; // number of grid levels required
            int nx_fine, ny_fine, nz_fine; // coarse resolution
    };


}


#endif