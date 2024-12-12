#ifndef __STACK__
#define __STACK__

#include "transfert_operators.hpp"
#include "boundary_condition.hpp"

namespace multigrid {

/****************** Stack class *************************
*  This is a stack of all the grid hierarchy from the coarsest level to the finest one.
*/
    /** @brief a stack to track grids at all level */
    template <typename T, int Ndim>
    class Stack : public std::vector<ndArray<T, Ndim>>
    {
        public:
            Stack(const Settings<T> &s)
            {
                nb_levels = s.nb_levels;
                nx_coarse = s.nx_coarse;
                ny_coarse = s.ny_coarse;
                nz_coarse = s.nz_coarse;
                finestLevel = (nb_levels - 1);
                this->resize(nb_levels);
                // set coarsest grid size
                int nx{nx_coarse}, ny{ny_coarse}, nz{nz_coarse};
                if (Ndim == 2)
                    nz = 1;
                nd::index_t coarse_size[3] = {nx,ny,nz};
                nd::index_t coarse_numel   = nx*ny*nz;
                (*this)[0] = ndArray<T,Ndim>(new T[coarse_numel], coarse_size, true);

                // Data set for all levels
                for (auto level = 1; level <= finestLevel; level++)
                {
                    /* THIS IS GOOD FOR vertex-based multigrid-grid*/
                    // nx = 2 * (nx - 1) + 1;
                    // ny = 2 * (ny - 1) + 1;
                    // nz = (Ndim == 3) ? (2 * (nz - 1) + 1) : 0;
                    nx = 2 * nx;
                    ny = 2 * ny;
                    nz = (Ndim == 3) ? (2 * nz) : 0;
                    auto lev_numel = (Ndim == 2) ? (nx * ny) : (nx * ny * nz);
                    if(Ndim == 2)
                        (*this)[0] = ndArray<T,Ndim>(new T [lev_numel], {nx,ny}, true);
                    else if(Ndim == 3)
                        (*this)[0] = ndArray<T,Ndim>(new T [lev_numel], {nx,ny,nz}, true);
                    else
                    {
                        std::cout << "###### Warning in Stack<T,Ndim>::Stack(const Settings &s) : " << "\n";
                        std::cout << "Only dimension = 2/3 can be supported." << "\n";
                    }
                }

                // TODO Boundary conditions for all levels
                // Initialize Boundaries conditions
            };
            virtual ~Stack() {};

            size_t finest_level;                    // finest level
            static const size_t coarsest_level = 0; // coarsest level. for future extension coarsestLevel can be different to zero.

            /** @brief function to coarse a grid from level l to level l-1 */
            template<TransfertOperator Op>
            inline void coarsen(size_t level);

            /** @brief function to coarse a grid from level l to level l-1 and store result in coarseLevelData */
            template<TransfertOperator Op>
            inline void coarsen(size_t level, ndArray<T, Ndim> &coarse_level_data);

            /** @brief function to refine a grid from level l to level l+1 */
            template<TransfertOperator Op>
            inline void refine(size_t level);
    
            /** @brief function to refine a grid from level l to level l+1 and store result in fineLevelData */
            template<TransfertOperator Op>
            inline void refine(size_t level, ndArray<T, Ndim> &fine_level_data);
            // BoundaryConditions bc; // boundary condition

        private:
            const size_t nb_levels; // number of grid levels required
            int nx_coarse, ny_coarse, nz_coarse; // coarse resolution
    };


    template <typename T,int Ndim>
    template<TransfertOperator Op>
    inline void Stack<T, Ndim>::coarsen(size_t level)
    {
        restriction_operator<T, Ndim>((*this)[level - 1], (*this)[level]);
    }

    template <typename T, int Ndim>
     template<TransfertOperator Op>
    inline void Stack<T, Ndim>::coarsen(size_t level, ndArray<T, Ndim> &coarse_level_data)
    {
        restriction_operator(coarse_level_data, (*this)[level]);
    }

    template <typename T, int Ndim>
     template<TransfertOperator Op>
    inline void Stack<T, Ndim>::refine(size_t level)
    {
        prolongation_operator((*this)[level], (*this)[level + 1]);
    }

    template <typename T, int Ndim>
    template<TransfertOperator Op>
    inline void Stack<T, Ndim>::refine(size_t level, ndArray<T, Ndim> &fine_level_data)
    {
        prolongation_operator((*this)[level], fine_level_data);
    }
}


#endif