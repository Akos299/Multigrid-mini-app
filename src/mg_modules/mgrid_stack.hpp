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
                nGrids = s.numberOfGrids;
                x_mini_resolution = s.x_minimumResolution;
                y_mini_resolution = s.y_minimumResolution;
                z_mini_resolution = s.z_minimumResolution;
                finestLevel = (nGrids - 1);
                this->resize(nGrids);
                // set coarsest grid size
                int nx{x_mini_resolution}, ny{y_mini_resolution}, nz{z_mini_resolution};
                if (Ndim == 2)
                    nz = 0;

                // Data set for all levels
                for (auto level = 1; level <= finestLevel; level++)
                {
                    nx = 2 * (nx - 1) + 1;
                    ny = 2 * (ny - 1) + 1;
                    nz = (Ndim == 3) ? (2 * (nz - 1) + 1) : 0;
                    auto lin_size = (Ndim == 1) ? (nx * ny) : (nx * ny * nz);
                    std::vector<T> tmp_init_data(0, lin_size);
                    (*this)[level] = ndArray<T, Ndim>(tmp_init_data.data(), {nx, ny, 1}, true);
                }

                // TODO Boundary conditions for all levels
            };
            virtual ~Stack() {};

            size_t finestLevel;                    // finest level
            static const size_t coarsestLevel = 0; // coarsest level

            /** @brief function to coarse a grid from level l to level l-1 */
            template<TransfertOperator Op>
            inline void coarsen(size_t level);

            /** @brief function to coarse a grid from level l to level l-1 and store result in coarseLevelData */
            template<TransfertOperator Op>
            inline void coarsen(size_t level, ndArray<T, Ndim> &coarseLevelData);

            /** @brief function to refine a grid from level l to level l+1 */
            template<TransfertOperator Op>
            inline void refine(size_t level);
    
            /** @brief function to coarse a grid from level l to level l-1 and store result in fineLevelData */
            template<TransfertOperator Op>
            inline void refine(size_t level, ndArray<T, Ndim> &fineLevelData);
            BoundaryConditions bc; // boundary condition

        private:
            const size_t nGrids; // number of grid levels required
            int x_mini_resolution, y_mini_resolution, z_mini_resolution; // minimum resolution
    };


    template <typename T,int Ndim>
    template<TransfertOperator Op>
    inline void Stack<T, Ndim>::coarsen(size_t level)
    {
        restriction_operator<T, Ndim>((*this)[level - 1], (*this)[level]);
    }

    template <typename T, int Ndim>
     template<TransfertOperator Op>
    inline void Stack<T, Ndim>::coarsen(size_t level, ndArray<T, Ndim> &coarseLevelData)
    {
        restriction_operator(coarseLevelData, (*this)[level]);
    }

    template <typename T, int Ndim>
     template<TransfertOperator Op>
    inline void Stack<T, Ndim>::refine(size_t level)
    {
        prolongation_operator((*this)[level], (*this)[level + 1]);
    }

    template <typename T, int Ndim>
    template<TransfertOperator Op>
    inline void Stack<T, Ndim>::refine(size_t level, ndArray<T, Ndim> &fineLevelData)
    {
        prolongation_operator((*this)[level], fineLevelData);
    }
}


#endif