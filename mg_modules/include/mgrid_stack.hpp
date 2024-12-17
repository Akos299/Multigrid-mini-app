#ifndef __STACK_HPP__
#define __STACK_HPP__

#include "../include/transfert_operators.hpp"
#include "../include/boundary_condition.hpp"
#include "utilities.hpp"
#include "../include/ndArray.h"
#include <vector>

namespace multigrid {

/****************** Stack class *************************
*  This is a stack of all the grid hierarchy from the coarsest level to the finest one.
*/
    /** @brief a stack to track grids at all level */
    template <typename T, int Ndim>
    class Stack : public std::vector<nd::ndArray<T, Ndim>>
    {
        public:
            Stack(const my_settings &s)
            {
                std::cout << "stack_constructor = "  << "\n";
                nb_levels    = s.nb_levels;
                nx_fine      = s.nx_fine;
                ny_fine      = s.ny_fine;
                nz_fine      = s.nz_fine;
                finest_level = (nb_levels - 1);
                // coarsest_level = 0;
                this->resize(nb_levels);
                // set finest grid size
                nd::index_t nx{nx_fine + 2*s.Nghost}, ny{ny_fine + 2*s.Nghost}, nz{nz_fine + 2*s.Nghost};
                if (Ndim == 2)
                    nz = 1;
                nd::index_t finest_size_3[3] = {nx,ny,nz}, finest_size_2[2]={nx,ny};
                nd::index_t finest_numel   = nx*ny*nz;
                std::cout << "stack_size = " << this->size() << "\n";
                // (*this)[finest_level] =     nd::ndArray<T,Ndim>(new T[finest_numel], finest_size_3, true);
                if(Ndim == 2)
                    (*this)[finest_level] =     nd::ndArray<T,Ndim>(new T[finest_numel], finest_size_2, true);
                else if(Ndim == 3)
                        (*this)[finest_level] = nd::ndArray<T,Ndim>(new T [finest_numel], finest_size_3, true);
                else
                {
                    std::cout << "###### Warning in Stack<T,Ndim>::Stack(const Settings &s) : " << "\n";
                    std::cout << "Only dimension = 2/3 can be supported." << "\n";
                }
                // Data set for all levels
                std::cout << "finest_level = " << this->finest_level;
                std::cout << "corsest_level = " << this->coarsest_level;
                for(int lev_ = finest_level-1; lev_ >= 0; --lev_)
                {
                    std::cout << " I'm in this loop lev = " << lev_ <<  "\n";
                     nx = 2 * (nx - 1) + 1;
                    ny = 2 * (ny - 1) + 1;
                    nz = (Ndim == 3) ? (2 * (nz - 1) + 1) : 0;
                    nx = (nx/2) +1;
                    ny = (ny/2) +1;
                    nz = (Ndim == 3) ? (nz/2)+1 : 0;
                    auto lev_numel = (Ndim == 2) ? (nx * ny) : (nx * ny * nz);
                    nd::index_t finest_size_2[2] = {nx, ny}, finest_size_3[3] = {nx,ny,nz};
                    if(Ndim == 2)
                        (*this)[lev_] = nd::ndArray<T,Ndim>(new T [lev_numel],finest_size_2, true);
                    else if(Ndim == 3)
                        (*this)[lev_] = nd::ndArray<T,Ndim>(new T [lev_numel], finest_size_3, true);
                    else
                    {
                        std::cout << "###### Warning in Stack<T,Ndim>::Stack(const Settings &s) : " << "\n";
                        std::cout << "Only dimension = 2/3 can be supported." << "\n";
                    }
                }
                // for (int level = this->finest_level-1; level >= this->coarsest_level; --level)
                // {
                    // std::cout << " I'm in this loop lev = " << level <<  "\n";
                    /* THIS IS GOOD FOR vertex-based multigrid-grid*/
                    // nx = 2 * (nx - 1) + 1;
                    // ny = 2 * (ny - 1) + 1;
                    // nz = (Ndim == 3) ? (2 * (nz - 1) + 1) : 0;
                    // nx = nd::index_t(nx/2) +1;
                    // ny = nd::index_t(ny/2) +1;
                    // nz = (Ndim == 3) ? (nd::index_t(nz/2)+1) : 0;
                    // auto lev_numel = (Ndim == 2) ? (nx * ny) : (nx * ny * nz);
                    // nd::index_t finest_size_2[2] = {nx, ny}, finest_size_3[3] = {nx,ny,nz};
                    // if(Ndim == 2)
                    //     (*this)[level] = nd::ndArray<T,Ndim>(new T [lev_numel],finest_size_2, true);
                    // else if(Ndim == 3)
                    //     (*this)[level] = nd::ndArray<T,Ndim>(new T [lev_numel], finest_size_3, true);
                    // else
                    // {
                    //     std::cout << "###### Warning in Stack<T,Ndim>::Stack(const Settings &s) : " << "\n";
                    //     std::cout << "Only dimension = 2/3 can be supported." << "\n";
                    // }
                // }

                // TODO Boundary conditions for all levels
                // Initialize Boundaries conditions
                // std::cout << "nb_level = " << this->finest_level;
                // for(int i = 0; i < this->nb_levels; i++)
                //     (*this)[i].set_zero();
            };
            virtual ~Stack() {};



            /** @brief function to coarse a grid from level l to level l-1 */
            // template<typename  Op>
            inline void coarsen(size_t level, TransfertOperator& Op);

            /** @brief function to coarse a grid from level l to level l-1 and store result in coarseLevelData */
            // template<typename  Op>
            inline void coarsen(size_t level, nd::ndArray<T, Ndim> &coarse_level_data, TransfertOperator& Op);

            /** @brief function to refine a grid from level l to level l+1 */
            // template<typename  Op>
            inline void refine(size_t level, TransfertOperator& Op);
    
            /** @brief function to refine a grid from level l to level l+1 and store result in fineLevelData */
            // template<typename  Op>
            inline void refine(size_t level, nd::ndArray<T, Ndim> &fine_level_data, TransfertOperator& Op);
            // BoundaryConditions bc; // boundary condition

            inline int  get_coarser_level() const {return coarsest_level;}
            inline int  get_finer_level() const {return finest_level;}

        private:
            size_t nb_levels; // number of grid levels required
            nd::index_t nx_fine, ny_fine, nz_fine; // coarse resolution
            int  finest_level;                    // finest level
           int  coarsest_level = 0; // coarsest level. for future extension coarsestLevel can be different to zero.
            
    };


    template <typename T,int Ndim>
    // template<typename Op>
    inline void Stack<T, Ndim>::coarsen(size_t level, TransfertOperator& Op)
    {
        restriction_operator<T,Ndim>((*this)[level - 1], (*this)[level],Op);
    }

    template <typename T, int Ndim>
    inline void Stack<T, Ndim>::coarsen(size_t level, nd::ndArray<T, Ndim> &coarse_level_data, TransfertOperator& Op)
    {
        restriction_operator<T,Ndim>(coarse_level_data, (*this)[level], Op);
    }

    template <typename T, int Ndim>
    //  template<typename Op>
    inline void Stack<T, Ndim>::refine(size_t level, TransfertOperator& Op)
    {
        prolongation_operator<T,Ndim>((*this)[level], (*this)[level + 1],Op);
    }

    template <typename T, int Ndim>
    // template<typename Op>
    inline void Stack<T, Ndim>::refine(size_t level, nd::ndArray<T, Ndim> &fine_level_data, TransfertOperator& Op)
    {
        prolongation_operator<T,Ndim>((*this)[level], fine_level_data, Op);
    }
}


#endif