#ifndef __MULTIGRID__
#define __MULTIGRID__

#include <vector>
#include <array>
#include "utilities.hpp"
#include "Bcondition.hpp"

namespace multigrid
{

    template <typename T, int dim>
    inline void restriction_operator(std::vector<std::vector<T>> &corseData, std::vector<std::vector<T>> &fineData)
    {
    }

    template <typename T, int dim>
    inline void interpolation_operator(std::vector<std::vector<T>> &corseData, std::vector<std::vector<T>> &fineData)
    {
    }

    template <typename T>
    struct Settings
    {
        // T aspectRatio{1};
        size_t numberOfGrids{8};
        size_t minimumResolution{4};
        T residualTolerance{1e-10};
        size_t maxIterations{400};
        CycleType cycleType{multigrid::wCycle};
        size_t Npre{1}, Npost{1}, maxLevel{21}, Nvcyle_cvg{20};
    };

    template <typename T, int dim>
    class Stack : public std::vector<std::vector<std::vector<T>>>
    {
    public:
        Stack(const Settings<T> &s)
        {
            nGrids = s.numberOfGrids;
            mini_resolution = s.minimumResolution;
            finestLevel(nGrids - 1);
            resize(nGrids);
            int nx{mini_resolution}, ny{mini_resolution}, nz{mini_resolution};
            if (dim == 2)
                nz = 0;

            // Data set for all levels
            for (auto level = 1; level <= finestLevel; level++)
            {
                nx = 2 * (nx - 1) + 1;
                ny = 2 * (ny - 1) + 1;
                nz = (dim == 3) ? (2 * (nz - 1) + 1) : 0;
                (*this)[level] = std::vector<std::vector<T>>(std::vector<T>(0, nx), ny);
            }

            // TODO Boundary conditions for all levels
        };
        virtual ~Stack() {};

        const size_t finestLevel;
        static const size_t coarsestLevel = 0;

        inline void coarsen(size_t level);
        inline void coarsen(size_t level, std::vector<std::vector<T>> &coarseLevelData);
        inline void refine(size_t level);
        inline void refine(size_t level, std::vector<std::vector<T>> &fineLevelData);

        BoundaryConditions bc;

    private:
        const size_t nGrids; // number of grid levels required
        int mini_resolution;
    };

    template <typename T, int dim>
    inline void Stack<T, dim>::coarsen(size_t level)
    {
        restriction_operator<T, dim>((*this)[level - 1], (*this)[level]);
    }

    template <typename T, int dim>
    inline void Stack<T, dim>::coarsen(size_t level, std::vector<std::vector<T>> &coarseLevelData)
    {
        restriction_operator(coarseLevelData, (*this)[level]);
    }

    template <typename T, int dim>
    inline void Stack<T, dim>::refine(size_t level)
    {
        interpolation_operator((*this)[level], (*this)[level + 1]);
    }

    template <typename T, int dim>
    inline void Stack<T, dim>::refine(size_t level, std::vector<std::vector<T>> &fineLevelData)
    {
        interpolation_operator((*this)[level], fineLevelData);
    }

    template <typename T, int dim>
    class MultigridBase
    {
    public:
        MultigridBase(const Settings<T> &settings);
        virtual ~MultigridBase() {};

        inline std::vector<std::vector<T>> &get_result();
        template <typename U>
        inline void initial_guess(U arg);
        inline std::vector<std::vector<T>> &source_term();
        template <typename U>
        inline void source_term(U arg);

        inline void evaluate_operator(size_t level, std::vector<std::vector<T>> &result);
        inline void evaluate_residual(size_t level, std::vector<std::vector<T>> &result);

        void relax(const size_t level, size_t N);
        void relax(const size_t level, const double tolerance);

        virtual inline void multigrid() {}
        virtual inline void solve() { multigrid(); }

        template <typename T>
        virtual T differential_operator(size_t, int, int) = 0;
        template <typename T>
        virtual T differential_operator(size_t, int, int, int) = 0;
        template <typename T>
        virtual void relaxation_updater(size_t, int, int) = 0;
        template <typename T>
        virtual void relaxation_updater(size_t, int, int, int) = 0;

    protected:
        Stack<T, dim> solution, source; // Grids for solution and source term
        Stack temp;                     // Extra storage for muligrid solver
        CycleType cycle_type;
        size_t Npre;
        size_t Npost;
        const size_t maxIter;
        const T residualTolerance;
        size_t finestLevel, coarsestLevel;
        size_t nxfine, nyfine, nzfine;
        bool isSourceSet;
        bool isInitialValueSet;

    private:
        T residualSum, normSum;
    };

    template <typename T, int dim>
    MultigridBase<T, dim>::MultigridBase(const Settings<T> &settings) :

                                                                        solution(settings),
                                                                        source(settings),
                                                                        temp(settings),
                                                                        cycle_type(settings.cycleType),
                                                                        Npre(settings.Npre),
                                                                        Npost(settings.Npost),
                                                                        residualTolerance(settings.residualTolerance),
                                                                        maxIter(settings.residualTolerance),
                                                                        isSourceSet(false),
                                                                        isInitialValueSet(false)
    {
        finestLevel = solution.finestLevel;
        coarsestLevel = solution.coarsestLevel;
        if (dim == 2)
        {
            nxfine = solution[finestLevel][0];
            nyfine = solution[finestLevel][0][0];
        }
        // else if (dim == 3)
        // {
        //     nxfine = solution[finestLevel][0];
        //     nyfine = solution[finestLevel][0][0];
        //     nzfine = solution[finestLevel][0][0][0];
        // }
    }

    template <typename T, int dim>
    inline std::vector<std::vector<T>> &MultigridBase<T, dim>::get_result()
    {
        return solution[finestLevel];
    }

    template <typename T, int dim>
    template <typename U>
    inline void MultigridBase<T, dim>::initial_guess(U arg)
    {
        solution[finestLevel] = arg;
        isInitialValueSet = true;
    }

    template <typename T, int dim>
    inline std::vector<std::vector<T>> &MultigridBase<T, dim>::source_term()
    {
        return source[finestLevel];
    }

    template <typename T, int dim>
    template <typename U>
    inline void MultigridBase<T, dim>::source_term(U arg)
    {
        source[finestLevel] = arg;
        isSourceSet = true;
    }

    template <typename T, int dim>
    inline void MultigridBase<T, dim>::evaluate_operator(size_t level, std::vector<std::vector<T>> &result)
    {
        size_t nx = result.size();
        size_t ny = result[0].size();
        for (auto i = 0; i < nx; i++)
        {
            for (auto j = 0; j < ny; j++)
            {
                result[i][j] = differential_operator(level, i, j);
            }
        }
    }

    template <typename T, int dim>
    inline void MultigridBase<T, dim>::evaluate_residual(size_t level, std::vector<std::vector<T>> &result){
        size_t nx = result.size();
        size_t ny = result[0].size();
        for (auto i = 0; i < nx; i++)
        {
            for (auto j = 0; j < ny; j++)
            {
                result[i][j] = source[level][i][j] - differential_operator(level, i, j);
            }
        }
    }

}

#endif