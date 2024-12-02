#ifndef _BCONDITIONS_
#define _BCONDITIONS_

#include "utilities.hpp"


namespace multigrid {

    class BoundaryConditions {

        public:
            BoundaryConditions();
            virtual ~BoundaryConditions();
            BoundaryConditions(const BoundaryConditions& bc);
            const BoundaryConditions& operator=(const BoundaryConditions& bc);

            // inline BoundaryCondition;
            
    };
}

#endif