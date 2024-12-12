#include "ndArray.h"
#include <iostream>


int main()
{
    std::cout << "begin of declaration and initialization " << "\n\n";
    const nd::index_t size[1] = {5};
    const nd::index_t numel = 5;
    nd::ndArray<double,1> A(new double[numel], size, true);
    nd::ndArray<double,1> B(new double[numel], size, true);
    nd::ndArray<double,1> C(new double[numel], size, true);

    for(unsigned int i = 0; i < size[0]; i++)
    {
        A({i}) = i+1;
        B({i}) = i*2 +1;
    }

    std::cout << "end of declaration and initialization " << "\n";

    for(nd::index_t i = 0; i < size[0]; i++)
    {
        std::cout << " A[" << i << "] = " << A({i}) << "\n"; 
       
    }

    std::cout << "\n\n";
    for(nd::index_t i = 0; i < size[0]; i++)
    {
      
        std::cout << " B[" << i << "] = " << B({i}) << "\n";
       
    }
    std::cout << "\n\n\n\n";
    A.swap(B);
    for(unsigned int i = 0; i < size[0]; i++)
    {
        std::cout << " A[" << i << "] = " << A({i}) << "\n"; 
        
    }
    std::cout << "\n\n"; 
    for(nd::index_t i = 0; i < size[0]; i++)
    {
       
        std::cout << " B[" << i << "] = " << B({i}) << "\n";
       
    }
    B+=A;
    C.copy(B);
    C*=0.68;
    std::cout << "\n\n"; 
    for(nd::index_t i = 0; i < size[0]; i++)
    {
       
        std::cout << " C[" << i << "] = " << C({i}) << "\n";
       
    }

    if (!A.empty())
            std::cout << "A is not empty" << "\n"; /*test the emptyness of the array*/
    A.info();

    A.set_zero();

    for(unsigned int i = 0; i < size[0]; i++)
    {
        std::cout << " A[" << i << "] = " << A({i}) << "\n"; 
       
    }

    // A.reset();
    A.clear();
    
    if (A.empty())
        std::cout << "A is not empty" << A.ndims() << " "  << A.numel() << "\n"; /*test the emptyness of the array*/
    A.info();




}