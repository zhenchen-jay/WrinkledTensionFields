// CppNumericalSolver
#ifndef STRONGWOLFE_H
#define STRONGWOLFE_H

#include "../meta.h"

namespace cppoptlib {

/**
 * @brief Implement the strong Wolfe line search.
 * @details: This is not embeded in the library. Implemented by Zhen Chen. All the details are based on the Algorithm 3.2 && 3.3 mentioned in "Numerical Optimization" Jorge Nocedal Stephen J. Wright,  Chapter 3. Line Search Methods
 *
 * @tparam T scalar type
 * @tparam P problem type
 * @tparam Ord order of solver
 */
template<typename ProblemType>
class StrongWolfe {
    
public:
    using Scalar = typename ProblemType::Scalar;
    using TVector = typename ProblemType::TVector;
    
    static Scalar linesearch( const TVector &x, const TVector &searchDir, ProblemType &objFunc,
                             Scalar alpha_init = 1, int iterationSteps = 1000,
                             Scalar c1 = 1e-4, Scalar c2 = 0.9, // for newton type problem set to 0.9 and for CG set to 0.1
                             Scalar alpha_min = 0, Scalar alpha_max = 1e16
                             )
    // Algorithm 3.2 mentioned in the book Numerical Optimization
    {
        if (alpha_init <= 0)
        {
            std::cout<<"init step = "<<alpha_init<<", which is less than 0."<<std::endl;
            return -1;
        }
        Scalar alpha1 = alpha_min, alpha2 = alpha_init;
        
        TVector grad;
        objFunc.gradient(x, grad);
        const Scalar phi0 = objFunc.value(x);
        const Scalar dphi0 = grad.dot(searchDir);
        
        Scalar phi1 = phi0, phi2, dphi1 = dphi0, dphi2;
        
        if(dphi0 > 0)
        {
            std::cout<<"Not a descent direction, maybe numerical inaccurate gradient!!"<<std::endl;
            return -1;
        }
        
        int i = 0;
        while(true)
        {
            TVector x2 = x + alpha2 * searchDir, grad2;
            phi2  = objFunc.value(x2);
            if(phi2 > phi0 + c1 * alpha2 * dphi0 || (phi2 >= phi1 && i > 0))
                return zoom(x, searchDir, objFunc, iterationSteps, c1, c2,
                            alpha1, alpha2,
                            phi1, phi2, phi0,
                            dphi1, dphi2, dphi0);
            
            objFunc.gradient(x2, grad2);
            dphi2 = grad2.dot(searchDir);
            
            if(std::abs(dphi2) <= - c2 * dphi0)
                return alpha2;
            
            if(dphi2 >= 0)
                return zoom(x, searchDir, objFunc, iterationSteps, c1, c2,
                            alpha2, alpha1,
                            phi2, phi1, phi0,
                            dphi2, dphi1, dphi0);
            
            alpha1 = alpha2;
            alpha2 = 2 * alpha2;    // increase alpha2;
            
            i++;
            
            if(i > iterationSteps)
            {
                std::cout<<"Reach the maximum iteration step during the line search: "<<iterationSteps<<", maybe numerical inaccurate gradient."<<std::endl;
                return -1;
            }
            
            if(alpha2 > alpha_max)
            {
                std::cout<<"Reach the largest time step: "<<alpha_max<<", maybe numerical inaccurate gradient."<<std::endl;
                return -1;
            }
        }
    }
    
    static Scalar zoom(const TVector &x, const TVector &searchDir, ProblemType &objFunc,
                       const int iterationSteps,
                       const Scalar c1, const Scalar c2, // for newton type problem set to 0.9 and for CG set to 0.1
                       Scalar alpha_low, Scalar alpha_high,
                       Scalar phi_low, Scalar phi_high, const Scalar phi0,
                       Scalar dphi_low, Scalar dphi_high, const Scalar dphi0)
    {
        int i = 0;
        while(true)
        {
            // we use cubic interpolation to obtain a trial between alpha_low and alpha_high
            Scalar d1 = dphi_low + dphi_high - 3 * (phi_low - phi_high) / (alpha_low - alpha_high);
            Scalar d2 = (d1 * d1 - dphi_low * dphi_high);
            assert(d2 >= 0);
            d2 = std::sqrt(d2);
            Scalar alpha = alpha_high - (alpha_high - alpha_low) * (dphi_high + d2 - d1) / (dphi_high - dphi_low + 2 * d2);
            
            if( alpha <= std::min(alpha_low, alpha_high) ||
               alpha >= std::max(alpha_low, alpha_high))
                alpha = (alpha_low + alpha_high) / 2.0;
            
            TVector x1 = x + alpha * searchDir;
            TVector grad1;
            Scalar phi = objFunc.value(x1);
            objFunc.gradient(x1, grad1);
            Scalar dphi = grad1.dot(searchDir);
            
            if(phi > phi0 + c1 * alpha * dphi0 || phi >= phi_low)
            {
                if(alpha == alpha_high)
                {
                    if(phi <= phi0 + c1 * alpha * dphi0)
                        return alpha;
                    else
                    {
                        std::cout<<"line search failed, maybe numerical inaccurate gradient."<<std::endl;
                        return -1;
                    }
                }
                alpha_high = alpha;
                phi_high = phi;
                dphi_high = dphi;
            }
            else
            {
                if(std::abs(dphi) <= -c2 * dphi0)
                    return alpha;
                if(dphi * (alpha_high - alpha_low) >= 0)
                {
                    alpha_high = alpha_low;
                    phi_high = phi_low;
                    dphi_high = dphi_low;
                }
                
                 if(alpha == alpha_low)
                {
                    if(phi <= phi0 + c1 * alpha * dphi0)
                        return alpha;
                    else
                    {
                        std::cout<<"line search failed, maybe numerical inaccurate gradient."<<std::endl;
                        return -1;
                    }
                }
                
                alpha_low = alpha;
                phi_low = phi;
                dphi_low = dphi;
            }
            i++;
            if(i > iterationSteps)
            {
                std::cout<<"Reach the maximum iteration step during the zoom process in line search: "<<iterationSteps<<", maybe numerical inaccurate gradient."<<std::endl;
                return -1;
            }
            
            if(std::abs(alpha_low - alpha_high) < 1e-15)
            {
                if(phi <= phi0 + c1 * alpha * dphi0)
                    return alpha;
                else
                {
                    std::cout<<"line search failed, maybe numerical inaccurate gradient."<<std::endl;
                    return -1;
                }
            }
            
        }
        
    }
    
};
}

#endif /* STRONGWOLFE_H */
