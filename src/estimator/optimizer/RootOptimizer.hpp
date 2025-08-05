#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_OPTIMIZER_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_OPTIMIZER_HPP_

#include <Math/Minimizer.h>

#include "estimator/optimizer/Optimizer.hpp"

/**
 * @class RootOptimizer
 * 
 * @brief Optimizer implementation using ROOT's Minuit2 Migrad algorithm.
 * 
 * This class wraps the ROOT Minuit2 minimizer interface, providing an
 * implementation of the abstract `Optimizer` base class.
 */
class RootOptimizer : public Optimizer {

public:

    RootOptimizer(unsigned max_calls, unsigned int max_iter, double tol);

    ~RootOptimizer() override = default;

private:

    ROOT::Math::Minimizer* m_opti;

    /**
     * @brief Implementation-specific optimization logic.
     *
     * @param func Cost function.
     * 
     * @return true if optimization succeeds, false otherwise.
     */
    bool optimize_impl(const std::function<double(const double*)>& func) override;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_OPTIMIZER_HPP_