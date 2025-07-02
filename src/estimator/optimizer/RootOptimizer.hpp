#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_OPTIMIZER_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_OPTIMIZER_HPP_

#include "estimator/optimizer/Optimizer.hpp"

#include <Math/Minimizer.h>

class RootOptimizer : public Optimizer {

public:

    RootOptimizer(unsigned max_calls, unsigned int max_iter, double tol);

    ~RootOptimizer() override = default;

private:

    ROOT::Math::Minimizer* m_opti;

    bool optimize_impl(const std::function<double(const double*)>& func) override;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_OPTIMIZER_HPP_