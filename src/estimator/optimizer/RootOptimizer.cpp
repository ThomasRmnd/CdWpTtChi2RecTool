#include "estimator/optimizer/RootOptimizer.hpp"

#include <Math/Factory.h>
#include <Math/Functor.h>

#include "SniperKernel/SniperLog.h"

RootOptimizer::RootOptimizer(unsigned int max_calls, unsigned int max_iter, double tol) :
    Optimizer(max_calls, max_iter, tol),
    m_opti(ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"))
{}

bool RootOptimizer::optimize_impl(const std::function<double(const double*)>& func) {
    m_opti->Clear();
    if (!m_size) {
        LogWarn << "No variable set, cannot optimize. Skipping\n";
        return false;
    }

    ROOT::Math::Functor f(func, m_size);

    m_opti->SetFunction(f);
    for (std::size_t k = 0ul; k < m_size; ++k) {
        if (m_steps[k] == 0.0) m_opti->SetFixedVariable(k, m_names[k].c_str(), m_params[k]);
        else m_opti->SetVariable(k, m_names[k].c_str(), m_params[k], m_steps[k]);
    }
    m_opti->SetMaxFunctionCalls(c_max_calls);
    m_opti->SetMaxIterations(c_max_iter);
    m_opti->SetTolerance(c_tol);
    m_opti->SetPrintLevel(0);

    m_opti->Minimize();

    const double* params = m_opti->X();
    m_params = std::vector<double>(params, params + m_size);
    m_cost = m_opti->MinValue();

    return true;
}