#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_OPTIMIZER_NELDERMEADOPTIMIZER_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_OPTIMIZER_NELDERMEADOPTIMIZER_HPP_

#include "estimator/optimizer/Optimizer.hpp"

class NelderMeadOptimizer : public Optimizer {

public:

    NelderMeadOptimizer(unsigned int max_calls, unsigned int max_iter, double tol);
    NelderMeadOptimizer(unsigned int max_calls, unsigned int max_iter, double tol, double alpha, double gamma, double rho, double sigma);

    ~NelderMeadOptimizer() override = default;

private:

    const double c_alpha, c_gamma, c_rho, c_sigma;

    std::vector<double(const double*)> m_func;
    
    std::vector<std::vector<double>> m_simplex;
    std::vector<double> m_costs;
    std::vector<double> m_x_0;

    std::vector<double> m_x_r;
    double m_cost_r;
    std::vector<double> m_x_e;
    double m_cost_e;
    std::vector<double> m_x_c;
    double m_cost_c; 

    void initSimplex(const std::function<double(const double*)>& func);
    void order();
    void centroid();
    void reflection(const std::function<double(const double*)>& func);
    void expansion(const std::function<double(const double*)>& func);
    void contraction(const std::vector<double>& x, const std::function<double(const double*)>& func);
    void homothety(const std::function<double(const double*)>& func);
    bool optimize_impl(const std::function<double(const double*)>& func) override;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_OPTIMIZER_NELDERMEADOPTIMIZER_HPP_