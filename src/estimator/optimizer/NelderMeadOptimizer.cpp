#include "estimator/optimizer/NelderMeadOptimizer.hpp"

#include <numeric>

#include "SniperKernel/SniperLog.h"

NelderMeadOptimizer::NelderMeadOptimizer(unsigned int max_calls, unsigned int max_iter, double tol) :
    Optimizer(max_calls, max_iter, tol),
    c_alpha(1.0),
    c_gamma(2.0),
    c_rho(0.5),
    c_sigma(0.5)
{}

NelderMeadOptimizer::NelderMeadOptimizer(unsigned int max_calls, unsigned int max_iter, double tol, double alpha, double gamma, double rho, double sigma) :
    Optimizer(max_calls, max_iter, tol),
    c_alpha(alpha),
    c_gamma(gamma),
    c_rho(rho),
    c_sigma(sigma)
{}

bool NelderMeadOptimizer::optimize_impl(const std::function<double(const double*)>& func) {
    if (!m_size) {
        LogError << "No variable set, cannot optimize. Skipping\n";
        return false;
    }

    initSimplex(func);
    order();

        unsigned int iter = 0;
        unsigned int calls = m_size + 1;
        double tol = m_costs[m_size] - m_costs[0];

        // LogDebug << "Iteration " << iter << ", function calls " << calls << ", cost " << m_costs[0] << ", tol " << tol << '\n';
        while (iter < c_max_iter && calls < c_max_calls && c_tol < tol) {
            ++iter;
            centroid();
            reflection(func);
            ++calls;
            if (m_costs[0] <= m_cost_r && m_cost_r < m_costs[m_size - 1]) {
                m_simplex[m_size] = m_x_r;
                m_costs[m_size] = m_cost_r;
            }
            else if (m_cost_r < m_costs[0]) {
                expansion(func);
                ++calls;
                if (m_cost_e < m_cost_r) {
                    m_simplex[m_size] = m_x_e;
                    m_costs[m_size] = m_cost_e;
                }
                else {
                    m_simplex[m_size] = m_x_r;
                    m_costs[m_size] = m_cost_r;
                }
            }
            else if (m_cost_r < m_costs[m_size]) {
                contraction(m_x_r, func);
                ++calls;
                if (m_cost_c < m_cost_r) {
                    m_simplex[m_size] = m_x_c;
                    m_costs[m_size] = m_cost_c;
                }
                else {
                    homothety(func);
                    calls += m_size;
                }
            }
            else {
                contraction(m_simplex[m_size], func);
                ++calls;
                if (m_cost_c < m_costs[m_size]) {
                    m_simplex[m_size] = m_x_c;
                    m_costs[m_size] = m_cost_c;
                }
                else {
                    homothety(func);
                    calls += m_size;
                }
            }
            order();
            tol = m_costs[m_size] - m_costs[0];
            /* LogDebug << "Simplex: ";
            for (std::size_t i = 0; i < m_size + 1; ++i) {
                for (std::size_t j = 0; j < m_size; ++j) {
                    std::cout << m_simplex[i][j] << ' ';
                }
                std::cout << "cost = " << m_costs[i] << '\n';
            }
            LogDebug << "Iteration " << iter << ", function calls " << calls << ", tol " << tol << '\n'; */
        }

    m_params = m_simplex[0];
    m_cost = m_costs[0];

    return true;
}

void NelderMeadOptimizer::initSimplex(const std::function<double(const double*)>& func) {
    m_simplex.assign(m_size + 1, m_params);
    for (std::size_t k = 0; k < m_size; ++k) {
        m_simplex[k + 1][k] += m_steps[k];
    }
    m_costs.assign(m_size + 1, 0.);
    for (std::size_t k = 0; k < m_size + 1; ++k) {
        m_costs[k] = func(m_simplex[k].data());
    }
    m_x_0.assign(m_size, 0.);
    m_x_r.assign(m_size, 0.);
    m_cost_r = 0.;
    m_x_e.assign(m_size, 0.);
    m_cost_e = 0.;
    m_x_c.assign(m_size, 0.);
    m_cost_c = 0.;
}

void NelderMeadOptimizer::order() {
    std::vector<std::size_t> idx(m_size + 1);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](std::size_t i, std::size_t j) { return m_costs[i] < m_costs[j]; });
    std::vector<std::vector<double>> tmp(m_size + 1, m_params);
    std::vector<double> tmp_costs(m_size + 1, 0.);
    for (std::size_t k = 0; k < m_size + 1; ++k) {
        tmp[k] = m_simplex[idx[k]];
        tmp_costs[k] = m_costs[idx[k]];
    }
    m_simplex = tmp;
    m_costs = tmp_costs;
}

void NelderMeadOptimizer::centroid() {
    m_x_0.assign(m_size, 0.);
    for (std::size_t i = 0; i < m_size; ++i) {
        for (std::size_t j = 0; j < m_size; ++j) {
            m_x_0[i] += m_simplex[j][i];
        }
        m_x_0[i] /= m_size;
    }
}

void NelderMeadOptimizer::reflection(const std::function<double(const double*)>& func) {
    for (std::size_t k = 0; k < m_size; ++k) {
        m_x_r[k] = m_x_0[k] + c_alpha * (m_x_0[k] - m_simplex[m_size][k]);
    }
    m_cost_r = func(m_x_r.data());
}

void NelderMeadOptimizer::expansion(const std::function<double(const double*)>& func) {
    for (std::size_t k = 0; k < m_size; ++k) {
        m_x_e[k] = m_x_0[k] + c_gamma * (m_x_r[k] - m_x_0[k]);
    }
    m_cost_e = func(m_x_e.data());
}

void NelderMeadOptimizer::contraction(const std::vector<double>& x, const std::function<double(const double*)>& func) {
    for (std::size_t k = 0; k < m_size; ++k) {
        m_x_c[k] = m_x_0[k] + c_rho * (x[k] - m_x_0[k]);
    }
    m_cost_c = func(m_x_c.data());
}

void NelderMeadOptimizer::homothety(const std::function<double(const double*)>& func) {
    for (std::size_t i = 1; i < m_size + 1; ++i) {
        for (std::size_t j = 0; j < m_size; ++j) {
            m_simplex[i][j] = m_simplex[0][j] + c_sigma * (m_simplex[i][j] - m_simplex[0][j]);
        }
        m_costs[i] = func(m_simplex[i].data());
    }
}