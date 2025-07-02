#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_OPTIMIZER_OPTIMIZER_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_OPTIMIZER_OPTIMIZER_HPP_

#include <functional>
#include <string>
#include <vector>

class Optimizer {

public:

    Optimizer(unsigned int max_calls, unsigned int max_iter, double tol);

    virtual ~Optimizer() = default;

    template<typename Func> 
    bool optimize(Func&& func) {
        std::function<double(const double*)> std_func = [&func](const double* x) { return func(x); };
        return optimize_impl(std_func);
    }

    bool setParams(const std::vector<double>& ivars, const std::vector<double>& steps, const std::vector<std::string>& names);

    const std::vector<double>& getParams() const;
    double getCost() const;

    void printParams();

protected:

    const unsigned int c_max_calls;
    const unsigned int c_max_iter;
    const double c_tol;

    std::size_t m_size;
    std::vector<double> m_params;
    std::vector<double> m_steps;
    std::vector<std::string> m_names;
    double m_cost;

    virtual bool optimize_impl(const std::function<double(const double*)>& func) = 0;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_OPTIMIZER_OPTIMIZER_HPP_