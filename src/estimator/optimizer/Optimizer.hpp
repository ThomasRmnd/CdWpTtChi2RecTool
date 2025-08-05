#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_OPTIMIZER_OPTIMIZER_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_OPTIMIZER_OPTIMIZER_HPP_

#include <functional>
#include <string>
#include <vector>

/**
 * @class Optimizer
 * 
 * @brief Abstract base class for numerical optimizers.
 */
class Optimizer {

public:

    Optimizer(unsigned int max_calls, unsigned int max_iter, double tol);

    virtual ~Optimizer() = default;

    /**
     * @brief Run the optimization using the provided cost function.
     *
     * @tparam Func Callable type: must be convertible to `double(const double*)`.
     * 
     * @param func The cost function to minimize.
     * 
     * @return true if optimization succeeds, false otherwise.
     */
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

    /**
     * @brief Implementation-specific optimization logic.
     *
     * @param func Cost function.
     * 
     * @return true if optimization succeeds, false otherwise.
     */
    virtual bool optimize_impl(const std::function<double(const double*)>& func) = 0;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_OPTIMIZER_OPTIMIZER_HPP_