#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_ESTIMATOR_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_ESTIMATOR_HPP_

#include <vector>

#include "utils/Configurable.hpp"
#include "utils/RecPmtProp.hpp"
#include "utils/TrackParams.hpp"

class Estimator : public ToolBase, public Configurable, public IParamsHandler, public OParamsHandler {

public:

    using Configurable::Configurable;

    virtual ~Estimator() = default;

    virtual bool estimate(RecPmtTable& table) = 0;

    bool setParams(const std::vector<double>& params, const std::vector<double>& steps, const std::vector<std::string>& names);

    const std::vector<double>& getParams() const;
    double getCost() const;

protected:

    std::size_t m_size;
    std::vector<double> m_params;
    std::vector<double> m_steps;
    std::vector<std::string> m_names;
    double m_cost;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_ESTIMATOR_HPP_