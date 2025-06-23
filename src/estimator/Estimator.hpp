#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_ESTIMATOR_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_ESTIMATOR_HPP_

#include "SniperKernel/ToolBase.h"
#include "utils/TrackParams.hpp"

#include <vector>

#include "utils/RecPmtProp.hpp"

class Estimator : public ToolBase, public IParamsHandler, public OParamsHandler {

public:

    using ToolBase::ToolBase;

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