#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_PIPELINE_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_PIPELINE_HPP_

#include "estimator/Estimator.hpp"

#include <memory>
#include <variant>

#include "initializer/Initializer.hpp"
#include "transformer/Transformer.hpp"

class Pipeline : public Estimator {

public:

    using Step = std::variant<std::shared_ptr<Estimator>, std::shared_ptr<Initializer<FhtMethodTag>>, std::shared_ptr<Transformer>>;

    using Estimator::Estimator;

    ~Pipeline() override = default;

    bool initialize() override;
    bool finalize() override;

    bool estimate(RecPmtTable& table) override;

    void addStep(const Step& step);
    
    ParamsType getIParamsType() override;
    ParamsType getOParamsType() override;

protected:

    std::vector<Step> m_pipe;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_PIPELINE_HPP_