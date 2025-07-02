#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_FHTMINIMIZERESTIMATOR_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_FHTMINIMIZERESTIMATOR_HPP_

#include "estimator/MinimizerEstimator.hpp"

class FhtMinimizerEstimator : public MinimizerEstimator<FhtMethodTag> {

public:

    using MinimizerEstimator<FhtMethodTag>::MinimizerEstimator;

    virtual ~FhtMinimizerEstimator() = default;

    bool operator()(RecPmtTable& table) override;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_FHTMINIMIZERESTIMATOR_HPP_