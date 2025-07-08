#ifndef CDWPTTCHI2RECTOOL_STRATEGY_CLASSIFIER_CLASSIFIER_HPP_
#define CDWPTTCHI2RECTOOL_STRATEGY_CLASSIFIER_CLASSIFIER_HPP_

#include "strategy/Strategy.hpp"
#include "utils/RecPmtProp.hpp"

class Classifier {

public:

    virtual StrategyType classify(const RecPmtTable& table) = 0;

};

#endif // CDWPTTCHI2RECTOOL_STRATEGY_CLASSIFIER_CLASSIFIER_HPP_