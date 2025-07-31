#ifndef CDWPTTCHI2RECTOOL_STRATEGY_CLASSIFIER_CLASSIFIER_HPP_
#define CDWPTTCHI2RECTOOL_STRATEGY_CLASSIFIER_CLASSIFIER_HPP_

#include "utils/RecPmtProp.hpp"
#include "utils/TrackParams.hpp"

class Classifier {

public:

    virtual ParamsType classify(const RecPmtTable& table) = 0;

};

#endif // CDWPTTCHI2RECTOOL_STRATEGY_CLASSIFIER_CLASSIFIER_HPP_