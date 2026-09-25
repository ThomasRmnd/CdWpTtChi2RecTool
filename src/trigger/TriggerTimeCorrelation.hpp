#ifndef CDWPTTCHI2RECTOOL_UTILS_TRIGGERTIMECORRELATION_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_TRIGGERTIMECORRELATION_HPP_

#include "EvtNavigator/NavBuffer.h"

#include "utils/Configurable.hpp"
#include "utils/RecPmtProp.hpp"

class TriggerTimeCorrelation : public Configurable {

public:

    using Configurable::Configurable;

    ~TriggerTimeCorrelation() override = default;

    void configure(const SniperJSON& config) override;

    bool correct(RecPmtTable& table, JM::NavBuffer* buf);

private:

    std::pair<double, double> m_time_window;

};

#endif // CDWPTTCHI2RECTOOL_UTILS_TRIGGERTIMECORRELATION_HPP_