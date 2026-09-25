#ifndef CDWPTTCHI2RECTOOL_UTILS_TRIGGERTIMECORRELATION_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_TRIGGERTIMECORRELATION_HPP_

#include "Context/TimeStamp.h"
#include "EvtNavigator/NavBuffer.h"

#include "utils/Configurable.hpp"
#include "utils/DetectorType.hpp"
#include "utils/RecPmtProp.hpp"

class TriggerTimeCorrelation : public Configurable {

public:

    using Configurable::Configurable;

    ~TriggerTimeCorrelation() override = default;

    void configure(const SniperJSON& config) override;

    bool correct(RecPmtTable& table, JM::NavBuffer* buf);

private:

    std::pair<double, double> m_time_window;

    TimeStamp m_cur_ts;
    TimeStamp m_other_ts;

    DetectorType m_cur_type;
    DetectorType m_other_type;

    DetectorType getDetectorType(JM::EvtNavigator* nav);
    void changeRefTime(RecPmtTable& table, double diff_ts);
    void changeRefTimeInRange(RecPmtTable& table, double diff_ts, const RecPmtType& type);

};

#endif // CDWPTTCHI2RECTOOL_UTILS_TRIGGERTIMECORRELATION_HPP_