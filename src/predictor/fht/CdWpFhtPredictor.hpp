#ifndef CDWPTTCHI2RECTOOL_PREDICTOR_FHT_CDWPFHTPREDICTOR_HPP_
#define CDWPTTCHI2RECTOOL_PREDICTOR_FHT_CDWPFHTPREDICTOR_HPP_

#include "predictor/Predictor.hpp"

#include "predictor/fht/fht/CdFht.hpp"
#include "predictor/fht/fht/WpFht.hpp"

template<ParamsType _Pt>
class CdWpFhtPredictor : public Predictor<FhtMethodTag> {

public:

    CdWpFhtPredictor(const std::shared_ptr<CdFht<_Pt>>& cd_fht, const std::shared_ptr<WpFht<_Pt>>& wp_fht) :
        m_cd_fht(cd_fht),
        m_wp_fht(wp_fht)
    {}

    ~CdWpFhtPredictor() = default;

    void predict(const_iterator first, const_iterator last, theo_iterator theo, const double* params) override {
        m_cd_fht->setTrack(params);
        m_wp_fht->setTrack(params);
        const_iterator mid = std::find_if(std::make_reverse_iterator(last), std::make_reverse_iterator(first), [&](const RecPmtProp& pmt) { return hasPmtType(pmt, RecPmtType::PMT_CD); }).base();
        std::transform(first, mid, theo, [&](const RecPmtProp& pmt) { 
            return m_cd_fht->calculate(pmt);
        });
        std::transform(mid, last, theo + std::distance(first, mid), [&](const RecPmtProp& pmt) { 
            return m_wp_fht->calculate(pmt);
        });
    }

    ParamsType getIParamsType() override {
        return _Pt;
    }

protected:

    std::shared_ptr<CdFht<_Pt>> m_cd_fht;
    std::shared_ptr<WpFht<_Pt>> m_wp_fht;

};

#endif // CDWPTTCHI2RECTOOL_PREDICTOR_FHT_CDWPFHTPREDICTOR_HPP_