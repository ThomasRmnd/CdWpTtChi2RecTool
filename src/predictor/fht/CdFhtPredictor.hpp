#ifndef CDWPTTCHI2RECTOOL_PREDICTOR_FHT_CDFHTPREDICTOR_HPP_
#define CDWPTTCHI2RECTOOL_PREDICTOR_FHT_CDFHTPREDICTOR_HPP_

#include "predictor/Predictor.hpp"

#include "predictor/fht/fht/CdFht.hpp"

template<typename _Pt>
class CdFhtPredictor : public Predictor<FhtMethodTag> {

public:

    CdFhtPredictor(const std::shared_ptr<CdFht<_Pt>>& fht) : m_fht(fht) {}

    ~CdFhtPredictor() = default;

    void predict(const_iterator first, const_iterator last, theo_iterator theo, const double* params) override {
        m_fht->setTrack(params);
        std::transform(first, last, theo, [&](const RecPmtProp& pmt) { 
            return m_fht->calculate(pmt); // assumes pmt is in CD
        });
    }

    ParamsType getIParamsType() const override {
        return _Pt;
    }

protected:

    std::shared_ptr<CdFht<_Pt>> m_fht;

};

#endif // CDWPTTCHI2RECTOOL_PREDICTOR_FHT_CDFHTPREDICTOR_HPP_