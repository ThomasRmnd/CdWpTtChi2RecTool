#ifndef CDWPTTCHI2RECTOOL_PREDICTOR_FHT_WPFHTPREDICTOR_HPP_
#define CDWPTTCHI2RECTOOL_PREDICTOR_FHT_WPFHTPREDICTOR_HPP_

#include "predictor/Predictor.hpp"

#include "predictor/fht/fht/WpFht.hpp"

template<typename _Pt>
class WpFhtPredictor : public Predictor<FhtMethodTag> {

public:

    WpFhtPredictor(const std::shared_ptr<WpFht<_Pt>>& fht) : m_fht(fht) {}

    ~WpFhtPredictor() = default;

    void predict(const_iterator first, const_iterator last, theo_iterator theo, const double* params) override {
        m_fht->setTrack(params);
        std::transform(first, last, theo, [&](const RecPmtProp& pmt) { 
            return m_fht->calculate(pmt); // assumes pmt is in WP
        });
    }

    ParamsType getIParamsType() const override {
        return _Pt;
    }

protected:

    std::shared_ptr<WpFht<_Pt>> m_fht;

};

#endif // CDWPTTCHI2RECTOOL_PREDICTOR_FHT_WPFHTPREDICTOR_HPP_