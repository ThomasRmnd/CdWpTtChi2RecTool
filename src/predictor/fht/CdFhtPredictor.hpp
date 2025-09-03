#ifndef CDWPTTCHI2RECTOOL_PREDICTOR_FHT_CDFHTPREDICTOR_HPP_
#define CDWPTTCHI2RECTOOL_PREDICTOR_FHT_CDFHTPREDICTOR_HPP_

#include "predictor/Predictor.hpp"
#include "predictor/fht/fht/CdFht.hpp"

/**
 * @class CdFhtPredictor
 * @tparam _ParamsTag 
 * 
 * @brief Derived class to calculate CD expected data for FHT method
 */
template<typename _ParamsTag>
class CdFhtPredictor : public Predictor<FhtMethodTag> {

    static_assert(std::is_base_of<ParamsTag, _ParamsTag>::value, "Tag must derive from ParamsTag");

public:

    CdFhtPredictor(const std::shared_ptr<CdFht<_ParamsTag>>& fht) : m_fht(fht) {}

    ~CdFhtPredictor() = default;

    ParamsType getIParamsType() override {
        return ParamsTraits<_ParamsTag>::type;
    }

    /**
     * @brief Calculate the expected data
     * 
     * @param first the beginning iterator of the experimental data 
     * @param last the end iterator of the experimental data
     * @param theo the beginning iterator of the expected data
     * @param params track parameters
     */
    void predict(const_iterator first, const_iterator last, theo_iterator theo, const double* params) override {
        m_fht->setTrack(params);
        std::transform(first, last, theo, [&](const RecPmtProp& pmt) { 
            return m_fht->calculate(pmt); // assumes pmt is in CD
        });
    }

protected:

    std::shared_ptr<CdFht<_ParamsTag>> m_fht;

};

#endif // CDWPTTCHI2RECTOOL_PREDICTOR_FHT_CDFHTPREDICTOR_HPP_