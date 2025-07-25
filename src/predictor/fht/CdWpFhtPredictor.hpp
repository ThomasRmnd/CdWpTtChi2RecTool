#ifndef CDWPTTCHI2RECTOOL_PREDICTOR_FHT_CDWPFHTPREDICTOR_HPP_
#define CDWPTTCHI2RECTOOL_PREDICTOR_FHT_CDWPFHTPREDICTOR_HPP_

#include "predictor/Predictor.hpp"
#include "predictor/fht/fht/CdFht.hpp"
#include "predictor/fht/fht/WpFht.hpp"

/**
 * @class CdWpFhtPredictor
 * @tparam _ParamsTag 
 * 
 * @brief Derived class to calculate CD+WP expected data for FHT method
 */
template<typename _ParamsTag>
class CdWpFhtPredictor : public Predictor<FhtMethodTag> {

    static_assert(std::is_base_of<ParamsTag, _ParamsTag>::value, "Tag must derive from ParamsTag");

public:

    CdWpFhtPredictor(const std::shared_ptr<CdFht<_ParamsTag>>& cd_fht, const std::shared_ptr<WpFht<_ParamsTag>>& wp_fht) :
        m_cd_fht(cd_fht),
        m_wp_fht(wp_fht)
    {}

    ~CdWpFhtPredictor() = default;

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
     * 
     * @return Raw \f$ \chi^2 \f$
     */
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

protected:

    std::shared_ptr<CdFht<_ParamsTag>> m_cd_fht;
    std::shared_ptr<WpFht<_ParamsTag>> m_wp_fht;

};

#endif // CDWPTTCHI2RECTOOL_PREDICTOR_FHT_CDWPFHTPREDICTOR_HPP_