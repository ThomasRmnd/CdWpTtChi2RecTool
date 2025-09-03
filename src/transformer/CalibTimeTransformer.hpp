#ifndef CDWPTTCHI2RECTOOL_TRANSFORMER_CALIBTIMETRANSFORMER_HPP_
#define CDWPTTCHI2RECTOOL_TRANSFORMER_CALIBTIMETRANSFORMER_HPP_

#include "transformer/Transformer.hpp"

/**
 * @class CalibTimeTransformer
 * 
 * @brief Derived class for transformer.
 * This transformer is used to shift/offset the FHT of PMTs (as a calibration)
 */
class CalibTimeTransformer : public Transformer {

public:

    CalibTimeTransformer(const std::string& name, const RecPmtType& type, double offset);

    ~CalibTimeTransformer() override = default;

    void configure(const SniperJSON& config) override;

protected:

    double m_offset;

    void transformPmt(RecPmtProp& pmt) override;

};

#endif // CDWPTTCHI2RECTOOL_TRANSFORMER_CALIBTIMETRANSFORMER_HPP_