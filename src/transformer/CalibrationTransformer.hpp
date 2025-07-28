#ifndef CDWPTTCHI2RECTOOL_TRANSFORMER_CALIBRATIONTRANSFORMER_HPP_
#define CDWPTTCHI2RECTOOL_TRANSFORMER_CALIBRATIONTRANSFORMER_HPP_

#include "transformer/Transformer.hpp"

/**
 * @class CalibrationTransformer
 * 
 * @brief Derived class for transformer.
 * This transformer is used to shift/offset the FHT of PMTs (as a calibration)
 */
class CalibrationTransformer : public Transformer {

public:

    CalibrationTransformer(const std::string& name, const RecPmtType& type, double offset);

    ~CalibrationTransformer() override = default;

    void configure(const SniperJSON& config) override;

protected:

    double m_offset;

    void transformPmt(RecPmtProp& pmt) override;

    friend DLElement* SniperCreateDLE_T<CalibrationTransformer>(const std::string&);
    CalibrationTransformer(const std::string& name);

};

#endif // CDWPTTCHI2RECTOOL_TRANSFORMER_CALIBRATIONTRANSFORMER_HPP_