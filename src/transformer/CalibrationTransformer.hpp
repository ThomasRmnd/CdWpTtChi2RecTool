#ifndef CDWPTTCHI2RECTOOL_TRANSFORMER_CALIBRATIONTRANSFORMER_HPP_
#define CDWPTTCHI2RECTOOL_TRANSFORMER_CALIBRATIONTRANSFORMER_HPP_

#include "transformer/Transformer.hpp"

class CalibrationTransformer : public Transformer {

public:

    CalibrationTransformer(const std::string& name);
    CalibrationTransformer(const std::string& name, const RecPmtType& type, double offset);

    ~CalibrationTransformer() override = default;

protected:

    double m_offset;

    void transform(RecPmtProp& pmt) override;

};

#endif // CDWPTTCHI2RECTOOL_TRANSFORMER_CALIBRATIONTRANSFORMER_HPP_