#ifndef CDWPTTCHI2RECTOOL_TRANSFORMER_FHTCHARGETHOLDTRANSFORMER_HPP_
#define CDWPTTCHI2RECTOOL_TRANSFORMER_FHTCHARGETHOLDTRANSFORMER_HPP_

#include "transformer/Transformer.hpp"

/**
 * @class FhtChargeTholdTransformer
 * 
 * @brief Derived class for transformer.
 * This transformer is used to remove low FHT charge PMTs  
 */
class FhtChargeTholdTransformer : public Transformer {

public:

    FhtChargeTholdTransformer(const std::string& name);
    FhtChargeTholdTransformer(const std::string& name, const RecPmtType& type, double q_thold);

    ~FhtChargeTholdTransformer() override = default;

    void configure(const SniperJSON& config);

protected:

    double m_q_thold;

    void transformPmt(RecPmtProp& pmt) override;

};

#endif // CDWPTTCHI2RECTOOL_TRANSFORMER_FHTCHARGETHOLDTRANSFORMER_HPP_