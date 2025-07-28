#ifndef CDWPTTCHI2RECTOOL_TRANSFORMER_PMTLOCTYPETRANSFORMER_HPP_
#define CDWPTTCHI2RECTOOL_TRANSFORMER_PMTLOCTYPETRANSFORMER_HPP_

#include "transformer/Transformer.hpp"

/**
 * @class PmtLocTypeTransformer
 * 
 * @brief Derived class for transformer.
 * This transformer is used to remove PMTs of a certain type
 */
class PmtLocTypeTransformer : public Transformer {

public:

    using Transformer::Transformer;
    
    ~PmtLocTypeTransformer() override = default;

protected:

    void transformPmt(RecPmtProp& pmt) override;

};

#endif // CDWPTTCHI2RECTOOL_TRANSFORMER_PMTLOCTYPETRANSFORMER_HPP_