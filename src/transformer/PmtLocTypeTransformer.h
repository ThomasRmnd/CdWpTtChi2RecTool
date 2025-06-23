#ifndef CDWPTTCHI2RECTOOL_TRANSFORMER_PMTLOCTYPETRANSFORMER_HPP_
#define CDWPTTCHI2RECTOOL_TRANSFORMER_PMTLOCTYPETRANSFORMER_HPP_

#include "transformer/Transformer.hpp"

class PmtLocTypeTransformer : public Transformer {

public:

    using Transformer::Transformer;
    
    ~PmtLocTypeTransformer() override = default;

protected:

    void transform(RecPmtProp& pmt) override;

};

#endif // CDWPTTCHI2RECTOOL_TRANSFORMER_PMTLOCTYPETRANSFORMER_HPP_