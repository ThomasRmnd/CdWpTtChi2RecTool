#ifndef CDWPTTCHI2RECTOOL_TRANSFORMER_PMTTYPETRANSFORMER_HPP_
#define CDWPTTCHI2RECTOOL_TRANSFORMER_PMTTYPETRANSFORMER_HPP_

#include "transformer/Transformer.hpp"

/**
 * @class PmtTypeTransformer
 * 
 * @brief Derived class for transformer.
 * This transformer is used to remove PMTs of a certain type
 */
class PmtTypeTransformer : public Transformer {

public:

    PmtTypeTransformer(const std::string& name, const RecPmtType& pmt_type);
    
    ~PmtTypeTransformer() override = default;

protected:

    void transformPmt(RecPmtProp& pmt) override;

};

#endif // CDWPTTCHI2RECTOOL_TRANSFORMER_PMTTYPETRANSFORMER_HPP_