#ifndef CDWPTTCHI2RECTOOL_TRANSFORMER_TRANSFORMER_H_
#define CDWPTTCHI2RECTOOL_TRANSFORMER_TRANSFORMER_H_

#include <algorithm>

#include "utils/Configurable.hpp"
#include "utils/RecPmtProp.hpp"

class Transformer : public Configurable, public PmtTypeChecker {

public:

    Transformer(const std::string& name);
    Transformer(const std::string& name, const RecPmtType& pmt_type);

    virtual ~Transformer() = default;

    virtual void transform(RecPmtTable& table);

protected:

    RecPmtTable::iterator m_ftable; // iterator to the first target PMT
    RecPmtTable::iterator m_ltable; // iterator to the last target PMT

    bool findRange(RecPmtTable& table);
    
    virtual void transformPmt(RecPmtProp& pmt) = 0;

};

#endif // CDWPTTCHI2RECTOOL_TRANSFORMER_TRANSFORMER_H_