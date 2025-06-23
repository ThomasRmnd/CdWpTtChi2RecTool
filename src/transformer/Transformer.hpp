#ifndef CDWPTTCHI2RECTOOL_TRANSFORMER_TRANSFORMER_H_
#define CDWPTTCHI2RECTOOL_TRANSFORMER_TRANSFORMER_H_

#include "SniperKernel/ToolBase.h"

#include <algorithm>

#include "utils/RecPmtProp.hpp"

class Transformer : public ToolBase {

public:

    Transformer(const std::string& name);
    Transformer(const std::string& name, const RecPmtType& type);

    virtual ~Transformer() = default;

    virtual void operator()(RecPmtTable& table);

protected:

    const RecPmtType c_type; // target type
    RecPmtTable::iterator m_ftable; // iterator to the first target PMT
    RecPmtTable::iterator m_ltable; // iterator to the last target PMT

    bool checkType(const RecPmtProp& pmt);
    bool findRange(RecPmtTable& table);
    
    virtual void transform(RecPmtProp& pmt) = 0;

};

#endif // CDWPTTCHI2RECTOOL_TRANSFORMER_TRANSFORMER_H_