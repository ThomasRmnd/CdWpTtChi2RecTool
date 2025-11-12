#ifndef CDWPTTCHI2RECTOOL_TRANSFORMER_TRANSFORMER_HPP_
#define CDWPTTCHI2RECTOOL_TRANSFORMER_TRANSFORMER_HPP_

#include <algorithm>

#include "utils/Configurable.hpp"
#include "utils/RecPmtProp.hpp"

/**
 * @class Transformer
 * 
 * @brief Base class for transformer (filters, correctors, etc)
 */
class Transformer : public Configurable, public PmtTypeChecker {

public:

    Transformer(const std::string& name, const RecPmtType& pmt_type);

    virtual ~Transformer() = default;

    /**
     * @brief Transform (remove/change) all PMTs of the target PMT type
     * 
     * @param table experiemental data vector
     */
    virtual void transform(RecPmtTable& table);

protected:

    RecPmtTable::iterator m_ftable; // iterator to the first target PMT
    RecPmtTable::iterator m_ltable; // iterator to the last target PMT

    bool findRange(RecPmtTable& table);
    
    virtual void transformPmt(RecPmtProp& pmt) = 0;

};

#endif // CDWPTTCHI2RECTOOL_TRANSFORMER_TRANSFORMER_HPP_