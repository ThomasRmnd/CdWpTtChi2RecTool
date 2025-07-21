#ifndef CDWPTTCHI2RECTOOL_INITIALIZER_TT_MINMAXHEIGHTINITIALIZER_HPP_
#define CDWPTTCHI2RECTOOL_INITIALIZER_TT_MINMAXHEIGHTINITIALIZER_HPP_

#include "initializer/Initializer.hpp"

class MinMaxHeightInitializer : public Initializer<TtMethodTag> {

public:

    using Initializer<TtMethodTag>::Initializer;
    
    ~MinMaxHeightInitializer() override = default;

    bool initiate(const vector_type& data) override;

    ParamsType getOParamsType() override;

};

#endif // CDWPTTCHI2RECTOOL_INITIALIZER_TT_MINMAXHEIGHTINITIALIZER_HPP_