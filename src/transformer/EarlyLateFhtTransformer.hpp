#ifndef CDWPTTCHI2RECTOOL_TRANSFORMER_EARLYLATEFHTTRANSFORMER_HPP_
#define CDWPTTCHI2RECTOOL_TRANSFORMER_EARLYLATEFHTTRANSFORMER_HPP_

#include "transformer/Transformer.hpp"

#include <memory>

#include <TH1D.h>

class EarlyLateFhtTransformer : public Transformer {

public:

    EarlyLateFhtTransformer(const std::string& name);
    EarlyLateFhtTransformer(const std::string& name, const RecPmtType& type, int nb_bins, double xmin, double xmax, double pmt_thold, double shift, double relative_cut);

    ~EarlyLateFhtTransformer() override = default;

    void transform(RecPmtTable& table) override;

protected:

    std::unique_ptr<TH1D> m_hist;

    double m_pmt_thold;
    double m_shift;
    double m_relative_cut;
    double m_itime;

    double getITime(const RecPmtTable& table);
    void transformPmt(RecPmtProp& pmt) override;

};

#endif // CDWPTTCHI2RECTOOL_TRANSFORMER_EARLYLATEFHTTRANSFORMER_HPP_