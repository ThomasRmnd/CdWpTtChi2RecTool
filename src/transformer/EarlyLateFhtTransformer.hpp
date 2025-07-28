#ifndef CDWPTTCHI2RECTOOL_TRANSFORMER_EARLYLATEFHTTRANSFORMER_HPP_
#define CDWPTTCHI2RECTOOL_TRANSFORMER_EARLYLATEFHTTRANSFORMER_HPP_

#include <memory>

#include <TH1D.h>

#include "transformer/Transformer.hpp"

/**
 * @class EarlyLateFhtTransformer
 * 
 * @brief Derived class for transformer.
 * This transformer is used to remove PMTs hit before/after the entry/exit of the track 
 */
class EarlyLateFhtTransformer : public Transformer {

public:

    EarlyLateFhtTransformer(const std::string& name, const RecPmtType& type, int nb_bins, double xmin, double xmax, double pmt_thold, double shift, double relative_cut);

    ~EarlyLateFhtTransformer() override = default;

    void configure(const SniperJSON& config) override;

    /**
     * @brief Transform (remove/change) all PMTs of the target PMT type
     * 
     * @param table experiemental data vector
     */
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