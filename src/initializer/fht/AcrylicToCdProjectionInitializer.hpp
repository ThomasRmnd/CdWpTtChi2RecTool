#ifndef CDWPTTCHI2RECTOOL_INITIALIZER_FHT_ACRYLICTOCDPROJECTIONINITIALIZER_HPP_
#define CDWPTTCHI2RECTOOL_INITIALIZER_FHT_ACRYLICTOCDPROJECTIONINITIALIZER_HPP_

#include <memory>

#include "initializer/Initializer.hpp"

/**
 * @class AcrylicToCdProjectionInitializer
 * 
 * @brief Derived class to calculate track parameters initial guess for FHT method
 * This initial guess is designed to project the initial guess from acrylic surface on to the CD surface for a single track 
 */
class AcrylicToCdProjectionInitializer : public Initializer<FhtMethodTag> {

public:

    AcrylicToCdProjectionInitializer(const std::string& name, const std::shared_ptr<Initializer<FhtMethodTag>>& acrylic_init);

    ~AcrylicToCdProjectionInitializer() override = default;

    void configure(const SniperJSON& config);

    ParamsType getOParamsType() override;

    /**
     * @brief Calculate the track parameters initial guess
     * 
     * @param data experimental data vector
     * 
     * @return Boolean whether the calculation is successfull or not
     */
    bool initiate(const RecPmtTable& table) override;

private:

    std::shared_ptr<Initializer<FhtMethodTag>> m_acrylic_init;

};

#endif // CDWPTTCHI2RECTOOL_INITIALIZER_FHT_ACRYLICTOCDPROJECTIONINITIALIZER_HPP_