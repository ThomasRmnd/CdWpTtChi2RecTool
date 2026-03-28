#include "initializer/fht/AcrylicToCdProjectionInitializer.hpp"

#include "SniperKernel/SniperLog.h"

#include "utils/constants.hpp"
#include "utils/TrackParams.hpp"

AcrylicToCdProjectionInitializer::AcrylicToCdProjectionInitializer(const std::string& name, const std::shared_ptr<Initializer<FhtMethodTag>>& acrylic_init) :
    Initializer<FhtMethodTag>(name),
    m_acrylic_init(acrylic_init)
{}

void AcrylicToCdProjectionInitializer::configure(const SniperJSON& config) {
    if (!config.valid()) return;
    m_acrylic_init->configure(config);
}

ParamsType AcrylicToCdProjectionInitializer::getOParamsType() {
    if (m_acrylic_init->getOParamsType() != ParamsType::SingleAcrylic) {
        LogError << "Acrylic initializer has not parameter type SingleAcrylic\n";
        return ParamsType::Unknown;
    }
    return ParamsType::SingleCd;
}

bool AcrylicToCdProjectionInitializer::initiate(const RecPmtTable& table) {
    if (!m_acrylic_init->initiate(table)) return false;

    m_params = m_acrylic_init->getParams();
    double t_0;
    vec3 orig, dir;
    double length;
    TrackSetterHelper<SingleCdParamsTag>::set(m_params.data(), t_0, orig, dir, length);
    
    double b_half = dot(dir, orig);
    double c = mag2(orig) - constants::r_cd * constants::r_cd;
    double disc_fourth = b_half * b_half - c;
    if (disc_fourth >= 0.0) {
        double d = -b_half - std::sqrt(disc_fourth);
        orig = orig + d * dir;
    }
    m_params = std::vector<double>{t_0, theta(orig), phi(orig), theta(dir), phi(dir)};
    return true;
}