#include "estimator/fht/CorrectionMapLoopEstimator.hpp"

CorrectionMapLoopEstimator::CorrectionMapLoopEstimator(const std::string& name, const std::shared_ptr<Estimator>& esti, const std::vector<std::shared_ptr<CorrectionMap>>& maps, std::size_t nb_loop) :
    Estimator{name},
    m_esti{esti},
    m_maps(maps),
    m_nb_loop{nb_loop}
{}

bool CorrectionMapLoopEstimator::initialize() {
    if (!m_esti) {
        LogError << "Estimator is nullptr\n";
        return false;
    }
    if (!m_esti->initialize()) return false;
    if (m_nb_loop == 0) LogWarn << "The number of loops is set to 0\n";
    if (m_maps.empty()) LogWarn << "There are no correction maps\n";
    for (const std::shared_ptr<CorrectionMap>& map : m_maps) {
        if (!map) {
            LogError << "One of the correction map is nullptr\n";
            return false;
        }
        if (!map->initialize()) return false;
    }
    return true;
}

bool CorrectionMapLoopEstimator::finalize() {
    if (m_esti && !m_esti->finalize()) return false;
    for (const std::shared_ptr<CorrectionMap>& map : m_maps) {
        if (map && !map->finalize()) return false;
    }
    return true;
}

ParamsType CorrectionMapLoopEstimator::getIParamsType() {
    ParamsType type = m_esti->getIParamsType();
    for (const std::shared_ptr<CorrectionMap>& map : m_maps) {
        if (type == map->getIParamsType()) continue;
        LogError << "Cost function and correction maps have different parameters type\n";
        return ParamsType::Unknown;
    }
    return type;
}

ParamsType CorrectionMapLoopEstimator::getOParamsType() {
    return getIParamsType();
}

bool CorrectionMapLoopEstimator::estimate(RecPmtTable& table) {
    if (m_params.empty()) {
        LogError << "No parameters were set as input. Abording\n";
        return false;
    }
    copyTable(table);
    for (std::size_t k = 0; k < m_nb_loop; ++k) {
        if (!applyCorrMap(table)) return false;
        m_esti->setParams(m_params, m_steps, m_names);
        if (!m_esti->estimate(table)) return false;
        m_params = m_esti->getParams();
        m_cost = m_esti->getCost();
        resetFht(table);
    }
    return true;
}

void CorrectionMapLoopEstimator::copyTable(const RecPmtTable& table) {
    m_fhts.assign(table.size(), 0.);
    for (std::size_t k = 0; k < table.size(); ++k) {
        m_fhts[k] = table[k].fht;
    }
}

bool CorrectionMapLoopEstimator::applyCorrMap(RecPmtTable& table) {
    for (const std::shared_ptr<CorrectionMap>& map : m_maps) {
        map->correct(table.begin(), table.end(), m_params.data());
    }
    return true;
}

void CorrectionMapLoopEstimator::resetFht(RecPmtTable& table) {
    for (std::size_t k = 0; k < table.size(); ++k) {
        table[k].fht = m_fhts[k];
    }
}