#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_CORRECTIONMAPLOOPESTIMATOR_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_CORRECTIONMAPLOOPESTIMATOR_HPP_

#include "estimator/MinimizerEstimator.hpp"

#include <memory>

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolFactory.h"

#include "estimator/fht/map/CorrectionMap.hpp"

template<typename _Pt>
class CorrectionMapLoopEstimator : public MinimizerEstimator<FhtMethodTag> {

public:

    CorrectionMapLoopEstimator(const std::string& name) :
        MinimizerEstimator<FhtMethodTag>(name)
    {};

    CorrectionMapLoopEstimator(const std::string& name, const std::shared_ptr<Optimizer>& opti, const std::shared_ptr<CostFunction<FhtMethodTag>>& func, const std::vector<std::shared_ptr<CorrectionMap<_Pt>>>& maps, std::size_t nb_loop) :
        MinimizerEstimator<FhtMethodTag>(name, opti, func),
        m_maps(maps),
        c_nb_loop(nb_loop)
    {};

    virtual ~CorrectionMapLoopEstimator() override = default;

    virtual bool initialize() override {
        if (c_nb_loop == 0) LogWarn << "The number of loops is set to 0\n";
        if (m_maps.empty()) LogWarn << "There are no correction maps\n";
        for (const std::shared_ptr<CorrectionMap<_Pt>>& map : m_maps) {
            if (!map->initialize()) return false;
        }
        return MinimizerEstimator<FhtMethodTag>::initialize();
    };

    virtual bool finalize() override {
        for (const std::shared_ptr<CorrectionMap<_Pt>>& map : m_maps) {
            if (!map->finalize()) return false;
        }
        return MinimizerEstimator<FhtMethodTag>::finalize();
    };

    ParamsType getIParamsType() override {
        if (this->m_func->getIParamsType() != _Pt::type) {
            LogError << "Cost function and correction maps have different parameters type\n";
            return ParamsType::Unknown;
        }
        return _Pt::type;
    }

    ParamsType getOParamsType() override {
        return getIParamsType();
    }

    bool operator()(RecPmtTable& table) override {
        if (m_params.empty()) {
            LogError << "No parameters were set as input. Abording\n";
            return false;
        }
        
        copyTable(table);
        
        for (std::size_t k = 0; k < c_nb_loop; ++k) {
            if (!applyCorrMap(this->m_ivars, table)) return false;

            this->m_func->set(table);
            if (!this->m_opti->setParams(this->m_ivars, this->m_steps, this->m_names)) return false;
            if (!this->m_opti->operator()(*this->m_func)) return false;

            this->m_ivars = this->m_opti->getFVars();
            resetFht(table);
        }

        this->m_fvars = this->m_ivars;
        this->m_cost = this->m_opti->getCost();

        this->m_opti->printFVars();
        return true;
    }

protected:

    std::vector<std::shared_ptr<CorrectionMap<_Pt>>> m_maps;
    const std::size_t c_nb_loop;

    std::vector<double> m_fhts;

    void copyTable(const RecPmtTable& table) {
        m_fhts.assign(table.size(), 0.);
        for (std::size_t k = 0; k < table.size(); ++k) {
            m_fhts[k] = table[k].fht;
        }
    }

    bool applyCorrMap(const std::vector<double>& vars, RecPmtTable& table) {
        for (const std::shared_ptr<CorrectionMap<_Pt>>& map : m_maps) {
            map->setTrack(vars.data());
            map->operator()(table);
        }
        return true;
    }

    void resetFht(RecPmtTable& table) {
        for (std::size_t k = 0; k < table.size(); ++k) {
            table[k].fht = m_fhts[k];
        }
    };

    template<typename T>
    friend class FhtTtCorrMapEstimator;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_FHT_CORRECTIONMAPLOOPESTIMATOR_HPP_