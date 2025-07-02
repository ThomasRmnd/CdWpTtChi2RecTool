#include "estimator/Pipeline.hpp"

#include "SniperKernel/SniperLog.h"
#include "SniperKernel/ToolFactory.h"

DECLARE_TOOL(Pipeline);

bool Pipeline::initialize() {
    if (m_pipe.empty()) {
        LogError << "No steps in the pipeline\n";
        return false;
    }
    for (const Step& step : m_pipe) {
        if (std::holds_alternative<std::shared_ptr<Estimator>>(step)) {
            std::shared_ptr<Estimator> esti = std::get<std::shared_ptr<Estimator>>(step);
            if (!esti) {
                LogError << "Estimator not set\n";
                return false;
            }
            if (!esti->initialize()) return false;
        } 
        else if (std::holds_alternative<std::shared_ptr<Initializer<FhtMethodTag>>>(step)) {
            std::shared_ptr<Initializer<FhtMethodTag>> init = std::get<std::shared_ptr<Initializer<FhtMethodTag>>>(step);
            if (!init) {
                LogError << "Initializer not set\n";
                return false;
            }
            if (!init->initialize()) return false;
        } 
        else if (std::holds_alternative<std::shared_ptr<Transformer>>(step)) {
            std::shared_ptr<Transformer> transfo = std::get<std::shared_ptr<Transformer>>(step);
            if (!transfo) {
                LogError << "Transformer not set\n";
                return false;
            }
            if (!transfo->initialize()) return false;
        }
    }
    if (getIParamsType() == ParamsType::Unknown) return false;
    if (getOParamsType() == ParamsType::Unknown) return false;
    return true;
}

ParamsType Pipeline::getIParamsType() {
    for (const Step& step : m_pipe) {
        if (std::holds_alternative<std::shared_ptr<Estimator>>(step)) {
            std::shared_ptr<Estimator> esti = std::get<std::shared_ptr<Estimator>>(step);
            return esti->getIParamsType();
        }
        else if (std::holds_alternative<std::shared_ptr<Initializer<FhtMethodTag>>>(step)) {
            std::shared_ptr<Initializer<FhtMethodTag>> init = std::get<std::shared_ptr<Initializer<FhtMethodTag>>>(step);
            return init->getOParamsType();
        }
    }
    return ParamsType::Unknown;
}

ParamsType Pipeline::getOParamsType() {
    ParamsType type = ParamsType::Unknown;
    bool was_set = false;
    for (const Step& step : m_pipe) {
        if (std::holds_alternative<std::shared_ptr<Estimator>>(step)) {
            std::shared_ptr<Estimator> esti = std::get<std::shared_ptr<Estimator>>(step);
            if (!was_set) was_set = true;
            else if (type != esti->getIParamsType()) {
                LogError << "Different parameters type with one of the step in the pipeline\n";
                return ParamsType::Unknown;
            }
            type = esti->getOParamsType();
        }
        else if (std::holds_alternative<std::shared_ptr<Initializer<FhtMethodTag>>>(step)) {
            std::shared_ptr<Initializer<FhtMethodTag>> init = std::get<std::shared_ptr<Initializer<FhtMethodTag>>>(step);
            if (!was_set) was_set = true;
            else if (type != init->getOParamsType()) {
                LogError << "Different parameters type with one of the step in the pipeline\n";
                return ParamsType::Unknown; 
            }
            type = init->getOParamsType();
        }
    }
    return type;
}

bool Pipeline::estimate(RecPmtTable& table) {
    for (Step& step : m_pipe) {
        if (std::holds_alternative<std::shared_ptr<Estimator>>(step)) {
            std::shared_ptr<Estimator> esti = std::get<std::shared_ptr<Estimator>>(step);
            esti->setParams(m_params, m_steps, m_names);
            if (!esti->estimate(table)) return false;
            m_cost = esti->getCost();
            m_params = esti->getParams();
        } 
        else if (std::holds_alternative<std::shared_ptr<Initializer<FhtMethodTag>>>(step)) {
            std::shared_ptr<Initializer<FhtMethodTag>> init = std::get<std::shared_ptr<Initializer<FhtMethodTag>>>(step);
            if (!init->operator()(table)) return false;
            m_params = init->getParams();
        } 
        else if (std::holds_alternative<std::shared_ptr<Transformer>>(step)) {
            std::shared_ptr<Transformer> transfo = std::get<std::shared_ptr<Transformer>>(step);
            transfo->operator()(table);
        }
    }

    m_params = m_params;

    return true;
}

bool Pipeline::finalize() {
    for (const Step& step : m_pipe) {
        if (std::holds_alternative<std::shared_ptr<Estimator>>(step)) {
            std::shared_ptr<Estimator> esti = std::get<std::shared_ptr<Estimator>>(step);
            if (!esti->finalize()) return false;
        } 
        else if (std::holds_alternative<std::shared_ptr<Initializer<FhtMethodTag>>>(step)) {
            std::shared_ptr<Initializer<FhtMethodTag>> init = std::get<std::shared_ptr<Initializer<FhtMethodTag>>>(step);
            if (!init->finalize()) return false;
        } 
        else if (std::holds_alternative<std::shared_ptr<Transformer>>(step)) {
            std::shared_ptr<Transformer> transfo = std::get<std::shared_ptr<Transformer>>(step);
            if (!transfo->finalize()) return false;
        }
    }   
    return true;
}