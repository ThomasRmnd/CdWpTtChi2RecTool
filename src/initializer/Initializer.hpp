#ifndef CDWPTTCHI2RECTOOL_INITIALIZER_INITIALIZER_HPP_
#define CDWPTTCHI2RECTOOL_INITIALIZER_INITIALIZER_HPP_

#include "utils/Configurable.hpp"
#include "utils/Method.hpp"
#include "utils/TrackParams.hpp"

/**
 * @class Initializer
 * @tparam _MethodTag Method tag
 * 
 * @brief Base class to calculate track parameters initial guess
 */
template<typename _MethodTag>
class Initializer : public OParamsHandler, public Configurable {

    static_assert(std::is_base_of<MethodTag, _MethodTag>::value, "Tag must be derived from MethodTag");

public:

    typedef typename MethodTraits<_MethodTag>::vector_type vector_type;

    Initializer(const std::string& name) : Configurable{name} {}

    virtual ~Initializer() = default;

    const std::vector<double>& getParams() const {
        return m_params;
    }

    void printParams() {
        LogInfo << "Initial parameters: ";
        for (std::size_t k = 0; k < m_params.size() - 1; ++k) {
            std::cout << m_params[k] << ", ";
        }
        std::cout << m_params.back() << '\n';
    }

    /**
     * @brief Calculate the track parameters initial guess
     * 
     * @param data experimental data vector
     * 
     * @return Boolean whether the calculation is successfull or not
     */
    virtual bool initiate(const vector_type& data) = 0;

protected:

    std::vector<double> m_params;

};

#endif // CDWPTTCHI2RECTOOL_INITIALIZER_INITIALIZER_HPP_