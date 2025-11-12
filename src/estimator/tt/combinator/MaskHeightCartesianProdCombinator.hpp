#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_TT_COMBINATOR_MASKHEIGHTCARTESIANPRODCOMBINATOR_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_TT_COMBINATOR_MASKHEIGHTCARTESIANPRODCOMBINATOR_HPP_

#include <unordered_map>

#include "estimator/tt/combinator/Combinator.hpp"

class MaskHeightCartesianProdCombinator : public Combinator {

public:

    MaskHeightCartesianProdCombinator(const std::string& name, std::size_t max_nb_heights);
    ~MaskHeightCartesianProdCombinator() override = default;

    void configure(const SniperJSON& config) override;

    bool combine(const std::vector<vec3>& hits) override;

private:

    bool getHeightMap(const std::vector<vec3>& hits);
    bool cartesianProduct(std::vector<vec3>& curr, std::unordered_map<double, std::vector<vec3>>::const_iterator it, std::vector<bool>::const_iterator mask_it);

    std::size_t m_max_nb_heights;
    std::unordered_map<double, std::vector<vec3>> m_height_map;
    std::vector<bool> m_mask;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_TT_COMBINATOR_MASKHEIGHTCARTESIANPRODCOMBINATOR_HPP_