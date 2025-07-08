#ifndef CDWPTTCHI2RECTOOL_ESTIMATOR_TT_CONVERTER_FUZENEIGHBORCONVERTER_HPP_
#define CDWPTTCHI2RECTOOL_ESTIMATOR_TT_CONVERTER_FUZENEIGHBORCONVERTER_HPP_

#include "estimator/tt/converter/Converter.hpp"

#include <unordered_map>
#include <vector>

class FuzeNeighborConverter : public Converter {

    struct WallHits {

        std::vector<double> x;
        std::vector<double> y;
        double z_x;
        double z_y;

    };

public:

    using Converter::Converter;

    ~FuzeNeighborConverter() override = default;

    bool operator()(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable) override;

protected:

    bool getWallHits(RecPmtTable::const_iterator ftable, RecPmtTable::const_iterator ltable);
    bool cleanWallHits();
    std::vector<double> cleanAxisHits(std::vector<double>& axis_hits);
    bool convertToHits();

    std::unordered_map<int, WallHits> m_wall_hits;
    std::vector<WallHits> m_cleaned_wall_hits;

};

#endif // CDWPTTCHI2RECTOOL_ESTIMATOR_TT_CONVERTER_FUZENEIGHBORCONVERTER_HPP_