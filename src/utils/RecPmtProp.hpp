#ifndef CDWPTTCHI2RECTOOL_UTILS_RECPMTPROP_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_RECPMTPROP_HPP_

#include <algorithm>
#include <unordered_map>
#include <utility>
#include <vector>

#include "SniperKernel/SniperLog.h"

#include "Identifier/CdID.h"
#include "Identifier/WpID.h"
#include "RecTools/PmtProp.h"

#include "utils/vec3.hpp"

enum class RecPmtType : unsigned char {

    PMT_UNKNOWN           = 0x00,
    PMT_20INCH_HAMAMATSU  = 0x01,
    PMT_20INCH_NNVT       = 0x02,
    PMT_20INCH_HIGHQENNVT = 0x04,
    PMT_20INCH            = PMT_20INCH_HAMAMATSU | PMT_20INCH_NNVT | PMT_20INCH_HIGHQENNVT,
    PMT_3INCH             = 0x08,
    PMT_CD                = PMT_20INCH | PMT_3INCH,
    PMT_WP                = 0x10,
    PMT_TT                = 0x20,

};

inline constexpr RecPmtType operator|(RecPmtType lhs, RecPmtType rhs) noexcept {
	return static_cast<RecPmtType>(static_cast<unsigned char>(lhs) | static_cast<unsigned char>(rhs));
}

inline constexpr RecPmtType operator&(RecPmtType lhs, RecPmtType rhs) noexcept {
	return static_cast<RecPmtType>(static_cast<unsigned char>(lhs) & static_cast<unsigned char>(rhs));
}

inline constexpr RecPmtType operator^(RecPmtType lhs, RecPmtType rhs) noexcept {
	return static_cast<RecPmtType>(static_cast<unsigned char>(lhs) ^ static_cast<unsigned char>(rhs));
}

inline constexpr RecPmtType operator~(RecPmtType rhs) noexcept {
	return static_cast<RecPmtType>(~static_cast<unsigned char>(rhs));
}

inline constexpr RecPmtType& operator|=(RecPmtType& lhs, RecPmtType rhs) noexcept {
	return (lhs = lhs | rhs);
}

inline constexpr RecPmtType& operator&=(RecPmtType& lhs, RecPmtType rhs) noexcept {
	return (lhs = lhs & rhs);
}

inline constexpr RecPmtType& operator^=(RecPmtType& lhs, RecPmtType rhs) noexcept {
	return (lhs = lhs ^ rhs);
}

namespace std {

template <>
struct hash<RecPmtType> {
    std::size_t operator()(const RecPmtType& t) const noexcept {
        return static_cast<std::size_t>(t);
    }
};

} // namespace std

inline std::unordered_map<RecPmtType, double> pmt_type_to_res = {
    {RecPmtType::PMT_20INCH_NNVT, 7.0}, // 5.5},
    {RecPmtType::PMT_20INCH_HIGHQENNVT, 7.0}, // 5.5},
    {RecPmtType::PMT_20INCH_HAMAMATSU, 5.0}, // 2.0},
    {RecPmtType::PMT_3INCH, 2.3},
    {RecPmtType::PMT_WP, 13.0}
};

struct RecPmtProp { // sizeof = 64 (PmtType : unsigned char)

    vec3 pos;
    double inv_res;
    double totq;
    double q; // fht q
    double fht;
    unsigned int id;
    RecPmtType type;
    bool used;

    // double inv_res() const { return pmt_type_to_invres[type]; }; // <-- this would allow to remove inv_res member and thus reduce sizeof(RecPmtProp) to 56

};

inline constexpr bool operator==(const RecPmtProp& lhs, const RecPmtProp& rhs) noexcept {
    return lhs.id == rhs.id;
}

inline constexpr bool operator!=(const RecPmtProp& lhs, const RecPmtProp& rhs) noexcept {
    return lhs.id != rhs.id;
}

inline constexpr bool hasPmtType(const RecPmtProp& pmt, const RecPmtType& type) noexcept {
    return ( (pmt.type & type ) == pmt.type );
}

typedef std::vector<RecPmtProp> RecPmtTable;

class TableConverter {

public:

    static void convert(const PmtTable* src, RecPmtTable& dst) {
        getTotUsedAndPE(src);
        dst.resize(m_tot_used);
        RecPmtTable::iterator it_dst = dst.begin();
        for (PmtTable::const_iterator it_src = src->begin(); it_src != src->end(); ++it_src) {
            if (!it_src->used) continue;
            convertToRecPmtProp(*it_src, *it_dst++);
        }
    }

    static double getTotPE() { return m_tot_pe; }

private:

    inline static std::size_t m_tot_used = 0ul;
    inline static double m_tot_pe = 0.0;

    static void getTotUsedAndPE(const PmtTable* src) {
        m_tot_used = 0ul;
        m_tot_pe = 0.0;
        std::for_each(src->begin(), src->end(), [](const PmtProp& pmt) {
            if (pmt.used) {
                ++m_tot_used;
                m_tot_pe += pmt.q;
            }
        });
    }

    static RecPmtType getRecPmtType(const PmtProp& src) {
        RecPmtType type = RecPmtType::PMT_UNKNOWN;
        Identifier id = Identifier{src.pmtid};
        if (src.loc == 3) {
            type = RecPmtType::PMT_TT;
        }
        else if (src.loc == 2 && WpID::is20inch(id)) {
            type = RecPmtType::PMT_WP;
        }
        else if (src.loc == 1 && CdID::is20inch(id) && CdID::pmtType(id) == 1) {
            type = RecPmtType::PMT_20INCH_HAMAMATSU;
        }
        else if (src.loc == 1 && CdID::is20inch(id) && CdID::pmtType(id) == 2) {
            type = RecPmtType::PMT_20INCH_NNVT;
        }
        else if (src.loc == 1 && CdID::is3inch(id) && CdID::pmtType(id) == 3) {
            type = RecPmtType::PMT_3INCH;
        }
        return type;
    }

    static void convertToRecPmtProp(const PmtProp& src, RecPmtProp& dst) {
        dst.pos = src.pos;
        dst.type = getRecPmtType(src);
        dst.inv_res = 1.0 / pmt_type_to_res[dst.type];
        dst.totq = src.q;
        dst.q = src.hitq[0];
        dst.fht = src.fht;
        dst.id = src.pmtid;
        dst.used = true;
    }

};

#endif // CDWPTTCHI2RECTOOL_UTILS_RECPMTPROP_HPP_