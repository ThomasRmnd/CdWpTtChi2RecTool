#ifndef CDWPTTCHI2RECTOOL_UTILS_RECPMTPROP_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_RECPMTPROP_HPP_

#include <algorithm>
#include <unordered_map>
#include <utility>
#include <vector>

#include "SniperKernel/SniperLog.h"

#include "Geometry/IPMTParamSvc.h"
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

inline std::unordered_map<RecPmtType, double> g_pmt_type_to_res = {
    {RecPmtType::PMT_20INCH_NNVT, 8.0},
    {RecPmtType::PMT_20INCH_HIGHQENNVT, 8.0},
    {RecPmtType::PMT_20INCH_HAMAMATSU, 6.0},
    {RecPmtType::PMT_3INCH, 2.3},
    {RecPmtType::PMT_WP, 13.0},
    {RecPmtType::PMT_TT, 2.0}
};

struct RecPmtProp { // sizeof = 64

    vec3 pos;
    double inv_res;
    double totq;
    double q; // fht q
    double fht;
    unsigned int id;
    RecPmtType type;
    bool used = false;

    // double inv_res() const { return g_pmt_type_to_res[type]; }; // <-- this would allow to remove inv_res member, reducing sizeof(RecPmtProp) to 56

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

class PmtTypeChecker {

public:

    PmtTypeChecker(const RecPmtType& pmt_type) : c_pmt_type(pmt_type) {}

    virtual ~PmtTypeChecker() = default;

    bool checkPmtType(const RecPmtProp& pmt) {
        return hasPmtType(pmt, c_pmt_type);
    }

protected:

    const RecPmtType c_pmt_type;

};

typedef std::vector<RecPmtProp> RecPmtTable;

// TODO: Maybe use PMTParamSvc instead of CdID::pmtType
class TableConverter {

public:

    static void convert(const PmtTable* src, RecPmtTable& dst) {
        getTotUsedAndPE(src);
        dst.resize(m_tot_used, RecPmtProp{});
        RecPmtTable::iterator it_dst = dst.begin();
        for (PmtTable::const_iterator it_src = src->begin(); it_src != src->end(); ++it_src) {
            if (!it_src->used) continue;
            convertToRecPmtProp(*it_src, *it_dst++);
        }
    }

    static void setPMTSvc(IPMTParamSvc* pmt_svc) { m_pmt_svc = pmt_svc; }
    static double getTotPE() { return m_tot_pe; }

private:

    inline static std::size_t m_tot_used = 0ul;
    inline static double m_tot_pe = 0.0;
    inline static IPMTParamSvc* m_pmt_svc = nullptr;

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
        // Identifier id = Identifier{src.pmtid};
        Identifier id = Identifier{static_cast<unsigned int>(src.value)};
        if (src.loc == 3) {
            type = RecPmtType::PMT_TT;
        }
        else if (src.loc == 2) {
            type = RecPmtType::PMT_WP;
        }
        else if (src.loc ==1) {
            type = resolveCdPmtType(id);
        }
        else {
            throw std::runtime_error("Unknown PMT location: " + std::to_string(src.loc));
        }
        return type;
    }

    static RecPmtType resolveCdPmtType(const Identifier& id) {
        RecPmtType type = RecPmtType::PMT_CD;
        if (CdID::is20inch(id)) {
            if (m_pmt_svc) {
                type = RecPmtType::PMT_20INCH;
                int id_module = CdID::module(id);
                if (m_pmt_svc->isHamamatsu(id_module)) {
                    type = RecPmtType::PMT_20INCH_HAMAMATSU;
                }
                else if (m_pmt_svc->isNormalNNVT(id_module)) {
                    type = RecPmtType::PMT_20INCH_NNVT;
                }
                else if (m_pmt_svc->isHighQENNVT(id_module)) {
                    type = RecPmtType::PMT_20INCH_HIGHQENNVT;
                }
                else {
                    throw std::runtime_error("Unknown CD PMT type from PMTParamSvc: " + std::to_string(static_cast<unsigned int>(id)));
                }
            }
            else {
                type = RecPmtType::PMT_20INCH;
                if (CdID::pmtType(id) == 1) {
                    type = RecPmtType::PMT_20INCH_HAMAMATSU;
                }
                else if (CdID::pmtType(id) == 2) {
                    type = RecPmtType::PMT_20INCH_NNVT;
                }
                else {
                    throw std::runtime_error("Unknown CD 20inch PMT type: " + std::to_string(static_cast<unsigned int>(id)));
                }
            }
        }
        else if (CdID::is3inch(id)) {
            type = RecPmtType::PMT_3INCH;
        }
        else {
            throw std::runtime_error("Unknown CD PMT type: " + std::to_string(static_cast<unsigned int>(id)));
        }
        return type;
    }

    static void convertToRecPmtProp(const PmtProp& src, RecPmtProp& dst) {
        dst.pos = src.pos;
        dst.type = getRecPmtType(src);
        dst.inv_res = 1.0 / g_pmt_type_to_res[dst.type];
        dst.totq = src.q;
        dst.q = src.hitq[0];
        dst.fht = src.fht;
        // dst.id = src.pmtid;
        dst.id = static_cast<unsigned int>(src.value);
        dst.used = true;
    }

};

#endif // CDWPTTCHI2RECTOOL_UTILS_RECPMTPROP_HPP_