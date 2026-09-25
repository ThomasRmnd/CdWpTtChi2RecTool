#ifndef CDWPTTCHI2RECTOOL_UTILS_DETECTORTYPE_HPP_
#define CDWPTTCHI2RECTOOL_UTILS_DETECTORTYPE_HPP_

enum class DetectorType {
    NONE = 0,
    CD = 1 << 0,
    WP = 1 << 1,
    TT = 1 << 2
};

inline DetectorType operator|(const DetectorType& a, const DetectorType& b) {
    return static_cast<DetectorType>(static_cast<int>(a) | static_cast<int>(b));
}

inline DetectorType operator&(const DetectorType& a, const DetectorType& b) {
    return static_cast<DetectorType>(static_cast<int>(a) & static_cast<int>(b));
}

inline DetectorType operator^(const DetectorType& a, const DetectorType& b) {
    return static_cast<DetectorType>(static_cast<int>(a) ^ static_cast<int>(b));
}

inline DetectorType operator~(const DetectorType& a) {
    return static_cast<DetectorType>(~static_cast<int>(a));
}

inline DetectorType& operator|=(DetectorType& a, const DetectorType& b) {
    return a = a | b;
}

inline DetectorType& operator&=(DetectorType& a, const DetectorType& b) {
    return a = a & b;
}

inline DetectorType& operator^=(DetectorType& a, const DetectorType& b) {
    return a = a ^ b;
}

#endif // CDWPTTCHI2RECTOOL_UTILS_DETECTORTYPE_HPP_