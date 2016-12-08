#ifndef OCCMAP_CONSTANTS_H_K73HUCVF
#define OCCMAP_CONSTANTS_H_K73HUCVF

namespace ca { namespace scrollgrid {

template<class T>
struct OccMapConstants;

template<>
struct OccMapConstants<float> {
  static constexpr float UNKNOWN = 0.5f;
  static constexpr float OCCUPIED = 0.9f;
  static constexpr float FREE = 0.0f;
  static constexpr float UPDATE_POS = 0.08f;
  static constexpr float UPDATE_NEG = 0.01f;

  static
  float update_pos(float x) {
    return (std::min(x + UPDATE_POS, OCCUPIED));
  }

  static
  float update_neg(float x) {
    return (std::max(x - UPDATE_NEG, FREE));
  }
};

template<>
struct OccMapConstants<uint8_t> {
  static constexpr uint8_t UNKNOWN = 128;
  static constexpr uint8_t OCCUPIED = 250;
  static constexpr uint8_t FREE = 0;
  static constexpr int32_t UPDATE_POS = 20;
  static constexpr int32_t UPDATE_NEG = 2;

  static
  uint8_t update_pos(uint8_t x) {
    // TODO avoid casts
    int32_t new_val = static_cast<int32_t>(x) + UPDATE_POS;
    new_val = std::min(static_cast<int32_t>(OCCUPIED), new_val);
    return static_cast<uint8_t>(new_val);
  }

  static
  uint8_t update_neg(uint8_t x) {
    // TODO avoid casts
    int32_t new_val = static_cast<int32_t>(x) - UPDATE_NEG;
    new_val = std::max(static_cast<int32_t>(FREE), new_val);
    return static_cast<uint8_t>(new_val);
  }

};

} }

#endif /* end of include guard: OCCMAP_CONSTANTS_H_K73HUCVF */
