#ifndef RecoTracker_MkFitCore_src_MkRZLimits_h
#define RecoTracker_MkFitCore_src_MkRZLimits_h

namespace mkfit {
 
  class LayerInfo;

  struct MkRZLimits {
    const LayerInfo *m_layer_info_1 = nullptr;
    const LayerInfo *m_layer_info_2 = nullptr;
    float m_rin = 0, m_rout = 0, m_zmin = 0, m_zmax = 0;
    bool m_is_barrel = true;
    bool m_is_outward = true;
    bool m_is_double = false;
    bool m_is_initialized = false;

    MkRZLimits() = default;
    MkRZLimits(const LayerInfo &li, bool is_outward);
    MkRZLimits(const LayerInfo &li1, const LayerInfo &li2, bool is_outward);

    void setup(const LayerInfo &li, bool is_outward);
    void setup(const LayerInfo &li1, const LayerInfo &li2, bool is_outward);

    void reset();

    const LayerInfo &layer_info_1() const { return *m_layer_info_1; }
    const LayerInfo &layer_info_2() const { return *m_layer_info_2; }

    bool is_outward() const { return m_is_outward; }
    bool is_inward() const { return !m_is_outward; }

    // Extremely crude track / layer intersection pre-check.
    // Returns true if there is a chance of intersection.
    // Can also add r, pr as args -- beware of loopers and near sagita cases.
    bool rz_quadrant_check(float z, float pz) {
      if (is_inward()) {
        pz = -pz;
      }
      if (pz < 0.0f) {
        if (z < m_zmin)
          return false;
      } else {
        if (z > m_zmax)
          return false;
      }
      return true;
    }
  };

}

#endif
