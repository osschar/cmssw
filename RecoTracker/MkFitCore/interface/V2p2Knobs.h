#ifndef RecoTracker_MkFitCore_interface_V2p2Knobs_h
#define RecoTracker_MkFitCore_interface_V2p2Knobs_h

//==============================================================================
// The table of the V2p2 knobs of src/V2p2Config.h, in a CMSSW build with
// MKFIT_V2P2_MUTABLE only.
//
// The knobs are constexpr in the normal CMSSW build. MKFIT_V2P2_MUTABLE makes
// them mutable globals, as in the standalone build, and fills this table: one
// entry per knob, named "<group>.<knob>" (e.g. "Window.dphi_trk_fac"), with
// its type, its address and its default. The entries are made by the knob
// declarations themselves (MKFIT_V2P2_KNOB), so the table cannot miss a knob.
//
// A debugging build: the knobs are process-wide and read without
// synchronisation, so change them only while no event is processed (e.g. in an
// EDLooper's endOfLoop, or before the first event). Never the production build.
//==============================================================================

#if defined(MKFIT_V2P2_MUTABLE) && !defined(MKFIT_STANDALONE)

#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace mkfit::Config::V2p2 {

  struct Knob {
    // The variant index is the type: 0 bool, 1 int, 2 float.
    using Address = std::variant<bool*, int*, float*>;
    using Value = std::variant<bool, int, float>;

    std::string name;  // "<group>.<knob>"
    Address address;
    Value default_value;

    const char* type_name() const;
    Value get() const;
    // Sets the knob; throws std::invalid_argument if the value's type is not the knob's.
    void set(const Value& v) const;
  };

  // All knobs, in the order of their declaration in src/V2p2Config.h.
  const std::vector<Knob>& knobs();
  // The knob of that name, or nullptr.
  const Knob* find_knob(std::string_view name);
  // Every knob back to its default.
  void reset_knobs();

  const char* knob_type_name(const Knob::Value& v);
  std::string knob_value_string(const Knob::Value& v);

  // Made by MKFIT_V2P2_KNOB in V2p2Config.cc: adds the knob to the table, with
  // its value at static initialisation as the default.
  struct KnobRegistration {
    KnobRegistration(const char* name, bool* address);
    KnobRegistration(const char* name, int* address);
    KnobRegistration(const char* name, float* address);
  };

}  // namespace mkfit::Config::V2p2

#endif

#endif
