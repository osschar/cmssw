// Standalone definitions of the knobs declared in V2p2Config.h. In the CMSSW
// build they are constexpr and this file defines nothing; in the CMSSW build
// with MKFIT_V2P2_MUTABLE it defines them and the table of interface/V2p2Knobs.h.
//
// Must be the first include, so the header is expanded with definitions.
#define MKFIT_V2P2_CONFIG_DEFINE
#include "RecoTracker/MkFitCore/src/V2p2Config.h"

#if defined(MKFIT_STANDALONE)

namespace mkfit::Config::V2p2 {

  void set_extra_dq(float f) {
    Window::dq_trk_fac = f;
    Window::dq_hit_fac = f * 1.2f;  // the old EXTRA_DQ ratio
  }

}  // namespace mkfit::Config::V2p2

#endif

#if defined(MKFIT_V2P2_MUTABLE) && !defined(MKFIT_STANDALONE)

#include <charconv>
#include <stdexcept>
#include <type_traits>

namespace mkfit::Config::V2p2 {

  namespace {
    std::vector<Knob> &knob_table() {
      static std::vector<Knob> table;
      return table;
    }

    template <typename T>
    void add_knob(const char *name, T *address) {
      knob_table().push_back(Knob{name, address, *address});
    }
  }  // namespace

  KnobRegistration::KnobRegistration(const char *name, bool *address) { add_knob(name, address); }
  KnobRegistration::KnobRegistration(const char *name, int *address) { add_knob(name, address); }
  KnobRegistration::KnobRegistration(const char *name, float *address) { add_knob(name, address); }

  const char *knob_type_name(const Knob::Value &v) {
    static constexpr const char *names[] = {"bool", "int", "float"};
    return names[v.index()];
  }

  std::string knob_value_string(const Knob::Value &v) {
    return std::visit(
        [](auto x) -> std::string {
          if constexpr (std::is_same_v<decltype(x), bool>) {
            return x ? "true" : "false";
          } else {
            char buf[32];
            auto r = std::to_chars(buf, buf + sizeof(buf), x);  // shortest form that reads back exactly
            return std::string(buf, r.ptr);
          }
        },
        v);
  }

  const char *Knob::type_name() const { return knob_type_name(default_value); }

  Knob::Value Knob::get() const {
    return std::visit([](auto *p) -> Value { return *p; }, address);
  }

  void Knob::set(const Value &v) const {
    if (v.index() != address.index())
      throw std::invalid_argument("V2p2 knob " + name + " is a " + type_name() + ", not a " + knob_type_name(v));
    std::visit([&v](auto *p) { *p = std::get<std::remove_pointer_t<decltype(p)>>(v); }, address);
  }

  const std::vector<Knob> &knobs() { return knob_table(); }

  const Knob *find_knob(std::string_view name) {
    for (const Knob &k : knob_table()) {
      if (k.name == name)
        return &k;
    }
    return nullptr;
  }

  void reset_knobs() {
    for (const Knob &k : knob_table())
      k.set(k.default_value);
  }

}  // namespace mkfit::Config::V2p2

#endif
