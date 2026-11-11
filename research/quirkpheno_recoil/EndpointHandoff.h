#pragma once
#include <cstdint>
#include <stdexcept>

// One-use ownership transfer. Independent of Geant4 for native testing.
struct EndpointHandoff {
  bool pending{false};
  int parent{0}, daughter{0};
  std::uint64_t own{0}, partner{0};
  void begin(int parentID, int daughterPDG, std::uint64_t a, std::uint64_t b) {
    if (pending || parentID <= 0 || daughterPDG == 0)
      throw std::logic_error("invalid or duplicate endpoint handoff");
    pending=true; parent=parentID; daughter=daughterPDG; own=a; partner=b;
  }
  bool isParent(int trackID) const { return pending && parent==trackID; }
  void consume(int parentID, int daughterPDG, std::uint64_t a, std::uint64_t b) {
    if (!pending || parent!=parentID || daughter!=daughterPDG || own!=a || partner!=b)
      throw std::logic_error("endpoint handoff identity or string-state mismatch");
    pending=false;
  }
};
