#pragma once
#include "element/element.hpp"

namespace abq {
std::unique_ptr<ElementGroup> make_solid_group();
std::unique_ptr<ElementGroup> make_shell_group();
std::unique_ptr<ElementGroup> make_beam_group();
std::unique_ptr<ElementGroup> make_connector_group();
std::unique_ptr<ElementGroup> make_truss_group();
}  // namespace abq
