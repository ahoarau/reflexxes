// Copyright 2026 Intrinsic Innovation LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// This file holds a set of functional programming utilities to work with
// collections and arrays.

#ifndef INTRINSIC_ICON_REFLEXXES_INTERNAL_FUNCTIONAL_H_
#define INTRINSIC_ICON_REFLEXXES_INTERNAL_FUNCTIONAL_H_

#include <algorithm>
#include <array>

namespace intrinsic {
namespace reflexxes {
namespace internal {

namespace zip_details {
// Utility class used to implement Zip below.
template <typename C1Iterator, typename C2Iterator>
class ZipHolder {
 public:
  using C1Reference = decltype(*std::declval<C1Iterator>());
  using C2Reference = decltype(*std::declval<C2Iterator>());

  using ResultPair = std::pair<C1Reference, C2Reference>;

  class ZipIterator {
   public:
    ZipIterator(C1Iterator first, C2Iterator second)
        : first_(first), second_(second) {}

    ResultPair operator*() { return ResultPair(*first_, *second_); }

    void operator++() {
      first_++;
      second_++;
    }

    bool operator!=(const ZipIterator& other) const {
      // Note the '&&' as we want it to be false if either are equal.  This
      // allows iteration checking for "end" to finish as soon as one iterator
      // has reached its final point.
      return first_ != other.first_ && second_ != other.second_;
    }

   private:
    C1Iterator first_;
    C2Iterator second_;
  };

  ZipHolder(C1Iterator c1_begin, C1Iterator c1_end, C2Iterator c2_begin,
            C2Iterator c2_end)
      : begin_(c1_begin, c2_begin), end_(c1_end, c2_end) {}

  auto begin() { return begin_; }
  auto end() { return end_; }

 private:
  ZipIterator begin_;
  ZipIterator end_;
};
template <typename C1Iterator, typename C2Iterator>
ZipHolder(C1Iterator, C1Iterator, C2Iterator, C2Iterator)
    -> ZipHolder<C1Iterator, C2Iterator>;

}  // namespace zip_details

// Utility function for zipping together two containers for an iteration.
// Example:
// std::vector<Foo> foo;
// std::vector<Bar> bar;
// for (auto [x, y] = Zip(foo, bar)) {
//  // do something
// }
//
// where foo and bar are references to the corresponding values at each
// iteration of foo and bar.
//
// The return value of this call must not outlive the result of calling begin()
// or end() on the constituent containers.  So for the above example, one cannot
// call "push_back" on foo or bar in the for loop, or the result will be
// undefined.
template <typename C1, typename C2>
auto Zip(C1&& c1, C2&& c2) {
  return zip_details::ZipHolder(c1.begin(), c1.end(), c2.begin(), c2.end());
}

// Returns the first ToSize members of an array as an array<T, ToSize>.
template <size_t ToSize, typename T, size_t FromSize>
std::array<T, ToSize> Head(const std::array<T, FromSize>& from) {
  static_assert(ToSize <= FromSize);
  std::array<T, ToSize> to;
  std::copy_n(from.begin(), ToSize, to.begin());
  return to;
}

}  // namespace internal
}  // namespace reflexxes
}  // namespace intrinsic

#endif  // INTRINSIC_ICON_REFLEXXES_INTERNAL_FUNCTIONAL_H_
