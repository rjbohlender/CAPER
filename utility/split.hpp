// -*- mode: C++ -*-
/**
 * @file split.hpp
 * @brief Implements a simple (read amateur) Splitter.
 *
 */
#ifndef PGEN_SPLIT_HPP
#define PGEN_SPLIT_HPP 1

#include <algorithm>
#include <string>
#include <string_view>
#include <vector>

namespace RJBUtil {
/**
 * @class Splitter split.hpp "src/util/split.hpp"
 * @brief A utility for quickly splitting strings into readable substring
 * segments.
 *
 * The implementation uses string_view to reference substrings without
 * copying, providing a light-weight view over the original data.
 */
template <typename String_t, typename Integral_t = size_t> class Splitter {
public:
  /** Aliases */
  using string_type = String_t;
  using int_type = Integral_t;
  using view_type = std::basic_string_view<typename string_type::value_type>;
  using container_type = std::vector<view_type>;
  using iterator = typename container_type::iterator;
  using const_iterator = typename container_type::const_iterator;

  /** Constructors */
  Splitter() = default;

  // View-based constructor (non-owning)
  Splitter(view_type str, view_type delim, int_type n = 0,
           bool split_adjacent = false)
      : n_(n), data_(str), delim_(delim), split_adjacent_(split_adjacent) {
    split_tokens();
  }

  // Lvalue string constructor (non-owning)
  Splitter(const string_type &str, view_type delim, int_type n = 0,
           bool split_adjacent = false)
      : n_(n), data_(str), delim_(delim), split_adjacent_(split_adjacent) {
    split_tokens();
  }

  // Rvalue string constructor (owning)
  Splitter(string_type &&str, view_type delim, int_type n = 0,
           bool split_adjacent = false)
      : n_(n), storage_(std::move(str)), data_(storage_), delim_(delim),
        split_adjacent_(split_adjacent) {
    split_tokens();
  }

  Splitter(const Splitter &rhs)
      : n_(rhs.n_), storage_(rhs.storage_), delim_(rhs.delim_),
        split_adjacent_(rhs.split_adjacent_) {
    data_ = storage_.empty() ? rhs.data_ : view_type(storage_);
    split_tokens();
  }

  Splitter(Splitter &&) noexcept = default;

  Splitter &operator=(const Splitter &rhs) {
    if (this != &rhs) {
      n_ = rhs.n_;
      storage_ = rhs.storage_;
      delim_ = rhs.delim_;
      split_adjacent_ = rhs.split_adjacent_;
      data_ = storage_.empty() ? rhs.data_ : view_type(storage_);
      split_tokens();
    }
    return *this;
  }

  Splitter &operator=(Splitter &&rhs) noexcept = default;

  iterator begin() { return tokens_.begin(); }

  const_iterator cbegin() const { return tokens_.cbegin(); }

  iterator end() { return tokens_.end(); }

  const_iterator cend() const { return tokens_.cend(); }

  size_t size() const { return tokens_.size(); }

  bool empty() const { return tokens_.empty(); }

  iterator erase(iterator pos) { return tokens_.erase(pos); }

  iterator erase(const_iterator pos) { return tokens_.erase(pos); }

  view_type &front() { return tokens_.front(); }

  view_type &back() { return tokens_.back(); }

  view_type &operator[](int_type i) { return tokens_[i]; }

  view_type &at(int_type i) { return tokens_.at(i); }

private:
  /**
   * @brief A split member function to simplify multiple constructors should be
   * need them.
   */
  void split() {
#if 0
	const_iter cit = data_.cbegin();
	size_type cur = 0;
	size_type next = 0;
	size_type npos = std::string::npos;
	while ((cit + cur) != data_.cend()) {
	  next = data_.find_first_of(delim_, cur);
	  // Skip adjacent.
	  if (next - cur == 0) {
		cur++;
		next = cur;
		continue;
	  }
	  if (next == npos) {
		if (data_.cend() - (cit + cur) > 0) {
		  tokens_.emplace_back(std::string(cit + cur, data_.cend()));
		}
		break;
	  }

	  tokens_.emplace_back(std::string(cit + cur, cit + next));

	  // Prepare for next loop
	  next++;
	  cur = next;

	}
#else
    // Faster version described at:
    // https://www.bfilipek.com/2018/07/string-view-perf-followup.html
    for (auto first = data_.data(), second = data_.data(),
              last = first + data_.size();
         second != last && first != last; first = second + 1) {
      second = std::find_first_of(first, last, std::cbegin(delim_),
                                  std::cend(delim_));

      if (split_adjacent_) {
        tokens_.emplace_back(first, second - first);
      } else if (first != second)
        tokens_.emplace_back(first, second - first);
    }
#endif
  }
  void split_up_to_n() {
    int cur = 0;
    tokens_.reserve(n_ + 1);
    // Faster version described at:
    // https://www.bfilipek.com/2018/07/string-view-perf-followup.html
    for (auto first = data_.data(), second = data_.data(),
              last = first + data_.size();
         second != last && first != last; first = second + 1) {
      second = std::find_first_of(first, last, std::cbegin(delim_),
                                  std::cend(delim_));

      if (split_adjacent_) {
        tokens_.emplace_back(first, second - first);
      } else if (first != second)
        tokens_.emplace_back(first, second - first);
      cur++;
      if (cur >= n_) {
        if (second != last) {
          tokens_.emplace_back(second + 1, last - (second + 1));
        }
        break;
      }
    }
  }

  /** Private members */
  void split_tokens() {
    tokens_.clear();
    if (n_ > 0)
      split_up_to_n();
    else
      split();
  }

  size_t n_{};
  string_type storage_{};
  view_type data_{};
  view_type delim_{};
  container_type tokens_{};
  bool split_adjacent_{};
};
} // namespace RJBUtil
#endif // PGEN_SPLIT_HPP