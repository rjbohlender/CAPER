// -*- mode: C++ -*-
/**
 * @file split.hpp
 * @brief Implements a simple (read amateur) Splitter.
 *
 */
#ifndef PGEN_SPLIT_HPP
#define PGEN_SPLIT_HPP 1

#include <string>
#include <vector>

namespace RJBUtil {
	/**
	 * @class Splitter split.hpp "src/util/split.hpp"
	 * @brief A utility for quickly splitting strings into readable substring segments.
	 *
	 * @remark  Removed string_view, as string_views can't be built from iterators.
	 */
	template<typename _String_t>
	class Splitter
	{
	public:
			/** Aliases */
			using string_type = _String_t;
			using const_iter = typename string_type::const_iterator;
			using size_type = typename _String_t::size_type;

			/** Constructors */
			Splitter(const string_type &str, const string_type &delim)
							: __data(str), delim(delim) {
				split();
			}

			Splitter(const string_type &&str, const string_type &&delim)
							: __data(str), delim(delim) {
				split();
			}

			Splitter(const string_type &str, const string_type &&delim)
							: __data(str), delim(delim) {
				split();
			}

			Splitter(const string_type &&str, const string_type &delim)
							: __data(str), delim(delim) {
				split();
			}

			auto begin() {
				return __tokens.begin();
			}

			auto cbegin() {
				return __tokens.cbegin();
			}

			auto end() {
				return __tokens.end();
			}

			auto cend() {
				return __tokens.cend();
			}

			auto size() {
				return __tokens.size();
			}

			template<typename _Integer>
			auto operator[](_Integer i) {
				return __tokens[i];
			}

			template<typename _Integer>
			auto at(_Integer i) {
				return __tokens.at(i);
			}

	private:
			/**
			 * @brief A split member function to simplify multiple constructors should be need them.
			 */
			void split() {
				const_iter cit = __data.cbegin();
				size_type cur = 0;
				size_type next = 0;
				size_type npos = std::string::npos;
				while ((cit + cur) != __data.cend()) {
					next = __data.find_first_of(delim, cur);
					// Skip adjacent.
					if (next - cur == 0) {
						cur++;
						next = cur;
						continue;
					}
					if (next == npos) {
					  	if(__data.cend() - (cit + cur) > 0) {
							__tokens.emplace_back(std::string(cit + cur, __data.cend()));
					  	}
					  	break;
					}

					__tokens.emplace_back(std::string(cit + cur, cit + next));

					// Prepare for next loop
					next++;
					cur = next;

				}
			}

			/** Private members */
			const string_type __data;
			const string_type delim;
			std::vector<std::string> __tokens;
	};
}
#endif // PGEN_SPLIT_HPP