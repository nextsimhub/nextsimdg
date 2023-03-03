/*!
 * @file TextTag.hpp
 *
 * @date 27 Feb 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef TEXTTAG_HPP
#define TEXTTAG_HPP

#include <string>

struct TextTag {
    constexpr TextTag(const char* textIn)
        : text(textIn)
    {
    }
    template <std::size_t N>
    constexpr TextTag(const char (&a)[N])
        : text(a)
    {
    }
    operator std::string() const { return std::string(text); };
    const char* text;
};

#endif /* TEXTTAG_HPP */
