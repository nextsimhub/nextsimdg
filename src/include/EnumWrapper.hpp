/*!
 * @file EnumWrapper.hpp
 *
 * @date Oct 15, 2021
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef ENUMWRAPPER_HPP
#define ENUMWRAPPER_HPP

#include <string>
#include <map>
#include <boost/program_options.hpp>

// A helper class to make configuring enums easier with boost::program_options.

namespace EnumWrap {
template <typename E> class EnumWrapper;
template <typename E> std::istream& operator>>(std::istream&, EnumWrapper<E>&);

/*!
 * A helper class to make configuring enums easier with boost::program_options.
 *
 * This class provides a templated implementation of the instream
 * operator necessary for enums to be read by
 * boost::program_options. Also provided is a macro which will format
 * the map that is required to go from the string representation in
 * the config file to the values of the enum.
 *
 * @section <usage> (Usage)
 * For an enum @c Enum with values @c x and @c y
 * @code
 * enum class Enum {x, y};
 * @endcode
 *
 * the mapping can be defined using the provided @c MAP_ENUM macro
 * where the first argument is the name of the enum type and the
 * following arguments are the pairs of string label and enum values
 *
 * @code
 * MAP_ENUM(Enum,
 *          {"x", Enum::x},
 *          {"wye", Enum::y});
 * @endcode
 *
 * For long @c enum type names, a typedef or even macro variable
 * substitution can be used without affecting the assignment of the
 * map. The map can also be set by manually writing the code to assign
 * to the static std::map element named @c map.
 *
 * When defining the @c boost::program_options options description,
 * the new option is defined with type @c EnumWrapper<Enum>. After
 * parsing the value can be retrieved using @c
 * .as<EnumWrapper<Enum>>(). The value can then be stored in a
 * variable of type @c Enum, as the conversion is handled
 * automatically.
 */
template<typename E>
class EnumWrapper {
public:

  typedef std::map<std::string, E> MapType;

  /*!
   * Set and return the value of the wrapped enum, based on the passed
   * string and the set map.
   *
   * @param key a key corresponding to one in the templated map which
   * will set the value of the wrapped enum.
   */
  E operator()(const std::string& key)
  {
    value = map.at(key);
    return value;
  }
  //! Cast the wrapped enum back to a plain one.
  operator E() const
  {
    return value;
  }

  /*!
   * Set the mapping between strings and enum types
   *
   * @param inMap the mapping to be used.
   */  
  static void setMap(const MapType& inMap)
  {
    map.clear();
    map = inMap;
  }

  friend std::istream& operator>><E>(std::istream& is, EnumWrapper<E>&);
private:

  //! The public, static mapping between input strings and enum values.
  static MapType map;
  
  E value;
};

template <typename E>
typename EnumWrapper<E>::MapType EnumWrapper<E>::map;

//! A templated input operator that uses the defined map to set the
//! value of the wrapped enum.
template <typename E>
std::istream& operator>>(std::istream& is, EnumWrapper<E>& e)
{
  std::string tok;
  is >> tok;
  try {
    e(tok);
  } catch (const std::out_of_range& oor) {
    throw boost::program_options::validation_error(
						   boost::program_options::validation_error::invalid_option);
  }
  return is;
}
}

#endif //ndef ENUMWRAPPER_HPP
