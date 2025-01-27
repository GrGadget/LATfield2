#ifndef INT2STRING_HPP
#define INT2STRING_HPP

/*! \file int2string.hpp
   \brief int2string a small function to convert integer to string
   \author Neil Bevis
 
 */

#include "LATfield2/macros.hpp"
#include <string>

namespace LATfield2
{

/*!  
 \brief integer to string method
 \param number  : int, input integer
 \param max     : int, max number (to specify number of digit in the string.
 \param zeropad : bool, true will pad zeros, fals will not. default is true.
 
 \return string.
*/
std::string int2string(int number, int max = 999, bool zeropad = true);

}
#endif
