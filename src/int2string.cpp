#include "config.h"
#include "int2string.hpp"

namespace LATfield2
{
/*!  
 \brief integer to string method
 \param number  : int, input integer
 \param max     : int, max number (to specify number of digit in the string.
 \param zeropad : bool, true will pad zeros, fals will not. default is true.
 
 \return string.
*/
std::string int2string(int number, int max, bool zeropad)
{
  std::string output;
  char c;
  int i;

  //Get number of digits needed for max
  int digits=1;
  for(i=10; (i-1)<max; i*=10) digits++;

  //Create string of length digits (padded with zeros)
  for(i=0; i<digits; i++)
    {
      c = '0' + number%10;
      if(c!='0' || zeropad)
	{
	  output = c + output;
	}
      number /= 10;
    }

  return output;
}

}
