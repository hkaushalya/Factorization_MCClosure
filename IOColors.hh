#ifndef IOCOLORS_HH
#define IOCOLORS_HH

/*******************************************************************
* These methods will colorize the terminal outputs.
* usage:
* 1. #include "IOColors.hh
* 2. std::cout << red  << "this text is red" << clearatt << std::endl;
* CONS:
* 1. These methods adds hiddend characters to the stream. If you save 
* the output to a file you'll see them.
* 2. Make sure to add 'reset' the default color scheme by adding
* 'clearatt' to the end of your stream.
* 3. The underline and bold does not seem to work.
*  Author: Sam Hewamanage <samantha@fnal.gov>
********************************************************************
*	$Id: IOColors.hh,v 1.1 2010/04/16 21:00:46 samantha Exp $
*	$Log: IOColors.hh,v $
*	Revision 1.1  2010/04/16 21:00:46  samantha
*	Creted this to get color outputs to my jobs. These were inside the
*	CommonTools.cc/hh. But it became circular and decided to move these to a new
*	file.
*	
********************************************************************/

#include <ostream>
//for color outputs

//foreground colors
inline std::ostream& red      (std::ostream &s) { s << "\033[31m"; return s; }
inline std::ostream& green    (std::ostream &s) { s << "\033[32m"; return s; }
inline std::ostream& yellow   (std::ostream &s) { s << "\033[33m"; return s; }
inline std::ostream& blue     (std::ostream &s) { s << "\033[34m"; return s; }
inline std::ostream& magenta  (std::ostream &s) { s << "\033[35m"; return s; }
inline std::ostream& cyan     (std::ostream &s) { s << "\033[36m"; return s; }
inline std::ostream& white    (std::ostream &s) { s << "\033[37m"; return s; }
inline std::ostream& black    (std::ostream &s) { s << "\033[0m" ; return s; }
inline std::ostream& clearatt (std::ostream &s) { s << "\033[0m" ; return s; }
//special
inline std::ostream& bold (std::ostream &s) { s << "\033[1m" ; return s; }   //does not seem to  work
inline std::ostream& undeline (std::ostream &s) { s << "\033[4m" ; return s; }
inline std::ostream& blink (std::ostream &s) { s << "\033[5m" ; return s; }   //does not seem to  work
//background colors
inline std::ostream& redbg      (std::ostream &s) { s << "\033[41m"; return s; }
inline std::ostream& greenbg    (std::ostream &s) { s << "\033[42m"; return s; }
inline std::ostream& yellowbg   (std::ostream &s) { s << "\033[43m"; return s; }
inline std::ostream& bluebg     (std::ostream &s) { s << "\033[44m"; return s; }
inline std::ostream& magentabg  (std::ostream &s) { s << "\033[45m"; return s; }
inline std::ostream& cyanbg     (std::ostream &s) { s << "\033[46m"; return s; }
inline std::ostream& whitebg    (std::ostream &s) { s << "\033[47m"; return s; }
inline std::ostream& blackwhite    (std::ostream &s) { s << "\033[7m" ; return s; }; //blk text on white background

#endif
