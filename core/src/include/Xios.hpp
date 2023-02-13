/*!
 * @file Xios.hpp
 * @date 7th February 2023
 * @author Dr Alexander Smith <as3402@cam.ac.uk>
 */

#ifndef SRC_INCLUDE_XIOS_HPP
#define SRC_INCLUDE_XIOS_HPP

namespace Nextsim {

//! Class to handle interfacing with the XIOS library
class Xios {
public:
    static void writeState();
    static void configure(int argc, char* argv[]);
protected:
    bool isConfigured;
private:


};

} /* end namespace Nextsim */

#endif