/*!
 * @file Xios.hpp
 * @date 7th February 2023
 * @author Dr Alexander Smith <as3402@cam.ac.uk>
 */

#ifndef SRC_INCLUDE_XIOS_HPP
#define SRC_INCLUDE_XIOS_HPP


#include "include/Configured.hpp"

namespace Nextsim {

//! Class to handle interfacing with the XIOS library
class Xios : public Configured<Xios> {
public:
    Xios();
    ~Xios();
    void Finalise();

    static void writeState();
    //void configure() override;
    void configure() override; 
    void initialise(int argc, char* argv[]);

    //Arguments TBC
    void setState();
    void getState();
    void writeStateData();
    void readStateData();

    enum {
        ENABLED_KEY,
    };

protected:
    bool m_isConfigured;
private:
    std::string m_enabledStr;
};

} /* end namespace Nextsim */

#endif