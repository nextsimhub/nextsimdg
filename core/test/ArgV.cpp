/*!
 * @file ArgV.cpp
 *
 * @date Oct 11, 2021
 * @author Tim Spain
 */

#include "ArgV.hpp"

namespace Nextsim {

ArgV::ArgV(std::vector<std::string> vs)
{
    nStrings = vs.size();
    ppc = new char*[nStrings + 1]; // Add 1 for the conventional final \0
    for (int i = 0; i < nStrings; ++i) {
        int len = vs[i].size();
        char* pc = new char[len + 1];
        for (int j = 0; j < len; ++j) {
            pc[j] = vs[i][j];
        }
        pc[len] = '\0';
        ppc[i] = pc;
    }
    // Ensure argv[argc] = "\0", as required by the standard.
    char* pnull = new char[1];
    pnull[0] = '\0';
    ppc[nStrings] = pnull;
}

ArgV::~ArgV()
{
    for (int i = 0; i < nStrings; ++i) {
        delete[] ppc[i];
    }
    delete[] ppc;
}

char** ArgV::operator()() { return ppc; }

int ArgV::argc() { return nStrings; }
}
