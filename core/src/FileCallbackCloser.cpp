/*!
 * @file FileCallbackCloser.cpp
 *
 * @date 15 May 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#include "include/FileCallbackCloser.hpp"

namespace Nextsim {

std::list<FileCallbackCloser::ClosingFn> FileCallbackCloser::closingFns;

void FileCallbackCloser::onClose(FileCallbackCloser::ClosingFn fn) {
    closingFns.push_back(fn);
}

void FileCallbackCloser::close(const std::string& filename) {
    for (auto& fn: closingFns) {
        fn(filename);
    }
}

void FileCallbackCloser::clearAllClose()
{
    closingFns.clear();
}
} /* namespace Nextsim */
