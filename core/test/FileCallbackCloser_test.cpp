/*!
 * @file FileCallbackCloser_test.cpp
 *
 * @date 15 May 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

#include "include/FileCallbackCloser.hpp"

#include <filesystem>
#include <functional>
#include <ostream>

namespace Nextsim {

const std::string testFileName = "testFileName.nc";

class FileHolder {
public:
    static void open(const std::string& fileName)
    {
        openFiles.emplace(fileName, fileName);
    }
    static void close(const std::string& fileName)
    {
        openFiles[fileName].close();
        openFiles.erase(fileName);
    }
    static bool isOpen(const std::string& fileName)
    {
        return openFiles.count(fileName) > 0;
    }
private:
    static std::map<std::string, std::ofstream> openFiles;
};

std::map<std::string, std::ofstream> FileHolder::openFiles;

TEST_SUITE_BEGIN("FileCallbackCloser");
TEST_CASE("Test the FileHolder class")
{
    FileHolder::open(testFileName);
    REQUIRE(FileHolder::isOpen(testFileName));
    FileHolder::close(testFileName);
    REQUIRE(!FileHolder::isOpen(testFileName));
    std::filesystem::remove(testFileName);
}

TEST_CASE("Close a file")
{
    FileCallbackCloser::onClose(FileHolder::close);

    FileHolder::open(testFileName);
    REQUIRE(FileHolder::isOpen(testFileName));
    FileCallbackCloser::close(testFileName);
    REQUIRE(!FileHolder::isOpen(testFileName));
    std::filesystem::remove(testFileName);
}
TEST_SUITE_END();
} /* namespace Nextsim */
