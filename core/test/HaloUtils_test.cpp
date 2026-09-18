/*!
 * @author  Mikolaj Adam Kowalski <mak60@cam.ac.uk>
 * @brief   Unit tests for the Halo and other MPI utilites that not require MPI to run
 * @details
 *
 * Ordinary doctest test (single rank, no MPI asserts over the test) covering
 * the pure utilities in Halo.hpp / Halo.cpp which do not need the MPI
 * machinery:
 *
 */

#include "include/Halo.hpp"
#include <doctest/doctest.h>

#include <vector>

namespace Nextsim {
namespace HaloExchange {

    TEST_SUITE("Halo utilities")
    {
        TEST_CASE("rotate180InPlace")
        {
            SUBCASE("Row pointers")
            {
                std::vector<std::vector<int>> matrix = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };

                // Create vectors of pointers to the start and end of each row
                std::vector<int*> starts;
                std::vector<int*> ends;
                for (auto& row : matrix) {
                    starts.push_back(row.data());
                    ends.push_back(row.data() + row.size());
                }

                // Call the rotate function
                rotate180InPlace(starts, ends);

                // Verify the result against the expected matrix
                std::vector<std::vector<int>> expected = { { 9, 8, 7 }, { 6, 5, 4 }, { 3, 2, 1 } };
                REQUIRE(matrix == expected);
            }

            SUBCASE("Iterators")
            {
                std::vector<std::vector<int>> matrix = { { 1, 2, 3 }, { 4, 5, 6 }, { 7, 8, 9 } };

                // Create vectors of iterators to the start and end of each row
                std::vector<std::vector<int>::iterator> starts;
                std::vector<std::vector<int>::iterator> ends;
                for (auto& row : matrix) {
                    starts.push_back(row.begin());
                    ends.push_back(row.end());
                }

                // Call the rotate function
                rotate180InPlace(starts, ends);

                // Verify the result against the expected matrix
                std::vector<std::vector<int>> expected = { { 9, 8, 7 }, { 6, 5, 4 }, { 3, 2, 1 } };
                REQUIRE(matrix == expected);
            }

            SUBCASE("Single row")
            {
                std::vector<std::vector<int>> matrix = { { 1, 2, 3 } };

                // Create vectors of pointers to the start and end of each row
                std::vector<int*> starts;
                std::vector<int*> ends;
                for (auto& row : matrix) {
                    starts.push_back(row.data());
                    ends.push_back(row.data() + row.size());
                }

                // Call the rotate function
                rotate180InPlace(starts, ends);

                // Verify the result against the expected matrix
                std::vector<std::vector<int>> expected = { { 3, 2, 1 } };
                REQUIRE(matrix == expected);
            }

            SUBCASE("Empty matrix")
            {
                std::vector<std::vector<int>> matrix = {};

                // Create vectors of pointers to the start and end of each row
                std::vector<int*> starts;
                std::vector<int*> ends;

                // Call the rotate function
                rotate180InPlace(starts, ends);

                // Verify the result against the expected (empty) matrix
                std::vector<std::vector<int>> expected = {};
                REQUIRE(matrix == expected);
            }

            SUBCASE("Even number of rows")
            {
                std::vector<std::vector<int>> matrix = { { 1, 2 }, { 3, 4 }, { 5, 6 }, { 7, 8 } };

                // Create vectors of pointers to the start and end of each row
                std::vector<int*> starts;
                std::vector<int*> ends;
                for (auto& row : matrix) {
                    starts.push_back(row.data());
                    ends.push_back(row.data() + row.size());
                }

                // Call the rotate function
                rotate180InPlace(starts, ends);

                // Verify the result against the expected matrix
                std::vector<std::vector<int>> expected = { { 8, 7 }, { 6, 5 }, { 4, 3 }, { 2, 1 } };
                REQUIRE(matrix == expected);
            }
        }

        TEST_CASE("rotate180InPlace with invalid input")
        {
            std::vector<std::vector<int>> matrix = { { 1, 2, 3 }, { 4, 5, 6 } };

            SUBCASE("Mismatched sizes")
            {
                std::vector<int*> starts = { matrix[0].data() };
                std::vector<int*> ends
                    = { matrix[0].data() + matrix[0].size(), matrix[1].data() + matrix[1].size() };

                REQUIRE_THROWS_AS(rotate180InPlace(starts, ends), std::invalid_argument);
            }

            SUBCASE("Different row sizes")
            {
                std::vector<int*> starts = { matrix[0].data(), matrix[1].data() };
                std::vector<int*> ends = { matrix[0].data() + matrix[0].size(),
                    matrix[1].data() + matrix[1].size() - 1 };

                REQUIRE_THROWS_AS(rotate180InPlace(starts, ends), std::invalid_argument);
            }
        }
    }
} // end of HaloExchange namespace
} // end of nextsim namespace
