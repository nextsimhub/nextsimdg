/*!
 * @file MainMPI.cpp
 *
 * @date 13 Sep 2024
 * @author Einar Ólason <einar.olason@nersc.no>
 *
 * This file is required so that the MPI enabled doc-test OASISCoupledOcean_test.cpp (and any
 * others) can run.
 */

#define DOCTEST_CONFIG_IMPLEMENT

#include <doctest/extensions/doctest_mpi.h>

int main(int argc, char** argv)
{
    doctest::mpi_init_thread(argc, argv, MPI_THREAD_MULTIPLE);

    doctest::Context ctx;
    ctx.setOption("reporters", "MpiConsoleReporter");
    ctx.setOption("reporters", "MpiFileReporter");
    ctx.setOption("force-colors", true);
    ctx.applyCommandLine(argc, argv);

    int test_result = ctx.run();

    doctest::mpi_finalize();

    return test_result;
}
