/*!
 * @file XiosInit_test.cpp
 * @brief Initialisation Testing for XIOS
 * @date Feb 28, 2023
 * @author Dr Tom Meltzer <tdm39@cam.ac.uk>
 * @author Dr Alexander Smith <as3402@cam.ac.uk>
 *
 * Test file intended to validate the initialization of the server process of
 * the XIOS class, including syncronising the state of the Nextsim grid. This
 * file serves as a set of unit tests, the system tests are elsewhere.
 */

#define DOCTEST_CONFIG_IMPLEMENT
#include <doctest/extensions/doctest_mpi.h>


int main(int argc, char** argv) {
  doctest::mpi_init_thread(argc,argv,MPI_THREAD_MULTIPLE);

  doctest::Context ctx;
  ctx.setOption("reporters", "MpiConsoleReporter");
  ctx.setOption("reporters", "MpiFileReporter");
  ctx.setOption("force-colors", true);
  ctx.applyCommandLine(argc, argv);

  int test_result = ctx.run();

  doctest::mpi_finalize();

  return test_result;
}
