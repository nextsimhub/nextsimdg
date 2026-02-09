#define DOCTEST_CONFIG_IMPLEMENT

#include "include/Finalizer.hpp"
#include <doctest/doctest.h>
#ifdef USE_KOKKOS
#include <Kokkos_Core.hpp>
#endif

// Main routine for tests that use Kokkos in some way.
int main(int argc, char** argv)
{
#ifdef USE_KOKKOS
    Kokkos::initialize(argc, argv);
#endif

    doctest::Context ctx;
    ctx.applyCommandLine(argc, argv);

    int test_result = ctx.run();

    // cleanup resources including Kokkos allocated buffers before finalizing Kokkos to prevent
    // errors
    Nextsim::Finalizer::finalize();
#ifdef USE_KOKKOS
    Kokkos::finalize();
#endif

    return test_result;
}
