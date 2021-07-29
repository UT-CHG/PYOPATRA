//
// Created by Georgia Stuart on 2/2/21.
//

#define CATCH_CONFIG_RUNNER
#include <catch2/catch.hpp>
#include "mpi.h"

int main( int argc, char* argv[] ) {
    MPI_Init(&argc, &argv);
    int result = Catch::Session().run( argc, argv );
    MPI_Finalize();
    return result;
}