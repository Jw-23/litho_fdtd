#include "test_support.hpp"

#include <cstdio>

void run_geometry_tests();
void run_material_tests();
void run_solver_tests();

int main() {
    run_geometry_tests();
    run_material_tests();
    run_solver_tests();
    if (test::failures == 0) {
        std::printf("litho_fdtd tests passed\n");
        return 0;
    }
    std::fprintf(stderr, "litho_fdtd tests failed: %d check(s)\n", test::failures);
    return 1;
}
