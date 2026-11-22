add_executable(litho_fdtd_tests
    geometry_test.cpp
    main.cpp
    material_test.cpp
    solver_test.cpp
)
target_link_libraries(litho_fdtd_tests PRIVATE litho_fdtd)
if(MSVC)
    target_compile_options(litho_fdtd_tests PRIVATE /W4)
else()
    target_compile_options(litho_fdtd_tests PRIVATE -Wall -Wextra -Wpedantic)
endif()
if(LITHO_FDTD_WARNINGS_AS_ERRORS)
    if(MSVC)
        target_compile_options(litho_fdtd_tests PRIVATE /WX)
    else()
        target_compile_options(litho_fdtd_tests PRIVATE -Werror)
    endif()
endif()
add_test(NAME litho_fdtd_tests COMMAND litho_fdtd_tests)
