
FIND_PACKAGE(GSL REQUIRED)

add_executable(test_integration test.cpp)
target_link_libraries(test_integration "${GSL_LIBRARIES}")
target_link_libraries(test_integration "${PROJECT_SOURCE_DIR}/lib/libcuba.a")
install(TARGETS test_integration DESTINATION "${PROJECT_SOURCE_DIR}/bin")

