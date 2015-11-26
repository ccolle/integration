

add_executable(test_integration test.cpp)
add_executable(test_genz test_genz.cpp)

set(LIBS "${GSL_LIBRARIES} ${PROJECT_SOURCE_DIR}/lib/libcuba.a" cubature)

target_link_libraries(test_integration "${LIBS}")
target_link_libraries(test_genz "${LIBS}")

set_target_properties(test_integration test_genz PROPERTIES
                        INSTALL_RPATH "${PROJECT_SOURCE_DIR}/lib"
                        INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS test_integration test_genz
    RUNTIME DESTINATION "${PROJECT_SOURCE_DIR}/bin")

