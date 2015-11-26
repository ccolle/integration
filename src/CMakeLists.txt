

add_executable(test_integration test.cpp)
target_link_libraries(test_integration "${GSL_LIBRARIES}")
target_link_libraries(test_integration "${PROJECT_SOURCE_DIR}/lib/libcuba.a" cubature)
set_target_properties(test_integration PROPERTIES
                        INSTALL_RPATH "${PROJECT_SOURCE_DIR}/lib"
                        INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS test_integration RUNTIME DESTINATION "${PROJECT_SOURCE_DIR}/bin")

