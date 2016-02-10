

add_executable(test_integration test.cpp)
add_executable(test_genz test_genz.cpp)
add_executable(test_cuba cubatest.cpp)

set(LIBS "${GSL_LIBRARIES} ${PROJECT_SOURCE_DIR}/lib/libcuba.a" cubature)

target_link_libraries(test_integration "${LIBS}")
target_link_libraries(test_genz "${LIBS}")
target_link_libraries(test_cuba "${LIBS}")

set_target_properties(test_integration test_genz test_cuba PROPERTIES
                        INSTALL_RPATH "${PROJECT_SOURCE_DIR}/lib"
                        INSTALL_RPATH_USE_LINK_PATH TRUE)

install(TARGETS test_integration test_genz test_cuba
    RUNTIME DESTINATION "${PROJECT_SOURCE_DIR}/bin")

