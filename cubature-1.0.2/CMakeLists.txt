
add_library(cubature SHARED pcubature.c hcubature.c) # clencurt_gen.c)

install(TARGETS cubature DESTINATION "${PROJECT_SOURCE_DIR}/lib")
