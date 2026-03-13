# add_library(GridManager
#     Grid3D.cpp
#     Field.cpp
#     ReMap3D.cpp
# )

# target_include_directories(GridManager PUBLIC
#     ${CMAKE_CURRENT_SOURCE_DIR}
# )

# target_include_directories(GridManager SYSTEM PUBLIC ${PETSC_INCLUDES})

# target_link_libraries(GridManager PUBLIC
#     PkgConfig::PETSC
#     MPI::MPI_CXX
# )

# target_compile_definitions(GridManager PUBLIC ${PETSC_DEFINITIONS})

add_library(GridManager)

target_sources(GridManager PRIVATE
    Grid3D.cpp
    Field.cpp
    ReMap3D.cpp
)

target_include_directories(GridManager PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_compile_features(GridManager PUBLIC cxx_std_17)

target_link_libraries(GridManager PUBLIC
    MPI::MPI_CXX
    PkgConfig::PETSC
)

target_compile_definitions(GridManager PUBLIC
    ${PETSC_DEFINITIONS}
)

target_include_directories(GridManager SYSTEM PUBLIC ${PETSC_INCLUDES})
