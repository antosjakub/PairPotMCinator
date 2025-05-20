
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")
add_executable(PairPotMCinator
    PairPotMCinator.cpp
)
target_sources(PairPotMCinator
    PRIVATE
    objects.h
    utility.h
    energy.h
    finder.h
    buckingham_pot_data.h
)

find_package(GSL REQUIRED)
find_package(nlohmann_json REQUIRED)
target_link_libraries(PairPotMCinator GSL::gsl GSL::gslcblas nlohmann_json::nlohmann_json)


install(TARGETS PairPotMCinator DESTINATION bin)