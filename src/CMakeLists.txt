set(NDHIST_SOURCES wiremod_ndhist.cc $ENV{SCECORR_PATH}/include/SCECorr.h)
add_executable(wiremod_ndhist ${NDHIST_SOURCES})
target_link_libraries(wiremod_ndhist
    ROOT::Core
    ROOT::Tree
    ROOT::Hist
    ROOT::RIO
    ROOT::Net
    ROOT::TreePlayer
    ROOT::Physics
)
target_include_directories(wiremod_ndhist
    PRIVATE $ENV{SCECORR_PATH}/include/
)

add_executable(ndfit ndfit.cc)
target_link_libraries(ndfit
    ROOT::Core
    ROOT::Tree
    ROOT::Hist
    ROOT::RIO
    ROOT::Net
    ROOT::TreePlayer
    ROOT::Physics
    ROOT::Minuit2
)
