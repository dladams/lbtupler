
art_make(BASENAME_ONLY
         LIB_LIBRARIES
           Geometry_service
           Geometry
           LArProperties_service
           DetectorProperties_service
           RawData
           RecoBase
           ${SIMULATIONBASE}
           ${ART_FRAMEWORK_CORE}
           ${ART_FRAMEWORK_PRINCIPAL}
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_FRAMEWORK_SERVICES_OPTIONAL}
           ${ART_FRAMEWORK_SERVICES_OPTIONAL_TFILESERVICE_SERVICE}
           ${ART_PERSISTENCY_COMMON}
           ${ART_PERSISTENCY_PROVENANCE}
           ${ART_UTILITIES}
           ${MF_MESSAGELOGGER}
           ${MF_UTILITIES}
           ${CETLIB}
           ${ROOT_BASIC_LIB_LIST}
           ${Boost_SYSTEM_LIBRARY}
         MODULE_LIBRARIES
           LbTupler
)

install_headers()
install_fhicl()
install_source( README.txt )
install_scripts()
