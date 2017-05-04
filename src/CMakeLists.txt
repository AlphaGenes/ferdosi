# List all src files below inside set
SET(FERDOSI_src

Ferdosi.f90
FerdosiTesterMain.f90
ModuleAccuracy.f90
ModuleGet.f90
ModuleParameters.f90
# ${SRCFERDOSI}ModuleReadFiles.f90
ModuleRunFerdosi.f90
ModuleSireType.f90
ModuleWriteFiles.f90
)


if(FERLIB)
    SET(FERDOSILIB_src  ${FERDOSI_src}) 
    ADD_LIBRARY(${FERLIB} STATIC ${FERDOSILIB_src})
    TARGET_LINK_LIBRARIES(${FERLIB} ${AHLIB})
else()
    add_definitions(-D COMMIT=${VERSION})
    ADD_EXECUTABLE(${PROGRAMEXE} ${FERDOSI_src})


    # add external links here
    TARGET_LINK_LIBRARIES(${PROGRAMEXE} ${AHLIB})
    INSTALL(TARGETS ${PROJECTNAMETOREPLACEEXE} RUNTIME DESTINATION bin)
endif(FERLIB)
