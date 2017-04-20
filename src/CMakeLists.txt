# List all src files below inside set
SET(FERDOSI_src

${SRCFERDOSI}Ferdosi.f90
${SRCFERDOSI}FerdosiTesterMain.f90
${SRCFERDOSI}ModuleAccuracy.f90
${SRCFERDOSI}ModuleGet.f90
${SRCFERDOSI}ModuleParameters.f90
${SRCFERDOSI}ModuleReadFiles.f90
${SRCFERDOSI}ModuleRunFerdosi.f90
${SRCFERDOSI}ModuleSireType.f90
${SRCFERDOSI}ModuleTLC.f90
${SRCFERDOSI}ModuleWriteFiles.f90
)

add_definitions(-D COMMIT=${VERSION})
ADD_EXECUTABLE(${PROGRAMEXE} ${FERDOSI_src})


# add external links here
TARGET_LINK_LIBRARIES(${PROGRAMEXE} ${AHLIB})
INSTALL(TARGETS ${PROJECTNAMETOREPLACEEXE} RUNTIME DESTINATION bin)