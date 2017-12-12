add_executable(${ExampleName} ${ExampleName}.F90)

target_link_libraries(${ExampleName} PRIVATE opencmisslibs)

install(TARGETS ${ExampleName} DESTINATION bin)

configure_file(zero-rate.xml zero-rate.xml)

configure_file(visualise.cmgui visualise.cmgui)
