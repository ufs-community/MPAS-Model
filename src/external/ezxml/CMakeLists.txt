
if(MPAS_UFS_DYCORE)
  add_library(ezxml ezxml_shim.c ezxml.c)
else()
  add_library(ezxml ezxml.c)
endif()
add_library(${PROJECT_NAME}::external::ezxml ALIAS ezxml)
target_include_directories(ezxml PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

install(TARGETS ezxml EXPORT ${PROJECT_NAME}ExportsExternal
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
