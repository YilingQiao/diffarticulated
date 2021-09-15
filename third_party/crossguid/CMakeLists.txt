
ADD_LIBRARY(crossguid STATIC guid.cpp)

IF(WIN32)
	target_compile_definitions(crossguid PUBLIC GUID_WINDOWS PRIVATE GUID_WINDOWS)
ELSE()
IF(APPLE)
	target_compile_definitions(crossguid PUBLIC GUID_LIBUUID PRIVATE GUID_LIBUUID)
	target_link_libraries(crossguid PUBLIC uuid PRIVATE uuid)
ENDIF()
ENDIF()

target_compile_definitions(crossguid PUBLIC GUID_LIBUUID PRIVATE GUID_LIBUUID)
target_include_directories(crossguid PUBLIC . PRIVATE .)


