#
# Copyright 2020 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


include_directories(. ..)

set(DL ${CMAKE_DL_LIBS})

find_package(PkgConfig)

add_definitions(-DBT_USE_DOUBLE_PRECISION)


SUBDIRS(../third_party/crossguid)
set(crossguid_INCLUDE_DIRS ../third_party/crossguid)
include_directories(${crossguid_INCLUDE_DIRS})



message("USING third_party/zmq")
SET(cppzmq_FOUND TRUE)
SUBDIRS(../third_party/zeromq)
set(cppzmq_INCLUDE_DIRS ../third_party/zeromq/include )
include_directories(${cppzmq_INCLUDE_DIRS})
add_definitions(-DZMQ_STATIC)



message ("USING third_party/nlohmann")
SET(nlohmann_json_FOUND TRUE)
set(nlohmann_json_INCLUDE_DIRS ../third_party/nlohmann/include)
include_directories(${nlohmann_json_INCLUDE_DIRS})



message ("USING third_party/tinyxml2")
set(tinyxml2_FOUND TRUE)
set(tinyxml2_INCLUDE_DIRS ../third_party/tinyxml2/include)
SUBDIRS(../third_party/tinyxml2)



set(MESHCAT_LIBRARIES zmq crossguid uuid  cpp_base64  LinearMath
            OpenGLWindow  gwen BussIK  pthread)
IF(WIN32)
    set(MESHCAT_LIBRARIES ${MESHCAT_LIBRARIES} ws2_32 Winmm iphlpapi)
ENDIF()
set(HAS_MESHCAT TRUE)
