# CMakeList.txt: RaMA-G 的 CMake 项目，在此处包括源代码并定义
# 项目特定的逻辑。
#
cmake_minimum_required (VERSION 3.8)

# 如果支持，请为 MSVC 编译器启用热重载。
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project ("RaMA-G")

set(SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(UTILS_DIR ${SOURCE_DIR}/utils)

add_subdirectory(${UTILS_DIR}/spdlog)
add_subdirectory(${UTILS_DIR}/CLI11)

# 将源代码添加到此项目的可执行文件。
add_executable(RaMA-G ${SOURCE_DIR}/RaMA-G.cpp ${UTILS_DIR}/config.hpp)


# Include directories for header files
target_include_directories(RaMA-G PRIVATE
    ${UTILS_DIR}
    ${UTILS_DIR}/spdlog/include/spdlog
    ${UTILS_DIR}/CLI11/include/CLI

)

target_link_libraries(RaMA-G PRIVATE
    spdlog::spdlog
    CLI11::CLI11
)

if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_target_properties(RaMA-G PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

endif()

# TODO: 如有需要，请添加测试并安装目标。
