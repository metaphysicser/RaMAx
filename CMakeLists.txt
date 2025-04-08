cmake_minimum_required(VERSION 3.8)

project("RaMA-G")

# MSVC 热重载支持
if(POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT
    "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,"
    "$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

# ============ 路径设置 ============
set(SOURCE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(UTILS_DIR         "${SOURCE_DIR}/utils")
set(PREPROCESS_DIR    "${SOURCE_DIR}/preprocess")
set(INDEX_DIR         "${SOURCE_DIR}/index")
set(BIGBWT_DIR        "${INDEX_DIR}/bigBWT")
set(SUFFIX_ARRAY_DIR  "${INDEX_DIR}/suffix_array")
set(R_INDEX_DIR       "${INDEX_DIR}/r-index")

# 假设 parlaylib 被下载到 utils/parlaylib 目录
set(PARLAYLIB_DIR "${UTILS_DIR}/parlaylib")
set(PARLAYLIB_INCLUDE_DIR "${PARLAYLIB_DIR}/include")

# 添加 parlaylib 头文件路径
include_directories(${PARLAYLIB_INCLUDE_DIR})

# 如果 parlaylib 是一个源代码库而不是编译后的库，你可以直接将其添加到项目中进行编译
add_subdirectory(${PARLAYLIB_DIR})

# ============ 添加子模块 ============
add_subdirectory("${UTILS_DIR}/spdlog")
add_subdirectory("${UTILS_DIR}/CLI11")
add_subdirectory("${UTILS_DIR}/sdsl-lite")

# ============ 安装 libdivsufsort（使用 install.sh） ============
# 检查是否需要安装 libdivsufsort
if (NOT EXISTS "${UTILS_DIR}/sdsl-lite/include/divsufsort.h")
    message(STATUS "divsufsort.h not found, running install.sh")

    add_custom_target(sdsl-install ALL
        COMMAND bash install.sh ${UTILS_DIR}/sdsl-lite/
        WORKING_DIRECTORY ${UTILS_DIR}/sdsl-lite
        COMMENT "Running install.sh from inside sdsl-lite directory"
    )
else()
    message(STATUS "divsufsort.h found, skipping install.sh")
endif()

# 如果安装任务存在，则添加依赖
if (TARGET sdsl-install)
    add_dependencies(RaMA-G sdsl-install)
endif()


# ============ 查找系统库 ============
find_package(CURL REQUIRED)
find_package(ZLIB REQUIRED)

# ============ 创建主可执行程序 ============
add_executable(RaMA-G
    "${SOURCE_DIR}/RaMA-G.cpp"
    "${PREPROCESS_DIR}/raw_data_process.cpp"
    "${PREPROCESS_DIR}/split_data.cpp"
    "${PREPROCESS_DIR}/fasta_manager.cpp"
    "${INDEX_DIR}/index.cpp"
    "${BIGBWT_DIR}/newscan.cpp"
    "${BIGBWT_DIR}/utils.c"
    "${BIGBWT_DIR}/xerrors.c"
    "${SUFFIX_ARRAY_DIR}/gsa/gsacak.c"
    "${SUFFIX_ARRAY_DIR}/CaPS-SA/Suffix_Array.cpp"
    )

# ============ 添加头文件路径 ============
target_include_directories(RaMA-G PRIVATE
    "${UTILS_DIR}"
    "${UTILS_DIR}/spdlog/include"
    "${UTILS_DIR}/CLI11/include"
    "${UTILS_DIR}/cereal/include"
    "${PREPROCESS_DIR}"
    "${INDEX_DIR}"
    "${BIGBWT_DIR}"
    "${SUFFIX_ARRAY_DIR}"
    "${R_INDEX_DIR}"
    "${UTILS_DIR}/sdsl-lite/include"
    "${UTILS_DIR}/sdsl-lite/external/libdivsufsort/include"
)

# ============ 链接库 ============
target_link_libraries(RaMA-G PRIVATE
    spdlog::spdlog
    CLI11::CLI11
    ZLIB::ZLIB
    CURL::libcurl
    sdsl
    divsufsort
    divsufsort64
    pthread
)


# ============ C++ 标准设置 ============
if(CMAKE_VERSION VERSION_GREATER 3.12)
  set_target_properties(RaMA-G PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
  )
endif()

# 启用 AVX2 支持
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2 -march=native")


# ============ Debug 宏 ============
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(RaMA-G PRIVATE _DEBUG_)
endif()

# ============ 警告和调试信息 ============
add_definitions("-Wall -g")

# 可选安装目标
# install(TARGETS RaMA-G DESTINATION bin)
