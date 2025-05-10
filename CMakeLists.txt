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
set(SOURCE_DIR         "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(INCLUDE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SUBMODULE_DIR      "${CMAKE_CURRENT_SOURCE_DIR}/submodule")
set(PREPROCESS_DIR     "${SOURCE_DIR}/preprocess")
set(INDEX_DIR          "${SOURCE_DIR}/index")
set(ALIGN_DIR          "${SOURCE_DIR}/align")
set(ANCHOR_DIR         "${SOURCE_DIR}/anchor")
set(BIGBWT_DIR         "${INDEX_DIR}/bigBWT")
set(BIGBWT_INCLUDE_DIR "${INCLUDE_DIR}/bigBWT")
set(SUFFIX_ARRAY_DIR   "${INDEX_DIR}/suffix_array")
set(SUFFIX_ARRAY_INCLUDE_DIR "${INCLUDE_DIR}/suffix_array")

# ============ 编译优化选项 ============
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -mavx2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

# ============ 添加子模块 ============
add_subdirectory("${SUBMODULE_DIR}/spdlog")
add_subdirectory("${SUBMODULE_DIR}/CLI11")
add_subdirectory("${SUBMODULE_DIR}/sdsl-lite")
add_subdirectory("${SUBMODULE_DIR}/WFA2-lib")

# ============ 安装 libdivsufsort（使用 install.sh） ============
if (NOT EXISTS "${SUBMODULE_DIR}/sdsl-lite/include/divsufsort.h")
    message(STATUS "divsufsort.h not found, running install.sh")

    add_custom_target(sdsl-install ALL
        # ① 先确保 sdsl-lite/build 目录存在
        # COMMAND ${CMAKE_COMMAND} -E make_directory build
        # ② 再执行 install.sh
        COMMAND bash install.sh ${SUBMODULE_DIR}/sdsl-lite/
        WORKING_DIRECTORY ${SUBMODULE_DIR}/sdsl-lite/
        COMMENT "Running install.sh from inside sdsl-lite directory"
        VERBATIM
    )
else()
    message(STATUS "divsufsort.h found, skipping install.sh")
endif()


# ============ 查找系统库 ============
find_package(CURL REQUIRED)
find_package(ZLIB REQUIRED)
find_package(OpenMP)
find_package(Threads REQUIRED)

# ============ 创建主可执行程序 ============
add_executable(RaMA-G
    "${SOURCE_DIR}/RaMA-G.cpp"
    "${SOURCE_DIR}/rare_aligner.cpp"
    "${PREPROCESS_DIR}/raw_data_process.cpp"
    "${PREPROCESS_DIR}/split_data.cpp"
    "${PREPROCESS_DIR}/fasta_manager.cpp"
    "${INDEX_DIR}/fm_index.cpp"
    "${BIGBWT_DIR}/newscan.cpp"
    "${BIGBWT_DIR}/utils.c"
    "${BIGBWT_DIR}/xerrors.c"
    "${SUFFIX_ARRAY_DIR}/gsa/gsacak.c"
    "${SUFFIX_ARRAY_DIR}/CaPS-SA/Suffix_Array.cpp"
    "${ANCHOR_DIR}/anchor.cpp"
    "${ALIGN_DIR}/cigar.cpp"
)

# ============ 添加头文件路径 ============
target_include_directories(RaMA-G PRIVATE
    "${INCLUDE_DIR}"
    "${SUBMODULE_DIR}/spdlog/include"
    "${SUBMODULE_DIR}/CLI11/include"
    "${SUBMODULE_DIR}/cereal/include"
    "${BIGBWT_INCLUDE_DIR}"
    "${SUFFIX_ARRAY_INCLUDE_DIR}"
    "${SUBMODULE_DIR}/sdsl-lite/include"
    "${SUBMODULE_DIR}/sdsl-lite/external/libdivsufsort/include"
    "${SUBMODULE_DIR}/WFA2-lib"
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
    Threads::Threads
    wfa2cpp_static
)

# ============ C++ 标准设置 ============
if(CMAKE_VERSION VERSION_GREATER 3.12)
  set_target_properties(RaMA-G PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
  )
endif()

# ============ 宏定义 ============
target_compile_definitions(RaMA-G PRIVATE M64)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(RaMA-G PRIVATE _DEBUG_)
endif()

# ============ 警告和调试信息 ============
add_definitions("-Wall -g")

# ============ 添加依赖关系（必须在定义 RaMA-G 后） ============
if (TARGET sdsl-install)
    add_dependencies(RaMA-G sdsl-install)
endif()

# 可选安装目标
# install(TARGETS RaMA-G DESTINATION bin)
