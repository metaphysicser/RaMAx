cmake_minimum_required(VERSION 3.10)

project("RaMAx")

# MSVC 热重载支持
if(POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT
    "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,"
    "$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

# 设置 C++20 标准
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ================= 路径设置 =================
set(SOURCE_DIR         "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(INCLUDE_DIR        "${CMAKE_CURRENT_SOURCE_DIR}/include")
set(SUBMODULE_DIR      "${CMAKE_CURRENT_SOURCE_DIR}/submodule")
set(PREPROCESS_DIR     "${SOURCE_DIR}/preprocess")
set(INDEX_DIR          "${SOURCE_DIR}/index")
set(ALIGN_DIR          "${SOURCE_DIR}/align")
set(ANCHOR_DIR         "${SOURCE_DIR}/anchor")
set(GRAPH_DIR          "${SOURCE_DIR}/graph")
set(BIGBWT_DIR         "${INDEX_DIR}/bigBWT")
set(BIGBWT_INCLUDE_DIR "${INCLUDE_DIR}/bigBWT")
set(SUFFIX_ARRAY_DIR   "${INDEX_DIR}/suffix_array")
set(SUFFIX_ARRAY_INCLUDE_DIR "${INCLUDE_DIR}/suffix_array")

# ============ 编译优化选项 ============
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -mavx2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")

# ================= 子模块依赖 =================
add_subdirectory("${SUBMODULE_DIR}/spdlog")
# add_subdirectory("${SUBMODULE_DIR}/CLI11")
# add_subdirectory("${SUBMODULE_DIR}/sdsl-lite")
add_subdirectory("${SUBMODULE_DIR}/WFA2-lib")
# 已经 add_subdirectory(...) 产生了目标 "divsufsort"

add_subdirectory("${SUBMODULE_DIR}/libdivsufsort")

# ================= 外部库依赖 =================
find_package(CURL REQUIRED)
find_package(ZLIB REQUIRED)
find_package(OpenMP)
find_package(Threads REQUIRED)
find_package(TBB REQUIRED)

# ================= 公共源文件 =================
set(RAMA_COMMON_SOURCES
    "${PREPROCESS_DIR}/raw_data_process.cpp"
    "${PREPROCESS_DIR}/fasta_manager.cpp"
    "${PREPROCESS_DIR}/newick_tree.cpp"
    "${PREPROCESS_DIR}/sequence_utils.cpp"
    "${PREPROCESS_DIR}/SeqPro.cpp"
    "${INDEX_DIR}/fm_index.cpp"
    "${INDEX_DIR}/bigBWT_build.cpp"
    "${INDEX_DIR}/bigBWT/newscan.cpp"
    "${INDEX_DIR}/bigBWT/utils.c"
    "${INDEX_DIR}/bigBWT/xerrors.c"
    "${SUFFIX_ARRAY_DIR}/gsa/gsacak.c"
    "${SUFFIX_ARRAY_DIR}/CaPS-SA/Suffix_Array.cpp"
    "${ANCHOR_DIR}/anchor.cpp"
    "${ANCHOR_DIR}/anchor_io.cpp"
    "${ANCHOR_DIR}/anchor_chunk.cpp"
    "${ANCHOR_DIR}/cluster.cpp"
    "${ALIGN_DIR}/align.cpp"
    "${ALIGN_DIR}/cigar.cpp"
    "${ALIGN_DIR}/ksw2_extz2_sse.c"
    "${ALIGN_DIR}/pair_rare_aligner.cpp"
    "${ALIGN_DIR}/multiple_rare_aligner.cpp"
    "${GRAPH_DIR}/ramesh_node.cpp"
    "${GRAPH_DIR}/ramesh_graph.cpp"
    "${GRAPH_DIR}/ramesh_export.cpp"
)

# ================= 主程序入口文件 =================
set(MAIN_G "${SOURCE_DIR}/RaMA-G.cpp")
set(MAIN_X "${SOURCE_DIR}/RaMAx.cpp")

# ================= 公共头文件路径 =================
set(RAMA_INCLUDE_DIRS
    "${INCLUDE_DIR}"
    "${SUBMODULE_DIR}/spdlog/include"
    # "${SUBMODULE_DIR}/CLI11/include"
#   "${SUBMODULE_DIR}/cereal/include"
    "${BIGBWT_INCLUDE_DIR}"
    "${SUFFIX_ARRAY_INCLUDE_DIR}"
    # "${SUBMODULE_DIR}/sdsl-lite/include"
    # "${SUBMODULE_DIR}/sdsl-lite/external/libdivsufsort/include"
    "${CMAKE_BINARY_DIR}/submodule/libdivsufsort/include"
    "${SUBMODULE_DIR}/WFA2-lib"
    "${SUBMODULE_DIR}/parlaylib/include"
)

# ================= 编译配置函数 =================
function(configure_rama_target target)
    target_include_directories(${target} PRIVATE ${RAMA_INCLUDE_DIRS})

    target_link_libraries(${target} PRIVATE
        spdlog::spdlog
        # CLI11::CLI11
        ZLIB::ZLIB
        CURL::libcurl
        # sdsl
        divsufsort
        divsufsort64
        Threads::Threads
        wfa2cpp_static
        TBB::tbb
        OpenMP::OpenMP_CXX
    )

    # 确保在构建目标前先构建libdivsufsort（生成头文件）
    add_dependencies(${target} divsufsort divsufsort64)

    target_compile_definitions(${target} PRIVATE M64)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
        target_compile_definitions(${target} PRIVATE _DEBUG_)
    endif()

    target_compile_options(${target} PRIVATE -Wall -g)
endfunction()

# ================= 创建两个可执行文件 =================
add_executable(RaMA-G ${MAIN_G} ${RAMA_COMMON_SOURCES})
configure_rama_target(RaMA-G)

add_executable(RaMAx ${MAIN_X} ${RAMA_COMMON_SOURCES})
configure_rama_target(RaMAx)

# ================= 安装规则（可选） =================
# install(TARGETS RaMAG RaMAx DESTINATION bin)
