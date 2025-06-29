#include "data_process.h" 
#include <cereal/types/vector.hpp> // cereal 序列化支持 vector
#include <cereal/types/list.hpp>   // cereal 序列化支持 list
#include <cereal/types/memory.hpp> // cereal 序列化支持智能指针
#include "cereal/types/unordered_map.hpp"
#include "cereal/archives/binary.hpp"
#include <cereal/types/string.hpp>

bool saveMatchVec3D(const std::string& filename, const MatchVec3DPtr& data) {
    if (!data) return false;

    std::ofstream os(filename, std::ios::binary);
    if (!os) return false;

    cereal::BinaryOutputArchive oar(os);

    uint64_t dim0 = data->size();  // 最外层
    oar(dim0);

    for (const auto& layer : *data) {
        uint64_t dim1 = layer.size();  // 第二层
        oar(dim1);

        for (const auto& vec : layer) {
            uint64_t dim2 = vec.size();  // 最内层
            oar(dim2);

            for (const auto& match : vec) {
                oar(match);  // 直接序列化 Match
            }
        }
    }

    return static_cast<bool>(os);
}

bool loadMatchVec3D(const std::string& filename, MatchVec3DPtr& data) {
    std::ifstream is(filename, std::ios::binary);
    if (!is) return false;

    cereal::BinaryInputArchive iar(is);

    uint64_t dim0;
    iar(dim0);

    data = std::make_shared<MatchVec3D>();
    data->resize(dim0);

    for (uint64_t i = 0; i < dim0; ++i) {
        uint64_t dim1;
        iar(dim1);

        (*data)[i].resize(dim1);

        for (uint64_t j = 0; j < dim1; ++j) {
            uint64_t dim2;
            iar(dim2);

            MatchVec vec;
            vec.reserve(dim2);

            for (uint64_t k = 0; k < dim2; ++k) {
                Match m;
                iar(m);
                vec.push_back(std::move(m));
            }

            (*data)[i][j] = std::move(vec);
        }
    }

    return static_cast<bool>(is);
}

bool saveSpeciesMatchMap(const std::string& filename, const SpeciesMatchVec3DPtrMapPtr& map_ptr) {
    if (!map_ptr) return false;

    std::ofstream os(filename, std::ios::binary);
    if (!os) return false;

    cereal::BinaryOutputArchive oar(os);

    uint64_t map_size = map_ptr->size();
    oar(map_size);

    for (const auto& [species_name, match_ptr] : *map_ptr) {
        oar(species_name);

        if (!match_ptr) {
            oar(uint64_t(0));  // 空指针表示为0层
            continue;
        }

        uint64_t dim0 = match_ptr->size();
        oar(dim0);

        for (const auto& layer : *match_ptr) {
            uint64_t dim1 = layer.size();
            oar(dim1);

            for (const auto& vec : layer) {
                uint64_t dim2 = vec.size();
                oar(dim2);

                for (const auto& m : vec) {
                    oar(m);
                }
            }
        }
    }

    return static_cast<bool>(os);
}

bool loadSpeciesMatchMap(const std::string& filename, SpeciesMatchVec3DPtrMapPtr& map_ptr) {
    std::ifstream is(filename, std::ios::binary);
    if (!is) return false;

    cereal::BinaryInputArchive iar(is);

    uint64_t map_size;
    iar(map_size);

    map_ptr = std::make_shared<SpeciesMatchVec3DPtrMap>();

    for (uint64_t i = 0; i < map_size; ++i) {
        std::string species;
        iar(species);

        uint64_t dim0;
        iar(dim0);

        MatchVec3DPtr data = std::make_shared<MatchVec3D>();
        data->resize(dim0);

        for (uint64_t i0 = 0; i0 < dim0; ++i0) {
            uint64_t dim1;
            iar(dim1);
            (*data)[i0].resize(dim1);

            for (uint64_t i1 = 0; i1 < dim1; ++i1) {
                uint64_t dim2;
                iar(dim2);

                MatchVec vec;
                vec.reserve(dim2);

                for (uint64_t i2 = 0; i2 < dim2; ++i2) {
                    Match m;
                    iar(m);
                    vec.push_back(std::move(m));
                }

                (*data)[i0][i1] = std::move(vec);
            }
        }

        (*map_ptr)[species] = std::move(data);
    }

    return static_cast<bool>(is);
}


