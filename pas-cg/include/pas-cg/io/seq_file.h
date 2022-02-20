#pragma once
#include "xmd/model/model.h"
#include <filesystem>

namespace xmd {
    class seq_file {
    public:
        explicit seq_file(std::filesystem::path const& seq_file_path);
        model const& to_model() const;

    private:
        model xmd_model;
    };
}