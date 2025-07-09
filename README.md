
#  RaMAx

##  Build Instructions

### ✅ Prerequisites (Ubuntu/Debian)

Before building the project, make sure the following packages are installed:

```bash
sudo apt update
sudo apt install -y \
    build-essential \
    cmake \
    libcurl4-openssl-dev \
    zlib1g-dev \
    libboost-dev \
    pkg-config \
    doxygen \
    libhdf5-dev
```

> **Note:**  
> - `libboost-dev` is required for Boost.Geometry, which is header-only — no need to link Boost libraries.
> - Other third-party libraries (e.g. `spdlog`, `CLI11`, `sdsl-lite`, and `libdivsufsort`) are included in the repository and will be built automatically.

---

### ⚙️ Build Steps

```bash
# Clone the repository
git clone https://github.com/metaphysicser/RaMAx --recursive
cd RaMAx

# Create a build directory
mkdir build && cd build

# Configure the project with CMake
cmake ..

# Build using all available CPU cores
make -j$(nproc)
```

After successful build, the executable `RaMAx` will be located in the `build/` directory.

