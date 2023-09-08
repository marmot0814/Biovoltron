Installation {#installation}
============

## Requirements
Other dependency will be locally downloaded by our submodule design.
- gcc >= 11.2
- cmake >= 3.16
- git
- zlib

```bash
sudo apt install git \
                 cmake \
                 make \
                 zlib1g-dev
```

## Recommended project layout
For this project, we recommend following layout.
```
$ tree app
app
├── CMakeLists.t
├── Biovoltron
│   ├── CMakeLists.txt
│   ├── include
│   ├── submodules
│   └── ...
├── build
└── src
    └── main.cpp
```

To set up this project layout, you can follow the script below.
```bash
mkdir app
cd app
mkdir src build
git clone --recurse-submodules http://github.com/jhhung/Biovoltron
```

> You can also install the dependency yourself if you are not satisfied with current submodules design.
> However, you will need to modify the current cmake file to make it work!


To test whether `Biovoltron` is properly set up, run the unit test.
```bash
cd Biovoltron
mkdir build && cd build
cmake ..
make -j
./tests/biovoltron-test
```


## How to build your first example with Biovoltron
After we passed the unit test, we can now compile and run a small example.
First, we create our `main.cpp` under the `src` directory with the following content.
```cpp
#include <iostream>
#include <sstream>
#include <biovoltron/file_io/fasta.hpp>

auto fasta = R"(>hello-world
ACTG)";

int main() {
    auto rec = biovoltron::FastaRecord<true>{};
    std::istringstream(fasta) >> rec;
    std::cout << rec.name << '\n';
    std::cout << rec.seq << '\n';
    return 0;
}
```

To compile the source code, create the following `CMakeLists.txt`:
```cmake
cmake_minimum_required (VERSION 3.16)

# your project name
project(app)

# include biovoltron CMakeLists.txt
include(${PROJECT_SOURCE_DIR}/Biovoltron/CMakeLists.txt)

# link your source code to your project
add_executable(app src/main.cpp)
target_link_libraries(app biovoltron)
```

The project layout should look like below by now:
```bash
$ tree app
app
├── CMakeLists.txt
├── Biovoltron
│   ├── CMakeLists.txt
│   ├── include
│   ├── submodules
│   └── ...
├── build
└── src
    └── main.cpp
```

Now we can compile and run the code.
```bash
cd build
cmake ..
make app -j
./app
```

The output should look like this:
```
hello-world
ACTG
```

Please follow our [Tutorial](tutorial.md) for more advanced topics!
