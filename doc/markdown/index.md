Biovoltron {#mainpage}
=========

Bioinformatic tools library using C++20.

Features
--------

* Use C++20.
* Fast (see [Benchmark](benchmark.md)).

Biovoltron at a glance
---------------------
- File io
```cpp
#include <biovoltron/file_io/fasta.hpp>
#include <fstream>
#include <iostream>
#include <ranges>

using namespace biovoltron;

int main() {
  auto fa = std::ifstream{"hs37d5.fa"};
  auto refs = std::ranges::istream_view<FastaRecord<>>(fa);
  for (const auto& ref : refs)
    std::cout << ref << "\n";
}
```
- Applications
```cpp
#include <biovoltron/pipe/algo_pipe.hpp>
#include <iostream>
#include <range/v3/all.hpp>

using namespace biovoltron;

int main() {
  auto ref = FastaRecord<>{};
  std::ifstream{"chr22.fa"} >> ref;
  auto index = ref | pipe::build{.LOOKUP_LEN = 8};

  auto fq1 = std::ifstream{"chr22.1.fq"};
  auto fq2 = std::ifstream{"chr22.2.fq"};
  auto reads = ranges::view::zip(ranges::istream_range<FastqRecord<>>(fq1),
                                 ranges::istream_range<FastqRecord<>>(fq2));

  const auto vars
    = reads
      | pipe::align{.ref = ref, .index = std::move(index)}
      | ranges::action::sort
      | pipe::call{.ref = ref};
  for (const auto& var : vars)
    std::cout << var << "\n";
}
```

Dependencies
------------

* g++ 11.1.0+
* boost 1.74.0+
* zlib 1.2.11+
* tbb 2020.3+
* range-v3 0.11.0+
* cmake 3.16+ (optional)

Installation
------------

We recommend that you use CMake to build your project.
See [Installation](installation.md).

Contribute
----------

See [Contribute](contribute.md).

License
-------

Biovoltron is licensed under the MIT license.
