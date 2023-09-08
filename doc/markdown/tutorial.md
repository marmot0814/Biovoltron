Tutorial {#tutorial}
========

This document provides example recipes on how to carry out particular tasks using **Biovoltron**.

## Read sequence files
- Combined with `std::ranges`, biovoltron enable users to easily input common bioinformation format.

```cpp=
#include <iostream>
#include <sstream>
#include <ranges>
#include <biovoltron/file_io/fasta.hpp>

using namespace biovoltron;

auto fasta = R"(>chr1
ACGT
>chr1
AGGCTGA
>chr1
GGAGTATAATATATATATATATAT)";

int main() {
    auto input = std::istringstream(fasta);
    for (auto& r : std::ranges::istream_view<FastaRecord<true>>(input)) {
        std::cout << "name: " << r.name << '\n';
        std::cout << "seq : " << r.seq << '\n';
    }
    return 0;
}
```

## Construction and assignment of biovoltron::istring symbols
- The design of `biovoltron::istring` makes dna/rna string convert to numeric/bit representation easily.

```cpp=
#include <cassert>
#include <iostream>
#include <biovoltron/utility/istring.hpp>

int main() {
    using namespace std::string_literals;
    using namespace biovoltron;

    const auto s = 12031_s;
    assert(Codec::hash(s) == 0b01'10'00'11'01);
    assert(Codec::rhash(0b01'10'00'11'01, 5) == s);
    assert(Codec::to_istring("ACNGGTT") == 0142233_s);
    assert(Codec::rev_comp("ACNGGTT") == "AACCNGT"s);
    assert(Codec::rev_comp(s) == 20312_s);
    std::cout << s << '\n';
    return 0;
}
```

## Compressed vector(DibitVector) along with std::ranges
- Not only `biovoltron::istring`, `biovoltron::DibitVector` can also paired with std::ranges easily.

```cpp=
#include <algorithm>
#include <iostream>
#include <ranges>
#include <biovoltron/container/xbit_vector.hpp>

int main() {
    auto v = biovoltron::DibitVector<>{3, 2, 1, 2, 3, 0, 0, 1, 2};
    assert(v.size() == 9);
    assert(v.num_blocks() == 3);

    auto base_view = std::views::transform([](auto c){ return "ACGT"[c]; });
    auto comp_view = std::views::transform([](auto c){ return 0b11u - c; });

    std::cout << "Reverse complement of ";
    std::ranges::copy(v | base_view, std::ostream_iterator<char>{std::cout, ""});
    std::cout << " is ";
    for (auto c : v | std::views::reverse | comp_view | base_view)
        std::cout << c;

    // result:
    // Reverse complement of TGCGTAACG is CGTTACGCA
}
```

## Serializer
- `biovoltron::Serializer` can archieve all kinds of random access container(satisfied by the std::ranges::random_access_range concept), including `biovoltron::xbit_vector`.

```cpp=
#include <string>
#include <biovoltron/container/xbit_vector.hpp>
#include <biovoltron/utility/archive/serializer.hpp>

template <typename R>
void test(const R& r, const char* path) {
    using biovoltron::Serializer;
    {
        auto fout = std::ofstream{path, std::ios::binary};
        Serializer::save(fout, r);
    }
    auto r2 = R{};
    {
        auto fin = std::ifstream{path, std::ios::binary};
        Serializer::load(fin, r2);
    }
    assert(r == r2);
}

struct Point {
    double x;
    double y;
    auto operator<=>(const Point&) const = default;
};

int main() {
    test(std::vector{true, false, false}, "vector_bool.bin");
    test(std::vector{7, 5, 16, 8}, "vector_int.bin");
    test(std::string{"Exemplar"}, "string.bin");
    test(std::vector{Point{0.0, 0.0}, {1.0, 2.0}, {3.0, 4.0}}, "vector_point.bin");
    test(biovoltron::DibitVector<>{1, 0, 2, 3, 3, 0, 2}, "DibitVector.bin");
    return 0;
}
```

## FMIndex
The memory usage for FMIndex in run time is affect by above parameters. Take a 3.1Gb human genome for example, the memory occupation can be calculate as follwing:
- bwt string:`3.1Gb / 4 = 0.775` Gb.
- hierarchical occurrence table:
    - L1 occ occupy fixed `3.1Gb * 16 / 256 = 0.194Gb` plus
    - L2 occ occupy `3.1Gb * 4 / occ_intv(16) = 0.775 Gb`.
- suffix array: `3.1Gb * 4 / sa_intv(1) = 12Gb`.
    - Noticed that in default mode we dont sampling suffix value since this can reduce frequently memory allocation and intense computation when occurs massive query.
- lookup table: `4^lookup_len(14) * 4 / 1024^3 = 1Gb`.

```cpp=
#include <iostream>
#include <ranges>
#include <biovoltron/algo/align/exact_match/fm_index.hpp>
#include <biovoltron/file_io/fasta.hpp>

using namespace biovoltron;

int main(int argc, char** argv) {
  if (argc != 2) {
    std::cout << "argv[1] should be your fasta file.\n";
    return 0;
  }

  auto fname = std::string(argv[1]);
  {
    // input
    auto ref = istring{};
    auto fin = std::ifstream(fname);
    for (const auto& r : std::ranges::istream_view<FastaRecord<true>>(fin))
      ref += r.seq;

    // fm-index only support 'ACTG', so we need to change 'N' to 'ACGT'
    std::ranges::transform(ref, ref.begin(), [](auto& c) { return c % 4; });

    // build
    auto fmi = FMIndex{};
    fmi.build(ref);

    // save
    auto fout = std::ofstream(fname + ".fmi", std::ios::binary);
    fmi.save(fout);
    std::cout << fname + ".fmi saved.\n";
  }

  // load
  auto fmi = FMIndex{};
  {
    auto fin = std::ifstream(fname + ".fmi", std::ios::binary);
    fmi.load(fin);
    std::cout << fname + ".fmi loaded.\n";
  }

  // query
  for (auto seed = istring{}; std::cin >> seed;) {
    const auto [beg, end, offset] = fmi.get_range(seed, 0);
    std::cout << "seed: " << seed << "\n";
    std::cout << "seed offset: " << offset << "\n";
    std::cout << "occurrence: " << end - beg << "\n";
    const auto offsets = fmi.get_offsets(beg, end);
    std::cout << "ref offsets: ";
    for (const auto offset : offsets) std::cout << offset << " ";
    std::cout << "\n";
  }
}
```

## SmithWaterman with simd acceleration
- `Biovoltron::SmithWaterman` is modern C++ Smith-Waterman wrapper which using SSW Library as a underlying code base.

```cpp=
#include <iostream>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/algo/align/inexact_match/smithwaterman_sse.hpp>

int main() {
    using namespace biovoltron;
    const auto read = 1102111000031323330032232203332323_s;
    const auto ref  = 2131100321101000010313231313001322323213_s;
    const auto prof = SseSmithWaterman::get_profile(read);
    const auto [score, score2, ref_beg, ref_end, read_beg, read_end, ref_end2, cigar] = SseSmithWaterman::align(prof, ref);
    std::cout << "score: " << score << "\n";
    std::cout << "ref_beg: " << ref_beg << "\n";
    std::cout << "ref_end: " << ref_end << "\n";
    std::cout << "read_beg: " << read_beg << "\n";
    std::cout << "read_end: " << read_end << "\n";
    std::cout << "cigar: " << cigar << "\n";
    return 0;
}
```

## Annotator
```cpp=
#include <cassert>
#include <biovoltron/algo/annotate/annotator.hpp>

using namespace biovoltron;

int main() {
    auto genes = Annotator<std::string>{};
    genes.insert_at("gene1", Interval{"chr1:5-15"});
    genes.insert_at("gene2", Interval{"chr1:2-10"});
    genes.insert_at("gene3", Interval{"chr1:20-30"});
    genes.insert_at("gene4", Interval{"-chr2:2-10"});
    genes.index();

    // result returns a vector of object, sorted by begin
    auto results = genes.find(Interval{"chr1:6-12"});
    assert(results.size() == 2);
    assert(results[0] == "gene2");
    assert(results[1] == "gene1");

    results = genes.find(Interval{"-chr2:6-9"});
    assert(results.size() == 1);
    assert(results[0] == "gene4");

    results = genes.find(Interval{"chr2:6-9"});
    assert(results.empty());

    results = genes.find(Interval{"chr2:0-2"});
    assert(results.empty());

    results = genes.find(Interval{"chr999:5-6"});
    assert(results.empty());
    return 0;
}
```

## Variant calling
```cpp=
#include <biovoltron/pipe/algo_pipe.hpp>
#include <iostream>
#include <range/v3/all.hpp>
#include <biovoltron/algo/sort/stable_sorter.hpp>

using namespace biovoltron;

int main() {
    auto ref = FastaRecord<>{};
    std::ifstream{"chr22.fa"} >> ref;
    auto index = ref | pipe::build<uint32_t, StableSorter<uint32_t>>{.LOOKUP_LEN = 8};

    auto fq1 = std::ifstream{"chr22.1.fq"};
    auto fq2 = std::ifstream{"chr22.2.fq"};
    auto reads = ranges::view::zip(ranges::istream_range<FastqRecord<>>(fq1),
        ranges::istream_range<FastqRecord<>>(fq2));

    const auto vars = reads
        | pipe::align{.ref = ref, .index = std::move(index)}
        | ranges::action::sort
        | pipe::call{.ref = ref};
    for (const auto& var : vars)
        std::cout << var << "\n";
    return 0;
}
```
