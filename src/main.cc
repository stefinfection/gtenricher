#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include "../contrib/htslibpp/htslibpp.h"
#include "../contrib/htslibpp/htslibpp_variant.h"

using namespace YiCppLib::HTSLibpp;

int main(int argc, const char* argv[]) {

    if(argc < 2) {
        std::cout<<"no subset defined on the command line. exiting...";
        exit(1);
    }

    // input and output file handles
    auto htsInHandle = htsOpen("-", "r");
    auto htsOutHandle = htsOpen("-", "w");

    auto header = htsHeader<bcfHeader>::read(htsInHandle);

    // if header is not present, fail and quit
    if(header.get() == nullptr) {
        std::cerr<<"unable to read header from input stream"<<std::endl;
        exit(1);
    }

    // get sample names
    std::vector<const std::string> sampleNames;
    sampleNames.reserve(bcf_hdr_nsamples(header.get()));
    std::transform( 
        htsHeader<bcfHeader>::dictBegin(header, htsHeader<bcfHeader>::DictType::SAMPLE), 
        htsHeader<bcfHeader>::dictEnd(header, htsHeader<bcfHeader>::DictType::SAMPLE),
        std::back_inserter(sampleNames),
        [](const auto &sampleRec){ const auto p = htsProxy(sampleRec); return p.key(); }
    );

    std::cerr<<"proband: ";
    for_each(sampleNames.cbegin(), sampleNames.cend(), [](auto& sample) { std::cerr<<sample<<" "; });
    std::cerr<<std::endl;

    const std::set<const std::string> subset(argv+1, argv+argc);

    std::cerr<<"subset: ";
    for_each(subset.cbegin(), subset.cend(), [](auto& sample) { std::cerr<<sample<<" "; });
    std::cerr<<std::endl;

    auto outHeader = bcfHeader{bcf_hdr_subset(header.get(), 0, nullptr, nullptr)};
    bcf_hdr_append(outHeader.get(), "##INFO=<ID=ENRICH,Number=8,Type=Integer,Description=\"Enrichment of affected genotypes in subset\">");
    bcf_hdr_write(htsOutHandle.get(), outHeader.get());

    
    std::for_each(begin(htsInHandle, header), end(htsInHandle, header), [&](auto &rec) {
        int ngt, *gt_arr = NULL, ngt_arr = 0;
        ngt = bcf_get_genotypes(header.get(), rec.get(), &gt_arr, &ngt_arr);

        std::vector<int32_t> counts(8, 0); // IDX 0-3: proband; IDX4-7: subset

        for(int sample=0; sample<sampleNames.size(); sample++) {

            int key = ((gt_arr[sample * 2] >> 1) - 1) + ((gt_arr[sample * 2 + 1] >> 1) - 1);
            if(key < 0) key = 3;

            ++counts[key];
            if(subset.find(sampleNames[sample]) != subset.end()) ++counts[key+4];
        }

        bcf_update_info_int32(outHeader.get(), rec.get(), "ENRICH", &counts[0], 8);
        bcf_subset(outHeader.get(), rec.get(), 0, nullptr);
        bcf_write(htsOutHandle.get(), outHeader.get(), rec.get());

        free(gt_arr);
    });

    return 0;
}
