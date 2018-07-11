#include "config.h"

#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <tuple>
#include "log.h"
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
    std::vector<std::string> samples;
    samples.reserve(bcf_hdr_nsamples(header.get()));
    std::transform( 
        htsHeader<bcfHeader>::dictBegin(header, htsHeader<bcfHeader>::DictType::SAMPLE), 
        htsHeader<bcfHeader>::dictEnd(header, htsHeader<bcfHeader>::DictType::SAMPLE),
        std::back_inserter(samples),
        [](const auto &sampleRec){ const auto p = htsProxy(sampleRec); return p.key(); }
    );

    // get proband names from command line
    const std::set<std::string> subset(argv+1, argv+argc);
    
    logger(LOGLV_INFO)<<samples.size()<<std::endl;
    { std::stringstream strbuf; for(auto &s : samples) strbuf<<s<<" "; logger(LOGLV_INFO)<<"proband: "<<strbuf.str()<<std::endl; }
    { std::stringstream strbuf; for(auto &s : subset)  strbuf<<s<<" "; logger(LOGLV_INFO)<<"subset:  "<<strbuf.str()<<std::endl; }

    auto outHeader = bcfHeader{bcf_hdr_subset(header.get(), 0, nullptr, nullptr)};
    bcf_hdr_append(outHeader.get(), "##INFO=<ID=ENRICH,Number=8,Type=Integer,Description=\"Enrichment of affected genotypes in subset\">");
    htsOutHandle << outHeader;
    
    std::for_each(begin(htsInHandle, header), end(htsInHandle, header), [&](auto &rec) {

        if(rec->n_allele > 2) {
            logger(LOGLV_WARN)<<"skipping multi-allelic record at "<<bcf_hdr_id2name(header.get(), rec->rid)<<":"<<rec->pos+1<<std::endl;
            return;
        }

        int32_t ngt, *gt_arr = nullptr, ngt_arr = 0;
        ngt = bcf_get_genotypes(header.get(), rec.get(), &gt_arr, &ngt_arr);
        std::unique_ptr<int32_t> uniq_gt_arr(gt_arr);

        std::vector<int32_t> counts(8, 0); // IDX 0-3: proband; IDX4-7: subset

        for(int sid=0; sid<samples.size(); sid++) {

            int key = ((gt_arr[sid * 2]) >> 1) + ((gt_arr[sid* 2 + 1]) >> 1) - 2;

            assert(key < 4);

            if(key < 0) key = 3;

            ++counts[key];
            if(subset.find(samples[sid]) != subset.end()) ++counts[key+4];
        }

        bcf_update_info_int32(outHeader.get(), rec.get(), "ENRICH", &counts[0], 8);
        bcf_subset(header.get(), rec.get(), 0, nullptr);
        htsOutHandle<<bcfHdrRecPair( outHeader, rec );
    });

    return 0;
}
