#include <iostream>
#include <vector>
#include <set>
#include <algorithm>
#include "../contrib/htslibpp/htslibpp.h"
#include "../contrib/htslibpp/htslibpp_variant.h"

using namespace YiCppLib::HTSLibpp;

int main(int argc, const char* argv[]) {

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
    std::vector<std::string> sampleNames;
    std::transform( 
        htsHeader<bcfHeader>::dictBegin(header, htsHeader<bcfHeader>::DictType::SAMPLE), 
        htsHeader<bcfHeader>::dictEnd(header, htsHeader<bcfHeader>::DictType::SAMPLE),
        std::back_inserter(sampleNames),
        [](auto &sampleRec){ auto p = htsProxy(sampleRec); return p.key(); }
    );

    std::set<std::string> subset{"11200X11", "11200X12", "11200X9"};

    bcf_hdr_write(htsOutHandle.get(), header.get());
    std::for_each(begin(htsInHandle, header), end(htsInHandle, header), [&htsOutHandle, &header, &sampleNames, &subset](auto &bcfRecord) {
        int ngt, *gt_arr = NULL, ngt_arr = 0;
        ngt = bcf_get_genotypes(header.get(), bcfRecord.get(), &gt_arr, &ngt_arr);

        std::vector<int> proband_counts(4, 0);
        std::vector<int> subset_counts(4, 0);

        for(int sample=0; sample<sampleNames.size(); sample++) {

            int key = ((gt_arr[sample * 2] >> 1) - 1) + ((gt_arr[sample * 2 + 1] >> 1) - 1);
            std::cout<<"key: "<<key<<std::endl;
            if(key < 0) key = 3;

            ++proband_counts[key];
            if(subset.find(sampleNames[sample]) != subset.end()) ++subset_counts[key];
        }

        for(int i=0; i<ngt; i+=2) {
            std::cout<<(gt_arr[i]>>1) - 1<<"/"<<(gt_arr[i+1] >> 1) - 1<<"\t";
        }
        std::cout<<std::endl;
        std::cout<<"proband: ";
        for(int i=0; i<4; i++) std::cout<<proband_counts[i]<<"\t";
        std::cout<<std::endl;
        std::cout<<"subset : ";
        for(int i=0; i<4; i++) std::cout<<subset_counts[i]<<"\t";
        std::cout<<std::endl;
    });

    

    return 0;
}
