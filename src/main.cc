#include <iostream>
#include "../contrib/htslibpp/htslibpp.h"
#include "../contrib/htslibpp/htslibpp_variant.h"

using namespace YiCppLib::HTSLibpp;

int main(int argc, const char* argv[]) {
    auto htsInHandle = htsOpen("-", "r");
    auto htsOutHandle = htsOpen("-", "w");

    auto header = htsHeader<bcfHeader>::read(htsInHandle);

    if(header.get() == nullptr) {
        std::cerr<<"unable to read header from input stream"<<std::endl;
        exit(1);
    }

    bcf_hdr_write(htsOutHandle.get(), header.get());

    for_each(begin(htsInHandle, header), end(htsInHandle, header), [&htsOutHandle, &header](auto &bcfRecord) {
        bcf_write(htsOutHandle.get(), header.get(), bcfRecord.get());
    });


    return 0;
}
