#include <cuda_runtime.h>
#include <cuda.h>

#include <iostream>
#include <assert.h>
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

__global__ void test_gen_detid(DetId* id) {
    DetId did;
    *id = did;
}

__global__ void test_gen_ecal_detid(EEDetId *id) {
    EEDetId did(30, 50, 1);
    *id = did;

    // trigger functions on the device

    did.ix();    
    did.iy();    
    did.zside(); 
    did.subdet();
//     did.sc(); //  ---> exception??????
    //     did.isc();   -->  error: identifier "EEDetId::isc const" is undefined in device code
    did.ic();
    did.iquadrant();
    did.positiveZ();
    did.iPhiOuterRing();
}

void test_detid() {
    // test det ids
    DetId h_id, h_id_test;
    DetId h_test0{1};
    DetId *d_id;

    cudaMalloc((void**)&d_id, sizeof(DetId));
    cudaMemcpy(d_id, &h_id, sizeof(DetId), cudaMemcpyHostToDevice);
    test_gen_detid<<<1,1>>>(d_id);
    cudaMemcpy(&h_id_test, d_id, sizeof(DetId), cudaMemcpyDeviceToHost);
    
    assert(h_id_test == h_id);
    assert(h_id != h_test0);
}

void test_ecal_detid() {
    EEDetId h_id;
    EEDetId h_id_test0{30, 50, 1};
    EEDetId *d_id;

    cudaMalloc((void**)&d_id, sizeof(EEDetId));
    cudaMemcpy(d_id, &h_id, sizeof(EEDetId), cudaMemcpyHostToDevice);
    test_gen_ecal_detid<<<1,1>>>(d_id);
    cudaMemcpy(&h_id, d_id, sizeof(EEDetId), cudaMemcpyDeviceToHost);

    std::cout << h_id_test0 << std::endl;
    std::cout << h_id << std::endl;
    assert(h_id_test0 == h_id);
}

int main(int argc, char** argv) {
    int nDevices;
    cudaGetDeviceCount(&nDevices);
    std::cout << "nDevices = " << nDevices << std::endl;

//     // test det id functionality
    if (nDevices>0)
        test_detid();

    // test ecal det ids
    if (nDevices>0)
        test_ecal_detid();

    return 0;
}
