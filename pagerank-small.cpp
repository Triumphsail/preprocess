#include <iostream>
#include <string.h>
#include <vector>
#include <x86intrin.h>
#include <chrono>
#include <unistd.h>
#include <sys/types.h>
#include "Graph/BigGraph-new.h"
#include "Rdma/IBController.h"
#include "Debug/Debug.h"
#include "Rdma/Buffer.h"
#include "Graph/atomic.h"
#include "Graph/BigVector.h"
#define IOSIZE 1048576 * 24
#define PAGESIZE 4096
uint64_t getBeginIdByPartitionId(BigGraph &graph, uint32_t id) {
    uint64_t sourceBeginId = 0;
    for (uint32_t i = 0; i < id; i++) {
        sourceBeginId += graph.getVertexNumbersByPartitionId(i);
    }
    return sourceBeginId;
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        printf("pagerank [IB Config Path] [Graph Path] [iterations] [memery budget in GB(default 8GB)]\n");
        return 0;
    }
    std::string ibConfigPath = argv[1];
    std::string GraphPath = argv[2];
    int iterations = atoi(argv[3]);
    long memeryBytes = (argc >= 5) ? atol(argv[4]) * 1024l * 1024l * 1024l : 8l * 1024l * 1024l * 1024l;
    IBController ibController(ibConfigPath);
    uint32_t localMachineId = ibController.getLocalMachineId();
    uint32_t machineNum = ibController.getMachineNum();
    bool inMemory = false;
    BigGraph graph(GraphPath, localMachineId);
    graph.setMemoryBytes(memeryBytes);
    uint32_t partitions = graph.getPartitions();
    if(partitions == machineNum) {
        inMemory = true;
    }
    uint32_t partitionSize = partitions / machineNum;
    uint64_t localVertexSize[partitionSize];
    uint64_t vertexSize = 0;
    uint32_t maxVertexSize = 0;
    
    VertexId vertices = graph.getVertices();
    if(vertices % partitions == 0) {
        maxVertexSize = vertices / partitions;
    } else {
        maxVertexSize = vertices / partitions + 1;
    }
    
    for (uint32_t i = 0; i < partitionSize; i++) {
        localVertexSize[i] = graph.getVertexNumbersByPartitionId(localMachineId * partitionSize + i);
        vertexSize += localVertexSize[i];
    }
    std::cout << "vertexSize: " << maxVertexSize << std::endl;
    // BigVector<VertexId> degree(GraphPath + "/degree", vertexSize);
    std::vector<VertexId> degree;
    degree.reserve(vertexSize);
    std::string path = GraphPath + "/vertex-" + std::to_string(localMachineId);
    int degreeFile = open(path.c_str(), O_RDONLY);
    char* buffer = (char *)memalign(PAGESIZE, IOSIZE);
    MY_ASSERT(buffer != NULL, "[PAGERANK] Cannot memalign memory.");
    int bytes = -1;
    uint32_t vId = 0;
    while(true) {
        bytes = read(degreeFile, buffer, IOSIZE);
        MY_ASSERT(bytes != -1, "[PAGERANK] Read Vertex File Failed.\n");
        if(bytes == 0) break;
        for (int pos = 0; pos < bytes; pos += 2 * sizeof(VertexId)) {
            // VertexId vertexId = *(VertexId*) (buffer + pos);
            VertexId d = *(VertexId*) (buffer + pos + sizeof(VertexId));
            if(vId < vertexSize){
                degree[vId] = d;
                vId++;
            }
        }
    }

    uint64_t sourceBeginId = getBeginIdByPartitionId(graph, localMachineId * partitionSize);
    std::vector<float*> sumByMachineId;
    std::vector<float*> recvBufferByMachineId;
    uint32_t sendNum = 0;
    for (uint32_t id = 0; id < machineNum; id++) {
        VertexId vertexNum = maxVertexSize;
        float *sum = (float*)malloc(vertexNum * sizeof(float));
        memset(sum, 0, sizeof(float) * (vertexNum));
        sumByMachineId.push_back(sum);
        float *recv = NULL;
        if (id != localMachineId) {
            recv = (float*)malloc((maxVertexSize) * sizeof(float));
            memset(recv, 0, sizeof(float) * (maxVertexSize));
        }
        recvBufferByMachineId.push_back(recv);
    }

    for (uint32_t id = 0; id < machineNum; id++) {
        if (localMachineId == id) {
            continue;
        } else {
            ibController.regionBuffer((void*)sumByMachineId[id], (maxVertexSize) * sizeof(float),
                                      (void*)recvBufferByMachineId[id], (maxVertexSize) * sizeof(float), id);
        }
    }

    ibController.useSocketToConnectIB();
    std::cout << "connect ib success" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();

    // BigVector<float> pagerank(GraphPath + "/pagerank", vertexSize);
    //pagerank.fill(0);
    float *pagerank = new float[localVertexSize];
    memset(pagerank, 0, sizeof(float) * (localVertexSize));

    for(uint32_t i = 0; i < partitionSize; i++) {
        graph.streamVertices(
            [&](){
                std::cout << "PageRank Init." << std::endl;
            },
            [&](VertexId id) {
                if (degree[id] != 0) {
                    pagerank[id] = 1.f / degree[id];
                }
                return 0;
            }, i, 0,
            [&](std::pair<VertexId, VertexId> range) {

            },
            [&](std::pair<VertexId, VertexId> range) {

            }
        );
    }
    auto degreeEnd = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = degreeEnd - start;
    std::cout << "Init Success, Time: " << diff.count() << std::endl;
    std::vector<bool> ifActive(machineNum, true);
    std::vector<int> iterByMachineId(machineNum, 0);
    std::vector<VertexId> vertexOffset(machineNum+1, 0);
    for (uint32_t i = 1; i <= partitionSize; i++) {
        vertexOffset[i] = vertexOffset[i - 1] + localVertexSize[i - 1];
    }
    for (int iter = 0; iter < iterations; iter++) {
        std::cout << "-----------" << "iterations: " << iter + 1 << "-----------" << std::endl;
        uint32_t vertexPatitionNum = partitionSize;
        uint32_t bigPartitionNum = vertexPatitionNum * vertexPatitionNum;
        for(uint32_t BPID = 0; BPID < bigPartitionNum; BPID++) {
            uint32_t vertexPartitionId = BPID / vertexPatitionNum;
            std::cout << "part: " << BPID << ", vertex part id: " << vertexPartitionId << std::endl;
            graph.streamEdges(
                [&](Edge &e, uint32_t machineId) {
                    uint64_t targetBeginId = getBeginIdByPartitionId(graph, machineId * partitionSize + (BPID % partitionSize));
                    VertexId target = e.target - targetBeginId;
                    VertexId source = e.source - sourceBeginId;
                    write_add(&sumByMachineId[machineId][target], pagerank[source]);
                    // std::cout << "source: " << e.source << ", target: " << e.target <<", sumByMachineId: " << sumByMachineId[machineId][target] << std::endl;
                    return 0;
                },
                [&](uint32_t remoteMachineId){
                    if (remoteMachineId != localMachineId) {
                        ibController.RdmaWriteWithImm(remoteMachineId, sizeof(float) * (maxVertexSize), localMachineId, 0, 0);
                    }
                }, ifActive, BPID, vertexPartitionId, 0,
                [&](std::pair<VertexId, VertexId> range){
                    
                },
                [&](std::pair<VertexId, VertexId> range){
                    
                });
            std::cout << "Start Vertex" << std::endl;
            auto startVertex = std::chrono::high_resolution_clock::now();
            if(iter == iterations -1) {
                graph.streamVertices(
                    [&]() {
                        iterByMachineId[localMachineId]++;
                        sendNum += ibController.ifRdmaReceive(iter * bigPartitionNum + BPID + 1, iterByMachineId, ifActive);
                        std::cout << "Calculation PageRank." << std::endl;
                    },
                    [&](VertexId id) {
                        VertexId localId = id - vertexOffset[vertexPartitionId];
                        for(uint32_t i = 0; i < machineNum; i++) {
                            if(i != localMachineId) {
                                sumByMachineId[localMachineId][localId] += recvBufferByMachineId[i][localId];
                            }
                        }
                        pagerank[id] = 0.15f + 0.85f * sumByMachineId[localMachineId][localMachineId];
                        return 0;
                    }, vertexPartitionId, 0,
                    [&](std::pair<VertexId, VertexId> range) {
                    },
                    [&](std::pair<VertexId, VertexId> range) {
                    }
                );
            } else {
                graph.streamVertices(
                    [&]() {
                        iterByMachineId[localMachineId]++;
                        sendNum += ibController.ifRdmaReceive(iter * bigPartitionNum + BPID + 1, iterByMachineId, ifActive);
                        std::cout << "Calculation PageRank." << std::endl;
                    },
                    [&](VertexId id) {
                        VertexId localId = id - vertexOffset[vertexPartitionId];
                        for(uint32_t i = 0; i < machineNum; i++) {
                            if(i != localMachineId) {
                                sumByMachineId[localMachineId][localId] += recvBufferByMachineId[i][localId];
                            }
                        }
                        pagerank[id] = (0.15f + 0.85f * sumByMachineId[localMachineId][localMachineId]) / degree[id];
                        return 0;
                    }, vertexPartitionId, 0,
                    [&](std::pair<VertexId, VertexId> range) {

                    },
                    [&](std::pair<VertexId, VertexId> range) {

                    }
                );
            }
            for (uint32_t i = 0; i < machineNum; i++) {
                memset(sumByMachineId[i], 0, sizeof(float) * (graph.getVertexNumbersByPartitionId(i)));
            }
            auto endVertex = std::chrono::high_resolution_clock::now();
            diff = endVertex - startVertex;
            std::cout <<"use "<< diff.count() <<" s\n";
        }
    }

    std::cout << "-----------result-----------" << std::endl;
    std::string pagerankFile = GraphPath + "/PageRank" + std::to_string(localMachineId);

    FILE * file = fopen(pagerankFile.c_str(), "wb");
	fclose(file);
    int fout = open(pagerankFile.c_str(), O_WRONLY);
    long bytes = write(fout, pagerank, localVertexSize * sizeof(float));
    std::cout << "PageRank Write Success, bytes: " << bytes << std::endl;
    // double endTime = getWallTime();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout <<"use "<< diff.count() <<" s\n";
    // std::cout << "time: " << endTime - startTime << "ms" << std::endl;
    delete pagerank;

    auto end = std::chrono::high_resolution_clock::now();
    diff = end - start;
    std::cout <<"use "<< diff.count() <<" s\n";
    sumByMachineId.clear();
    recvBufferByMachineId.clear();
    while (sendNum < iterations * (machineNum - 1)) {
        sendNum += ibController.ifRdmaReceive(0, iterByMachineId, ifActive);
    }
    return 0;
}