#include <iostream>
#include <string.h>
#include <x86intrin.h>
#include <unordered_map>
#include <vector>
#include <chrono>
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
        printf("wcc [IB Config Path] [Graph Path] [memery budget in GB(default 8GB)]\n");
        return 0;
    }
    std::string ibConfigPath = argv[1];
    std::string GraphPath = argv[2];
    long memeryBytes = (argc >= 4) ? atol(argv[3]) * 1024l * 1024l * 1024l : 8l * 1024l * 1024l * 1024l;
    IBController ibController(ibConfigPath);
    uint32_t localMachineId = ibController.getLocalMachineId();
    uint32_t machineNum = ibController.getMachineNum();
    
    BigGraph graph(GraphPath, localMachineId);
    graph.setMemoryBytes(memeryBytes);
    uint32_t partitions = graph.getPartitions();
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
    uint64_t sourceBeginId = getBeginIdByPartitionId(graph, localMachineId * partitionSize);

    std::vector<VertexId*> labelByMachineId;
    std::vector<VertexId*> recvBufferByMachineId;
    for (uint32_t id = 0; id < machineNum; id++) {
        VertexId *label = (VertexId*)malloc(maxVertexSize * sizeof(VertexId));
        for (uint32_t i = 0; i < maxVertexSize; i++) {
            label[i] = vertices;
        }
        labelByMachineId.push_back(label);
        VertexId *recv = NULL;
        if (id != localMachineId) {
            recv = (VertexId*)malloc((maxVertexSize) * sizeof(VertexId));
            memset(recv, 0, sizeof(VertexId) * (maxVertexSize));
        }
        recvBufferByMachineId.push_back(recv);
    }

    for (uint32_t id = 0; id < machineNum; id++) {
        if (localMachineId == id) {
            continue;
        } else {
            ibController.regionBuffer((void*)labelByMachineId[id], maxVertexSize * sizeof(VertexId),
                                      (void*)recvBufferByMachineId[id], maxVertexSize * sizeof(VertexId), id);
        }
    }
    
    ibController.useSocketToConnectIB();
    std::cout << "connect ib success" << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    // BigVector<VertexId> localLabel(GraphPath + "/wcc", vertexSize);
    // BigVector<bool> activeVertex(GraphPath + "/active", vertexSize);
    // activeVertex.fill(1);
    std::vector<VertexId> localLabel;
    degree.reserve(vertexSize);
    std::vector<bool> activeVertex(vertexSize, 0); 
    VertexId activeVertices = 0;
    for(uint32_t i = 0; i < partitionSize; i++) {
        activeVertices += graph.streamVertices(
            [&](){
                std::cout << "WCC Init." << std::endl;
            },
            [&](VertexId id) {
                localLabel[id] = id + sourceBeginId;
                // std::cout << "id: " << id << ", label: " << labelByMachineId[localMachineId][id] << std::endl;
                return 1;
        }, i, 0,
            [&](std::pair<VertexId, VertexId> range) {
                // localLabel.load(range.first, range.second);
            },
            [&](std::pair<VertexId, VertexId> range) {
                // localLabel.save();
            }
        );
    }

    std::cout << "Init sucess" << std::endl;
    std::vector<bool> ifActive(machineNum, true);
    std::vector<int> iterByMachineId(machineNum, 0);
    int iter = 0;
    uint32_t sendNum = 0;
    std::vector<VertexId> vertexOffset(machineNum+1, 0);
    for (uint32_t i = 1; i <= partitionSize; i++) {
        vertexOffset[i] = vertexOffset[i - 1] + localVertexSize[i - 1];
    }
    while(true) {
        std::cout << "--------------iter: " << iter << "--------------" << std::endl;
        activeVertices = 0;
        uint32_t vertexPatitionNum = partitionSize;
        uint32_t bigPartitionNum = vertexPatitionNum * vertexPatitionNum;
        for (uint32_t BPID = 0; BPID < bigPartitionNum; BPID++) {
            iter++;
            iterByMachineId[localMachineId]++;
            uint32_t vertexPartitionId = BPID / vertexPatitionNum;
            std::cout << "edge start" << std::endl;
            graph.streamEdges(
                [&](Edge &e, uint32_t machineId) {
                    VertexId source = e.source - sourceBeginId;
                    if(activeVertex[source]) {
                        uint64_t targetBeginId = getBeginIdByPartitionId(graph, machineId * partitionSize + (BPID % partitionSize));
                        VertexId target = e.target - targetBeginId;
                        write_min(&labelByMachineId[machineId][target], localLabel[source]);
                    }
                    // std::cout << "source: " << e.source << ", target: " << e.target <<", labelByMachineId: " << labelByMachineId[machineId][target] << std::endl;
                    return 0;
                },
                [&](uint32_t remoteMachineId){
                    if (remoteMachineId != localMachineId) {
                        if (ifActive[localMachineId]) {
                            ibController.RdmaWriteWithImm(remoteMachineId, sizeof(VertexId) * (maxVertexSize), localMachineId, 0, 0);
                        } else {
                            ibController.RdmaWriteWithImm(remoteMachineId, 0, localMachineId + machineNum, 0, 0);
                        }
                    }
                }, ifActive, BPID, vertexPartitionId, 0,
                [&](std::pair<VertexId, VertexId> range){
                    // std::cout << "Lock Active Vertex" << std::endl;
                    // localLabel.lock(range.first, range.second);
                    // activeVertex.lock(range.first, range.second);
                },
                [&](std::pair<VertexId, VertexId> range){
                    // std::cout << "Unlock Active Vertex" << std::endl;
                    // localLabel.unlock(range.first, range.second);
                    // activeVertex.unlock(range.first, range.second);
                });
            std::cout << "edge end" << std::endl;
            std::cout << "vertex start" << std::endl;

            auto startVertex = std::chrono::high_resolution_clock::now();
            activeVertices += graph.streamVertices(
            [&]() {
                sendNum += ibController.ifRdmaReceive(iter, iterByMachineId, ifActive);
                std::cout << sendNum << std::endl;
            },
            [&](VertexId id) {
                VertexId flag = 0;
                VertexId localId = id - vertexOffset[vertexPartitionId];
                for (uint32_t i = 0; i < machineNum; i++) {
                    if (i != localMachineId) {
                        if (localLabel[id] > recvBufferByMachineId[i][localId]) {
                            localLabel[id] = recvBufferByMachineId[i][localId];
                            flag = 1;
                            activeVertex[id] = true;
                        } else {
                            activeVertex[id] = false;
                        }
                    } else {
                        if (localLabel[id] > labelByMachineId[i][localId]) {
                            localLabel[id] = labelByMachineId[i][localId];
                            activeVertex[id] = true;
                            flag = 1;
                        } else {
                            activeVertex[id] = false;
                        }
                    }
                }
                return flag;
            }, vertexPartitionId, 0,
                [&](std::pair<VertexId, VertexId> range) {
                    // localLabel.load(range.first, range.second);
                    // activeVertex.load(range.first, range.second);
                },
                [&](std::pair<VertexId, VertexId> range) {
                    // localLabel.save();
                    // activeVertex.save();
                }
            );

            std::cout << "vertex end" << std::endl;
            auto endVertex = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> diff = endVertex - startVertex;
            std::cout <<"use "<< diff.count() <<" s\n";
        }

        if (activeVertices == 0) {
            ifActive[localMachineId] = false;
            for (uint32_t remoteMachineId = 0; remoteMachineId < machineNum; remoteMachineId++) {
                if(remoteMachineId == localMachineId) continue;
                ibController.RdmaWriteWithImm(remoteMachineId, 0, localMachineId + machineNum, 0, 0);
            }
        } else {
            ifActive[localMachineId] = true;
            for (uint32_t remoteMachineId = 0; remoteMachineId < machineNum; remoteMachineId++) {
                if(remoteMachineId == localMachineId) continue;
                ibController.RdmaWriteWithImm(remoteMachineId, 0, localMachineId, 0, 0);
            }
        }
        iter++;
        iterByMachineId[localMachineId]++;
        sendNum += ibController.ifRdmaReceive(iter, iterByMachineId, ifActive);
        std::cout << sendNum << std::endl;
        uint32_t numbers = 0;
        for (uint32_t id = 0; id < machineNum; id++) {
            if(ifActive[id]) {
                break;
            } else {
                numbers++;
            }
        }
        if (numbers == machineNum) {
            break;
        }
        
    }
    iterByMachineId[localMachineId]++;
    std::cout << "-----------result-----------" << std::endl;
    BigVector<VertexId> labelStat(GraphPath + "/labelStat", vertices);
    labelStat.fill(0);
    for(uint32_t i = 0; i < partitionSize; i++) {
        graph.streamVertices(
            [&](){
                std::cout << "WCC." << std::endl;
            },
            [&](VertexId id) {
                labelStat[localLabel[id]]++;
                return 0;
        }, i, 0,
            [&](std::pair<VertexId, VertexId> range) {
            },
            [&](std::pair<VertexId, VertexId> range) {
            }
        );
    }
    VertexId localComponents = 0;
    for(uint32_t i = 0; i < partitions; i++) {
        localComponents += graph.streamVertices(
            [&](){
                std::cout << "WCC." << std::endl;
            },
            [&](VertexId id) {
                return labelStat[id] != 0;
        }, i, 0,
            [&](std::pair<VertexId, VertexId> range) {
            },
            [&](std::pair<VertexId, VertexId> range) {
            }
        );
    }
    std::cout << "Local Components is " << localComponents << std::endl;
    auto end = std::chrono::high_resolution_clock::now();
    
    std::chrono::duration<double> diff = end - start;
    std::cout <<"use "<< diff.count() <<" s\n";
    // std::cout << "time: " << endTime - startTime << "ms" << std::endl;


    return 0;
}