#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <malloc.h>
#include <errno.h>
#include <assert.h>
#include <string.h>
#include<algorithm>
#include <string>
#include <string.h>
#include <vector>
#include <thread>
#include <tuple>
#include <mutex>
#include "Graph/Queue.h"
#include "PreProcess/fileSystem.h"
//#include "Debug.h"
typedef int VertexId;
typedef long EdgeId;
typedef float Weight;
#define IOSIZE 1048576 * 24
#define PAGESIZE 4096
//#define machineNum 4
void gridPartition(std::string input, std::string output, VertexId vertices, EdgeId edges,
				 int edgeUnit, int partitions, int edgeType, int machineNum, std::vector<int> &verticesOut);
int get_partition_id(const VertexId vertices, const int partitions, const VertexId vertex_id) {
    if (vertices % partitions==0) {
            const int partition_size = vertices / partitions;
            return vertex_id / partition_size;
    }
    const VertexId partition_size = vertices / partitions + 1;
    const VertexId split_point = vertices % partitions * partition_size;
    return (vertex_id < split_point) ? vertex_id / partition_size : (vertex_id - split_point) / (partition_size - 1) + (vertices % partitions);
}

void preProcess(std::string input, std::string output, VertexId vertices, int partitions, int edgeType, int machineNum) { 
    MY_ASSERT((partitions % machineNum == 0), "Partition = Machine Number * N.\n");
	int edgeUnit;
    EdgeId edges;
    if(edgeType == 0)
        edgeUnit = sizeof(VertexId) * 2;
    else if(edgeType == 1)
        edgeUnit = sizeof(VertexId) * 2 + sizeof(Weight);
    else{
        printf("Edge type is error.\n");
        return;
    }
    edges = file_size(input) / edgeUnit;
    std::vector<VertexId> verticesOut(vertices, 0);
	printf("Vertices: %d, Edges: %ld\n", vertices, edges);
    
	char *buffer =  (char *)memalign(PAGESIZE, IOSIZE);
    int fin = open(input.c_str(), O_RDONLY);
    MY_ASSERT(fin != -1, "[preprocess] Cannot open input file.\n");
    while(true) {
        long bytes = read(fin, buffer, IOSIZE);
        MY_ASSERT(bytes != -1, "[preprocess] Cannot read input file.\n");
        if(bytes == 0) break;
        VertexId source;
        for (long pos=0;pos<bytes;pos+=edgeUnit) {
            source = *(VertexId*)(buffer + pos);
            //target = *(VertexId*)(buffer + pos + sizeof(VertexId));
            if(source>=vertices) {
                printf("Error, source: %d\n", source);
                return;
            }
            verticesOut[source]++;
        }
    }
    close(fin);
    std::vector<std::pair<VertexId, VertexId> > tt(vertices);
    for(int i=0;i<vertices;i++) {
        tt[i] = std::make_pair(i, verticesOut[i]);
    }
    
    sort(tt.begin(), tt.end(),[](const std::pair<VertexId, VertexId> &x, const std::pair<VertexId, VertexId> &y){return x.second < y.second;});
    VertexId *gridOffset = new VertexId [partitions];
	VertexId *tempOffset = new VertexId [partitions];
    const int partitionSplit = vertices % partitions;
    const int partitionSize = vertices / partitions + 1;
    for(int i = 0; i < partitions; i++) {
		if (i < partitionSplit) {
				gridOffset[i] = i * partitionSize;
		} else {
				gridOffset[i] = partitionSplit * partitionSize + (i - partitionSplit) * (partitionSize - 1);
		}
    }
	for (int i = 0; i < partitions - 1; i++) {
		tempOffset[i] = gridOffset[i+1];
	}
	tempOffset[partitions - 1] = vertices;
    const int vertex_buffer_size = sizeof(VertexId) * 2;
    VertexId start = 0, end = vertices-1;
    char *vertexBuffer = (char *)memalign(PAGESIZE, vertex_buffer_size*vertices);
    while (start<end)
     {
        VertexId pID = start % partitions;
		while(gridOffset[pID] >= tempOffset[pID]) {
			pID = (pID + 1) % partitions;
		}
        verticesOut[tt[start].first] = gridOffset[pID];
		*(VertexId*) (vertexBuffer + gridOffset[pID] * vertex_buffer_size) = gridOffset[pID];
		*(VertexId*) (vertexBuffer + gridOffset[pID] * vertex_buffer_size + sizeof(VertexId)) = tt[start].second;
		// printf("pid:%d, %d, %d, %d\n",pID, gridOffset[pID], tt[start].second, start);
        gridOffset[pID]++;
		while(gridOffset[pID] >= tempOffset[pID]) {
			pID = (pID + 1) % partitions;
		}
        verticesOut[tt[end].first] = gridOffset[pID];
		*(VertexId*) (vertexBuffer + gridOffset[pID] * vertex_buffer_size) = gridOffset[pID];
		*(VertexId*) (vertexBuffer + gridOffset[pID] * vertex_buffer_size + sizeof(VertexId)) = tt[end].second;
        // printf("pid:%d, %d, %d, %d\n",pID, gridOffset[pID], tt[end].second, end);
		// printf("******************\n");
        gridOffset[pID]++;
        start++;
        end--;
    }
    if(start == end){
        VertexId pID = start % partitions;
		while(gridOffset[pID] >= tempOffset[pID]) {
			pID = (pID + 1) % partitions;
		}
        verticesOut[tt[start].first] = gridOffset[pID];
		*(VertexId*) (vertexBuffer + gridOffset[pID] * vertex_buffer_size) = gridOffset[pID];
		*(VertexId*) (vertexBuffer + gridOffset[pID] * vertex_buffer_size + sizeof(VertexId)) = tt[start].second;
		// printf("pid:%d, %d, %d, %d\n",pID, gridOffset[pID], tt[start].second, start);
    }

	// 创建输出文件夹
    if (file_exists(output)) {
		remove_directory(output);
	}
	create_directory(output);
	printf("create dir: %s\n", output.c_str());
	int fout;
	int pOffset = partitions / machineNum;
	for (int mid = 0; mid < machineNum; mid++) {
		char filename[4096];
		sprintf(filename, "%s/vertex-%d", output.c_str(), mid);
		fout = open(filename, O_WRONLY|O_APPEND|O_CREAT, 0644);
		for (int i = mid * pOffset; i < (mid + 1) * pOffset; i++) {
			if (i < partitionSplit) {
					gridOffset[i] = i * partitionSize;
			} else {
					gridOffset[i] = partitionSplit * partitionSize + (i - partitionSplit) * (partitionSize - 1);
			}
			
			if(i < partitionSplit) {
				write(fout, vertexBuffer + gridOffset[i] * vertex_buffer_size, partitionSize * vertex_buffer_size);
			} else {
				write(fout, vertexBuffer + gridOffset[i] * vertex_buffer_size, (partitionSize - 1) * vertex_buffer_size);
			}
		}
		close(fout);
	}

    /*for(int i=0;i<vertices;i++) {
        printf("vertice id: %d, out: %d; New vertice id: %d\n", tt[i].first, tt[i].second, verticesOut[i]);
    }*/
	tt.clear();
	tt.shrink_to_fit();
	gridPartition(input, output, vertices, edges, edgeUnit, partitions, edgeType, machineNum, verticesOut);
	verticesOut.clear();
	verticesOut.shrink_to_fit();
}

void gridPartition(std::string input, std::string output, VertexId vertices, EdgeId edges, int edgeUnit, int partitions, int edgeType, int machineNum, std::vector<int> &verticesOut) {
    int parallelism = std::thread::hardware_concurrency();
    char ** buffers = new char * [parallelism*2];
	bool * occupied = new bool [parallelism*2];
	for (int i=0;i<parallelism*2;i++) {
		buffers[i] = (char *)memalign(PAGESIZE, IOSIZE);
		occupied[i] = false;
	}

    Queue<std::tuple<int, long> > tasks(parallelism);
    /***************************创建输出文件夹*******************************/
	int ** fout; //输出文件的文件描述符
	std::mutex ** mutexes; //互斥锁
	fout = new int * [partitions];
	mutexes = new std::mutex * [partitions];
        if (!file_exists(output)) {
	    printf("output dir is not exist.\n");
	    return;
	} else {
	    printf("%s\n", output.c_str());
	}	
    /***************************创建输出文件夹*******************************/

    /************************创建缓存并计算偏移量****************************/
    const int grid_buffer_size = 768; // 12 * 8 * 8
	char * global_grid_buffer = (char *) memalign(PAGESIZE, grid_buffer_size * partitions * partitions);
	char *** grid_buffer = new char ** [partitions];
	int ** grid_buffer_offset = new int * [partitions];
	for (int i=0;i<partitions;i++) {
		mutexes[i] = new std::mutex [partitions];
		fout[i] = new int [partitions];
		grid_buffer[i] = new char * [partitions];
		grid_buffer_offset[i] = new int [partitions];
		for (int j=0;j<partitions;j++) {
			char filename[4096];
			sprintf(filename, "%s/block-%d-%d", output.c_str(), i, j);
			fout[i][j] = open(filename, O_WRONLY|O_APPEND|O_CREAT, 0644);
			grid_buffer[i][j] = global_grid_buffer + (i * partitions + j) * grid_buffer_size;
			grid_buffer_offset[i][j] = 0;
		}
	}
    /************************创建缓存并计算偏移量****************************/

    /******************创建多线程vector与多线程执行函数**********************/
    std::vector<std::thread> threads;
    for (int ti=0;ti<parallelism;ti++) {
		threads.emplace_back([&]() {
			char * local_buffer = (char *) memalign(PAGESIZE, IOSIZE);
			int * local_grid_offset = new int [partitions * partitions];
			int * local_grid_cursor = new int [partitions * partitions]; //光标、标记
			VertexId source, target;
			Weight weight;
			while (true) {
				int cursor;
				long bytes;
				std::tie(cursor, bytes) = tasks.pop();
				if (cursor==-1) break;
				memset(local_grid_offset, 0, sizeof(int) * partitions * partitions);
				memset(local_grid_cursor, 0, sizeof(int) * partitions * partitions);
				char * buffer = buffers[cursor];
				for (long pos=0;pos<bytes;pos+=edgeUnit) {
						source = *(VertexId*)(buffer+pos);
						target = *(VertexId*)(buffer+pos+sizeof(VertexId));
						int i = get_partition_id(vertices, partitions, verticesOut[source]);
						int j = get_partition_id(vertices, partitions, verticesOut[target]);
						local_grid_offset[i*partitions+j] += edgeUnit;
						//printf("source: %d, target: %d, i: %d, j: %d, local_grid_offset[%d]: %d\n", verticesOut[source], verticesOut[target], i, j, i*partitions+j, local_grid_offset[i*partitions+j]);
						//local_grid_offset[i*partitions+j] += edgeUnit;
				}
				local_grid_cursor[0] = 0;
				for (int ij=1;ij<partitions*partitions;ij++) {
						local_grid_cursor[ij] = local_grid_offset[ij - 1];
						local_grid_offset[ij] += local_grid_cursor[ij];
				}
				//printf("local_grid_offset[partitions*partitions-1]: %d, bytes: %d\n", local_grid_offset[partitions*partitions-1], bytes);
				assert(local_grid_offset[partitions*partitions-1]==bytes);
				for (long pos=0;pos<bytes;pos+=edgeUnit) {
						source = *(VertexId*)(buffer+pos);
						target = *(VertexId*)(buffer+pos+sizeof(VertexId));
						int i = get_partition_id(vertices, partitions, verticesOut[source]);
						int j = get_partition_id(vertices, partitions, verticesOut[target]);
						*(VertexId*)(local_buffer+local_grid_cursor[i*partitions+j]) = verticesOut[source];
						*(VertexId*)(local_buffer+local_grid_cursor[i*partitions+j]+sizeof(VertexId)) = verticesOut[target];
						if (edgeType==1) {
								weight = *(Weight*)(buffer+pos+sizeof(VertexId)*2);
								*(Weight*)(local_buffer+local_grid_cursor[i*partitions+j]+sizeof(VertexId)*2) = weight;
						}
						local_grid_cursor[i*partitions+j] += edgeUnit;
				}
				int start = 0;
				for (int ij=0;ij<partitions*partitions;ij++) {
						assert(local_grid_cursor[ij]==local_grid_offset[ij]);
						int i = ij / partitions;
						int j = ij % partitions;
						std::unique_lock<std::mutex> lock(mutexes[i][j]);
						if (local_grid_offset[ij] - start > edgeUnit) {
								write(fout[i][j], local_buffer+start, local_grid_offset[ij]-start);
						} else if (local_grid_offset[ij] - start == edgeUnit) {
								memcpy(grid_buffer[i][j]+grid_buffer_offset[i][j], local_buffer+start, edgeUnit);
								grid_buffer_offset[i][j] += edgeUnit;
								if (grid_buffer_offset[i][j]==grid_buffer_size) {
										write(fout[i][j], grid_buffer[i][j], grid_buffer_size);
										grid_buffer_offset[i][j] = 0;
								}
						}
						start = local_grid_offset[ij];
				}
				occupied[cursor] = false;
			}
		});
	}
    /******************创建多线程vector与多线程执行函数**********************/

    /**************打开输入文件，并以IOSIZE为大小读入边数据*******************/
    int fin = open(input.c_str(),O_RDONLY);
    MY_ASSERT(fin != -1, "[preprocess] Cannot open input file.\n");
    int cursor = 0;
	long total_bytes = file_size(input);
	long read_bytes = 0;
	//double start_time = get_time();
	while (true) {
		long bytes = read(fin, buffers[cursor], IOSIZE);
		assert(bytes!=-1);
		if (bytes==0) break;
		occupied[cursor] = true; //当前标记cursor被使用，故设为true
		tasks.push(std::make_tuple(cursor, bytes));
		read_bytes += bytes;
		printf("progress: %.2f%%\r", 100. * read_bytes / total_bytes);
		fflush(stdout);
		//当图数据大小过大存在死循环的可能
		while (occupied[cursor]) {
			cursor = (cursor + 1) % (parallelism * 2);
		}
	}
	close(fin);
    /**************打开输入文件，并以IOSIZE为大小读入边数据*******************/

    /**************************输入线程退出条件*****************************/
    MY_ASSERT(read_bytes==edges*edgeUnit, "[preprocess] Do not read all edges from input file.\n");
    for (int ti=0;ti<parallelism;ti++) {
		tasks.push(std::make_tuple(-1, 0));
	}
    /**************************输入线程退出条件*****************************/

    /*****************************多线程运行*******************************/
	for (int ti=0;ti<parallelism;ti++) {
		threads[ti].join();
	}
    /*****************************多线程运行*******************************/

    /**********************将分区结果写到多个输出文件***********************/
    long ts = 0;
	for (int i=0;i<partitions;i++) {
		for (int j=0;j<partitions;j++) {
			if (grid_buffer_offset[i][j]>0) {
				ts += grid_buffer_offset[i][j];
				write(fout[i][j], grid_buffer[i][j], grid_buffer_offset[i][j]);
			}
		}
	}
    for (int i=0;i<partitions;i++) {
		for (int j=0;j<partitions;j++) {
			close(fout[i][j]);
		}
	}
    /**********************将分区结果写到多个输出文件***********************/

    /***********将图分区文件按照行进行合并，并记录对应的偏移量************/
    long offset[machineNum];
	int fout_row[machineNum];
	int fout_row_offset[machineNum];
	for(int i=0;i<machineNum;i++) {
		char rowName[4096];
		sprintf(rowName, "%s/row-%d", output.c_str(), i);
		fout_row[i] = open(rowName, O_WRONLY|O_APPEND|O_CREAT, 0644);
		char rowOffsetName[4096];
		sprintf(rowOffsetName, "%s/row_offset-%d", output.c_str(), i);
		fout_row_offset[i] = open(rowOffsetName, O_WRONLY|O_APPEND|O_CREAT, 0644);
		offset[i] = 0;
	}
	int partitionsOffset = partitions / machineNum;
	int partitionsStart = 0;
	for (int i=0;i<machineNum;i++) {
		int nextMachine = (i + 1) % machineNum;
		for(int k= partitionsStart;k<partitionsStart+partitionsOffset;k++) {
			for (int x = 0; x < partitionsOffset; x++) {
				for (int y = 0; y < machineNum; y++) {
					// printf("progress: %.2f%%\r", 100. * offset[i] / total_bytes);
					// fflush(stdout);
					write(fout_row_offset[i], &offset[i], sizeof(offset[i]));
					char filename[4096];
					int next = (nextMachine + y) % machineNum;
					int nextPartition = (next * partitionsOffset + x) % partitions;
					sprintf(filename, "%s/block-%d-%d", output.c_str(), k, nextPartition);
					std::cout << k << " " << nextPartition << std::endl;
					offset[i] += file_size(filename);
					fin = open(filename, O_RDONLY);
					while (true) {
						long bytes = read(fin, buffers[0], IOSIZE);
						assert(bytes!=-1);
						if (bytes==0) break;
						write(fout_row[i], buffers[0], bytes);
					}
					close(fin);
				}
			}
		}
		partitionsStart += partitionsOffset;
	}
	for(int i=0;i<machineNum;i++) {
		write(fout_row_offset[i], &offset[i], sizeof(offset[i]));
		close(fout_row_offset[i]);
		close(fout_row[i]);
	}
	printf("row oriented grid generated\n");

	//printf("it takes %.2f seconds to generate edge grid\n", get_time() - start_time);

	FILE * fmeta = fopen((output+"/meta").c_str(), "w");
	fprintf(fmeta, "%d %d %ld %d %d", edgeType, vertices, edges, partitions, machineNum);
	fclose(fmeta);
    /********************将每个分区的属性写到对应文件中*********************/
}

int main(int argc, char ** argv) {
    //printf("./preGraphProcess [input path] [output path] [vertices number] [partitions] [edge type]\n");
    std::string input = "";
    std::string output = "";
    VertexId vertices = 0;
    int partitions = 0;
    int edgeType = 0;
	int machineNum = 0;
    if(argc != 7){
        printf("./preGraphProcess [input path] [output path] [vertices number] [partitions] [edge type] [machine number]\n");
        return 0;
    }
    else {
        input = argv[1];
        output = argv[2];
        vertices = atoi(argv[3]);
        partitions = atoi(argv[4]);
        edgeType = atoi(argv[5]);
		machineNum = atoi(argv[6]);
    }
    preProcess(input, output, vertices, partitions, edgeType, machineNum);
    //gridPartition(input, output, vertices, edges, edgeUnit, partitions, edgeType, machineNum, verticesOut);
    printf("Success!\n");


}
