cmake_minimum_required(VERSION 3.10)
project(m_core)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_DEBUG   "${CMAKE_CXX_FLAGS_DEBUG} -std=c++11 -g -Wall -Wno-unused-variable -pthread")
set(CMAKE_CXX_FLAGS_RELEASE   "${CMAKE_CXX_FLAGS_RELEASE} -std=c++11 -O3 -pthread -fopenmp")

add_executable(abSearch abSearch.cpp RStarTree.h skylineCores.h graph.h utilities/UnionFind.cpp)
add_executable(klSearch klSearch.cpp RStarTree.h kl-core/KlSkylineCore.h kl-core/DiGraph.h utilities/UnionFind.cpp kl-core/DiGraph.cpp)
add_executable(klUpdate klUpdate.cpp RStarTree.h kl-core/KlSkylineCore.h kl-core/DiGraph.h utilities/UnionFind.cpp kl-core/DiGraph.cpp)
add_executable(abUpdate abUpdate.cpp RStarTree.h skylineCores.h graph.h utilities/UnionFind.cpp)
add_executable(multilayerSearch multilayerSearch.cpp RStarTree.h multilayer/multilayer_graph.h multilayer/multilayer_graph.cpp utilities/UnionFind.cpp)
add_executable(multilayerUpdate multilayerUpdate.cpp RStarTree.h multilayer/multilayer_graph.h multilayer/multilayer_graph.cpp utilities/UnionFind.cpp)

add_executable(m_core main.cpp)
