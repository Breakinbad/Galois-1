/*
 * This file belongs to the Galois project, a C++ library for exploiting parallelism.
 * The code is being released under the terms of the 3-Clause BSD License (a
 * copy is located in LICENSE.txt at the top-level directory).
 *
 * Copyright (C) 2018, The University of Texas at Austin. All rights reserved.
 * UNIVERSITY EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES CONCERNING THIS
 * SOFTWARE AND DOCUMENTATION, INCLUDING ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR ANY PARTICULAR PURPOSE, NON-INFRINGEMENT AND WARRANTIES OF
 * PERFORMANCE, AND ANY WARRANTY THAT MIGHT OTHERWISE ARISE FROM COURSE OF
 * DEALING OR USAGE OF TRADE.  NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH
 * RESPECT TO THE USE OF THE SOFTWARE OR DOCUMENTATION. Under no circumstances
 * shall University be liable for incidental, special, indirect, direct or
 * consequential damages or loss of profits, interruption of business, or
 * related expenses which may arise from use of Software or Documentation,
 * including but not limited to those resulting from defects in Software and/or
 * Documentation, or loss or inaccuracy of data of any kind.
 */

#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <array>
#include <unordered_set>

#include "bipart.h"
#include "galois/graphs/Util.h"
#include "galois/Timer.h"
//#include "GraphReader.h"
#include "Lonestar/BoilerPlate.h"
#include "galois/graphs/FileGraph.h"
#include "galois/LargeArray.h"

namespace cll = llvm::cl;

static const char* name = "HYPAR";
static const char* desc =
    "Partitions a hypergraph into K parts and minimizing the graph cut";
static const char* url = "HyPar";

static cll::opt<scheduleMode> schedulingMode(
    cll::desc("Choose a inital scheduling mode:"),
    cll::values(clEnumVal(PLD, "PLD"), clEnumVal(PP, "PP"), clEnumVal(WD, "WD"),
                clEnumVal(RI, "RI"), clEnumVal(MRI, "MRI"),clEnumVal(MDEG, "MDEG"),clEnumVal(DEG, "DEG"),clEnumVal(MWD, "MWD"),clEnumVal(HIS, "HIS"),clEnumVal(RAND, "random"),clEnumValEnd),
    cll::init(PLD));

static cll::opt<bool>
    mtxInput("mtxinput",
             cll::desc("Use text mtx files instead of binary galois gr files"),
             cll::init(false));
static cll::opt<bool> weighted("weighted", cll::desc("weighted"),
                               cll::init(false));
static cll::opt<bool>
    verbose("verbose",
            cll::desc("verbose output (debugging mode, takes extra time)"),
            cll::init(false));
static cll::opt<std::string> outfile("output",
                                     cll::desc("output partition file name"));
static cll::opt<std::string>
    orderedfile("ordered", cll::desc("output ordered graph file name"));
static cll::opt<std::string>
    permutationfile("permutation", cll::desc("output permutation file name"));
static cll::opt<std::string> filename(cll::Positional,
                                      cll::desc("<input file>"), cll::Required);
static cll::opt<unsigned> numPartitions(cll::Positional,
                                   cll::desc("<number of partitions>"),
                                   cll::Required);
static cll::opt<unsigned> csize(cll::Positional,
                                   cll::desc("<size of coarsest graph>"),
                                   cll::Required);

static cll::opt<double> imbalance(
    "balance",
    cll::desc("Fraction deviated from mean partition size (default 0.01)"),
    cll::init(0.01));

// const double COARSEN_FRACTION = 0.9;

/*int cutsize(GGraph& g) { 
  unsigned size = std::distance(g.cellList().begin(), g.cellList().end());
  unsigned sizen = std::distance(g.getNets().begin(), g.getNets().end());
  int cutsize = 0;
  std::vector<int> cells;
  for (auto n : g.getNets()) { 
    bool cut_status = false;
    for (auto e : g.edges(n)) {
      auto cell1 = g.getEdgeDst(e);
    for (auto c : g.edges(n)) {
        auto cell2 = g.getEdgeDst(c);
        if(g.getData(cell1).getPart() != g.getData(cell2).getPart() && cell1 != cell2) {
          cutsize++;
          cut_status = true;
          break;
        }
      }
      if (cut_status == true)
        break;
    }
  }
  return cutsize;
}*/
/**
 * Partitioning 
 */
void Partition(MetisGraph* metisGraph, unsigned coarsenTo, unsigned K) {
  std::cout<<"in partition\n";
  galois::StatTimer TM;
  TM.start();

  galois::StatTimer T("CoarsenSEP");
  T.start();
  std::cout<<"in coarsening\n";
  //MetisGraph* mcg = coarsen(metisGraph, coarsenTo, schedulingMode);
  T.stop();

 galois::StatTimer T2("PartitionSEP");
  T2.start();
  partition(metisGraph, K);
  T2.stop();


  galois::StatTimer T3("Refine");
  T3.start();
  //refine(mcg, K);
  T3.stop();
  int one = 0;
  int zero = 0;
  std::cout << "coarsen:," << T.get() << "\n";
  std::cout << "clustering:," << T2.get() << '\n';
  std::cout << "Refinement:," << T3.get() << "\n";
  return;
}

int computingCut(GGraph& g) {

  GNodeBag bag;
  galois::GAccumulator<unsigned> edgecut;
  galois::do_all(
      galois::iterate((size_t)0, g.hedges),
      [&](GNode n) {
        std::set<unsigned> nump;
        for (auto cell : g.edges(n)) {
          auto c   = g.getEdgeDst(cell);
          int part = g.getData(c).getPart();
          nump.insert(part);
        }
        edgecut += (nump.size() - 1);
      },
      galois::loopname("cutsize"));
  return edgecut.reduce();
}

int computingBalance(GGraph& g) {
  int zero = 0, one = 0;
  for (int c = g.hedges; c < g.size(); c++) {
    int part = g.getData(c).getPart();
    if (part == 0) zero++;
    else one++;
  }
  return std::abs(zero - one);
}
// printGraphBeg(*graph)

int hash(unsigned val) {
  unsigned long int seed = val * 1103515245 + 12345;
  return((unsigned)(seed/65536) % 32768);
}

int main(int argc, char** argv) {
  galois::SharedMemSys G;
  LonestarStart(argc, argv, name, desc, url);

 // srand(-1);
  MetisGraph metisGraph;
  GGraph& graph = *metisGraph.getGraph();
  std::ifstream f(filename.c_str());
  //GGraph graph;// = *metisGraph.getGraph();
  std::string line;
  std::getline(f, line);
  std::stringstream ss(line);
  uint32_t i1;
  uint64_t i2;
  ss >> i1 >> i2;
  const int hedges = i1, nodes = i2;
  printf("hedges: %d\n", hedges);
  printf("nodes: %d\n\n", nodes);

  galois::StatTimer T("buildingG");
  T.start();
  // read rest of input and initialize hedges (build hgraph)
  galois::gstl::Vector<galois::PODResizeableArray<uint32_t> > edges_id(hedges+nodes);
  std::vector<std::vector<EdgeTy> > edges_data(hedges+nodes);
  std::vector<uint64_t> prefix_edges(nodes+hedges);
  int cnt = 0, edges = 0;
  while (std::getline(f, line)) {
    if (cnt >= hedges) {printf("ERROR: too many lines in input file\n"); exit(-1);}
    std::stringstream ss(line);
    int val;
    while (ss >> val) {
      //if ((val < 1) || (val > nodes)) {printf("ERROR: node value %d out of bounds--", val); std::cout<<val<<"SSS ";exit(-1);}
      unsigned newval = hedges + (val);// - 1);
      edges_id[cnt].push_back(newval);
      edges_id[newval].push_back(cnt);
      edges+=2;
    }
    cnt++;
  }
  f.close();
  graph.hedges = hedges;
  graph.hnodes = nodes;
  std::cout<<"number of edges "<<edges<<"\n";
  uint32_t sizes = hedges+nodes;
  galois::do_all(galois::iterate((uint32_t)0, sizes),
                [&](uint32_t c){
                  prefix_edges[c] = edges_id[c].size();
                });
  
  for (uint32_t c = 1; c < nodes+hedges; ++c) {
    prefix_edges[c] += prefix_edges[c - 1];
  }
  // edges = #edges, hedgecount = how many edges each node has, edges_id: for each node, which ndoes it is connected to
  // edges_data: data for each edge = 1
  graph.constructFrom(nodes+hedges, edges, prefix_edges, edges_id);
  galois::do_all(galois::iterate(graph),
                  [&](GNode n) {
                    if (n < hedges)
                      graph.getData(n).netnum = n;
                    else
                      graph.getData(n).netnum = INT_MAX;
                    graph.getData(n).netrand = INT_MAX;
                    graph.getData(n).netval = INT_MAX;
                    graph.getData(n).nodeid = n;
  
  });
  T.stop();
  std::cout<<"time to build a graph "<<T.get()<<"\n";
  graphStat(graph);
  std::cout<<"\n";
  galois::preAlloc(galois::runtime::numPagePoolAllocTotal() * 5);
  galois::reportPageAlloc("MeminfoPre");
  galois::do_all(
      galois::iterate(graph.hedges, graph.size()),
      [&](GNode item) {
        // accum += g->getData(item).getWeight();
        graph.getData(item, galois::MethodFlag::UNPROTECTED)
            .initRefine(0, true);
        graph.getData(item, galois::MethodFlag::UNPROTECTED).initPartition();
      },
      galois::loopname("initPart"));

  Partition(&metisGraph, csize, 2);
  /*std::ofstream out("amazonembbeding.txt");
  for (auto n = 0; n < graph.hedges; n++) {
    //out<<n<<"\t";
    for (int i = 0; i < dim; i++)
      out<<graph.getData(n).EMvec[i]<<"\t";
    out<<"\n";
  }*/
  //std::cout<<"Total Edge Cut: "<<computingCut(graph)<<"\n";
  //galois::runtime::reportStat_Single("HyPar", "Edge Cut", computingCut(graph));
  //galois::runtime::reportStat_Single("HyParzo", "zero-one", computingBalance(graph));
   galois::reportPageAlloc("MeminfoPost");

  return 0;
}

