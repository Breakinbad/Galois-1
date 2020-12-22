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
#include "galois/LargeArray.h"

namespace cll = llvm::cl;

static const char* name = "bipart";
static const char* desc =
    "Partitions a hypergraph into K parts and minimizing the graph cut";
static const char* url = "bipart";


static cll::opt<std::string> outfile("output",
                                     cll::desc("output partition file name"));
static cll::opt<std::string> filename(cll::Positional,
                                      cll::desc("<input input>"), cll::Required);
static cll::opt<unsigned> nodes(cll::Positional,
                                   cll::desc("<number of nodes>"),
                                   cll::Required);
static cll::opt<unsigned> hedges(cll::Positional,
                                   cll::desc("<number of hedges>"),
                                   cll::Required);

static cll::opt<double> imbalance(
    "balance",
    cll::desc("Fraction deviated from mean partition size (default 0.01)"),
    cll::init(0.01));
struct myNode {
  unsigned id;
  std::vector<double> EMvec = std::vector<double>(dim);

};

int main(int argc, char** argv) {
  galois::SharedMemSys G;
  LonestarStart(argc, argv, name, desc, url);

 // srand(-1);
  MetisGraph metisGraph;
  GGraph& graph = *metisGraph.getGraph();
  //GGraph graph;// = *metisGraph.getGraph();
  std::map<unsigned,myNode> node;// = std::vector<myNode>(hedges+nodes);
  std::map<unsigned, unsigned> newnodes;
  galois::StatTimer T("buildingG");
  T.start();
  // read rest of input and initialize hedges (build hgraph)
  

  unsigned hyperedges = 1195945;
  std::ifstream in(filename.c_str());
  
  std::string line1;
  unsigned lnum = 0;
  while(std::getline(in, line1)) {
    std::stringstream ss(line1);
    unsigned n;
    ss >> n;
   // std::cout<<n<<"\n";
    double val;
    unsigned iter = 0;
    while(ss >> val) {
      node[n].EMvec[iter] = val;
      iter++;
     // std::cout<<iter<<" ";
    }

  }


  std::ofstream out("LivejournalEmbbeding.txt");
  for (auto n = 0; n < hyperedges; n++) {
    //out<<n<<"\t";
    for (int i = 0; i < dim; i++)
      out<<node[n].EMvec[i]<<"\t";
    out<<"\n";
  }

  return 0;
}

