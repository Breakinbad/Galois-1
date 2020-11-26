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

#include "bipart.h"
#include "galois/Galois.h"
#include "galois/AtomicHelpers.h"
#include "galois/Reduction.h"
#include "galois/runtime/Profile.h"
#include "galois/substrate/PerThreadStorage.h"
#include "galois/gstl.h"

#include <iostream>
#include <unordered_set>

int TOTALW;
int LIMIT;
bool FLAG = false;
const int dim = 200;
namespace {

typedef galois::GAccumulator<unsigned> Pcounter;

double computeScore(GGraph& g, GNode n, GNode m) {
  double score = 0.0;
  //double normN = 0.0;
  //double normM = 0.0;
  for (int i = 0; i < dim; i++) {
    double val = g.getData(n).EMvec[i] * g.getData(m).EMvec[i];
    score += val;
   // normN += pow(g.getData(n).EMvec[i], 2);
   // normM += pow(g.getData(m).EMvec[i], 2);
  }
  
  return score ;/// (sqrt(normN) * sqrt(normM));
}

std::vector<double> mergeEmd(GGraph& g, std::vector<GNode> vec) {
  std::vector<double> tmp(dim);
  for (int i = 0; i < dim; i++)
    tmp[i] = 0.0;
  for (auto v : vec) {
    for (int j = 0; j < dim; j++)
      tmp[j] += g.getData(v).EMvec[j];
  }
  for (int i = 0; i < dim; i++)
    tmp[i] /= vec.size();
  
 return tmp;
}

int hash(unsigned val) {
  unsigned long int seed = val * 1103515245 + 12345;
  return((unsigned)(seed/65536) % 32768);
}

void parallelRand(MetisGraph* graph, int iter) {

  GGraph* fineGGraph  = graph->getFinerGraph()->getGraph();
  galois::do_all(
      galois::iterate((uint64_t)0, fineGGraph->hedges),
      [&](GNode item) {
          fineGGraph->getData(item).netrand = hash(fineGGraph->getData(item).netnum);
    },
    galois::loopname("rand"));
}

using MatchingPolicy  = void(GNode, GGraph*);

void PLD_f(GNode node, GGraph* fineGGraph) {
  int ss = std::distance(fineGGraph->edge_begin(node), fineGGraph->edge_end(node));
  fineGGraph->getData(node).netval = -ss;
}
void RAND_f(GNode node, GGraph* fineGGraph) {
  unsigned id = fineGGraph->getData(node).netrand;
  fineGGraph->getData(node).netval = -id;
  fineGGraph->getData(node).netrand = -fineGGraph->getData(node).netnum; 
}
void PP_f(GNode node, GGraph* fineGGraph) {
  int ss = std::distance(fineGGraph->edge_begin(node), fineGGraph->edge_end(node));
  fineGGraph->getData(node).netval = ss;
}
void WD_f(GNode node, GGraph* fineGGraph) {
  int w = 0;
  for (auto n : fineGGraph->edges(node)) {
    auto nn = fineGGraph->getEdgeDst(n);
    w += fineGGraph->getData(nn).getWeight();
  }
  fineGGraph->getData(node).netval = -w;
}
void MWD_f(GNode node, GGraph* fineGGraph) {
  int w = 0;
  for (auto n : fineGGraph->edges(node)) {
    auto nn = fineGGraph->getEdgeDst(n);
    w += fineGGraph->getData(nn).getWeight();
  }
  fineGGraph->getData(node).netval = w;
}
void RI_f(GNode node, GGraph* fineGGraph) {
  int ss = std::distance(fineGGraph->edge_begin(node), fineGGraph->edge_end(node));
  fineGGraph->getData(node).netval = ss;
}
void MRI_f(GNode node, GGraph* fineGGraph) {
  int ss = std::distance(fineGGraph->edge_begin(node), fineGGraph->edge_end(node));
  fineGGraph->getData(node).netval = ss;
}
void DEG_f(GNode node, GGraph* fineGGraph) {
  int w = 0;
  int ss = std::distance(fineGGraph->edge_begin(node), fineGGraph->edge_end(node));
  fineGGraph->getData(node).netval = ss;
  for (auto n : fineGGraph->edges(node)) {
    auto nn = fineGGraph->getEdgeDst(n);
    w += fineGGraph->getData(nn).getWeight();
  }
  fineGGraph->getData(node).netval = - (w / ss);
}
void MDEG_f(GNode node, GGraph* fineGGraph) {
  int w = 0;
  int ss = std::distance(fineGGraph->edge_begin(node), fineGGraph->edge_end(node));
  fineGGraph->getData(node).netval = ss;
  for (auto n : fineGGraph->edges(node)) {
    auto nn = fineGGraph->getEdgeDst(n);
    w += fineGGraph->getData(nn).getWeight();
  }
  fineGGraph->getData(node).netval = w/ss;
}



template <MatchingPolicy matcher> 
void parallelPrioRand(MetisGraph* graph, int iter) {

  GGraph* fineGGraph   = graph->getFinerGraph()->getGraph();
  parallelRand(graph, iter); 
  
  // Making deterministic
  galois::do_all(
      galois::iterate((uint64_t)0, fineGGraph->hedges),
      [&](GNode item) {
          matcher(item, fineGGraph);
          for (auto c : fineGGraph->edges(item)) {
            auto dst = fineGGraph->getEdgeDst(c);
            galois::atomicMin(fineGGraph->getData(dst).netval, fineGGraph->getData(item).netval.load());
          }
      },
      galois::steal(),  galois::loopname("atomicMin"));
  galois::do_all(
      galois::iterate((uint64_t)0, fineGGraph->hedges),
      [&](GNode item) {
            for (auto c : fineGGraph->edges(item)) {
                auto dst = fineGGraph->getEdgeDst(c);
                if (fineGGraph->getData(dst).netval == fineGGraph->getData(item).netval)
                galois::atomicMin(fineGGraph->getData(dst).netrand, fineGGraph->getData(item).netrand.load());
            }  
     },
      galois::steal(),  galois::loopname("secondMin2"));
  galois::do_all(
      galois::iterate((uint64_t)0, fineGGraph->hedges),
      [&](GNode item) {
        for (auto c : fineGGraph->edges(item)) {
          auto dst = fineGGraph->getEdgeDst(c);
          if (fineGGraph->getData(dst).netrand == fineGGraph->getData(item).netrand)
            galois::atomicMin(fineGGraph->getData(dst).netnum, fineGGraph->getData(item).netnum.load());
        }
      },
      galois::steal(),  galois::loopname("secondMin"));
}


// hyper edge matching
void parallelHMatchAndCreateNodes(MetisGraph* graph,
                                 int iter, GNodeBag& bag, std::vector<bool>& hedges, std::vector<unsigned>& weight) {
  //parallelPrioRand<matcher>(graph, iter);
  GGraph* fineGGraph   = graph->getFinerGraph()->getGraph();
  GGraph* coarseGGraph = graph->getGraph();
  assert(fineGGraph != coarseGGraph);
  galois::do_all(
      galois::iterate((uint64_t)0,fineGGraph->hedges),
      [&](GNode item) {
        for (auto c : fineGGraph->edges(item)) {
          auto dst = fineGGraph->getEdgeDst(c);
          double score = computeScore(*fineGGraph, item, dst);
          if (fineGGraph->getData(dst).netval < score) {
            fineGGraph->getData(dst).netnum.store(fineGGraph->getData(item).netnum);
            fineGGraph->getData(dst).netval.store(score);
          }
        }
      },
  galois::steal(),  galois::loopname("atomicScore"));
  typedef std::vector<GNode> VecTy;
  typedef galois::substrate::PerThreadStorage<VecTy> ThreadLocalData;
  ThreadLocalData edgesThreadLocal;
  std::string name = "phaseI";
  // hyperedge coarsening 
  galois::do_all(
      galois::iterate((uint64_t)0,fineGGraph->hedges),
      [&](GNode item) {
        unsigned id = fineGGraph->getData(item).netnum;
        bool flag = false;
        unsigned nodeid = INT_MAX;
        auto& edges = *edgesThreadLocal.getLocal();
        edges.clear();
        for (auto c : fineGGraph->edges(item)) {
          auto dst = fineGGraph->getEdgeDst(c);
          auto& data = fineGGraph->getData(dst);
          unsigned nid = data.nodeid;
          if (data.isMatched()) {
            flag = true;
            continue;
          }
          if (data.netnum == fineGGraph->getData(item).netnum) {
            edges.push_back(dst);
            nodeid = std::min(nodeid, dst);
          }
          else { 
             flag = true;
          }
        }

        if (!edges.empty()) {
          if (flag && edges.size() == 1) return; 
            fineGGraph->getData(item).setMatched();
          if (flag) hedges[item] = true;
          bag.push(nodeid);
          unsigned ww = 0;
          for (auto pp : edges) {
            ww += fineGGraph->getData(pp).getWeight();
            fineGGraph->getData(pp).setMatched();
            fineGGraph->getData(pp).setParent(nodeid);
            fineGGraph->getData(pp).netnum = fineGGraph->getData(item).netnum;
          }
          fineGGraph->getData(nodeid).EMvec = mergeEmd(*fineGGraph, edges);
          weight[nodeid-fineGGraph->hedges] = ww;
        }
      },
      galois::loopname("phaseI"));
}

void moreCoarse(MetisGraph* graph, int iter, std::vector<unsigned>& weight) {
  
  GGraph* fineGGraph   = graph->getFinerGraph()->getGraph();
  typedef std::set<int> SecTy;
  typedef std::vector<GNode> VecTy;
  GNodeBag bag;
  typedef galois::substrate::PerThreadStorage<VecTy> ThreadLocalData;
  ThreadLocalData edgesThreadLocal;
  galois::do_all(
      galois::iterate((uint64_t)0, fineGGraph->hedges),
      [&](GNode item) {
        if (fineGGraph->getData(item).isMatched()) return;
          for (auto c : fineGGraph->edges(item)) {
            auto dst = fineGGraph->getEdgeDst(c);
            if (fineGGraph->getData(dst).isMatched()) 
                fineGGraph->getData(dst).netval = INT_MIN;
          }
      },
      galois::steal(),  galois::loopname("atomicMin2"));
  galois::do_all( 
      galois::iterate((uint64_t)0, fineGGraph->hedges),
      [&](GNode item) {
          if (fineGGraph->getData(item).isMatched()) return;
          auto& cells = *edgesThreadLocal.getLocal();
          cells.clear();
          int best = INT_MAX;
          GNode b;
          //int w = 0;
          for (auto edge : fineGGraph->edges(item)) {
	      auto e = fineGGraph->getEdgeDst(edge);
              auto& data = fineGGraph->getData(e);
              if (!fineGGraph->getData(e).isMatched()) {
                  if (data.netnum == fineGGraph->getData(item).netnum) {
                      cells.push_back(e);
                  }
              }
              else if (fineGGraph->getData(e).netval == INT_MIN) {
                  auto nn = fineGGraph->getData(e).getParent();
                  if (fineGGraph->getData(e).getWeight() < best) {
                    best = fineGGraph->getData(e).getWeight();
                    b = e;
                  }
                  else if (fineGGraph->getData(e).getWeight() == best) {
                    if (e < b)
                      b = e;
                  }
              }

          }
          if (cells.size() > 0) {
              if (best < INT_MAX) {
                  auto nn = fineGGraph->getData(b).getParent();
                  for (auto e : cells) {
	            bag.push(e);
                    fineGGraph->getData(e).setMatched();
                    fineGGraph->getData(e).setParent(nn);
                    fineGGraph->getData(e).netnum = fineGGraph->getData(b).netnum;
                    for (int i = 0; i < dim; i++)
                       fineGGraph->getData(nn).EMvec[i] = 
                         (fineGGraph->getData(nn).EMvec[i] + fineGGraph->getData(e).EMvec[i])/2;
                  }
              }
                             
          }        
      },
         galois::steal(),galois::loopname("moreCoarse"));
      for (auto c : bag) {
        auto nn = fineGGraph->getData(c).getParent();
        int ww = weight[nn-fineGGraph->hedges];
        ww += fineGGraph->getData(c).getWeight();
        weight[nn-fineGGraph->hedges] = ww;
      }
}

// Coarsening phaseII
void coarsePhaseII(MetisGraph* graph,
                    int iter, std::vector<bool>& hedges, std::vector<unsigned> & weight) {

  GGraph* fineGGraph   = graph->getFinerGraph()->getGraph();
  typedef std::set<int> SecTy;
  typedef std::vector<GNode> VecTy;
  typedef galois::substrate::PerThreadStorage<SecTy> ThreadLocalData;
  ThreadLocalData edgesThreadLocal;
  typedef galois::substrate::PerThreadStorage<VecTy> ThreadLocalDataV;
  ThreadLocalDataV edgesThreadLocalV;
  std::string name = "CoarseningPhaseII";
  galois::GAccumulator<int> hhedges;
  galois::GAccumulator<int> hnode;
  moreCoarse(graph, iter, weight);

  galois::do_all( 
      galois::iterate((uint64_t)0, fineGGraph->hedges),
      [&](GNode item) {
        if (fineGGraph->getData(item).isMatched()) return;
        unsigned id = fineGGraph->getData(item).netnum;
        unsigned ids;
        int count = 0;
        for (auto c : fineGGraph->edges(item)) {
          auto dst = fineGGraph->getEdgeDst(c);
          auto& data = fineGGraph->getData(dst);
          if (data.isMatched()) {
            if (count == 0) {
              ids = data.getParent();
              count++;
            }
            else if (ids != data.getParent()) {
              count++;
              break;
            }
          }
          else { 
              count = 0;
              break;
          }
        }
        if (count == 1) {
            fineGGraph->getData(item).setMatched();
          
        }
        else {
           hedges[item] = true;
           fineGGraph->getData(item).setMatched();
         
        }

      },
      galois::loopname("count # Hyperedges"));

   //hedges += hhedges.reduce();
}

void parallelCreateEdges(MetisGraph* graph, GNodeBag& bag, std::vector<bool> hedges, std::vector<unsigned> weight) {

  /*GGraph* fineGGraph   = graph->getFinerGraph()->getGraph();
  GGraph* coarseGGraph = graph->getGraph();
  assert(fineGGraph != coarseGGraph);
  galois::GAccumulator<unsigned> hg;
  galois::do_all(
      galois::iterate((uint64_t)0, fineGGraph->hedges),
      [&](GNode n) {
          if (hedges[n])
              hg += 1;
      },
      galois::steal(), galois::loopname("number of hyperedges loop"));
  galois::do_all(
      galois::iterate((uint64_t)fineGGraph->hedges, fineGGraph->size()),
      [&](GNode ii) {
            if (!fineGGraph->getData(ii).isMatched()) { 
              bag.push(ii);
              fineGGraph->getData(ii).setMatched();
              fineGGraph->getData(ii).setParent(ii);
              weight[ii-fineGGraph->hedges] = fineGGraph->getData(ii).getWeight();
            }
      },
  galois::steal(), galois::loopname("noedgebag match"));
  unsigned nodes = std::distance(bag.begin(), bag.end());// + numnodes;
  int hnum = hg.reduce();
  unsigned newval = hnum;
  std::vector<unsigned> idmap(fineGGraph->hnodes);
  std::vector<unsigned> newrand(nodes);
  std::vector<unsigned> newWeight(nodes);
  galois::StatTimer Tloop("for loop");
  Tloop.start();
  std::vector<unsigned> v;
  for (auto n : bag) v.push_back(n);
  std::sort(v.begin(), v.end());
  for (auto n : v) {
    newrand[newval-hnum] = n;
    idmap[n-fineGGraph->hedges] = newval++;
    newWeight[idmap[n-fineGGraph->hedges]-hnum] = weight[n-fineGGraph->hedges];
  }
  std::vector<unsigned> EMvecNode(nodes); // EMveceddings for the nodes

  galois::do_all(
      galois::iterate((uint64_t)fineGGraph->hedges, fineGGraph->size()),
      [&](GNode n) {
        unsigned id = fineGGraph->getData(n).getParent();
           unsigned newi = idmap[id-fineGGraph->hedges];
           EMvecNode[newi] = id;
        fineGGraph->getData(n).setParent(idmap[id-fineGGraph->hedges]);
      },
      galois::steal(), galois::loopname("first loop"));
  Tloop.stop();
  //std::cout<<"total first loop "<<Tloop.get()<<"\n";

  uint32_t num_nodes_next = nodes + hnum;
  uint64_t num_edges_next; 
  std::vector<unsigned> EMvecHedge(hnum);
  galois::gstl::Vector<galois::PODResizeableArray<uint32_t>> edges_id(num_nodes_next);
  std::vector<unsigned> old_id(hnum);
  unsigned h_id = 0;
  //galois::StatTimer sloop("for loop II");
  //sloop.start();
  for (GNode n = 0; n < fineGGraph->hedges; n++) {
    if (hedges[n]) {
       old_id[h_id] = fineGGraph->getData(n).netnum;
       fineGGraph->getData(n).nodeid = h_id++;
    }
  }
  galois::do_all(galois::iterate((uint64_t)0, fineGGraph->hedges),
                [&](GNode n) {
                    if (!hedges[n]) return;
                        auto data = fineGGraph->getData(n, flag_no_lock);
                        unsigned id =  fineGGraph->getData(n).nodeid;
                        EMvecHedge[id] = n;
                    for (auto ii : fineGGraph->edges(n)) { 
                        GNode dst = fineGGraph->getEdgeDst(ii);
                        auto dst_data = fineGGraph->getData(dst, flag_no_lock);
                          unsigned pid = dst_data.getParent();
                          auto f = std::find(edges_id[id].begin(), edges_id[id].end(), pid);
                          if (f == edges_id[id].end()) {
                            edges_id[id].push_back(pid);
                          }
                    } // End edge loop
                }, galois::steal(),
                   galois::loopname("BuildGrah: Find edges"));


  std::vector<uint64_t> prefix_edges(num_nodes_next);
  galois::GAccumulator<uint64_t> num_edges_acc;
  galois::do_all(galois::iterate((uint32_t)0, num_nodes_next),
                [&](uint32_t c){
                  prefix_edges[c] = edges_id[c].size();
                  num_edges_acc += prefix_edges[c];
                }, galois::steal(),
                   galois::loopname("BuildGrah: Prefix sum"));

  num_edges_next = num_edges_acc.reduce();
  for (uint32_t c = 1; c < num_nodes_next; ++c) {
    prefix_edges[c] += prefix_edges[c - 1];
  }
  //galois::StatTimer TimerConstructFrom("Timer_Construct_From");
  //TimerConstructFrom.start();
  coarseGGraph->constructFrom(num_nodes_next, num_edges_next, prefix_edges, edges_id);
  //TimerConstructFrom.stop();
  //std::cout<<"graph cons time "<<TimerConstructFrom.get()<<"\n";
  coarseGGraph->hedges = hnum;
  coarseGGraph->hnodes = nodes;*/
  GGraph* fineGGraph   = graph->getFinerGraph()->getGraph();
  GGraph* coarseGGraph = graph->getGraph();
  assert(fineGGraph != coarseGGraph);
  galois::GAccumulator<unsigned> hg;
  galois::do_all(
      galois::iterate(size_t{0}, fineGGraph->hedges),
      [&](GNode n) {
        if (hedges[n])
          hg += 1;
      },
      galois::steal(), galois::loopname("number of hyperedges loop"));
  galois::do_all(
      galois::iterate(fineGGraph->hedges, fineGGraph->size()),
      [&](GNode ii) {
        if (!fineGGraph->getData(ii).isMatched()) {
          bag.push(ii);
          fineGGraph->getData(ii).setMatched();
          fineGGraph->getData(ii).setParent(ii);
          fineGGraph->getData(ii).netnum  = INT_MAX;
          weight[ii - fineGGraph->hedges] = fineGGraph->getData(ii).getWeight();
        }
      },
      galois::steal(), galois::loopname("noedgebag match"));
  unsigned hnum   = hg.reduce();
  unsigned nodes  = std::distance(bag.begin(), bag.end()); // + numnodes;
  unsigned newval = hnum;
  // std::map<unsigned, unsigned > idmap;
  std::vector<unsigned> idmap(fineGGraph->hnodes);
  std::vector<unsigned> newrand(nodes);
  std::vector<unsigned> newWeight(nodes);
  galois::StatTimer Tloop("for loop");
  Tloop.start();
  std::vector<unsigned> v;
  for (auto n : bag)
    v.push_back(n);
  // std::copy(bag.begin(), bag.end(), v.begin());
  std::sort(v.begin(), v.end());
  for (auto n : v) {
    newrand[newval - hnum]        = n;
    idmap[n - fineGGraph->hedges] = newval++;
    newWeight[idmap[n - fineGGraph->hedges] - hnum] =
        weight[n - fineGGraph->hedges];
  }
  // for (GNode n = fineGGraph->hedges; n < fineGGraph->size(); n++) {
  galois::do_all(
      galois::iterate(fineGGraph->hedges, fineGGraph->size()),
      [&](GNode n) {
        unsigned id = fineGGraph->getData(n).getParent();
        fineGGraph->getData(n).setParent(idmap[id - fineGGraph->hedges]);
      },
      galois::steal(), galois::loopname("first loop"));
  Tloop.stop();
  // std::cout<<"total first loop "<<Tloop.get()<<"\n";
  uint32_t num_nodes_next = nodes + hnum;
  uint64_t num_edges_next;
  galois::gstl::Vector<galois::PODResizeableArray<uint32_t>> edges_id(
      num_nodes_next);
  std::vector<std::vector<EdgeTy>> edges_data(num_nodes_next);
  std::vector<unsigned> old_id(hnum);
  unsigned h_id = 0;
  // galois::StatTimer sloop("for loop II");
  // sloop.start();
  std::vector<unsigned> EMvecHedge(hnum);
  for (GNode n = 0; n < fineGGraph->hedges; n++) {
    if (hedges[n]) {
      old_id[h_id]                  = fineGGraph->getData(n).netnum;
      EMvecHedge[h_id] = n;
      fineGGraph->getData(n).nodeid = h_id++;
    }
  }
  // sloop.stop();
  // std::cout<<"second for loop "<<sloop.get()<<"\n";
  galois::do_all(
      galois::iterate(size_t{0}, fineGGraph->hedges),
      [&](GNode n) {
        if (!hedges[n])
          return;
        auto data   = fineGGraph->getData(n, flag_no_lock);
        unsigned id = fineGGraph->getData(n).nodeid;

        for (auto ii : fineGGraph->edges(n)) {
          GNode dst     = fineGGraph->getEdgeDst(ii);
          auto dst_data = fineGGraph->getData(dst, flag_no_lock);
          unsigned pid  = dst_data.getParent();
          auto f = std::find(edges_id[id].begin(), edges_id[id].end(), pid);
          if (f == edges_id[id].end()) {
            edges_id[id].push_back(pid);
          }
        } // End edge loop
      },
      galois::steal(), galois::loopname("BuildGrah: Find edges"));

  std::vector<uint64_t> prefix_edges(num_nodes_next);
  galois::GAccumulator<uint64_t> num_edges_acc;
  galois::do_all(
      galois::iterate(uint32_t{0}, num_nodes_next),
      [&](uint32_t c) {
        prefix_edges[c] = edges_id[c].size();
        num_edges_acc += prefix_edges[c];
      },
      galois::steal(), galois::loopname("BuildGrah: Prefix sum"));

  num_edges_next = num_edges_acc.reduce();
  for (uint32_t c = 1; c < num_nodes_next; ++c) {
    prefix_edges[c] += prefix_edges[c - 1];
  }
  // galois::StatTimer TimerConstructFrom("Timer_Construct_From");
  // TimerConstructFrom.start();
  coarseGGraph->constructFrom(num_nodes_next, num_edges_next, prefix_edges,
                              edges_id, edges_data);
  // TimerConstructFrom.stop();
  // std::cout<<"graph cons time "<<TimerConstructFrom.get()<<"\n";
  coarseGGraph->hedges = hnum;
  coarseGGraph->hnodes = nodes;
  galois::do_all(
      galois::iterate(*coarseGGraph),
      [&](GNode ii) {
        if (ii < hnum) {
          coarseGGraph->getData(ii).netval = INT_MAX;
          coarseGGraph->getData(ii).netnum = old_id[ii]; 
          coarseGGraph->getData(ii).EMvec = fineGGraph->getData(EMvecHedge[ii]).EMvec;

        } 
        else {
            coarseGGraph->getData(ii).netval = INT_MIN;
            coarseGGraph->getData(ii).netnum = INT_MAX;
            coarseGGraph->getData(ii).netrand = INT_MAX;
            coarseGGraph->getData(ii).nodeid = ii;//fineGGraph->getData(id).nodeid;
            coarseGGraph->getData(ii).setWeight(newWeight[ii-coarseGGraph->hedges]);
            coarseGGraph->getData(ii).EMvec = fineGGraph->getData(newrand[ii-hnum]).EMvec;
        }
      },
      galois::steal(), galois::loopname("noedgebag match"));
 // sloop.stop();
 // std::cout<<"building CSR "<<sloop.get()<<"\n";
}


void findMatching(MetisGraph* coarseMetisGraph,
                      
                       int iter) {
  MetisGraph* fineMetisGraph = coarseMetisGraph->getFinerGraph();
  GNodeBag nodes;
  int sz = coarseMetisGraph->getFinerGraph()->getGraph()->hedges;
  std::vector<bool> hedges(sz, false);
  //for (int i = 0; i < sz; i++) hedges[i] = false; 
  std::vector<unsigned> weight(fineMetisGraph->getGraph()->hnodes);
  
        parallelHMatchAndCreateNodes(coarseMetisGraph,
                                            iter, nodes, hedges, weight);
       coarsePhaseII(coarseMetisGraph, iter, hedges, weight);
       parallelCreateEdges(coarseMetisGraph, nodes, hedges, weight);
}

MetisGraph* coarsenOnce(MetisGraph* fineMetisGraph, 
                         int iter) {
  MetisGraph* coarseMetisGraph = new MetisGraph(fineMetisGraph);
  findMatching(coarseMetisGraph, iter);
  return coarseMetisGraph;
}

} // namespace

MetisGraph* coarsen(MetisGraph* fineMetisGraph, unsigned coarsenTo) {

  MetisGraph* coarseGraph = fineMetisGraph;
  unsigned size           = fineMetisGraph->getGraph()->hnodes;//, fineMetisGraph->getGraph()->cellList().end());
  unsigned hedgeSize = 0;
  const float ratio = 55.0 / 45.0;  // change if needed
  const float tol = std::max(ratio, 1 - ratio) - 1;
  const int hi = (1 + tol) * size / (2 + tol);
  const int lo = size - hi;
  LIMIT = hi / 4;
  int totw = 0;
  
  //std::cout<<"inital weight is "<<totw<<"\n";
  unsigned Size = size;
  unsigned iterNum        = 0;
  unsigned newSize = size;
  while (Size > 100) { 
    //if (iterNum > 25) break;
    if (Size - newSize <= 20 && iterNum > 2) break; //final
     newSize = coarseGraph->getGraph()->hnodes;
     coarseGraph      = coarsenOnce(coarseGraph, iterNum);
     Size = coarseGraph->getGraph()->hnodes;
     hedgeSize = coarseGraph->getGraph()->hedges; 
     std::cout<<"SIZE IS "<<coarseGraph->getGraph()->hnodes<<" and net is "<<hedgeSize<<"\n";
     if (hedgeSize < 1000) return coarseGraph->getFinerGraph();
     
    ++iterNum;
    
  }
  return coarseGraph;
}
