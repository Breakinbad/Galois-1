/*
 * This file belongs to the Galois project, a C++ library for exploiting
 * parallelism. The code is being released under the terms of the 3-Clause BSD
 * License (a copy is located in LICENSE.txt at the top-level directory).
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

#include "galois/Galois.h"
#include "galois/Timer.h"
#include "bipart.h"
#include <set>
#include "galois/Galois.h"
#include "galois/AtomicHelpers.h"
#include <map>
#include <set>
#include <cstdlib>
#include <iostream>
#include <stack>
#include <climits>
#include <array>
#include <fstream>
#include <random>

namespace {
// final
int hash(unsigned val) {
  unsigned long int seed = val * 1103515245 + 12345;
  return ((unsigned)(seed / 65536) % 32768);
}

__attribute__((unused)) int cut(GGraph& g) {

  GNodeBag bag;
  galois::do_all(
      galois::iterate(g),
      [&](GNode n) {
        if (g.hedges <= n)
          return;
        for (auto cell : g.edges(n)) {
          auto c   = g.getEdgeDst(cell);
          int part = g.getData(c).getPart();
          for (auto x : g.edges(n)) {
            auto cc   = g.getEdgeDst(x);
            int partc = g.getData(cc).getPart();
            if (partc != part) {
              bag.push(n);
              return;
            }
          }
        }
      },
      galois::loopname("cutsize"));
  return std::distance(bag.begin(), bag.end());
}

void initGain(GGraph& g) {
  galois::do_all(
      galois::iterate(g),
      [&](GNode n) {
        if (n < g.hedges)
          return;
        g.getData(n).FS.store(0);
        g.getData(n).TE.store(0);
      },
      galois::loopname("firstinit"));

  typedef std::map<GNode, int> mapTy;
  typedef galois::substrate::PerThreadStorage<mapTy> ThreadLocalData;
  ThreadLocalData edgesThreadLocal;
  galois::do_all(
      galois::iterate(g),
      [&](GNode n) {
        if (g.hedges <= n)
          return;
        int p1 = 0;
        int p2 = 0;
        for (auto x : g.edges(n)) {
          auto cc  = g.getEdgeDst(x);
          int part = g.getData(cc).getPart();
          if (part == 0)
            p1++;
          else
            p2++;
          if (p1 > 1 && p2 > 1)
            break;
        }
        if (!(p1 > 1 && p2 > 1) && (p1 + p2 > 1)) {
          for (auto x : g.edges(n)) {
            auto cc  = g.getEdgeDst(x);
            int part = g.getData(cc).getPart();
            int nodep;
            if (part == 0)
              nodep = p1;
            else
              nodep = p2;
            if (nodep == 1) {
              galois::atomicAdd(g.getData(cc).FS, 1);
            }
            if (nodep == (p1 + p2)) {
              galois::atomicAdd(g.getData(cc).TE, 1);
            }
          }
        }
      },
      galois::steal(), galois::loopname("initGainsPart"));
}

} // namespace

// Final
void partition(MetisGraph* mcg, unsigned K) {
  GGraph* g = mcg->getGraph();
  galois::GAccumulator<unsigned int> accum;
  galois::do_all(
      galois::iterate((uint64_t)0, g->size()),
      [&](GNode item) {
        unsigned size = std::distance(g->edges(item).begin(), g->edges(item).end());
        accum += size;
      },
      galois::loopname("initPart"));

  unsigned avg = accum.reduce() / g->size();
  //waccum           = accum.reduce() - accumZ.reduce();
  // unsigned targetWeight = accum.reduce() / 2;
  std::cout<<"in Part\n";
  std::ofstream out("train.txt");
  std::ofstream testfile("test.txt");
  std::ofstream validate("validate.txt");
  std::map<int, std::set<int> > output;
  std::map<int, int> tout;
  galois::GAccumulator<unsigned> totalcount;
  unsigned counter = 0;
  for (auto v = 0; v < g->size(); v++ ) {
    bool flag = false;
    for (auto hnode : g->edges(v)) {
      if (output[v].size() > avg) break;
      auto h = g->getEdgeDst(hnode);
      if (!flag) { 
        tout[v] = h;
        flag = true;
      }
      else
        output[v].insert(h);
      for (auto node : g->edges(h)) {
      if (output[v].size() > avg) break;
        auto tmp = g->getEdgeDst(node);
        if (tmp != v)
          output[v].insert(tmp);
      }
    }
    counter += output[v].size();
  }
  unsigned ccs = counter + tout.size();
  std::cout<<"after first loop\n";


 std::cout<<"the real counter "<<ccs<<"\n";
/**************/

 for (auto v : tout) {

   out<<v.first<<"\t2\t"<<v.second<<"\n";
 }
 unsigned total = counter;
 unsigned testSeed = ccs * 0.1;
 std::map<unsigned, unsigned> test; 
 
 std::random_device rd; // obtain a random number from hardware
 std::mt19937 gen(rd()); // seed the generator
 std::uniform_int_distribution<> distr(0, total); // define the range
  unsigned validSeed = ccs * 0.09;
  std::cout<<testSeed+validSeed<<"  val+test\n";
  for (unsigned i = 0; i < testSeed + validSeed; i++) {
    unsigned t = distr(gen);
    while(test[t])
      t = distr(gen);
    test[t] = 1;
  }
  counter = 0;
  unsigned val = 0;
  std::cout<<"size of map "<<test.size()<<"\n";
  for (auto m : output) {
    auto res = m.second;
    unsigned cc = 0;
    for (auto a : res) {
      if (cc > avg) break;
      cc++;
      if (test[counter]) {
        if (val >= validSeed) {
          if (m.first < g->hedges) {
            if (a < g->hedges)
              testfile << m.first << "\t" <<"1\t"<< a <<"\n";
            else 
              testfile << m.first << "\t" <<"2\t"<< a <<"\n";
          }
          else {
            if (a >= g->hedges)
              testfile << m.first << "\t" <<"1\t"<< a <<"\n";
            else
              testfile << m.first << "\t" <<"2\t"<< a <<"\n";
          }
        }
        else {
          if (counter%3==0) {
            if (m.first < g->hedges) {
              if (a < g->hedges)
                testfile << m.first << "\t" <<"1\t"<< a <<"\n";
              else 
                testfile << m.first << "\t" <<"2\t"<< a <<"\n";
            }
            else {
              if (a >= g->hedges)
                testfile << m.first << "\t" <<"1\t"<< a <<"\n";
              else
                testfile << m.first << "\t" <<"2\t"<< a <<"\n";
            }
          }
          else {
            val++;
            if (m.first < g->hedges) {
              if (a < g->hedges)
                validate << m.first << "\t" <<"1\t"<< a <<"\n";
              else 
                validate << m.first << "\t" <<"2\t"<< a <<"\n";
            }
            else {
              if (a >= g->hedges)
                validate << m.first << "\t" <<"1\t"<< a <<"\n";
              else
                validate << m.first << "\t" <<"2\t"<< a <<"\n";
            }
          }
        }
      }
      else
      {
            if (m.first < g->hedges) {
              if (a < g->hedges)
                out << m.first << "\t" <<"1\t"<< a <<"\n";
              else 
                out << m.first << "\t" <<"2\t"<< a <<"\n";
            }
            else {
              if (a >= g->hedges)
                out << m.first << "\t" <<"1\t"<< a <<"\n";
              else
                out << m.first << "\t" <<"2\t"<< a <<"\n";
            }
      }
      counter++;
    }
  }


  /**************/

  out.close();
  validate.close();
  testfile.close();
}
