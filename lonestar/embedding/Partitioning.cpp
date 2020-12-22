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
  int waccum;
  galois::GAccumulator<unsigned int> accumZ;
  GNodeBag nodelist;
  galois::do_all(
      galois::iterate(g->hedges, g->size()),
      [&](GNode item) {
        accum += g->getData(item).getWeight();
        g->getData(item, galois::MethodFlag::UNPROTECTED).initRefine(1, true);
        g->getData(item, galois::MethodFlag::UNPROTECTED).initPartition();
      },
      galois::loopname("initPart"));

  galois::do_all(
      galois::iterate(size_t{0}, g->hedges),
      [&](GNode item) {
        for (auto c : g->edges(item)) {
          auto n = g->getEdgeDst(c);
          g->getData(n).setPart(0);
        }
      },
      galois::loopname("initones"));
  GNodeBag nodelistoz;
  galois::do_all(
      galois::iterate(g->hedges, g->size()),
      [&](GNode item) {
        if (g->getData(item).getPart() == 0) {
          accumZ += g->getData(item).getWeight();
          nodelist.push(item);
        } else
          nodelistoz.push(item);
      },
      galois::loopname("initones"));
  unsigned newSize = accum.reduce();
  waccum           = accum.reduce() - accumZ.reduce();
  // unsigned targetWeight = accum.reduce() / 2;
  unsigned kvalue        = (K + 1) / 2;
  unsigned targetWeight0 = (accum.reduce() * kvalue) / K;
  unsigned targetWeight1 = accum.reduce() - targetWeight0;

  if (static_cast<long>(accumZ.reduce()) > waccum) {
    int gain = waccum;
    // initGain(*g);
    while (1) {
      initGain(*g);
      std::vector<GNode> nodeListz;
      GNodeBag nodelistz;
      galois::do_all(
          galois::iterate(nodelist),
          [&](GNode node) {
            unsigned pp = g->getData(node).getPart();
            if (pp == 0) {
              nodelistz.push(node);
            }
          },
          galois::loopname("while"));

      for (auto c : nodelistz)
        nodeListz.push_back(c);
      std::sort(
          nodeListz.begin(), nodeListz.end(), [&g](GNode& lpw, GNode& rpw) {
            if (fabs((float)((g->getData(lpw).getGain()) *
                             (1.0f / g->getData(lpw).getWeight())) -
                     (float)((g->getData(rpw).getGain()) *
                             (1.0f / g->getData(rpw).getWeight()))) < 0.00001f)
              return (float)g->getData(lpw).nodeid <
                     (float)g->getData(rpw).nodeid;
            return (float)((g->getData(lpw).getGain()) *
                           (1.0f / g->getData(lpw).getWeight())) >
                   (float)((g->getData(rpw).getGain()) *
                           (1.0f / g->getData(rpw).getWeight()));
          });
      int i = 0;
      for (auto zz : nodeListz) {
        // auto zz = *nodeListz.begin();
        g->getData(zz).setPart(1);
        gain += g->getData(zz).getWeight();
        // std::cout<<" weight "<<g->getData(zz).getWeight()<<"\n";

        i++;
        if (gain >= static_cast<long>(targetWeight1))
          break;
        if (i > sqrt(newSize))
          break;
      }

      if (gain >= static_cast<long>(targetWeight1))
        break;
      // updateGain(*g,zz);
    }

  } else {

    int gain = accumZ.reduce();
    // std::cout<<"gain is "<<gain<<"\n";
    // initGain(*g);
    while (1) {
      initGain(*g);
      std::vector<GNode> nodeListz;
      GNodeBag nodelistz;
      galois::do_all(
          galois::iterate(nodelistoz),
          [&](GNode node) {
            // for (auto node : nodelist) {
            unsigned pp = g->getData(node).getPart();
            if (pp == 1) {
              nodelistz.push(node);
            }
          },
          galois::loopname("while"));
      for (auto c : nodelistz)
        nodeListz.push_back(c);

      std::sort(
          nodeListz.begin(), nodeListz.end(), [&g](GNode& lpw, GNode& rpw) {
            if (fabs((float)((g->getData(lpw).getGain()) *
                             (1.0f / g->getData(lpw).getWeight())) -
                     (float)((g->getData(rpw).getGain()) *
                             (1.0f / g->getData(rpw).getWeight()))) < 0.00001f)
              return (float)g->getData(lpw).nodeid <
                     (float)g->getData(rpw).nodeid;
            return (float)((g->getData(lpw).getGain()) *
                           (1.0f / g->getData(lpw).getWeight())) >
                   (float)((g->getData(rpw).getGain()) *
                           (1.0f / g->getData(rpw).getWeight()));
          });

      int i = 0;
      for (auto zz : nodeListz) {
        // auto zz = *nodeListz.begin();
        g->getData(zz).setPart(0);
        gain += g->getData(zz).getWeight();

        i++;
        if (gain >= static_cast<long>(targetWeight0))
          break;
        if (i > sqrt(newSize))
          break;
      }

      if (gain >= static_cast<long>(targetWeight0))
        break;

      // updateGain(*g,zz);
    }
  }
  std::ofstream out("train.txt");
  std::ofstream testfile("test.txt");
  std::ofstream validate("validate.txt");
  std::map<int, std::set<int> > output;
  galois::GAccumulator<unsigned> totalcount;
      /*galois::do_all(
          galois::iterate(size_t{0}, g->size()),
          [&](GNode node) {
            for (auto hnode : g->edges(node)) {
              totalcount += 1;
              if (node < g->hedges) return;
              
              auto h = g->getEdgeDst(hnode);
              for (auto mnode : g->edges(h)) {
                auto hv = g->getEdgeDst(mnode);
                if (hv != node)
                  totalcount += 1;
              }
            }
          },
          galois::loopname("while"));*/
  
  unsigned total = 99140;//totalcount.reduce();

 std::cout<<total<<"\n";

  for (auto v = g->hedges; v < g->size(); v++ ) {
    for (auto hnode : g->edges(v)) {
      auto h = g->getEdgeDst(hnode);
      for (auto node : g->edges(h)) {
        auto tmp = g->getEdgeDst(node);
        if (tmp != v)
          output[v].insert(tmp);
      }
    }
  }
 /*unsigned testSeed = total * 0.1;
 std::vector<unsigned> test(total,0); 

 std::random_device rd; // obtain a random number from hardware
 std::mt19937 gen(rd()); // seed the generator
 std::uniform_int_distribution<> distr(0, total); // define the range
  unsigned validSeed = total * 0.09;
  std::cout<<testSeed+validSeed<<"  val+test\n";
  for (unsigned i = 0; i < testSeed + validSeed; i++) {
    unsigned t = distr(gen);
    while(test[t])
      t = distr(gen);
    test[t] = 1;
  }
  unsigned blah = 0;
  for (int i = 0; i < test.size(); i++)
    if (test[i]) blah++;
  std::cout<<" real ones "<<blah<<"\n";
  unsigned counter = 0;
  unsigned val = 0;
  std::cout<<"size of map "<<test.size()<<"\n";
  for (auto m : output) {
    auto res = m.second;
    for (auto a : res) {
      if (test[counter]) {
        if (val >= validSeed)
          testfile << m.first << "\t" <<"1\t"<< a <<"\n";
        else {
          if (counter%3==0)
          testfile << m.first << "\t" <<"1\t"<< a <<"\n";
          else {
            val++;
            validate << m.first << "\t" <<"1\t"<< a <<"\n";
          }
        }
      }
      else
        out << m.first << "\t" <<"1\t"<< a <<"\n";
      counter++;
    }
  }

  for (auto v = 0; v < g->size(); v++ ) {
    for (auto hnode : g->edges(v)) {
      auto h = g->getEdgeDst(hnode);
      if (test[counter]) {
        if (val >= validSeed)
        testfile << v << "\t2\t"<< h <<"\n";
        else {
          val++;
          validate << v << "\t2\t"<< h <<"\n";
        }
      }
      else
        out << v << "\t2\t"<< h <<"\n";
      counter++;
    }
  }

  std::cout<<"the real counter "<<counter<<"\n";
  out.close();
  validate.close();
  testfile.close();*/
}
