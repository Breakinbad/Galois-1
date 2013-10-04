#ifndef APPS_CONNECTEDCOMPONENTS_LIGRAALGO_H
#define APPS_CONNECTEDCOMPONENTS_LIGRAALGO_H

#include "Galois/DomainSpecificExecutors.h"
#include "Galois/Graph/OCGraph.h"
#include "Galois/Graph/LCGraph.h"
#include "Galois/Graph/GraphNodeBag.h"
#include "llvm/Support/CommandLine.h"

#include <boost/mpl/if.hpp>

extern llvm::cl::opt<unsigned int> memoryLimit;

template<typename Graph>
void readInOutGraph(Graph& graph);

template<bool UseGraphChi>
struct LigraAlgo: public Galois::LigraGraphChi::ChooseExecutor<UseGraphChi>  {
  struct LNode {
    typedef unsigned int component_type;
    unsigned int id;
    unsigned int comp;
    unsigned int oldComp;
    
    component_type component() { return comp; }
    bool isRep() { return id == comp; }
  };

  typedef typename Galois::Graph::LC_CSR_Graph<LNode,void>
    ::template with_no_lockable<true>::type
    ::template with_numa_alloc<true>::type InnerGraph;
  typedef typename boost::mpl::if_c<UseGraphChi,
          Galois::Graph::OCImmutableEdgeGraph<LNode,void>,
          Galois::Graph::LC_InOut_Graph<InnerGraph>>::type
          Graph;
  typedef typename Graph::GraphNode GNode;

  template<typename Bag>
  struct Initialize {
    Graph& graph;
    Bag& bag;

    Initialize(Graph& g, Bag& b): graph(g), bag(b) { }
    void operator()(GNode n) {
      LNode& data = graph.getData(n, Galois::MethodFlag::NONE);
      data.comp = data.id;
      bag.push(n, graph.sizeEdges() / graph.size());
    }
  };

  struct Copy {
    Graph& graph;

    Copy(Graph& g): graph(g) { }
    void operator()(size_t n, Galois::UserContext<size_t>&) {
      (*this)(n);
    }
    void operator()(size_t id) {
      GNode n = graph.nodeFromId(id);
      LNode& data = graph.getData(n, Galois::MethodFlag::NONE);
      data.oldComp = data.comp;
    }
  };

  struct EdgeOperator {
    template<typename GTy>
    bool cond(GTy& graph, typename GTy::GraphNode) { return true; }
    
    template<typename GTy>
    bool operator()(GTy& graph, typename GTy::GraphNode src, typename GTy::GraphNode dst, typename GTy::edge_data_reference) {
      LNode& sdata = graph.getData(src, Galois::MethodFlag::NONE);
      LNode& ddata = graph.getData(dst, Galois::MethodFlag::NONE);

      typename LNode::component_type orig = ddata.comp;
      while (true) {
        typename LNode::component_type old = ddata.comp;
        if (old <= sdata.comp)
          return false;
        if (__sync_bool_compare_and_swap(&ddata.comp, old, sdata.comp)) {
          return orig == ddata.oldComp;
        }
      }
      return false;
    }
  };

  void readGraph(Graph& graph) {
    readInOutGraph(graph);
    this->checkIfInMemoryGraph(graph, memoryLimit);
  }

  void operator()(Graph& graph) {
    typedef Galois::WorkList::dChunkedFIFO<256> WL;
    typedef Galois::GraphNodeBagPair<> BagPair;
    BagPair bags(graph.size());

    Galois::do_all_local(graph, Initialize<typename BagPair::bag_type>(graph, bags.next()));
    while (!bags.next().empty()) {
      bags.swap();
      Galois::for_each_local<WL>(bags.cur(), Copy(graph));
      this->outEdgeMap(memoryLimit, graph, EdgeOperator(), bags.cur(), bags.next(), false);
    } 
  }
};

#endif