/** Per Thread workLists-*- C++ -*-
 * @file
 * @section License
 *
 * Galois, a framework to exploit amorphous data-parallelism in irregular
 * programs.
 *
 * Copyright (C) 2011, The University of Texas at Austin. All rights reserved.
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
 *
 * @section Description
 *
 * a thread local stl container for each thread
 *
 * @author <ahassaan@ices.utexas.edu>
 */
#ifndef GALOIS_RUNTIME_PER_THREAD_WORK_LIST_H_
#define GALOIS_RUNTIME_PER_THREAD_WORK_LIST_H_

#include <vector>
#include <deque>
#include <list>
#include <set>
#include <limits>
#include <iterator>
#include <iostream>

#include <cstdio>

#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include "Galois/Threads.h"
#include "Galois/TwoLevelIterator.h"
#include "Galois/Runtime/PerThreadStorage.h"
#include "Galois/Runtime/ThreadPool.h"
#include "Galois/Runtime/mm/Mem.h"
#include "Galois/Runtime/ll/gio.h"

namespace Galois {
namespace Runtime {

namespace HIDDEN {

enum GlobalPos {
  GLOBAL_BEGIN, GLOBAL_END
};

#define ADAPTOR_BASED_OUTER_ITER

// XXX: use a combination of boost::transform_iterator and
// boost::counting_iterator to implement the following OuterPerThreadWLIter
#ifdef ADAPTOR_BASED_OUTER_ITER

template <typename PerThrdWL>
struct WLindexer: 
  public std::unary_function<unsigned, typename PerThrdWL::Cont_ty&> 
{

  typedef typename PerThrdWL::Cont_ty Ret_ty;

  PerThrdWL* wl;

  WLindexer (): wl (NULL) {}

  WLindexer (PerThrdWL& _wl): wl(&_wl) {}

  Ret_ty& operator () (unsigned i) const {
    assert (wl != NULL);
    assert (i < wl->numRows ());
    return const_cast<Ret_ty&> (wl->get (i));
  }
};

template <typename PerThrdWL>
struct TypeFactory {
  typedef boost::transform_iterator<WLindexer<PerThrdWL>, boost::counting_iterator<unsigned> > OuterIter;

  typedef std::reverse_iterator<OuterIter> RvrsOuterIter;
};


template <typename PerThrdWL>
typename TypeFactory<PerThrdWL>::OuterIter make_outer_begin (PerThrdWL& wl) {

  typename TypeFactory<PerThrdWL>::OuterIter outer = boost::make_transform_iterator (boost::counting_iterator<unsigned> (0), WLindexer<PerThrdWL> (wl));

  outer->begin ();
  (*outer).begin ();

  return boost::make_transform_iterator (
      boost::counting_iterator<unsigned> (0), WLindexer<PerThrdWL> (wl));
}

template <typename PerThrdWL>
typename TypeFactory<PerThrdWL>::OuterIter make_outer_end (PerThrdWL& wl) {
  return boost::make_transform_iterator (
      boost::counting_iterator<unsigned> (wl.numRows ()), WLindexer<PerThrdWL> (wl));
}

template <typename PerThrdWL>
typename TypeFactory<PerThrdWL>::RvrsOuterIter make_outer_rbegin (PerThrdWL& wl) {
  return typename TypeFactory<PerThrdWL>::RvrsOuterIter (make_outer_end (wl));
}

template <typename PerThrdWL>
typename TypeFactory<PerThrdWL>::RvrsOuterIter make_outer_rend (PerThrdWL& wl) {
  return typename TypeFactory<PerThrdWL>::RvrsOuterIter (make_outer_begin (wl));
}

#else

template <typename PerThrdWL>
class OuterPerThreadWLIter: public std::iterator<std::random_access_iterator_tag, typename PerThrdWL::Cont_ty> {

  typedef typename PerThrdWL::Cont_ty Cont_ty;
  typedef std::iterator<std::random_access_iterator_tag, Cont_ty> Super_ty;
  typedef typename Super_ty::difference_type Diff_ty;

  PerThrdWL* workList;
  // using Diff_ty due to reverse iterator, whose 
  // end is -1, and,  begin is numRows - 1
  Diff_ty row;

  void assertInRange () const {
    assert ((row >= 0) && (row < workList->numRows ()));
  }

  Cont_ty& getWL () {
    assertInRange ();
    return (*workList)[row];
  }

  const Cont_ty& getWL () const {
    assertInRange ();
    return (*workList)[row];
  }


public:

  OuterPerThreadWLIter (): Super_ty (), workList (NULL), row (0) {}

  OuterPerThreadWLIter (PerThrdWL& wl, const GlobalPos& pos)
    : Super_ty (), workList (&wl), row (0) {

    switch (pos) {
      case GLOBAL_BEGIN:
        row = 0;
        break;
      case GLOBAL_END:
        row = wl.numRows ();
        break;
      default:
        std::abort ();
    }
  }

  typename Super_ty::reference operator * () { return getWL (); }

  const typename Super_ty::reference operator * () const { return getWL (); }

  typename Super_ty::pointer operator -> () { return &(getWL ()); }

  const typename Super_ty::value_type* operator -> () const { return &(getWL ()); }

  OuterPerThreadWLIter& operator ++ () {
    ++row;
    return *this;
  }

  OuterPerThreadWLIter operator ++ (int) {
    OuterPerThreadWLIter tmp (*this);
    operator ++ ();
    return tmp;
  }

  OuterPerThreadWLIter& operator -- () {
    --row;
    return *this;
  }

  OuterPerThreadWLIter operator -- (int) {
    OuterPerThreadWLIter tmp (*this);
    operator -- ();
    return tmp;
  }

  OuterPerThreadWLIter& operator += (Diff_ty d) {
    row = unsigned (Diff_ty (row) + d);
    return *this;
  }

  OuterPerThreadWLIter& operator -= (Diff_ty d) {
    row = unsigned (Diff_ty (row) - d);
    return *this;
  }

  friend OuterPerThreadWLIter operator + (const OuterPerThreadWLIter& it, Diff_ty d) {
    OuterPerThreadWLIter tmp (it);
    tmp += d;
    return tmp;
  }

  friend OuterPerThreadWLIter operator + (Diff_ty d, const OuterPerThreadWLIter& it) {
    return it + d;
  }

  friend OuterPerThreadWLIter operator - (const OuterPerThreadWLIter& it, Diff_ty d) {
    OuterPerThreadWLIter tmp (it);
    tmp -= d;
    return tmp;
  }

  friend Diff_ty operator - (const OuterPerThreadWLIter& left, const OuterPerThreadWLIter& right) {
    return Diff_ty (left.row) - Diff_ty (right.row);
  }

  typename Super_ty::reference operator [] (Diff_ty d) {
    return *((*this) + d);
  }

  friend bool operator == (const OuterPerThreadWLIter& left, const OuterPerThreadWLIter& right) {

    assert (left.workList == right.workList);
    return (left.row == right.row);
  }

  friend bool operator != (const OuterPerThreadWLIter& left, const OuterPerThreadWLIter& right) {
    return !(left == right);
  }

  friend bool operator < (const OuterPerThreadWLIter& left, const OuterPerThreadWLIter& right) {

    assert (left.workList == right.workList);

    return (left.row < right.row);
  }

  friend bool operator <= (const OuterPerThreadWLIter& left, const OuterPerThreadWLIter& right) {

    return (left == right) || (left < right);
  }

  friend bool operator > (const OuterPerThreadWLIter& left, const OuterPerThreadWLIter& right) {

    return !(left <= right);
  }

  friend bool operator >= (const OuterPerThreadWLIter& left, const OuterPerThreadWLIter& right) {

    return !(left < right);
  }

};


template <typename PerThrdWL>
OuterPerThreadWLIter<PerThrdWL> make_outer_begin (PerThrdWL& wl) {
  return OuterPerThreadWLIter<PerThrdWL> (wl, GLOBAL_BEGIN);
}

template <typename PerThrdWL>
OuterPerThreadWLIter<PerThrdWL> make_outer_end (PerThrdWL& wl) {
  return OuterPerThreadWLIter<PerThrdWL> (wl, GLOBAL_END);
}

template <typename PerThrdWL>
std::reverse_iterator<OuterPerThreadWLIter<PerThrdWL> > 
  make_outer_rbegin (PerThrdWL& wl) {
  typedef std::reverse_iterator<OuterPerThreadWLIter<PerThrdWL> > Ret_ty;
  return Ret_ty (make_outer_end (wl));
}

template <typename PerThrdWL>
std::reverse_iterator<OuterPerThreadWLIter<PerThrdWL> > 
  make_outer_rend (PerThrdWL& wl) {
  typedef std::reverse_iterator<OuterPerThreadWLIter<PerThrdWL> > Ret_ty;
  return Ret_ty (make_outer_begin (wl));
}

#endif

} // end namespace HIDDEN


template <typename Cont_tp> 
class PerThreadWorkList {

public:
  typedef Cont_tp Cont_ty;
  typedef typename Cont_ty::value_type value_type;
  typedef typename Cont_ty::reference reference;
  typedef typename Cont_ty::pointer pointer;
  typedef typename Cont_ty::size_type size_type;

  typedef typename Cont_ty::iterator local_iterator;
  typedef typename Cont_ty::const_iterator local_const_iterator;
  typedef typename Cont_ty::reverse_iterator local_reverse_iterator;
  typedef typename Cont_ty::const_reverse_iterator local_const_reverse_iterator;

  typedef PerThreadWorkList This_ty;

#ifdef ADAPTOR_BASED_OUTER_ITER

  typedef typename HIDDEN::TypeFactory<This_ty>::OuterIter OuterIter;
  typedef typename HIDDEN::TypeFactory<This_ty>::RvrsOuterIter RvrsOuterIter;

#else

  typedef HIDDEN::OuterPerThreadWLIter<This_ty> OuterIter;
  typedef std::reverse_iterator<OuterIter> RvrsOuterIter;

#endif

  typedef typename Galois::ChooseStlTwoLevelIterator<OuterIter, typename Cont_ty::iterator>::type global_iterator;

  typedef typename Galois::ChooseStlTwoLevelIterator<OuterIter, typename Cont_ty::const_iterator>::type global_const_iterator;

  typedef typename Galois::ChooseStlTwoLevelIterator<RvrsOuterIter, typename Cont_ty::reverse_iterator>::type global_reverse_iterator;

  typedef typename Galois::ChooseStlTwoLevelIterator<RvrsOuterIter, typename Cont_ty::const_reverse_iterator>::type global_const_reverse_iterator;



private:

  // XXX: for testing only

#if 0
  struct FakePTS {
    std::vector<Cont_ty*> v;

    FakePTS () { 
      v.resize (size ());
    }

    Cont_ty** getLocal () const {
      return getRemote (Galois::Runtime::LL::getTID ());
    }

    Cont_ty** getRemote (size_t i) const {
      assert (i < v.size ());
      return const_cast<Cont_ty**> (&v[i]);
    }

    size_t size () const { return Galois::Runtime::LL::getMaxThreads(); }

  };
#endif



  // typedef FakePTS PerThrdCont_ty;
  typedef Galois::Runtime::PerThreadStorage<Cont_ty*> PerThrdCont_ty;
  PerThrdCont_ty perThrdCont;

  void destroy () {
    for (unsigned i = 0; i < perThrdCont.size (); ++i) {
      delete *perThrdCont.getRemote (i);
      *perThrdCont.getRemote (i) = NULL;
    }
  }

protected:
  PerThreadWorkList (): perThrdCont () {
    for (unsigned i = 0; i < perThrdCont.size (); ++i) {
      *perThrdCont.getRemote (i) = NULL;
    }
  }

  void init (const Cont_ty& cont) {
    for (unsigned i = 0; i < perThrdCont.size (); ++i) {
      *perThrdCont.getRemote (i) = new Cont_ty (cont);
    }
  }


  ~PerThreadWorkList () { 
    destroy ();
  }


public:
  unsigned numRows () const { return perThrdCont.size (); }

  Cont_ty& get () { return **(perThrdCont.getLocal ()); }

  const Cont_ty& get () const { return **(perThrdCont.getLocal ()); }

  Cont_ty& get (unsigned i) { return **(perThrdCont.getRemote (i)); }

  const Cont_ty& get (unsigned i) const { return **(perThrdCont.getRemote (i)); }

  Cont_ty& operator [] (unsigned i) { return get (i); }

  const Cont_ty& operator [] (unsigned i) const { return get (i); }


  global_iterator begin_all () { 
    return Galois::stl_two_level_begin (
        HIDDEN::make_outer_begin (*this), HIDDEN::make_outer_end (*this)); 
  }

  global_iterator end_all () { 
    return Galois::stl_two_level_end (
        HIDDEN::make_outer_end (*this), HIDDEN::make_outer_end (*this)); 
  }

  global_const_iterator begin_all () const { 
    return Galois::stl_two_level_cbegin (
        HIDDEN::make_outer_begin (*this), HIDDEN::make_outer_end (*this));
  }

  global_const_iterator end_all () const { 
    return Galois::stl_two_level_cend (
        HIDDEN::make_outer_begin (*this), HIDDEN::make_outer_end (*this));
  }

  global_reverse_iterator rbegin_all () { 
    return Galois::stl_two_level_rbegin (
        HIDDEN::make_outer_rbegin (*this), HIDDEN::make_outer_rend (*this)); 
  }

  global_reverse_iterator rend_all () { 
    return Galois::stl_two_level_rend (
        HIDDEN::make_outer_rbegin (*this), HIDDEN::make_outer_rend (*this)); 
  }

  global_const_reverse_iterator rbegin_all () const { 
    return Galois::stl_two_level_crbegin (
        HIDDEN::make_outer_rbegin (*this), HIDDEN::make_outer_rend (*this));
  }

  global_const_reverse_iterator rend_all () const { 
    return Galois::stl_two_level_crend (
        HIDDEN::make_outer_rbegin (*this), HIDDEN::make_outer_rend (*this));
  }


  size_type size_all () const {
    size_type sz = 0;

    for (unsigned i = 0; i < perThrdCont.size (); ++i) {
      sz += get (i).size ();
    }

    return sz;
  }


  void clear_all () {
    for (unsigned i = 0; i < perThrdCont.size (); ++i) {
      get (i).clear ();
    }
  }

  bool empty_all () const {
    bool res = true;
    for (unsigned i = 0; i < perThrdCont.size (); ++i) {
      res = res && get (i).empty ();
    }

    return res;
  }


  // TODO: fill parallel

  template <typename Iter, typename R>
  void fill_serial (Iter begin, Iter end,
      R (Cont_ty::*pushFn) (const value_type&)=&Cont_ty::push_back) {

    const unsigned P = Galois::getActiveThreads ();

    typedef typename std::iterator_traits<Iter>::difference_type Diff_ty;

    // integer division, where we want to round up. So adding P-1
    Diff_ty block_size = (std::distance (begin, end) + (P-1) ) / P;

    assert (block_size >= 1);

    Iter block_begin = begin;

    for (unsigned i = 0; i < P; ++i) {

      Iter block_end = block_begin;

      if (std::distance (block_end, end) < block_size) {
        block_end = end;

      } else {
        std::advance (block_end, block_size);
      }

      for (; block_begin != block_end; ++block_begin) {
        // workList[i].push_back (Marked<Value_ty> (*block_begin));
        ((*this)[i].*pushFn) (value_type (*block_begin));
      }

      if (block_end == end) {
        break;
      }
    }
  }

};


// TODO: rename to per thread heap factory, move outside
template <typename T>
struct PerThreadAllocatorFactory {


  typedef MM::SimpleBumpPtrWithMallocFallback<MM::FreeListHeap<MM::SystemBaseAlloc> > BasicHeap;

  typedef MM::ThreadAwarePrivateHeap<BasicHeap> PerThreadHeap;

  typedef MM::ExternRefGaloisAllocator<T, PerThreadHeap> PerThreadAllocator;

};

// TODO: remove code reuse here.
template <typename T>
class PerThreadVector: 
  public PerThreadWorkList<std::vector<T, typename PerThreadAllocatorFactory<T>::PerThreadAllocator> > {

public:
  typedef typename PerThreadAllocatorFactory<T>::PerThreadHeap Heap_ty;
  typedef typename PerThreadAllocatorFactory<T>::PerThreadAllocator Alloc_ty;
  typedef std::vector<T, Alloc_ty> Cont_ty;

protected:
  typedef PerThreadWorkList<Cont_ty> Super_ty;

  Heap_ty heap;
  Alloc_ty alloc;

public:
  PerThreadVector (): Super_ty (), heap (), alloc (&heap) {

    Super_ty::init (Cont_ty (alloc));
  }

  

};


template <typename T>
class PerThreadDeque: 
  public PerThreadWorkList<std::deque<T, typename PerThreadAllocatorFactory<T>::PerThreadAllocator> > {

public:
  typedef typename PerThreadAllocatorFactory<T>::PerThreadHeap Heap_ty;
  typedef typename PerThreadAllocatorFactory<T>::PerThreadAllocator Alloc_ty;
  typedef std::deque<T, Alloc_ty> Cont_ty;

protected:
  typedef PerThreadWorkList<Cont_ty> Super_ty;

  Heap_ty heap;
  Alloc_ty alloc;

public:
  PerThreadDeque (): Super_ty (), heap (), alloc (&heap) {

    Super_ty::init (Cont_ty (alloc));

  }

};

template <typename T>
class PerThreadList:
  public PerThreadWorkList<std::list<T, typename PerThreadAllocatorFactory<T>::PerThreadAllocator> > {

public:
  typedef typename PerThreadAllocatorFactory<T>::PerThreadHeap Heap_ty;
  typedef typename PerThreadAllocatorFactory<T>::PerThreadAllocator Alloc_ty;
  typedef std::list<T, Alloc_ty> Cont_ty;

protected:
  typedef PerThreadWorkList<Cont_ty> Super_ty;

  Heap_ty heap;
  Alloc_ty alloc;

public:
  PerThreadList (): Super_ty (), heap (), alloc (&heap) {

    Super_ty::init (Cont_ty (alloc));

  }

};

template <typename T, typename C=std::less<T> >
class PerThreadSet: 
  public PerThreadWorkList<std::set<T, C, typename PerThreadAllocatorFactory<T>::PerThreadAllocator> > {

public:
  typedef typename PerThreadAllocatorFactory<T>::PerThreadHeap Heap_ty;
  typedef typename PerThreadAllocatorFactory<T>::PerThreadAllocator Alloc_ty;
  typedef std::set<T, C, Alloc_ty> Cont_ty;

protected:
  typedef PerThreadWorkList<Cont_ty> Super_ty;

  Heap_ty heap;
  Alloc_ty alloc;

public:
  explicit PerThreadSet (const C& cmp=C ()): Super_ty (), heap (), alloc (&heap) {

    Super_ty::init (Cont_ty (cmp, alloc));
  }

};

}
} // end namespace Galois

#endif // GALOIS_RUNTIME_PER_THREAD_WORK_LIST_H_
