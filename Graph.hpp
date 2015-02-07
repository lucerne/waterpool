#ifndef MODEL_GRAPH_HPP
#define MODEL_GRAPH_HPP

/**
 * @file Graph.hpp
 */

#include "MODEL/Util.hpp"
#include "Point.hpp"
#include "MortonCoder.hpp"
#include "BoundingBox.hpp"
#include <utility>
#include <algorithm>
#include <vector>
#include <assert.h>

#include <map>
#include <list>

using namespace std::rel_ops;

/**************************************************************************************
 *	 
 * 	 Graph
 *
 ************************************************************************************** 	 
 */
/** @class Graph
 * @brief
 */
template <typename V, typename E, int L=5>
class Graph {
 public:
 
  /** Type of this graph. */
  typedef Graph graph_type;

  /** Type of graph nodes. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  typedef Node node_type; 	
  typedef V node_value_type;
  
  /** Type of graph edges. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  typedef Edge edge_type;
  typedef E edge_value_type;    
  
  /** Type of indexes and sizes. Return type of Node::index() and
      Graph::num_nodes(), argument type of Graph::node. */
  typedef unsigned size_type;

	typedef typename MortonCoder<L>::code_type code_type;

  /** Accessing internal node 
    * 
    * Stores position @a np, node value @a nv, index @a idx */
  struct internal_node {
    Point np;   
    node_value_type nv;      
    size_type idx;
    code_type cd;
  };

  /** Accessing internal edges 
    * 
    * Stores position @a np, node value @a nv, index @a idx */
  struct internal_edge {   
  	size_type row_uid_;  
    size_type col_uid_;
    edge_value_type ev;
  };

  /** Construct an empty graph. */
  Graph() 
  	: total_edge_(0), v_node_(), idx_uid_(), adj_edge_(), cd_uid_(pow(2, 3*L)), mc_(MortonCoder<L>(BoundingBox(Point(-5,-5,-5), Point(5,5,5))) ) {
  	
  }

  /** Destroy the graph. */
  ~Graph() {
    clear();    
  }

/**************************************************************************************
 *	 
 * 	 Nodes
 *
 ************************************************************************************** 	 
 */
  class incident_iterator; 
  
  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
   
  class Node {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Node x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * @endcode
     */
    Node(const Graph* g, size_type uid) 
    : g_(const_cast<Graph*>(g)), uid_(uid) {    
    }

    /** Return the graph that contains this node. */
    graph_type* graph() const {
      return g_;
    }

    /** Return this node's position. 
      * 
      * Nodes are stored in a vector of internal_node struct */
    Point position() const {	 
      return g_->v_node_.at(uid_).np;
    }

    /** Return this node's index.
      * @pre 0 <= @a i < num_nodes()
	  *
      * Complexity: O(1).
       */
    size_type index() const { 
      return g_->v_node_.at(uid_).idx;
    }

	/** Inequality between two nodes. 
	  * Strict ordering on different graphs by pointer addresses.
	  */
    bool operator<(const Node& n) const {
    	if (g_ < n.g_) return true;
    	else if (g_ == n.g_) return uid_ < n.uid_;
    	else return false;
    }
    
	/** Equality between two nodes by pointer addresses. */
    bool operator==(const Node& n) const {
    	if (g_ == n.g_) return uid_ == n.uid_;
    	else return false;
    }    
    
    /** Return this node's value. */
    node_value_type& value(){
    	return g_->v_node_.at(uid_).nv;
    }
    
    /** Return this node's Morton Coder value. */
    code_type& code(){
    	return g_->v_node_.at(uid_).cd;
    }
      

	/** Return this node's number of adjacency vertices only if add_edge has been called 
		* for that node
		*/
	size_type degree() const{
		if (uid_ < g_->adj_edge_.size()) return g_->adj_edge_[uid_].size();
		else return 0;
	}

	/** Return an iterator to the beginning of this node's adjacency list. */
	incident_iterator edge_begin() const{
		incident_iterator Iter;
		Iter.p_ = g_->adj_edge_[uid_].begin();
		Iter.g_ = const_cast<Graph*>(g_);
  	return Iter;
	}
	
	/** Return an iterator to the past-the-end element of this node's adjacency list. */
	incident_iterator edge_end() const{
		incident_iterator Iter; 
		Iter.p_ = g_->adj_edge_[uid_].end();
		Iter.g_ = const_cast<Graph*>(g_);
  	return Iter;
	}
	
	/** Set this node to a new position. */
	void set_position(const Point& p){
	
		g_->v_node_.at(uid_).np = p;

		code_type m = g_->mc_.code(p);
		code_type old_m = g_->v_node_.at(uid_).cd;
		
		typename std::vector<code_type>::iterator it;
		typename std::vector<code_type>::iterator it_found = g_->cd_uid_.at(old_m).end();
		typename std::vector<std::vector<code_type> >::iterator rit;
		
		g_->v_node_.at(uid_).cd = m;
		g_->cd_uid_[m].push_back(uid_);
		
		for(it = g_->cd_uid_.at(old_m).begin(); it != g_->cd_uid_.at(old_m).end(); ++it){
			if ((*it) == uid_) it_found = it;	
		}
		
		g_->cd_uid_.at(old_m).erase(it_found);
		
			
	}
	
   private:
    friend class Graph;
    /** this graph */
  	Graph* g_;
  	/** this node's uid */
  	size_type uid_;
  };

/**************************************************************************************
 *	 
 * 	 Graph Node Method
 *
 ************************************************************************************** 	 
 */
  /** Return the number of nodes in the graph. */
  size_type size() const {
    return idx_uid_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return Node(this, idx_uid_.at(i));		
  }


  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new size() == old size() + 1
   * @post The returned node's index() == old size()
   *
   * Complexity: O(1) amortized time.  */
  Node add_node(const Point& position, const node_value_type& nv = node_value_type()) {
    size_type old_size = size(); 
    
		code_type m = mc_.code(position);
    internal_node e = {position, nv, old_size, m};
    v_node_.push_back(e);
    idx_uid_.push_back(old_size);
    
    cd_uid_[m].push_back(old_size);
    
    return Node(this, old_size);
  }

  /** Remove a node from the graph.
   * @param[in] n Node to be removed
   * @return true
   * @pre @a n == node(@a i) for some @a i with 0 <= @a i < size()
   * @post new size() == old size() - 1
   *
   * Can invalidate outstanding iterators. @a n becomes invalid, as do any
   * other Node objects equal to @a n. All other Node objects remain valid.
   *
   * Complexity: O( degree()) */
  bool remove_node(const Node& n) {
    
    /* Base Case : graph has 1 node, clear all */
    if ( idx_uid_.size() == 1 ) {
    	clear();     
    	return true;
    }
    
    size_type n_idx = n.index();
    size_type n_uid = n.uid_;
    size_type last_uid = idx_uid_.back();
    
    idx_uid_[n_idx] = last_uid;
    idx_uid_.erase ( idx_uid_.end()-1 );
    v_node_.at(last_uid).idx = n_idx;
        
    /** Removing all edges connected to the node 
     * 
     * LI : @a n_uid
     * Not LI : @a num_edges() decreases upon deletion 
     *
     * Base case: adj_edge_.empty() == true, no edge is added 
     * Remove the 0-th element until size() decreases to 0
     */
    if ( !adj_edge_.empty() ) { 
     	size_type k = 0; 
     	while ( k != adj_edge_.at(n_uid).size() ) 
     		remove_edge( Node(this, adj_edge_.at(n_uid).at(k).col_uid_), Node(this, n_uid) ); 
     }   
    return true;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects. */
  void clear() {
    v_node_.clear();
    adj_edge_.clear();
    idx_uid_.clear();
    total_edge_ = 0;
  }

/**************************************************************************************
 *	 
 * 	 Edges
 *
 ************************************************************************************** 	 
 */

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order. */
  class Edge {
   public:
  	/* Additional data structure for Edge */
  	typedef edge_value_type value_type;    
  	
    /** Construct an invalid Edge. */
    Edge(const Graph* g, std::pair<size_type,size_type> p) 
    : g_(const_cast<Graph*>(g)), uid1_(p.first), uid2_(p.second) {
    }
    /** Return a node of this Edge */
    Node node1() const {
      return Node(g_, uid1_);
    }    
    /** Return the other node of this Edge */
    Node node2() const {
      return Node(g_,uid2_);
    }
    
    
    /** 
     * Impose strict ordering on different graphs by pointer address, and uids
     * @pre Edge proxies are created with @a uid1 < @a uid2 
     */    
    bool operator<(const Edge& n) const {
    	if (g_ < n.g_) return true;
    	if (g_ == n.g_) {
    		if (uid1_ < n.uid1_) return true;
    		else if ( (uid1_ == n.uid1_) && (uid2_ < n.uid2_) ) return true;
    		else return false;
    	}
    }
    /** Equality: two edges are equal if they are the same graph with the 
      * same unordered pair of uid_'s */ 	
    bool operator==(const Edge& n) const {
    	if (g_ == n.g_) return ( (uid1_ == n.uid1_) && (uid2_ == n.uid2_) ) || ( (uid1_ == n.uid2_) && (uid2_ == n.uid1_) )  ;
    	else return false;
    }  
   
   	/** Calculate length of edge */ 	
    double length() const{
    	return (node1().position() - node2().position()).mag();
    }

	/* Return a sequence of edge values associated with node1's adjacency vertices */
	const value_type& value() const{
		typename std::vector<internal_edge>::const_iterator it;  
		for(it = g_->adj_edge_.at(uid1_).begin(); it != g_->adj_edge_.at(uid1_).end(); ++it){
			if ((*it).col_uid_ ==  uid2_) return (*it).ev;
			}
		assert(0);
		return (*g_->adj_edge_.at(0).begin()).ev;
	}

	/* Return an edge value. If no edge value exists, return 0 value
	 * @pre edge exists
	 */
	edge_value_type& value(){
		typename std::vector<internal_edge>::iterator it;  
		for(it = g_->adj_edge_.at(uid1_).begin(); it != g_->adj_edge_.at(uid1_).end(); ++it){
			if ((*it).col_uid_ ==  uid2_) return (*it).ev;
			}
		assert(0);
		return (*g_->adj_edge_.at(0).begin()).ev;
	}
		

   private:
    friend class Graph;
  	Graph* g_;
  	/* uids of node1, node2 */
  	size_type uid1_, uid2_;
  };
  
/**************************************************************************************
 *	 
 * 	 Edges Methods
 *
 ************************************************************************************** 	 
 */

  /** Return the total number of edges.
   *
   * Complexity: O(1) */ 
  size_type num_edges() const {
    return total_edge_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(num_edges())
   *
   * A very slow implementation by sequentially exam all adjacency lists */
  Edge edge(size_type i) const { 
    size_type sum = 0;
    for(size_type k = 0; k != adj_edge_.size(); k++)
		for(size_type j = 0; j != adj_edge_.at(k).size(); j++){
			if ( k < adj_edge_.at(k).at(j).col_uid_ ){
    			++sum;
    			if ( sum == (i+1) ) 
    				return Edge(this, std::make_pair(k, adj_edge_.at(k).at(j).col_uid_) );  
    		}
    	}
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return true if, for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(degree()) amortized time 
   *
   * find()'s worst case is O(degree()). 
   * Optimized by searching the smaller of the two adjacency lists */
  bool has_edge(const Node& a, const Node& b) const {
    
	if ( (adj_edge_.size() <= a.uid_) || (adj_edge_.size() <= b.uid_) ) return false;
	
	typename std::vector<internal_edge>::const_iterator it;
	
	if ( adj_edge_.at(a.uid_).size() < adj_edge_.at(b.uid_).size() ){
		for(it = adj_edge_.at(a.uid_).begin(); it != adj_edge_.at(a.uid_).end(); ++it)
			if ((*it).col_uid_ ==  b.uid_) return true;			
	} else{
		for(it = adj_edge_.at(b.uid_).begin(); it != adj_edge_.at(b.uid_).end(); ++it)
			if ((*it).col_uid_ ==  a.uid_) return true;		
	}
	
	return false;
    
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b) == true, then new num_edges() ==
   *   old num_edges(). Otherwise, new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: O(degree()), limited by has_edge() 
   * 
   * Number of entries of adjacency lists = 2*num_edges() 
   * Return one proxy, where row @a uid < column @a uid 
   */ 
    /** Delete both copies of the edge. Decrement size. 
    Impose an ordering on Edge proxies, such that the first uid of any Edge proxy
    is always the smaller uid */   
   
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& ev = edge_value_type() ) {    
    if (a.uid_ > b.uid_ ) {
	    if ( (a.uid_+1) > adj_edge_.size()) {
    		adj_edge_.resize(a.uid_+1);
    	} 
    } else {
 	    if ( (b.uid_+1) > adj_edge_.size()) {
    		adj_edge_.resize(b.uid_+1);
    	}    
    }     
    
    if (!adj_edge_.empty() && !has_edge(a,b)) {
    	++total_edge_;
    	
			adj_edge_.at(a.uid_).push_back( {a.uid_, b.uid_, ev} );
			adj_edge_.at(b.uid_).push_back( {b.uid_, a.uid_, ev} );
		
			return Edge(this, std::make_pair(a.uid_, b.uid_));	
    } 
    else return Edge(this, std::make_pair(a.uid_, b.uid_)); 
  }

  /** Remove an edge, if any, returning the number of edges removed.
   * @param[in] a,b Nodes composing an edge to be removed
   * @return false if old has_edge(@a a, @a b) == false, true otherwise
   * @pre @a a and @a b are valid nodes of this graph
   * @post !has_edge(@a a, @a b)
   * @post new num_edges() == old num_edges() - return value
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Can invalidate all edge and incident iterators.
   * Invalidates any edges equal to Edge(@a a, @a b). Must not invalidate
   * other outstanding Edge objects.
   *
   * Complexity: O(2 degree()), limited by has_edge() and find() 
   * 
   * Remove edges (a, b) and (b, a) */
  bool remove_edge(const Node& a, const Node& b) {
    if (has_edge(a,b)) { 	
     	
     	typename std::vector<internal_edge>::iterator it;  
    	typename std::vector<internal_edge>::iterator it_found;  
	
			for(it = adj_edge_.at(a.uid_).begin(); it != adj_edge_.at(a.uid_).end(); ++it)
				if ((*it).col_uid_ ==  b.uid_) it_found = it;
    	
	   	adj_edge_.at(a.uid_).erase(it_found);

			for(it = adj_edge_.at(b.uid_).begin(); it != adj_edge_.at(b.uid_).end(); ++it)
				if ((*it).col_uid_ ==  a.uid_) it_found = it;
    		
	   	adj_edge_.at(b.uid_).erase(it_found);
	   	
    	--total_edge_;
		
		return true;
	}
	return false;
  }

  /** Remove an edge, if any, returning the number of edges removed.
   * @param[in] e Edge to be removed
   * @return true
   * @pre @a e is a valid edge of this graph
   * @pre has_edge(@a e.node1(), @a e.node2()) == true
   * @post has_edge(@a e.node1(), @a e.node2()) == false
   * @post new num_edges() == old num_edges() - 1
   *
   * This is a synonym for remove_edge(@a e.node1(), @a e.node2()), but its
   * implementation can assume that @a e is definitely an edge of the graph.
   * This might allow a faster implementation.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Can invalidate all edge and incident iterators.
   * Invalidates any edges equal to Edge(@a a, @a b). Must not invalidate
   * other outstanding Edge objects.
   *
   * Complexity: O(degree()), limited by find() 
   * 
   * Remove edges (a, b) and (b, a)
   */
  bool remove_edge(const Edge& e) {
    	--total_edge_;
    	
    	size_type a_uid = e.uid1_;
    	size_type b_uid = e.uid2_;

     	typename std::vector<internal_edge>::iterator it;  
    	typename std::vector<internal_edge>::iterator it_found;  
	
			for(it = adj_edge_.at(a_uid).begin(); it != adj_edge_.at(a_uid).end(); ++it)
				if ((*it).col_uid_ ==  b_uid) it_found = it;
    	
	   	adj_edge_.at(a_uid).erase(it_found);

			for(it = adj_edge_.at(b_uid).begin(); it != adj_edge_.at(b_uid).end(); ++it)
				if ((*it).col_uid_ ==  a_uid) it_found = it;
    		
	   	adj_edge_.at(b_uid).erase(it_found);
    	
    return true;
  }

  /** @class Graph::node_iterator
   * @brief Iterator class for nodes. A forward iterator. */
  class node_iterator {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Node value_type;
    /** Type of pointers to elements. */
    typedef Node* pointer;
    /** Type of references to elements. */
    typedef Node& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;

    /** Construct an invalid node_iterator. */
    node_iterator() : p_(nullptr){
    
    /** Return the content of iterator */
    Node operator*() const{
    	return g_->node(p_->idx);  
    }
    
    /** Increment the iterator */
    node_iterator& operator++(){
    	++p_;
    	return *this;
    }
    
    /** Equality holds if the iterator points to the same address */
    bool operator==(const node_iterator& rhs) const{
    	return p_ == rhs.p_;
    }
    
    /** 
      * Less than holds if the iterators are not equal */
    bool operator<(const node_iterator& rhs) const{
    	return p_ != rhs.p_;
    }

   private:
    friend class Graph;
    friend class edge_iterator;
	Graph* g_;
    typename std::vector<internal_node>::const_iterator p_;
  };

  /** Return an iterator pointing to the first element of node sequence */
  node_iterator node_begin() const{ 
    node_iterator Iter; 
	Iter.p_ = v_node_.begin();
	Iter.g_ = const_cast<Graph*>(this);
  	 return Iter;
  }
  /** Return an iterator pointing to the past-the-end element of node sequence */  
  node_iterator node_end() const{
	node_iterator Iter; 
	Iter.p_ = v_node_.end();
	Iter.g_ = const_cast<Graph*>(this);
  	 return Iter;
  }

  /** @class Graph::edge_iterator
   * @brief Iterator class for edges. A forward iterator. */
  class edge_iterator {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;

    /** Construct an invalid edge_iterator. */
    edge_iterator() : pR_(nullptr), pC_(nullptr){
    }

	/** Return the content of iterator 
	 *  @pre row uid < column uid, *pC is an internal_edge */
    Edge operator*() const{
    	return Edge(g_, std::make_pair((*pC_).row_uid_, (*pC_).col_uid_));  
    }    
    
    
    /** Increment the row and column iterators of the adjacency list to the next valid 
    iterators. If no such element exists, return the past-the-end element of the 
    adjacency list, and the past-the-end element of the dereferenced row iterator.
     *    
   	 *  @pre @a pR_, @a pC_ are valid iterators in the range[adj_edge_.begin(), adj_edge_.end())   
     *  @post @a pC_ is valid such that row @a uid < column @a uid, or @a pR_ reaches 
     *  the end of an adjacency list and && @a pC_ is the past-the-end element
     *  
     *	LI : @a pC_ is either a valid iterator, or @a pC_ is the past-the-end element, or 
     *  pC_ is such that !(row @a uid < column @a uid)  
     * 
     *  Only entries in the adjacency list such that row @a uid < column @a uid are 
   	 *  iterated. User only see these entries.
     */
    edge_iterator& operator++(){	
    	pC_++;
    	
    	while ( (pC_ != (*pR_).end()) && ((*pC_).col_uid_ < (*pC_).row_uid_) ){
    		pC_++;	 		
    	}
    	
    	while( (pC_ == (*pR_).end()) && (pR_ != g_->adj_edge_.end()) ){
    	
    		++pR_;
    		if ( pR_ != g_->adj_edge_.end() ) {
    		
    			pC_ = (*pR_).begin();
    			while ( (pC_ != (*pR_).end()) && ((*pC_).col_uid_ < (*pC_).row_uid_) ){
    				pC_++;			
    			}
    		}
    	}
    	
    	/** if reached past-the-end element, set iterator to edge_iterator.end() */		
    	if ( pR_ == g_->adj_edge_.end() ) pC_ = g_->adj_edge_.back().end();
    	
    	return *this;
    }
   
    /** Equality holds if the iterator points to the same address 
     * @re(uid1, uid2) is not the same entry as (ui2, uid1) in adjacency list, 
     * hence they are not equal
     */
    bool operator==(const edge_iterator& rhs) const{
    	return (pR_ == rhs.pR_) && (pC_ == rhs.pC_) ;
    }
    
    /** Less than holds if the iterators are not equal 
     * @re(uid1, uid2) is not the same entry as (ui2, uid1) in adjacency list, 
     * hence if they are not equal, inequality holds
     */
    bool operator<(const edge_iterator& rhs) const{
    	return (pR_ != rhs.pR_) || (pC_ != rhs.pC_) ;
    }

   private:
    friend class Graph;
    Graph* g_;
    typename std::vector< std::vector<internal_edge> >::const_iterator pR_;
    typename std::vector<internal_edge>::const_iterator pC_;
  };

  /** Return a row iterator pointing to the first vector of the adjacency list, and 
   a column iterator pointing to the first uid of the dereferenced row iterator */
  edge_iterator edge_begin() const{
    edge_iterator Iter; 
		Iter.pR_ = adj_edge_.begin();
		Iter.pC_ = adj_edge_.front().begin();
		Iter.g_ = const_cast<Graph*>(this);
  	return Iter;
  }
  
  /** Return a row iterator pointing to the past-the-end element of the adjacency list,
  and a column iterator pointing to the past-the-end element of the dereferenced row 
  iterator */  
  edge_iterator edge_end() const{
    edge_iterator Iter; 
		Iter.pR_ = adj_edge_.end();
		Iter.pC_ = adj_edge_.back().end();
		Iter.g_ = const_cast<Graph*>(this);
  	return Iter;
  }


  /** @class Graph::incident_iterator
   * @brief Iterator class for edges incident to a given node. A forward
   * iterator. */
  class incident_iterator {
   public:
    // These type definitions help us use STL's iterator_traits.
    /** Element type. */
    typedef Edge value_type;
    /** Type of pointers to elements. */
    typedef Edge* pointer;
    /** Type of references to elements. */
    typedef Edge& reference;
    /** Iterator category. */
    typedef std::input_iterator_tag iterator_category;

    /** Construct an invalid incident_iterator. */
    incident_iterator() : p_(){
    }

    /** Return the content of iterator */
    Edge operator*() const{
    	return Edge(g_, std::make_pair((*p_).row_uid_, (*p_).col_uid_) );  
    }
    
    /** Increment the iterator 
     * @pre iterate down teh same adjacency list over a valid (non-null valued) range 
     */
    incident_iterator& operator++(){
    	++p_;
    	return *this;
    }
    
    /** Equality holds if the iterator points to the same address. 
     * (uid1, uid2) is not the same entry as (ui2, uid1) in adjacency lists, 
     * hence they are not equal
     * @pre Comparison only between edges in the same adjacency list 
     */
    bool operator==(const incident_iterator& rhs) const{
    	return p_ == rhs.p_;
    }
    
	/** Less than holds if the iterators are not equal. 
     * (uid1, uid2) is not the same entry as (ui2, uid1) in adjacency lists, 
     * hence if they are not equal, inequality holds
     * @pre Comparison only between edges in the same adjacency list 
     */
    bool operator<(const incident_iterator& rhs) const{
    	return p_ != rhs.p_;
	}

   private:
    friend class Graph;
    Graph* g_;
    typename std::vector<internal_edge>::const_iterator p_;
  };
  
  /** @class Graph::neighborhood_iterator
   * @brief Iterator class that visits nodes inside a bounding box. 
   */
	class neighborhood_iterator{
	
	public:
	  typedef Node value_type;
		typedef Node* pointer;
		typedef Node& reference;
		typedef std::input_iterator_tag iterator_category;
			
		neighborhood_iterator(const BoundingBox& bb) : pR_(nullptr), pC_(nullptr), bb_(bb) {
		}
		
    /** pC_ is not empty. return the node of iterator */			
		Node operator*() const{
			size_type uid = *pC_;
			return Node(g_, uid);
		}
		
    /** Increment the row and col iterator to the next element of morton code adjacency 
    	* list that is inside the constraint's bounding box and that is not empty. If such 
    	* element does not exist, set row iterator to the past-the-end element of the 
    	* adjacency list, and the col iterator to past-the-end element of the dereferenced 
    	* row iterator.
    	*
    	* @pre @a pR_, @a pC_ are non-empty iterators in range 
    	* [cd_uid_.begin(), cd_uid_.end)
    	* @pre min_ <= cd_ <= max_
    	*
    	* @post @a pR_, @a pC_ are non-empty iterators in range  
    	* [cd_uid_.begin(), cd_uid_.end) and @a min_ <= @a cd_ <= @a max_ OR 
    	* @a pR_ = cd_uid_.end(), @a pC_ = cd_uid_.back.end() and @a cd_ > @a max_
    	* 
    	*/
    neighborhood_iterator& operator++(){
    
  		++pC_;
  		
  		if ( pC_ == (*pR_).end() ){ 
  		
  			if ( pR_ != g_->cd_uid_.end()-1){
  			
  				cd_ += 1;  
  				cd_ = g_->mc_.advance_to_box(cd_, min_, max_);
  				pR_ = g_->cd_uid_.begin() + cd_;
  				
  				while ((*pR_).empty() && pR_ != g_->cd_uid_.end()-1 && cd_ <= max_){
  					
  					cd_ += 1;
	  				cd_ = g_->mc_.advance_to_box(cd_, min_, max_);
  				
  					pR_ = g_->cd_uid_.begin() + cd_;
  				}
  				
  				if (cd_ > max_ || ((*pR_).empty() && pR_ == g_->cd_uid_.end()-1)) 
  					pR_ = g_->cd_uid_.end();
  			}
  			
  			else{
  				pR_ = g_->cd_uid_.end(); 		
    		}
    		
    		if (pR_ == g_->cd_uid_.end()) {
    			pC_ = g_->cd_uid_.back().end();	 	  			
    		}
    		else{
    			pC_ = (*pR_).begin();
    		}
    		
    		if ( pC_ != (*pR_).end()){
					auto node = Node(g_, (*pC_));
  	  		if (node.index() == 0 || node.index() == 986)
  	  		}
  	  	}
    	return *this;
    }

    /** Equality holds if the iterator points to the same address 
     * @re(uid1, uid2) is not the same entry as (ui2, uid1) in adjacency list, 
     * hence they are not equal
     */
    bool operator==(const neighborhood_iterator& rhs) const{
    	return (pR_ == rhs.pR_) && (pC_ == rhs.pC_) ;
    }
    
    /** Less than holds if the iterators are not equal 
     * @re(uid1, uid2) is not the same entry as (ui2, uid1) in adjacency list, 
     * hence if they are not equal, inequality holds
     */
    bool operator<(const neighborhood_iterator& rhs) const{
    	return (pR_ != rhs.pR_) || (pC_ != rhs.pC_) ;
    }
			
		private:
			friend class Graph;
			Graph* g_;
			/** Going across uid adjacency list */
    	typename std::vector< std::vector<size_type> >::const_iterator pR_;
    	/** Going down uid adjacency list */
    	typename std::vector<size_type>::const_iterator pC_;
    	/** bounding box of constraint */
			BoundingBox bb_;
			/** mc associated with bounding box of constraint */
			code_type max_, min_, cd_;
  };
  
  /** Return a row iterator pointing to the first non-empty element of the morton code
  	* adjacency list that is also not empty, and a column iterator pointing to the first 
  	* uid of the dereferenced row iterator. Row and col iterator set to past-the-end 
  	* element if such element does not exist.  
  	*/
  neighborhood_iterator neighbor_begin(const BoundingBox& bb) const{
    neighborhood_iterator Iter(bb); 
    
    Point p_min = bb.min();
		Point p_max = bb.max();
			
		Iter.max_ = mc_.code(p_max);
		Iter.min_ = mc_.code(p_min);
			
		code_type t;
		if (Iter.max_ < Iter.min_){
			t = Iter.max_;
			Iter.max_ = Iter.min_;
			Iter.min_ = t;
		}
		
		Iter.cd_ = mc_.advance_to_box(0, Iter.min_, Iter.max_);
		Iter.pR_ = cd_uid_.begin() + Iter.cd_; 
		
		while((*Iter.pR_).empty() && Iter.pR_ != cd_uid_.end() ){
				
			Iter.cd_ += 1;
				
			if (!mc_.is_in_box(Iter.cd_, Iter.min_, Iter.max_)) 
			 	Iter.cd_ = mc_.advance_to_box(Iter.cd_, Iter.min_, Iter.max_);
			 	
			Iter.pR_ = cd_uid_.begin() + Iter.cd_; 
			
			if (Iter.cd_ > Iter.max_) Iter.pR_ = cd_uid_.end(); 
			else Iter.pR_ = cd_uid_.begin() + Iter.cd_; 
		}

		if (Iter.pR_ == cd_uid_.end()) Iter.pC_ = cd_uid_.back().end();
		else Iter.pC_ = (*Iter.pR_).begin();		
			
		Iter.g_ = const_cast<Graph*>(this);
		    
  	return Iter;
  }
  
  /** Return a row iterator pointing to the past-the-end element of morton code list,
  and a column iterator pointing to the past-the-end element of the dereferenced row 
  iterator */  
  neighborhood_iterator neighbor_end(const BoundingBox& bb) const{
    neighborhood_iterator Iter(bb); 
  
    Iter.pR_ = cd_uid_.end();
    Iter.pC_ = cd_uid_.back().end();
    Iter.g_ = const_cast<Graph*>(this);
    return Iter; 
  }

 private:
    /* total size of edge */
  	size_type total_edge_;
  	/* node container */
		std::vector<internal_node> v_node_;
		/* look up table for node index and uid */
		std::vector<size_type> idx_uid_;
		/* edge container */
		std::vector< std::vector<internal_edge> > adj_edge_;	
		/* Morton code to node uid lookup table */
		std::vector< std::vector<size_type> > cd_uid_;	
		/* Morton Code */
		MortonCoder<L> mc_;
		/* BoundingBox and Morton Code */
		//BoundingBox bb_;
		//typename MC_type mc_(bb_);
    
};

#endif
