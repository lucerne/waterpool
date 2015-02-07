#ifndef MODEL_FORCE_HPP
#define MODEL_FORCE_HPP

#include "MODEL/SDLViewer.hpp"
#include "MODEL/Util.hpp"
#include "MODEL/Color.hpp"

#include "Graph.hpp"
#include "Point.hpp"
#include "Simplex.hpp"

#include <fstream>
#include <cmath>
#include <assert.h>

/* Gravity in meters/sec */
static constexpr double grav = 9.81;
/* mass */
double m;
/* viscosity */
double c;

/* Node attributes */
template <typename M, typename V, typename R>
	struct node_value {
		M mass;
		V velocity;
		R radius;
};

/* Edge attributes */
template <typename K, typename L>
	struct edge_value {
		K spring_constant;
		L spring_length;
};

/* Graph Types definitions */
typedef Graph<node_value<double, Point, double>, edge_value<double, double>> GraphType;
typedef GraphType::Node Node;
typedef GraphType::Edge Edge;
typedef GraphType::node_value_type node_value_type;
typedef GraphType::edge_value_type edge_value_type;


/* Making a combined ForcePair using F1, F2 constructors */
template< typename F1, typename F2>
struct ForcePair{
	F1 force1;
	F2 force2;
	
	/* Making a combined ForcePair using F1, F2 constructors */
	ForcePair(F1 f1, F2 f2) : force1(f1), force2(f2) {
	}
		
	template <typename NODE>
  	Point operator()(NODE n, double t) {
  		return force1(n,t) + force2(n,t);
  	}	
};

/** Gravitational force */
struct GravityForce {

	template <typename NODE>
  	Point operator()(NODE n, double t) {
  		(void) t;
  		(void) n;
   		Point Fg = Point(0,0,-1) * m * grav;
   		return Fg;
  	}
};
  
/** Spring force */
struct MassSpringForce {
  
  template <typename NODE>
  Point operator()(NODE n, double t) {
  	(void) t;
    Point Fs = Point(0,0,0);

    // Calculate from positions of all adjacent vertices 
    auto pos1 = n.position();
    
  	for(auto it = n.edge_begin(); it != n.edge_end(); ++it){
  		auto edge = (*it);
  		auto pos2 = edge.node2().position();
  			
  		auto K = edge.value().spring_constant;
  		auto L = edge.value().spring_length;
  		
  		Fs = Fs - ( pos1 - pos2 ) * K * ((pos1 - pos2).mag() - L) / (pos1 - pos2).mag();
		
	    }
    	return Fs;
  	}
};


/** Damping force */
struct DampingForce {

	template <typename NODE>
  	Point operator()(NODE n, double t) {
  		(void) t;
  		/** Calculate from velocity */
  		auto v = n.value().velocity;
  		
   		Point Fd = - v * c;
   		return Fd;
  	}
};


  template< typename F1, typename F2>
  ForcePair<F1, F2> make_combined_force(F1 first, F2 second){ 
  	return ForcePair<F1, F2>(first, second);
  }

#endif
