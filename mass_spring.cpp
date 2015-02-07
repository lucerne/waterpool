/**
 * @file mass_spring.cpp
 * Implementation of mass-spring system using Graph
 *
 * @brief Reads in two files specified on the command line.
 * First file: 3D Points (one per line) defined by three doubles
 * Second file: Tetrahedra (one per line) defined by 4 indices into the point
 * list
 */

#include "MODEL/SDLViewer.hpp"
#include "MODEL/Util.hpp"
#include "MODEL/Color.hpp"

#include "Graph.hpp"
#include "Point.hpp"
#include "Simplex.hpp"
#include "Force.hpp"
#include "Constraint.hpp"

#include <fstream>
#include <cmath>
#include <assert.h>

/** Change a graph's nodes according to a step of the symplectic Euler
 *    method with the given node force.
 * @param[in,out] g Graph
 * @param[in] t The current time (useful for time-dependent forces)
 * @param[in] dt The time step
 * @param[in] force Function object defining the force per node
 * @return the next time step (usually @a t + @a dt)
 *
 * @a force is called as @a force(n, @a t), where n is a node of the graph
 * and @a t is the current time parameter. @a force must return a Point
 * representing the node's force at time @a t.
 */ 
template <typename G, typename F, typename C>
double symp_euler_step(G& g, double t, double dt, F force, C constr) {
  
  // Update position with positions at t = n
  for(auto nit = g.node_begin(); nit != g.node_end(); ++nit){
  	auto node = (*nit);
  	auto pos = node.position();
  
  	// Retrieve the velocity
  	auto vel = node.value().velocity;
  
	  // Calculate and update the next position. 
	  Point pos_next = pos + vel * dt;
  	node.set_position(pos_next);
  }
    
  // Impose constraint and reset node position and velocity at t = n+1
  constr(g, t);  
 
  // Update velocity with positions at t = n+1
  for(auto nit = g.node_begin(); nit != g.node_end(); ++nit){	
    auto node = (*nit);
  	auto vel = node.value().velocity;

  	// Calculate force
  	Point f_next = force(node,t);
  	 	
  	// Calculate and update the next velocity. 
  	Point v_next = vel + f_next * dt / m;
  		
  	node.value().velocity = v_next;
  }  
  return t + dt;
}

/** Initialize spring constant and spring length at t = 0
 * @param[in,out] g Graph
 *
 * Spring constant, spring length for each edge is initialized.
 *
 * Complexity: O(num_edges()).
 */ 
template <typename G>
void init(G& g) {		
  // Initialize node mass and velocity 
  for (auto nit = g.node_begin(); nit != g.node_end(); ++nit) {
  	auto node = (*nit);
  	
  	// initialize r with some big number
  	double r = 10;
  	
  	for(auto eit = node.edge_begin(); eit != node.edge_end(); ++eit){
  		auto edge = (*eit);
  		edge.value() = {100, edge.length()};
  		
  		if (r > edge.value().spring_length) r = edge.value().spring_length;
  	}
  	
  	node.value() = {m, Point(0,0,0),r};
	}
};

int main(int argc, char* argv[]) {

  if (argc < 2) {
    exit(1);
  }

  GraphType graph;

  // Read all Points and add them to the Graph
  std::ifstream nodes_file(argv[1]);
  Point p;
  while (MODEL::getline_parsed(nodes_file, p))
    graph.add_node(p);

  // Read all mesh squares and add their edges to the Graph
  std::ifstream tets_file(argv[2]);
  Tetrahedron t;		// Reuse this type
  while (MODEL::getline_parsed(tets_file, t)) {
    if (t.n[0] < graph.size() && t.n[1] < graph.size()
	&& t.n[2] < graph.size() && t.n[3] < graph.size()) {
      graph.add_edge(graph.node(t.n[0]), graph.node(t.n[1]));
      graph.add_edge(graph.node(t.n[0]), graph.node(t.n[2]));
//#if 0
      // Diagonal edges: 
      graph.add_edge(graph.node(t.n[0]), graph.node(t.n[3]));
      graph.add_edge(graph.node(t.n[1]), graph.node(t.n[2]));
//#endif
      graph.add_edge(graph.node(t.n[1]), graph.node(t.n[3]));
      graph.add_edge(graph.node(t.n[2]), graph.node(t.n[3]));
    }
  }
  
  // Set initial conditions for your nodes if necessary.
  m = 1.0/graph.num_nodes();
  c = 1.0/graph.num_nodes();
  
  // Initialize graph
  init(graph);
  
  // Make force
  auto f_g = make_combined_force(GravityForce(), DampingForce());
  auto f   = make_combined_force(f_g, MassSpringForce());


  // Make constraint
  std::vector<std::pair<int, Point>> constrained_nodes;
  for (auto nit = graph.node_begin(); nit != graph.node_end(); ++nit) {
  	if ((*nit).position() == Point(0,0,0) || (*nit).position() == Point(1,0,0))
  		constrained_nodes.push_back(std::make_pair((*nit).index(), (*nit).position()));
  }
  
  auto c_sphere = make_combined_constraint(NodeConstraint(constrained_nodes), Sphere());								
  auto c = make_combined_constraint(c_sphere, Neighbor());													
 
  // Print out the stats
  std::cout << graph.num_nodes() << " " << graph.num_edges() << std::endl;

  // Launch the SDLViewer
  MODEL::SDLViewer viewer;
  viewer.launch();

  auto node_map = viewer.empty_node_map(graph);
  viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
  viewer.add_edges(graph.edge_begin(), graph.edge_end(), node_map);

  viewer.center_view();

  // Begin the mass-spring simulation
  double dt = 0.001;
  double t_start = 0;
  double t_end = 3.0;
  
  
  // Testing Morton Coder using a simple L=2 level 
  // Change Morton Coder constructor to BoundingBox(Point(0,0,0), Point(2,2,2))
  auto bb = BoundingBox(Point(-0,-0,-0), Point(2,2,2));  
  auto mc = MortonCoder<2>(bb);
  Point p1 = Point(0,0,0);
  auto c1 = mc.code(p1);
  auto bb1 = mc.cell(c1);
  
  // Testing Morton Coder in the graph functions add node, set position 
  auto inc_p = Point(0.3125, 0.3125, 0.3125);
  for (auto nit = graph.node_begin(); nit != graph.node_end(); ++nit) {
  	auto node = (*nit);
  	
	}
  
  auto dd =  BoundingBox(Point(0,0,0), Point(0.5,0.5,1.5)); 

  for (auto hit = graph.neighbor_begin(dd); hit != graph.neighbor_end(dd); ++hit){
  	auto node = (*hit);
  }
  
  MODEL::Clock start;
  start.start();
    
  for (double t = t_start; t < t_end; t += dt) {
		symp_euler_step(graph, t, dt, f, c);

    // Update viewer with nodes' new positions
    viewer.add_nodes(graph.node_begin(), graph.node_end(), node_map);
    viewer.set_label(t);

    // These lines slow down the animation for small graphs, like grid0_*.
    // Feel free to remove them or tweak the constants.
    if (graph.size() < 100)
      MODEL::sleep(0.001);
  }

  return 0;
}
