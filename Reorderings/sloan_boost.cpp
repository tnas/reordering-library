/*
 * Copyright 2016 Thiago Nascimento nascimenthiago@gmail.com
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#include "sloan_boost.h"
#include <boost/config.hpp>
#include <vector>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/sloan_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/graph/profile.hpp>
#include <boost/graph/wavefront.hpp>

using namespace boost;
using namespace std;
using std::cout;
using std::endl;

typedef adjacency_list<
    setS, 
    vecS, 
    undirectedS, 
    property<
    vertex_color_t, 
    default_color_type,
    property<
    vertex_degree_t,
    int,
    property<
    vertex_priority_t,
    double > > > > SloanGraph;
typedef graph_traits<SloanGraph>::vertex_descriptor SloanVertex;
typedef graph_traits<SloanGraph>::vertices_size_type sloan_size_type;
typedef std::pair<std::size_t, std::size_t> SloanPair;


SloanGraph build_boost_sloan_graph(const METAGRAPH* mgraph)
{
	int node, neigh, degree, num_nodes;
	int* neighboors;
	SloanPair edge;
	
	num_nodes = mgraph->size; 
	SloanGraph bgraph(num_nodes);
	
	for (node = 0; node < num_nodes; ++node)
	{
		degree     = mgraph->graph[node].degree;
		neighboors = mgraph->graph[node].neighboors;
		
		for (neigh = 0; neigh < degree; ++neigh)
		{
			edge = SloanPair(node, neighboors[neigh]);
			add_edge(edge.first, edge.second, bgraph);
		}
	}
	
	return bgraph;
}


void Boost_Sloan(const METAGRAPH* mgraph, int** permutation, int start_node, int end_node)
{
	SloanGraph bgraph = build_boost_sloan_graph(mgraph);
  
	*permutation  = new int[mgraph->size];
	
	//Creating two iterators over the vertices
	graph_traits<SloanGraph>::vertex_iterator ui, ui_end;

	//Creating a property_map with the degrees of the degrees of each vertex
	property_map<SloanGraph,vertex_degree_t>::type deg = get(vertex_degree, bgraph);
	for (boost::tie(ui, ui_end) = vertices(bgraph); ui != ui_end; ++ui)
		deg[*ui] = degree(*ui, bgraph);

	//Creating a property_map for the indices of a vertex
// 	property_map<SloanGraph, vertex_index_t>::type index_map = get(vertex_index, bgraph);

// 	std::cout << "original bandwidth: " << bandwidth(bgraph) << std::endl;
// 	std::cout << "original profile: " << profile(bgraph) << std::endl;
// 	std::cout << "original max_wavefront: " << max_wavefront(bgraph) << std::endl;
// 	std::cout << "original aver_wavefront: " << aver_wavefront(bgraph) << std::endl;
// 	std::cout << "original rms_wavefront: " << rms_wavefront(bgraph) << std::endl;
  
	//Creating a vector of vertices  
	std::vector<SloanVertex> sloan_order(num_vertices(bgraph));
	//Creating a vector of size_type  
	std::vector<sloan_size_type> perm(num_vertices(bgraph));

//     //Setting the start node
//     SloanVertex s = vertex(0, bgraph);
//     int ecc;   //defining a variable for the pseudoperipheral radius
//     
//     //Calculating the pseudoeperipheral node and radius
//     SloanVertex e = pseudo_peripheral_pair(bgraph, s, ecc, get(vertex_color, bgraph), get(vertex_degree, bgraph) );
// 
//     cout << endl;
//     cout << "Starting vertex: " << s << endl;
//     cout << "Pseudoperipheral vertex: " << e << endl;
//     cout << "Pseudoperipheral radius: " << ecc << endl << endl;

	// Sloan ordering
	sloan_ordering(bgraph, start_node, end_node, sloan_order.begin(), get(vertex_color, bgraph), 
				get(vertex_degree, bgraph), get(vertex_priority, bgraph));
    
	std::copy(sloan_order.begin(), sloan_order.end(), *permutation);
	
// 	for (std::vector<SloanVertex>::const_iterator i = sloan_order.begin(); i != sloan_order.end(); ++i)
// 		cout << index_map[*i] << " ";
// 	cout << endl;

//     for (sloan_size_type c = 0; c != sloan_order.size(); ++c)
//       perm[index_map[sloan_order[c]]] = c;
//     std::cout << "  bandwidth: " 
//               << bandwidth(bgraph, make_iterator_property_map(&perm[0], index_map, perm[0]))
//               << std::endl;
//     std::cout << "  profile: " 
//               << profile(bgraph, make_iterator_property_map(&perm[0], index_map, perm[0]))
//               << std::endl;
//     std::cout << "  max_wavefront: " 
//               << max_wavefront(bgraph, make_iterator_property_map(&perm[0], index_map, perm[0]))
//               << std::endl;
//     std::cout << "  aver_wavefront: " 
//               << aver_wavefront(bgraph, make_iterator_property_map(&perm[0], index_map, perm[0]))
//               << std::endl;
//     std::cout << "  rms_wavefront: " 
//               << rms_wavefront(bgraph, make_iterator_property_map(&perm[0], index_map, perm[0]))
//               << std::endl;
    
	
}
