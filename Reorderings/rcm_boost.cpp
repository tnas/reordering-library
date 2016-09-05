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
#include <boost/config.hpp>
#include <vector>
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

#include "rcm_boost.h"
#include <boost/config.hpp>
#include <vector>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

using namespace boost;
using namespace std;

typedef adjacency_list<vecS, vecS, undirectedS, 
	property<vertex_color_t, default_color_type,
        property<vertex_degree_t,int> > > Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;
typedef graph_traits<Graph>::vertices_size_type size_type;
typedef std::pair<std::size_t, std::size_t> Pair;


Graph build_boost_graph(const METAGRAPH* mgraph)
{
	int node, neigh, degree, num_nodes;
	int* neighboors;
	Pair edge;
	
	num_nodes = mgraph->size; 
	Graph bgraph(num_nodes);
	
	for (node = 0; node < num_nodes; ++node)
	{
		degree     = mgraph->graph[node].degree;
		neighboors = mgraph->graph[node].neighboors;
		
		for (neigh = 0; neigh < degree; ++neigh)
		{
			edge = Pair(node, neighboors[neigh]);
			add_edge(edge.first, edge.second, bgraph);
		}
	}
	
	return bgraph;
}



void Boost_RCM(const METAGRAPH* mgraph, int** permut, const int root)
{
	Graph bgraph = build_boost_graph(mgraph);
  
	*permut  = new int[mgraph->size];
	
	graph_traits<Graph>::vertex_iterator ui, ui_end;

	property_map<Graph,vertex_degree_t>::type node_degree = get(vertex_degree, bgraph);
	for (boost::tie(ui, ui_end) = vertices(bgraph); ui != ui_end; ++ui)
		node_degree[*ui] = degree(*ui, bgraph);

// 	property_map<Graph, vertex_index_t>::type index_map = get(vertex_index, bgraph);

// 	std::cout << "original bandwidth: " << bandwidth(bgraph) << std::endl;

	std::vector<Vertex> inv_perm(num_vertices(bgraph));
	std::vector<size_type> perm(num_vertices(bgraph));
	
	Vertex root_vertex = vertex(root, bgraph);
	
	// Reverse Cuthill-Mckee ordering
	cuthill_mckee_ordering(bgraph, root_vertex, inv_perm.rbegin(), 
		get(vertex_color, bgraph), get(vertex_degree, bgraph));
	
// 	cout << "Reverse Cuthill-McKee ordering starting at: " << root_vertex << endl;
// 	cout << "  ";    
	
// 	for (std::vector<Vertex>::const_iterator i = inv_perm.begin(); i != inv_perm.end(); ++i)
// 		cout << index_map[*i] << " ";
// 	cout << endl;

	std::copy(inv_perm.begin(), inv_perm.end(), *permut);
	
// 	for (size_type c = 0; c != inv_perm.size(); ++c)
// 	{
// 		perm[index_map[inv_perm[c]]] = c;
// 	}

// 	cout << "  ";
// 	for (int no = 0; no < mgraph->size; ++no)
// 		cout << (*permut)[no] << " ";
// 	cout <<endl;
	
// 	std::cout << "  bandwidth: " 
// 		  << bandwidth(bgraph, make_iterator_property_map(&perm[0], index_map, perm[0]))
// 		  << std::endl;
}

