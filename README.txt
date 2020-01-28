— This is the implementation of Girvan Newman Algorithm for Detecting Communities in Networks

— The basic Algorithm is as follows
	* Calculate Betweenness for all the edges in the network
	* Remove the edge with the highest betweenness
	* Recalculate the betweenness for all the edges affected by the removal
	* Repeat from step2 until no edges remain
	
— Edge betweenness is the no.of the shortest paths that go through an edge in a network. It can be calculated using BFS(undirected) and Dijkstra(directed). Time Complexity : O(m*n).
	M — no.of edges
	n — no.of nodes
	
— The total time complexity is O(m^2*n). Since you run edge betweenness for each time you remove an edge.

— In this implementation we use a numerical index called modularity to evaluate how good a division is. 0 modularity indicates that it is no better than randomly choosing communities. The final answer would be the division with the maximum modularity.



Observation:

The results show that this algorithm is an accurate method with a high degree of success for extracting communities but disadvantage of O(m^2*n) is clearly visible. Speed can be improved by a bit using multi-threading in the BFS step but the difference is not much.
