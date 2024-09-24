/******************************************************************************
 * File: Trailblazer.cpp
 *
 * Implementation of the graph algorithms that comprise the Trailblazer
 * assignment.
 */

#include "Trailblazer.h"
#include "TrailblazerGraphics.h"
#include "TrailblazerTypes.h"
#include "TrailblazerPQueue.h"
#include "Grid.h"
#include "Vector.h"
#include "Stack.h"
#include "error.h"
#include <limits>
#include "map.h"
#include "set.h"
#include "random.h"
#include "foreach.h"
#include <iostream>

using namespace std;



void connectNodes(Map<Loc, Set<Loc> > &map, Loc &start, Loc &end);

/* Function: shortestPath
 * 
 * Finds the shortest path between the locations given by start and end in the
 * specified world.	 The cost of moving from one edge to the next is specified
 * by the given cost function.	The resulting path is then returned as a
 * Vector<Loc> containing the locations to visit in the order in which they
 * would be visited.	If no path is found, this function should report an
 * error.
 *
 * In Part Two of this assignment, you will need to add an additional parameter
 * to this function that represents the heuristic to use while performing the
 * search.  Make sure to update both this implementation prototype and the
 * function prototype in Trailblazer.h.
 */


//creates graph, where every dot is painted grey, its cost is maximum int and does not have previos dot;
void createGraph(Grid<Node*>& graph){
	for(int r = 0; r < graph.numRows(); r++){
		for(int c = 0; c < graph.numCols(); c++){
			graph[r][c] = new Node;
			graph[r][c]->color = GRAY;
			graph[r][c]->cost = numeric_limits<double>::infinity();
		} 
	}
}

//this function paints neighbors and changes its cost and adds in the queue, if needed
void paintNeighbors(Loc end, Loc curr, Loc neighbor, Grid<Node*>& graph, 
					double costFn(Loc from, Loc to, Grid<double>& world),
					TrailblazerPQueue<Loc>& queue, Grid<double>& world, Node* currNode){
					//double heuristic(Loc start, Loc end, Grid<double>& world)){

				if(graph.inBounds(neighbor.row, neighbor.col)){
				Node* neighNode = graph[neighbor.row][neighbor.col];
				double currCost = costFn(curr,neighbor, world) + currNode->cost;
				if(neighNode->color == GRAY){
				//paint neighbor in yellow
				neighNode->color = YELLOW;
				colorCell(world, neighbor, YELLOW);

				neighNode->cost = currCost;
				neighNode->previos = curr;
				queue.enqueue(neighbor, neighNode->cost);
				//queue.enqueue(neighbor, neighNode->cost + heuristic(neighbor, end, world));
				}else if(neighNode->color == YELLOW && currCost < neighNode->cost){
				//give neighbor appropriate cost
					neighNode->cost = currCost;
					neighNode->previos = curr;
					queue.decreaseKey(neighbor, currCost);
					//queue.decreaseKey(neighbor, currCost + heuristic(neighbor, end, world));
				}
			}
		}
	
void shortestPathInternal(Loc start,
             Loc end,
             Grid<double>& world,
			 double costFn(Loc from, Loc to, Grid<double>& world), Grid<Node*>& graph){
			// double heuristic(Loc start, Loc end, Grid<double>& world)) {

	//makes  every dot painted grey,
	//its cost is maximum int and does not have previos dot;
	createGraph(graph);

	TrailblazerPQueue<Loc> queue;
	graph[start.row][start.col]->color = YELLOW;
	graph[start.row][start.col]->cost = 0;
	queue.enqueue(start, 0);
	//queue.enqueue(start, heuristic(start, end, world));

	while(!queue.isEmpty()){
		Loc curr = queue.dequeueMin();
		Node* currNode = graph[curr.row][curr.col];
		if(currNode->color != GREEN){
			colorCell(world, curr, GREEN);
			currNode->color = GREEN;
				if(curr == end){
					return;
			}
		//now paint current point's neighbor points; it has 8 neighbor;
		for(int r = -1; r <= 1; r++){
			for(int c = -1; c <= 1; c++){
				Loc neighbor;
				neighbor.row = curr.row + r;
				neighbor.col = curr.col + c;
				paintNeighbors(end, curr, neighbor, graph, costFn,queue, world, currNode);
				//paintNeighbors(end, curr, neighbor, graph, costFn,queue, world, currNode, heuristic);
				}
			}
		}
	}
	error("No path found");
}

Vector<Loc> shortestPath(Loc start,
             Loc end,
             Grid<double>& world,
			 double costFn(Loc from, Loc to, Grid<double>& world),
			 double heuristic(Loc start, Loc end, Grid<double>& world)) {
	
	Grid<Node*> graph(world.numRows(), world.numCols());
	
	shortestPathInternal(start, end, world, costFn, graph);
	//shortestPathInternal(start, end, world, costFn, graph, heuristic);

	//go backwards
	Node* first = graph[end.row][end.col];
	Node* curr = first;
	Stack<Loc> path;
	path.push(end);
	while(curr != graph[start.row][start.col]){
		Loc prev = curr->previos;
		path.push(prev);
		curr = graph[prev.row][prev.col];
		}
	Vector<Loc> result;
	while(!path.isEmpty()){
		result.add(path.pop());
	}
    return result;
}


bool isOneClusterRemaining(Map<Loc, Set<Loc> > &map, 
						   Grid<Loc>& grid, int NodeCount){
	if(map.size() == 0)
		error("map should have elements");
	//int numbOfNodes = grid.numCols()*grid.numRows();
	Loc oneLoc = grid.get(0,0);
		if(map[oneLoc].size() == NodeCount)
		return true;
		else return false;
}


Set<Edge> createMaze(int numRows, int numCols) {

	Set<Edge> spanningTree;
	Grid<Loc> graph(numRows, numCols);
	TrailblazerPQueue<Edge> queue;

	//create clusters with map, where key is node, and key is its cluster
	Map<Loc, Set<Loc> > map;

	int NodeCount = 0;

	for(int r = 0; r < graph.numRows(); r++){
		for(int c = 0; c <graph.numCols(); c++){
			Loc curr;
			curr.row = r;
			curr.col = c;
			Set<Loc> cluster;
			cluster.add(curr);
			map.put(curr, cluster); // put nodes in map
			if(graph.inBounds(r, c + 1)){
			Loc first = makeLoc(r, c);
			Loc second = makeLoc(r, c + 1);
			Edge currEdge;

			currEdge = makeEdge(first, second);
			queue.enqueue(currEdge, randomReal(0, 1));
			NodeCount++;
			}
			if(graph.inBounds(r + 1, c)){
			//	cout << "node r = " << r << " c = " << c<< endl;
				//cout << "node second.r = " << r + 1 << " c = " << c << endl;
				Loc first = makeLoc(r, c);
			Loc second = makeLoc(r + 1, c);
			Edge currEdge = makeEdge(first, second);
			queue.enqueue(currEdge, randomReal(0, 1));
			NodeCount++;
			}
		}
	}
	//While there are two or more clusters remaining
	while(!isOneClusterRemaining(map, graph, NodeCount) && !queue.isEmpty()){
		Edge e = queue.dequeueMin();
		//If the endpoints of e are not in the same cluster
			if(!map[e.start].contains(e.end)){
			//Merge the clusters containing the endpoints of e
			connectNodes(map, e.start, e.end);
			//Add e to the resulting spanning tree
			spanningTree.add(e);
		}
	 }
    return spanningTree;
}


void connectNodes(Map<Loc, Set<Loc> > &map, Loc &start, Loc &end){
	Set<Loc> newCluster = map[start] + map[end];
	foreach(Loc curr in map[start]){
		map[curr] = newCluster;
	}
	foreach(Loc curr in map[end]){
		map[curr] = newCluster;
	}
	map[start] = newCluster;
	map[end] = newCluster;
}