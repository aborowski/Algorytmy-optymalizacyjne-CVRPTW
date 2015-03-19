//Andrzej Borowski 109719
#include <stdlib.h>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <ctime>

using namespace std;
//--------------------------------------------VERTEX - CUSTOMER
class vert {
public:
	int vertNumber;
	double xCoord, yCoord, demand, readyTime, dueDate, unloadTime;
	vector<double> distances;
	vert(int, double, double, double, double, double, double);
};

vert::vert(int id, double x, double y, double demand, double ready, double due, double unload) {
	this->vertNumber = id;
	this->xCoord = x;
	this->yCoord = y;
	this->demand = demand;
	this->readyTime = ready;
	this->dueDate = due;
	this->unloadTime = unload;
}

//-----------------------------------------------TRUCK
class truck {
public:
	//Basic route data
	double usedCap, distanceTravelled;
	vector<int> vertsVisited;
	truck(double, double, unsigned int);
	//Data calculated and used for merging purposes
	unsigned int bestMergeOption;
	double bestMergeSavings;
};

truck::truck(double capacity, double distance, unsigned int vertex) {
	this->usedCap = capacity;
	this->distanceTravelled = distance;
	this->bestMergeOption = -1;
	this->bestMergeSavings = -1;
	this->vertsVisited.push_back(vertex);
}

//-------------------------------------GLOBAL VARIABLES
int vertexCount, routes, vehicles;
double truckCapacity;
vector<truck*> truckVector;
vector<vert*> vertexVector;

//-------------------------------------SET DISTANCE VERT TO VERT
//Calculates distances between every set of two vertices
void vertexDistanceCalculator() {
	int index = 0;
	for (; index < vertexCount; index++) {
		vert *temp = vertexVector.at(index);
		int other = 0;
		for (; other < vertexCount; other++) {
			if (other != index) {
				double calDistance = (double) sqrt(pow(vertexVector.at(other)->xCoord - temp->xCoord, 2.0) + pow(vertexVector.at(other)->yCoord - temp->yCoord, 2.0));
				temp->distances.push_back(calDistance);
			}
			else
				temp->distances.push_back(0);
		}
	}
}

//-------------------------------------DEBUG DISTANCES PRINTER
void printDistances() {
	int index = 0;
	for (; index < vertexCount; index++) {
		vert *temp = vertexVector.at(index);
		cout << "Index: " << index << " distances: ";
		int other = 0;
		for (; other < vertexCount; other++) {
			cout << temp->distances.at(other) << " ";
		}
		cout << endl;
	}
}

//-------------------------------------INITIAL ROUTES CONSTRUCTOR
//Creates a single route to every vertex on the map
//Assumes that each vertex's demand CAN be covered with a single truck
//Checks if the initial solution is permitted
//Returns true if permitted, else false
bool initializeRoutes() {
	unsigned int index = 1;
	vert *start = vertexVector.at(0);
	vert *destination;
	double distance, capacityUsed;
	for (; index < (unsigned int) vertexCount; index++) {
		destination = vertexVector.at(index);
		distance = start->distances.at(index);
		//First test distance against due date
		if (distance <= destination->dueDate) {
			//Extends distance to be at least equal to ready time if needs to
			if (distance < destination->readyTime) {
				distance = destination->readyTime;
			}
			//Test against the due date of starting point
			if ((distance = distance + start->distances.at(index) + destination->unloadTime) <= start->dueDate) {
				capacityUsed = destination->demand;
				truckVector.push_back(new truck(capacityUsed, distance, index));
			}
			else
				return false;
		}
		else
			return false;
	}
	routes = truckVector.size();
	return true;
}

//-------------------------------------DEBUG ROUTE PRINTER
void printRoutes() {
	int index = 0;
	truck *current;
	for (; index < routes; index++) {
		current = truckVector.at(index);
		cout << "Route: " << index << " cap: " << current->usedCap << " distance: " << current->distanceTravelled << " verts:" << endl;
		unsigned int vert = 0;
		for (; vert < current->vertsVisited.size(); vert++) {
			cout << current->vertsVisited.at(vert) << " ";
		}
		cout << endl;
	}
}

//-------------------------------------CALCULATE BEST MERGE OPTION
//First and easiest test is to compare capacity; if merged routes exceed maximum capacity the routes cannot be merged
void mergeEvaluation(int index) {
	vector<int> possibleMergeVector;
	vector<double> savings;
	truck *merger = truckVector.at(index); //We're testing possible merges for this truck.
	truck *tested;
	unsigned int other = 0;
	double temporarySavings, temporaryDistance;
	for (; other < (unsigned int) routes; other++) { //Test for each route
		if (other == (unsigned int) index)
			continue; //The same truck index
		else {
			tested = truckVector.at(other);
			if (merger->usedCap + tested->usedCap > truckCapacity)
				continue; //Exceeds capacity
			else {
				//Remove the distance travelled to return home for the merger.
				temporaryDistance = merger->distanceTravelled - vertexVector.front()->distances.at(merger->vertsVisited.back());
				//Now we have to step by step try to add verts from tested route to the merger and see if we make it in the time window.
				unsigned int vertIterator = 0;
				for (; vertIterator < tested->vertsVisited.size(); vertIterator++) {
					if (vertIterator == 0) { //First vertex in the tested's vector.
						//Add to temporary distance the distance between last vertex and the next one attempted to merge
						temporaryDistance += vertexVector.at(merger->vertsVisited.back())->distances.at(tested->vertsVisited.front());
						if (temporaryDistance <= vertexVector.at(tested->vertsVisited.front())->dueDate) { //Arrived before due date.
							if (temporaryDistance < vertexVector.at(tested->vertsVisited.front())->readyTime)
								temporaryDistance = vertexVector.at(tested->vertsVisited.front())->readyTime;
							temporaryDistance += vertexVector.at(tested->vertsVisited.front())->unloadTime;
						}
						else {
							temporaryDistance = -1;
							break; //Out of time window.
						}
					}
					else { //All the other cases
						//Add to temporary distance the distance between last vertex and the next one attempted to merge
						temporaryDistance += vertexVector.at(tested->vertsVisited.at(vertIterator - 1))->distances.at(tested->vertsVisited.at(vertIterator));
						if (temporaryDistance <= vertexVector.at(tested->vertsVisited.at(vertIterator))->dueDate) { //Arrived before due date.
							if (temporaryDistance < vertexVector.at(tested->vertsVisited.at(vertIterator))->readyTime)
								temporaryDistance = vertexVector.at(tested->vertsVisited.at(vertIterator))->readyTime;
							temporaryDistance += vertexVector.at(tested->vertsVisited.at(vertIterator))->unloadTime;
						}
						else {
							temporaryDistance = -1;
							break; //Out of time window.
						}
					}
				}
				if (temporaryDistance == -1) //Failed to make it inside time window
					continue;
				//Attempt to return home in time window.
				temporaryDistance += vertexVector.front()->distances.at(tested->vertsVisited.back());
				if (temporaryDistance > vertexVector.front()->dueDate) //Out of window.
					continue;
				//If we made it this far, it means that this is indeed a possible merging option, now we have to calculate the savings.
				temporarySavings = merger->distanceTravelled + tested->distanceTravelled - temporaryDistance;
				if (temporarySavings > 0) { //If the difference between summed routes of two trucks and the merged one is positive, we have savings.
					possibleMergeVector.push_back(other);	//The position of the route in truckVector is added as possible merge.
					savings.push_back(temporarySavings);	//The savings are added at the same position in a different vector.
				}
			}
		}
	}
	//After the full search is done, we have to pick the winner and write it down for the route.
	if (possibleMergeVector.empty() == false) {
		other = 0;
		temporarySavings = 0;
		unsigned int index2;
		for (; other < possibleMergeVector.size(); other++) {
			if (temporarySavings < savings.at(other)) {
				temporarySavings = savings.at(other);
				index2 = possibleMergeVector.at(other); //index in truckVector will be the same as this.
			}
		}
		//We have picked the best option, now let's write it down to remember.
		merger->bestMergeOption = index2;
		merger->bestMergeSavings = temporarySavings;
	}
}

//-------------------------------------CLEAR BEST MERGE OPTIONS AFTER MERGE
void clearBestOptions() {
	unsigned int index = 0;
	for (; index < truckVector.size(); index++) {
		truckVector.at(index)->bestMergeOption = -1;
		truckVector.at(index)->bestMergeSavings = -1;
	}
}

//-------------------------------------PICK BEST MERGE OPTION AND MERGE
//Returns true if a merge was successful, false if not. Used to exit main loop.
bool makeMerge() {
	unsigned int index = 0;
	int id = -1, idto = -1;
	double savings = 0;
	//Picking the best option from all routes.
	for (; index < truckVector.size(); index++) {
		if (savings < truckVector.at(index)->bestMergeSavings) {
			id = index;
			idto = truckVector.at(index)->bestMergeOption;
			savings = truckVector.at(index)->bestMergeSavings;
		}
	}
	if (id == -1) { //In case no profitable merges can be done - either equal cost or none exist - exit
		return false;
	}
	truck *merger = truckVector.at(id);
	truck *target = truckVector.at(idto);
	//After we have picked the best option we merge it.
	//First the capacity, since it's always the same.
	merger->usedCap += target->usedCap;
	//Now we have to carefully change the distances, taking into consideration time window changes.
	//First, remove the distance it takes the merger to return home.
	double temporaryDistance = merger->distanceTravelled - vertexVector.front()->distances.at(merger->vertsVisited.back());
	//Calculate the distance taken to get to the first merged vertex and unloading there.
	temporaryDistance += vertexVector.at(merger->vertsVisited.back())->distances.at(target->vertsVisited.front());
	if (temporaryDistance < vertexVector.at(target->vertsVisited.front())->readyTime)
		temporaryDistance = vertexVector.at(target->vertsVisited.front())->readyTime;
	temporaryDistance += vertexVector.at(target->vertsVisited.front())->unloadTime;
	//Now add the vertex to merger's vector.
	merger->vertsVisited.push_back(target->vertsVisited.front());
	//Now let's repeat the same thing over for all of the remaining verts.
	index = 1;
	for (; index < target->vertsVisited.size(); index++) {
		temporaryDistance += vertexVector.at(target->vertsVisited.at(index - 1))->distances.at(target->vertsVisited.at(index));
		if (temporaryDistance < vertexVector.at(target->vertsVisited.at(index))->readyTime)
			temporaryDistance = vertexVector.at(target->vertsVisited.at(index))->readyTime;
		temporaryDistance += vertexVector.at(target->vertsVisited.at(index))->unloadTime;
		//Add the vertex index to the merger's vector.
		merger->vertsVisited.push_back(target->vertsVisited.at(index));
	}
	//Now let's add the distance it takes to return home.
	temporaryDistance += vertexVector.front()->distances.at(target->vertsVisited.back());
	//replace the old distance value with new one.
	merger->distanceTravelled = temporaryDistance;
	//Remove the merge target from truckVector.
	truckVector.erase(truckVector.begin() + idto);
	routes--;
	//Finally, clear the best merge options from every remaining route.
	clearBestOptions();
	return true;
}

//-------------------------------------SUM DISTANCES
double sumDistances() {
	int index = 0;
	double sum = 0;
	for (; index < routes; index++) {
		sum += truckVector.at(index)->distanceTravelled;
	}
	return sum;
}

//-------------------------------------SAVING TO FILE
//Depending on whether boolean parameter is true or false,
//prints the full information on a permitted solution or prints -1
//to indicate a solution not permitted.
void saveFile(string fName, bool good) {
	ofstream outputFile(fName);
	if (good == false) { //Solution not permitted
		outputFile << "-1" << endl;
		outputFile.close();
	}
	else { //Solution permitted
		outputFile << routes << " ";
		outputFile << fixed;
		outputFile.precision(5);
		outputFile << sumDistances() << endl;
		unsigned int index = 0;
		for (; index < routes; index++) { //Prints all vertices for each route
			int other = 0;
			for (; other < truckVector.at(index)->vertsVisited.size(); other++) {
				outputFile << truckVector.at(index)->vertsVisited.at(other) << " ";
			}
			outputFile << endl;
		}
		outputFile.close();
	}
}

//-------------------------------------LOADING FROM FILE
//Processes the input file to retrieve a set number of vertices and places them into the vertexVector.
int loadFile(string fName, int verts) {
	string buffer;
	ifstream inputFile;
	inputFile.open(fName);
	if (inputFile.good() == true) {
		//HEADER reading and writing
		//All the console outputs were here for debug purposes only
		inputFile >> buffer;
		//cout << buffer << endl;
		inputFile.ignore(100, '\n');
		inputFile >> buffer;
		//cout << buffer << endl;
		inputFile.ignore(100, '\n');
		getline(inputFile, buffer);
		//cout << buffer << endl;
		inputFile >> vehicles >> truckCapacity;
		//cout << vehicles << "\t" << truckCapacity << endl;

		inputFile.ignore(100, '\n');
		getline(inputFile, buffer);
		//cout << buffer << endl;
		getline(inputFile, buffer);
		//cout << buffer << endl;
		getline(inputFile, buffer);
		//cout << buffer << endl;
		//BODY
		int i, vertNum;
		double x, y, demand, ready, dueDate, unload;
		if (verts < 0) {
			i = 0;
			while (!inputFile.eof()) {
				inputFile >> vertNum >> x >> y >> demand >> ready >> dueDate >> unload;
				//cout << vertNum << "\t" << x << "\t\t" << y << "\t" << demand << "\t" << ready << "\t\t" << dueDate << "\t\t" << unload << endl;
				inputFile.ignore(1, '\n');
				vertexVector.push_back(new vert(vertNum, x, y, demand, ready, dueDate, unload));
				i++;
			}
			vertexVector.pop_back();
			i--;
			return i;
		}
		else {
			i = 0;
			while (i < verts && !inputFile.eof()) {
				inputFile >> vertNum >> x >> y >> demand >> ready >> dueDate >> unload;
				//cout << vertNum << "\t" << x << "\t\t" << y << "\t" << demand << "\t" << ready << "\t\t" << dueDate << "\t\t" << unload << endl;
				inputFile.ignore(1, '\n');
				vertexVector.push_back(new vert(vertNum, x, y, demand, ready, dueDate, unload));
				i++;
			}
			if (inputFile.eof()) {
				i--;
				vertexVector.pop_back();
			}
			return i;
		}
	}
	else {
		cout << "Something is wrong with the file you're trying to open.\n";
		return -1;
	}
}

//-------------------------------------MAIN FUNCTION
int main(int argc, char** argv)
{
	if (argc < 2) {
		cout << "Niewystarczajaca liczba parametrow.\n";
	}
	else {
		if (argc == 3) {
			vertexCount = atoi(argv[2]) + 1;
			vertexCount = loadFile(argv[1], atoi(argv[2]) + 1);
		}
		else {
			vertexCount = -1;
			vertexCount = loadFile(argv[1], -1);
		}
		vertexDistanceCalculator();
		if (initializeRoutes()) { //Initial solution is permitted.
			srand(time(NULL));
			cout << "Initial cost: " << fixed;
			cout.precision(5);
			cout <<sumDistances();
			clock_t begin = clock();
			unsigned int index;
			bool goOn = true;
			int counter = 0;
			while (goOn == true) {// Algorithm loop.		
				if (routes > 50) {
					for (; counter < 50; counter++) { // If there are more than 50 routes remaining algorithm is run for 50 of them chosen at random.
						index = rand() % routes; //Random selection.
						while (truckVector.at(index)->bestMergeOption != -1) { //If the random selection hit a route already used for evaluation, next one is used.
							index = (index + 1) % routes;
						}
						mergeEvaluation(index);
					}
				}
				else { //When there are only 50 or less routes remaining.
					index = 0;
					for (; index < truckVector.size(); index++) { //Calculating best merge options for each route.
						mergeEvaluation(index);
					}
				}
				goOn = makeMerge();
				counter = 0;
			}
			clock_t end = clock();
			double time = double(end - begin) / CLOCKS_PER_SEC;
			if (argc == 3)
				saveFile("res-" + (string)argv[2] + "-" + (string)argv[1], true);
			else
				saveFile("res-full-" + (string)argv[1], true);
			cout << "Seconds taken: " << time;
		}
		else { //Initial solution is not permitted.
			if (argc == 3)
				saveFile("res-" + (string)argv[2] + "-" + (string)argv[1], false);
			else
				saveFile("res-full-" + (string)argv[1], false);
			cout << "Rozwiazanie niemozliwe" << endl;
		}
	}
	return 0;
}

