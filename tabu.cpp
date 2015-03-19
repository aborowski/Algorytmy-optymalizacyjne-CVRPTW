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
	~truck();
};

truck::truck(double capacity, double distance, unsigned int vertex) {
	this->usedCap = capacity;
	this->distanceTravelled = distance;
	if (vertex != -1)
		this->vertsVisited.push_back(vertex);
}

truck::~truck() {
	this->vertsVisited.clear();
}

//-------------------------------------SOLUTION CLASS
class solutionTrucks {
public:
	truck* firstTruck;
	truck* secondTruck;
	solutionTrucks(truck*, truck*);
	~solutionTrucks();
};

solutionTrucks::solutionTrucks(truck* fRoute, truck* sRoute) {
	this->firstTruck = fRoute;
	this->secondTruck = sRoute;
}

solutionTrucks::~solutionTrucks() {
	delete this->firstTruck;
	delete this->secondTruck;
}

//-------------------------------------TABU ELEMENT
struct tabuElement {
	int vertOne;
	int vertTwo;
};

//-------------------------------------GLOBAL VARIABLES
int vertexCount, routes, vehicles;
double truckCapacity;
vector<truck*> truckVector;
vector<vert*> vertexVector;
vector<tabuElement*> tabuVector;


//-------------------------------------SET DISTANCE VERT TO VERT
//Calculates distances between every set of two vertices
void vertexDistanceCalculator() {
	int index = 0;
	for (; index < vertexCount; index++) {
		vert *temp = vertexVector.at(index);
		int other = 0;
		for (; other < vertexCount; other++) {
			if (other != index) {
				double calDistance = (double)sqrt(pow(vertexVector.at(other)->xCoord - temp->xCoord, 2.0) + pow(vertexVector.at(other)->yCoord - temp->yCoord, 2.0));
				temp->distances.push_back(calDistance);
			}
			else
				temp->distances.push_back(0);
		}
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
	for (; index < (unsigned int)vertexCount; index++) {
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

//-------------------------------------CLONE ROUTE
truck* cloneRoute(truck* route) {
	truck* newRoute;
	if (route->vertsVisited.empty())
		newRoute = new truck(route->usedCap, route->distanceTravelled, -1);
	else {
		newRoute = new truck(route->usedCap, route->distanceTravelled, route->vertsVisited.front());
		int index = 1;
		for (; index < route->vertsVisited.size(); index++) {
			newRoute->vertsVisited.push_back(route->vertsVisited.at(index));
		}
	}
	return newRoute;
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

//-------------------------------------PICK RANDOM VERT
int pickRandomVert(int routeIndex) {
	int random = rand() % truckVector.at(routeIndex)->vertsVisited.size();
	return truckVector.at(routeIndex)->vertsVisited.at(random);
}

//-------------------------------------CALCULATION OF DISTANCE
//Used to calculate distance for a route, or to test if it's viable
//Returns the distance or -1 when not viable.
double distanceCalculation(truck* route) {
	if (route->vertsVisited.empty()) {
		return 0;
	}
	else {
		double tempDist = vertexVector.front()->distances.at(route->vertsVisited.front());
		//Test against due date.
		if (tempDist > vertexVector.at(route->vertsVisited.front())->dueDate)
			return -1;
		else {
			//Extend to ready time if needed.
			if (tempDist < vertexVector.at(route->vertsVisited.front())->readyTime)
				tempDist = vertexVector.at(route->vertsVisited.front())->readyTime;
			tempDist += vertexVector.at(route->vertsVisited.front())->unloadTime;
			int index = 1;
			//Repeat above process for all vertices present in the routes vector.
			for (; index < route->vertsVisited.size(); index++) {
				tempDist += vertexVector.at(route->vertsVisited.at(index - 1))->distances.at(route->vertsVisited.at(index));
				if (tempDist > vertexVector.at(route->vertsVisited.at(index))->dueDate)
					return -1;
				else {
					if (tempDist < vertexVector.at(route->vertsVisited.at(index))->readyTime)
						tempDist = vertexVector.at(route->vertsVisited.at(index))->readyTime;
					tempDist += vertexVector.at(route->vertsVisited.at(index))->unloadTime;
				}
			}
			//Finally add a returning home distance
			tempDist += vertexVector.at(route->vertsVisited.back())->distances.front();
			if (tempDist > vertexVector.front()->dueDate)
				return -1;
			else
				return tempDist;
		}
	}
}

//-------------------------------------PICK BEST SOLUTION
//Input neighbourhood of solutions.
//Output best solution or a null pointer.
solutionTrucks* pickBestChoice(vector<solutionTrucks*> neighbourhood) {
	double tempDist, bestFound = -1;
	int solIndex = -1;
	int index = 0;
	for (; index < neighbourhood.size(); index++) {
		if (bestFound == -1) {
			bestFound = neighbourhood.at(index)->firstTruck->distanceTravelled + neighbourhood.at(index)->secondTruck->distanceTravelled;
			solIndex = index;
		}
		else {
			tempDist = neighbourhood.at(index)->firstTruck->distanceTravelled + neighbourhood.at(index)->secondTruck->distanceTravelled;
			if (tempDist < bestFound) {
				bestFound = tempDist;
				solIndex = index;
			}
		}
	}
	if (bestFound == -1 || solIndex == -1)
		return NULL;
	else
		return neighbourhood.at(solIndex);
}


//-------------------------------------NEIGHBOURHOOD GENERATION
//Generates a neighbourhood of viable vertex arrangements on the chosen 2 routes according
//to two available moves: insertion of a vertex from first route into the second or swapping
//of verts between the two routes.
//Returns the best solution from entire neighbourhood or a null pointer if no viable solutions
//were generated in the neighbourhood.
solutionTrucks* neighbourhoodSearch(unsigned int firstIndex, unsigned int secondIndex, int vertOne, int vertTwo) {
	vector<solutionTrucks*> solutionVector;
	solutionTrucks* solution;
	//Make copies of the chosen routes.
	//Clones contain quick clone ready copies
	truck* cloneFirst = cloneRoute(truckVector.at(firstIndex));
	truck* cloneSecond = cloneRoute(truckVector.at(secondIndex));
	truck* fRoute;
	truck* sRoute;
	truck* tempRoute;
	//First insertion neighbourhood.
	//Technically the first route only loses a vertex so can be calculated once and reused.
	int index = 0;
	int secIndex = 0;
	//Removal of the inserted vert
	for (; index < cloneFirst->vertsVisited.size(); index++) {
		if (cloneFirst->vertsVisited.at(index) == vertOne) {
			cloneFirst->vertsVisited.erase(cloneFirst->vertsVisited.begin() + index);
			break;
		}
	}
	//Prepare capacities so we don't have to alter them furhter for insertion.
	cloneFirst->usedCap -= vertexVector.at(vertOne)->demand;
	cloneSecond->usedCap += vertexVector.at(vertOne)->demand;
	//Calculate distances and see if it's still viable. If it's not, we can save a lot of time.
	cloneFirst->distanceTravelled = distanceCalculation(cloneFirst);
	if ((cloneFirst->distanceTravelled >= 0) && (cloneSecond->usedCap <= truckCapacity)) {
		index = 0;
		//All iterations except the last one.
		for (; index < cloneSecond->vertsVisited.size() - 1; index++) {
			fRoute = cloneRoute(cloneFirst);
			sRoute = cloneRoute(cloneSecond);
			sRoute->vertsVisited.insert(sRoute->vertsVisited.begin() + index, vertOne);
			sRoute->distanceTravelled = distanceCalculation(sRoute);
			if (sRoute->distanceTravelled < 0) {
				delete fRoute;
				delete sRoute;
			}
			else {
				solution = new solutionTrucks(fRoute, sRoute);
				solutionVector.push_back(solution);
			}
		}
		fRoute = cloneRoute(cloneFirst);
		sRoute = cloneRoute(cloneSecond);
		sRoute->vertsVisited.push_back(vertOne);
		sRoute->distanceTravelled = distanceCalculation(sRoute);
		if (sRoute->distanceTravelled < 0) {
			delete fRoute;
			delete sRoute;
		}
		else {
			solution = new solutionTrucks(fRoute, sRoute);
			solutionVector.push_back(solution);
		}
	}
	//Then swap neighbourhood.
	//Let's adapt the cloned routes for quick modification.
	index = 0;
	for (; index < cloneSecond->vertsVisited.size(); index++) {
		if (cloneSecond->vertsVisited.at(index) == vertTwo) {
			cloneSecond->vertsVisited.erase(cloneSecond->vertsVisited.begin() + index);
			break;
		}
	}
	///Adjust the capacities for the swapping move.
	cloneFirst->usedCap += vertexVector.at(vertTwo)->demand;
	cloneSecond->usedCap -= vertexVector.at(vertTwo)->demand;
	//Quick test if capacities are ok, if not then no need to waste time.
	if ((cloneFirst->usedCap <= truckCapacity) && (cloneSecond->usedCap <= truckCapacity)) {
		//We need two loops to generate all available pairs of placements of the chosen verts.
		index = 0;
		for (; index < cloneFirst->vertsVisited.size() + 1; index++) {
			tempRoute = cloneRoute(cloneFirst);
			//Last iteration, inserts at the end.
			if (index == cloneFirst->vertsVisited.size())
				tempRoute->vertsVisited.push_back(vertTwo);
			//All iterations except the last.
			else
				tempRoute->vertsVisited.insert(tempRoute->vertsVisited.begin() + index, vertTwo);
			tempRoute->distanceTravelled = distanceCalculation(tempRoute);
			if (tempRoute->distanceTravelled >= 0) {
				for (; secIndex < cloneSecond->vertsVisited.size() + 1; secIndex++) {
					fRoute = cloneRoute(tempRoute);
					sRoute = cloneRoute(cloneSecond);
					if (secIndex == cloneSecond->vertsVisited.size())
						sRoute->vertsVisited.push_back(vertOne);
					else
						sRoute->vertsVisited.insert(sRoute->vertsVisited.begin() + secIndex, vertOne);
					sRoute->distanceTravelled = distanceCalculation(sRoute);
					if (sRoute->distanceTravelled < 0) {
						delete fRoute;
						delete sRoute;
					}
					else {
						solution = new solutionTrucks(fRoute, sRoute);
						solutionVector.push_back(solution);
					}
				}
			}
			delete tempRoute;
			secIndex = 0;
		}
	}
	solutionTrucks* tempSol = pickBestChoice(solutionVector);
	while (!solutionVector.empty()) {
		if (solutionVector.back() == tempSol)
			solutionVector.pop_back();
		else {
			delete solutionVector.back();
			solutionVector.pop_back();
		}
	}
	delete cloneFirst;
	delete cloneSecond;
	return tempSol;
}

//-------------------------------------MAIN TABU LOOP
//Checks time before each iteration.
//Checks whether chosen vertices are not on tabu list.
//Checks if capacities are not exceeded for at least one move.
void tabuLoop(clock_t startTime) {
	bool endThroughTime = false;
	bool endThroughWorsening = false;
	bool permittedVertChoice;
	int tabuCounter = 0;
	int tabuLimit, index;
	int worseningCount = 0;
	int worseningMax = 25;
	int vertOne, vertTwo, routeOne, routeTwo;
	int tabuLoopBreaker = 0;
	vector<truck*> currentBest;
	tabuElement* element;
	solutionTrucks* solution;
	//Algorithm
	tabuLimit = int(vertexCount / 10);
	if (tabuLimit < 1)
		tabuLimit = 1;
	while (!endThroughTime && !endThroughWorsening) {
		tabuLoopBreaker = 0;
		if (double(clock() - startTime) / CLOCKS_PER_SEC > 180)
			endThroughTime = true;
		if (worseningCount >= worseningMax)
			endThroughWorsening = true;
		if (routes < 2)
			break;
		permittedVertChoice = false;
		while (!permittedVertChoice) {
			permittedVertChoice = true;
			//Pick routes.
			routeOne = rand() % routes;
			routeTwo = rand() % routes;
			while (routeTwo == routeOne)
				routeTwo = rand() % routes;
			//Pick verts.
			vertOne = pickRandomVert(routeOne);
			vertTwo = pickRandomVert(routeTwo);
			//Check on tabu list.
			index = 0;
			for (; index < tabuVector.size(); index++) {
				if ((tabuVector.at(index)->vertOne == vertOne && tabuVector.at(index)->vertTwo == vertTwo)
					||
					(tabuVector.at(index)->vertOne == vertTwo && tabuVector.at(index)->vertTwo == vertOne)) {
					permittedVertChoice = false;
					tabuLoopBreaker++;
				}
			}
			if (tabuLoopBreaker > 100)
				break;
		}
		//When picked verts are allowed.
		//Update tabu list.
		if (tabuCounter >= tabuLimit) {
			delete tabuVector.front();
			tabuVector.erase(tabuVector.begin());
			element = new tabuElement();
			element->vertOne = vertOne;
			element->vertTwo = vertTwo;
			tabuVector.push_back(element);
		}
		else {
			element = new tabuElement();
			element->vertOne = vertOne;
			element->vertTwo = vertTwo;
			tabuVector.push_back(element);
			tabuCounter++;
		}
		//Begin generation.
		solution = neighbourhoodSearch(routeOne, routeTwo, vertOne, vertTwo);
		if (solution != NULL) {
			//Check if better.
			if ((solution->firstTruck->distanceTravelled + solution->secondTruck->distanceTravelled)
				<
				(truckVector.at(routeOne)->distanceTravelled + truckVector.at(routeTwo)->distanceTravelled)) {
				worseningCount = 0;
				while (!currentBest.empty()) {
					delete currentBest.back();
					currentBest.pop_back();
				}
				//Check if we need to remove a route, then replace old routes with new.
				delete truckVector.at(routeTwo);
				truckVector.at(routeTwo) = cloneRoute(solution->secondTruck);
				if (solution->firstTruck->distanceTravelled == 0) {
					delete truckVector.at(routeOne);
					truckVector.erase(truckVector.begin() + routeOne);
					routes--;
				}
				else {
					delete truckVector.at(routeOne);
					truckVector.at(routeOne) = cloneRoute(solution->firstTruck);
				}
			}
			else {
				if (worseningCount == 0) {
					index = 0;
					for (; index < routes; index++) {
						currentBest.push_back(cloneRoute(truckVector.at(index)));
					}
					worseningCount++;
				}
				else {
					if (worseningCount == 20)
						endThroughWorsening;
					else
						worseningCount++;
				}
			}
		}
		delete solution;
	}
	if (endThroughWorsening || worseningCount > 0) {
		truckVector = currentBest;
	}
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
		cout << "B³¹d wczytywania pliku.\n";
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
			cout << sumDistances();
			clock_t begin = clock();
			tabuLoop(begin);
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

