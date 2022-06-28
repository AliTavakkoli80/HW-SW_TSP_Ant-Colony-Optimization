 //All systemc modules should include systemc.h header file
#include "systemc.h"
#include <vector>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <set>
#include <iostream>

using namespace std;
 //Hello_world is module name
 SC_MODULE(hardware) {
	SC_CTOR(hardware) {

		// Nothing in constructor
	}
	SC_METHOD (calculate) {

		//Print "Hello World" to the console.

		cout << "Hello World.\n";

		int age;

		cin >> age;
	}
};


 class Data
 {
 public:
     vector<pair<double, double>> cityGraphCoords;
     vector<vector<double>> cost; // Cost for going from i to j
     vector<vector<double>> visibility;
     vector<vector<double>> T; // Pheromone intensity cities i and j
     int N;
     // char s[30];

     Data()
     { 
         N = 5;

         cityGraphCoords.push_back(make_pair(1, 0));
         cityGraphCoords.push_back(make_pair(4, 1));
         cityGraphCoords.push_back(make_pair(4, 2));
         cityGraphCoords.push_back(make_pair(1, 3));
         cityGraphCoords.push_back(make_pair(0, 1.5));

         // exit(0);
         vector<double> L, M, U;
         cost.push_back(L);
         T.push_back(M);
         visibility.push_back(U);

         // for (int i = 0; i < N; i++)
         // {
         //     double x, y;
         //     // cout << "Enter city coordinates : ";
         //     // cin >> x >> y;
         //     // cityCoords.push_back(make_pair(x, y));
         // }

         for (int i = 0; i < N; i++)
         {
             vector<double> V(N + 1);
             vector<double> t(N + 1);
             vector<double> l(N + 1);
             // cout << "co : "
             //      << "(" << cityGraphCoords[i].first << "," << cityGraphCoords[i].second << ")\n";
             for (int j = 1; j <= N; j++)
             {
                 // cin >> V[j];
                 V[j] = sqrt(pow((cityGraphCoords[i].first - cityGraphCoords[j - 1].first), 2) +
                     pow((cityGraphCoords[i].second - cityGraphCoords[j - 1].second), 2));
                 // cout << "V: "
                 //      << "(" << cityGraphCoords[i].first << "," << cityGraphCoords[i].second << ") ->"
                 //      << "(" << cityGraphCoords[j - 1].first << "," << cityGraphCoords[j - 1].second << ") : " << V[j] << "\n";
                 t[j] = 1.0;
                 if (V[j] != 0)
                 {
                     l[j] = 1 / V[j];
                 }
             }
             cost.push_back(V);
             T.push_back(t);
             visibility.push_back(l);
         }
     }

     double tourDist(vector<int> inputPath)
     {
         int len = inputPath.size() - 1;
         double tourDist = 0.0;
         for (int i = 0; i < len; i++)
         {
             tourDist += cost[inputPath[i]][inputPath[i + 1]];
         }
         tourDist += cost[inputPath[len]][inputPath[0]];
         return tourDist;
     }

     void printVector(vector<int> C)
     {
         for (int i = 0; i < C.size(); i++)
             cout << C[i] << "->";
         cout << "\n";
     }
 };

 Data globalData; // global data construct

 class Ant
 {
 public:
     vector<int> trail;
     set<int, greater<int>> available;
     double alpha; // importance of the pheromone level T
     double beta;  // importance of the visibility
     Ant(double a, double b)
     {
         alpha = a;
         beta = b;
         trail.push_back(1); // always start from the city (1)
         for (int i = 2; i <= globalData.N; i++)
         {
             available.insert(i);
         }
     }
     void reset()
     {
         vector<int> tmp;
         tmp.push_back(1);
         trail = tmp; // reset to nest.
         for (int i = 2; i <= globalData.N; i++)
         {
             available.insert(i);
         }
     }

     void pheroRemain()
     {
         double tourDist = globalData.tourDist(trail);
         const int constNum = 1;
         double depositAmount = constNum / tourDist;
         int n = trail.size() - 1;
         for (int i = 0; i < n; i++)
         {
             globalData.T[trail[i]][trail[i + 1]] += depositAmount;
         }
         globalData.T[trail[n]][trail[0]] += depositAmount;
     }

     vector<int> stop()
     {
         pheroRemain();
         vector<int> temp = trail;
         reset();
         return temp;
     }

     void step()
     {
         int currentCity = trail[trail.size() - 1];
         double norm = probabilityNorm(currentCity);
         double p, gp;
         bool moved = false;
         double highestProb = 0;
         double cityHighest = 0;

         for (set<int>::iterator i = available.begin(); i != available.end(); i++)
         {
             p = moveProbability(currentCity, *i, norm);
             if (p > highestProb)
             {
                 cityHighest = *i;
                 highestProb = p;
             }
             gp = getRand();
             if (gp <= p)
             { // move
                 moved = true;
                 trail.push_back(*i);
                 available.erase(i);
                 break;
             }
         }
         if (!moved)
         {
             // make a move to the highest available prob city
             // move to cityHighest
             trail.push_back(cityHighest);
             available.erase(cityHighest);
         }
     }

     double getRand()
     {
         double p = (rand() / (RAND_MAX + 1.0));
         return p;
     }
     double moveProbability(int i, int j, double norm)
     {
         double p = (pow(globalData.T[i][j], alpha)) * (pow(globalData.visibility[i][j], beta));
         p /= norm;
         return p;
     }

     double probabilityNorm(int currentCity)
     {
         int size = available.size();
         double norm = 0.0;
         for (set<int>::iterator i = available.begin(); i != available.end(); i++)
         {
             norm += (pow(globalData.T[currentCity][*i], alpha)) * (pow(globalData.visibility[currentCity][*i], beta));
         }
         return norm;
     }
 };

 class ACO
 {
 public:
     int N;              // cities
     int M;              // no.of ants
     vector<Ant> ant;    // ants
     double evaporation; // evaporation rate
     double alpha;       // importance of the pheromone level
     double beta;        // importance of the visibility

     ACO(double a, double b, double e)
     {
         alpha = a;
         beta = b;
         evaporation = e;
         N = globalData.N;
         M = 5; // ants
         for (int i = 0; i < M; i++)
         {
             Ant a(alpha, beta);
             ant.push_back(a);
         }
     }

     void run()
     {
         vector<int> finalPath;
         double minTourCost, currenTourCost;
         for (int n = 0; n < 20; n++)
         {
             for (int p = 0; p < (N - 1); p++)
             {
                 for (int i = 0; i < M; i++)
                 {
                     ant[i].step();
                 }
             }
             for (int i = 0; i < M; i++)
             {
                 vector<int> p = ant[i].stop();
                 if (!finalPath.size())
                 {
                     finalPath = p;
                     minTourCost = globalData.tourDist(p);
                     continue;
                 }
                 currenTourCost = globalData.tourDist(p);
                 if (currenTourCost < minTourCost)
                 {
                     minTourCost = currenTourCost;
                     finalPath = p;
                 }
             }
             for (int i = 1; i <= N; i++)
             {
                 for (int j = 1; j <= N; j++)
                 {
                     globalData.T[i][j] *= evaporation;
                 }
             }
             // for (int m = 0; m < N; m++)
             // {
             //     for (int n = 0; n < N; n++)
             //     {
             //         cout << "global data.T(" << m << "," << n << ") =" << globalData.T[m + 1][n + 1] << endl;
             //     }
             // }
             // cout << "***********************************************************************************************" << endl;
         }
         printf("%lf\n", minTourCost);
         globalData.printVector(finalPath);
     }
 };

 int sc_main(int argc, char* argv[])
 {
     double alpha = 1;         // pheromone importance T
     double beta = 2;          // visibility importance
     double evaporation = 0.1; // evaporation rate
     ACO colony(alpha, beta, evaporation);
     colony.run();
     return 0;
 }