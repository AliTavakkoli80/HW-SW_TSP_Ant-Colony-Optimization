#include "systemc.h"
#include <cstdlib>
#include <vector>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <set>
#include <iostream>

using namespace std;

SC_MODULE(hardware) {
public:
    vector<pair<double, double>> cityGraphCoords;               // cities coordinations
    vector<vector<double>> cost;                                // Cost (distance) for going from i to j
    vector<vector<double>> visibility;                          // Visibility from city i to j
    vector<vector<double>> T;                                   // Pheromone intensity cities i and j
    sc_int<64> N;                                                      // Number of cities

    SC_CTOR(hardware)                                                      // Main class of program , general things ...
    {
               
    }

    void init() {
        vector<double> L, M, U;                                 // Putting Zero vectors in cost and T and visibility
        cost.push_back(L);                                      // to avoiding errors
        T.push_back(M);
        visibility.push_back(U);

        for (int i = 0; i < N; i++)                             // Calculating distance between cities
        {
            vector<double> C(N + 1);
            vector<double> t(N + 1);
            vector<double> l(N + 1);
            for (int j = 1; j <= N; j++)
            {
                // cin >> V[j];
                C[j] = sqrt(pow((cityGraphCoords[i].first - cityGraphCoords[j - 1].first), 2) +
                    pow((cityGraphCoords[i].second - cityGraphCoords[j - 1].second), 2));
                t[j] = 1.0;
                if (C[j] != 0)
                {
                    l[j] = 1 / C[j];
                }
            }
            cost.push_back(C);                                  // Setting distances
            T.push_back(t);                                     // Setting default evapuration
            visibility.push_back(l);                            // Setting visibility
        }
    }

    void resetHW() {
        //  Resetting harware module
        cityGraphCoords.clear();
        cost.clear();
        visibility.clear();
        T.clear();
    }

    double tourDist(vector<int> inputPath)                      // Calculating distance to the city where an ant has gone     
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

    void printVector(vector<int> C)                             // Printing a vector in route format
    {
        for (int i = 0; i < C.size(); i++)
            cout << C[i] << "->";
        cout << 1 << " : END" << endl;
    }
};

hardware globalData("globalData");                                 // Global data construct


class Ant                                                       // Ant implementation
{
public:
    vector<int> trail;                                          // Routes an ant has gone
    set<int, greater<int>> available;                           // Cities that are available for ant to go 
    double alpha;                                               // Importance of the pheromone level T
    double beta;                                                // Importance of the visibility
    Ant(double a, double b)
    {
        alpha = a;
        beta = b;
        trail.push_back(1);                                     // Always start from the city (1)
        for (int i = 2; i <= globalData.N; i++)
        {
            available.insert(i);                                // All cities are available except city 1
        }
    }
    void reset()                                                // Reseting an ant 
    {
        vector<int> tmp;
        tmp.push_back(1);
        trail = tmp;                                            // reset to city 1.
        for (int i = 2; i <= globalData.N; i++)
        {
            available.insert(i);
        }
    }

    void pheroRemain()                                          // Remaining event of phermone on routes
    {
        double tourDist = globalData.tourDist(trail);
        const int Q = 100;
        double depositAmount = Q / tourDist;
        int n = trail.size() - 1;
        for (int i = 0; i < n; i++)
        {
            globalData.T[trail[i]][trail[i + 1]] += depositAmount;
        }
        globalData.T[trail[n]][trail[0]] += depositAmount;
    }

    vector<int> stop()                                          // Stopping an ant from stepping
    {
        pheroRemain();
        vector<int> temp = trail;
        reset();
        return temp;
    }

    void step()                                                 // Stepping an ant from a city to another 
    {
        int currentCity = trail[trail.size() - 1];              // Calculating current city
        double norm = probabilityNorm(currentCity);             // Callculating sum of all probability of available cities
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
            // Probability  a random ant dont go to the highest probability city 
            //gp = getRand();
            //if (gp <= p)
            //{                                           // move
            //    cout << (gp) << endl;
            //    moved = true;
            //    trail.push_back(*i);
            //    available.erase(*i);
            //    break;
            //}
        }
        if (!moved)
        {
            // make a move to the highest available prob city
            // move to cityHighest
            trail.push_back(cityHighest);
            available.erase(cityHighest);           // moved city is not available anymore
        }
    }

    double getRand()
    {
        //srand(time(0));
        double p = ((double)rand() / ((double)RAND_MAX));
        return p;
    }
    double moveProbability(int i, int j, double norm)                       // Calculating the moving probability of an ant using the formula
    {
        double p = (pow(globalData.T[i][j], alpha)) * (pow(globalData.visibility[i][j], beta));
        p /= norm;
        return p;
    }

    double probabilityNorm(int currentCity)                     // Calculating all available distance probability for the formolla
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
        N = globalData.N;                           // Number of cities
        M = 5;                                      // Number of ants
        for (int i = 0; i < M; i++)                 // Creating ants
        {
            Ant a(alpha, beta);
            ant.push_back(a);
        }
    }

    void run()                                      // Simulation method
    {
        vector<int> finalPath;                      // Final path thett is optimum
        double minTourCost = DBL_MAX, currenTourCost;
        for (int i = 0; i < 20; i++)                // Loop For implementing time
        {
            for (int j = 0; j < (N - 1); j++)       // Number of total routes an ant must go
            {
                for (int k = 0; k < M; k++)         // All ants take one step
                {
                    ant[k].step();                  // going from one city to other
                }
            }
            for (int k = 0; k < M; k++)             // Stopping all ants who get into first city
            {
                vector<int> antsFinalTrial = ant[k].stop();     // Finding best (MIN) route
                //globalData.printVector(antsFinalTrial);
                if (!finalPath.size())
                {
                    finalPath = antsFinalTrial;
                    minTourCost = globalData.tourDist(antsFinalTrial);
                    continue;
                }

                currenTourCost = globalData.tourDist(antsFinalTrial);
                if (currenTourCost < minTourCost)
                {
                    minTourCost = currenTourCost;
                    finalPath = antsFinalTrial;
                }

            }
            for (int k = 1; k <= N; k++)          // Implementing evapurate
            {
                for (int l = 1; l <= N; l++)
                {
                    globalData.T[k][l] *= evaporation;
                }
            }
            
            // for (int m = 0; m < N; m++)
            // {
            //     for (int n = 0; n < N; n++)
            //     {
            //         cout << "global data.T(" << m << "," << n << ") =" << globalData.T[m + 1][n + 1] << endl;
            //     }
            // }
        }
        cout << "Minimum Tour Cost: " << minTourCost << endl;
        globalData.printVector(finalPath);
    }
};

void setCities(pair < double, double > cityGraph[], int n) {
    for (int i = 0; i < n; i++)
    {
        globalData.cityGraphCoords.push_back(cityGraph[i]);             // Setting cities coordinations
    }
}

/*---------------------------Test Bench--------------------------------*/

int sc_main(int argc, char* argv[])
{
    sc_start();

    /*------------------------INIT Vars------------------------------------*/

    double alpha = 1;         // pheromone importance T
    double beta = 2;          // visibility importance
    double evaporation = 0.1; // evaporation rate

    /*------------------------Test1---------------------------------------*/ cout << "\n+< test1 >+" << endl;

    int n1 = 5;
    globalData.N = n1;                                                  // Setting number of citeis

    pair < double, double > cityGraph1[] = { make_pair(1, 0),
                                        make_pair(4, 1),
                                        make_pair(4, 2),
                                        make_pair(1, 3),
                                        make_pair(0, 1.5) };
    setCities(cityGraph1, n1);
    globalData.init();

    //sc_start();

    ACO colony1(alpha, beta, evaporation);
    colony1.run();

    /*------------------------Reset Hardware-------------------*/

    globalData.resetHW();

    /*-------------------------Test2------------------------------------------*/ cout << "\n+< test2 >+" << endl;

    int n2 = 4;
    globalData.N = n2;                                                  // Setting number of citeis
    pair < double, double > cityGraph2[] = { make_pair(0, 0),
                                        make_pair(1, 0),
                                        make_pair(0, 1),
                                        make_pair(1, 1)};

    setCities(cityGraph2, n2);
    globalData.init();

    //sc_start();
    ACO colony2(alpha, beta, evaporation);
    colony2.run();


    return 0;
}

