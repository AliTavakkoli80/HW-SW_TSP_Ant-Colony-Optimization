#include "systemc.h"
#include "TSP_Ant-Colony-Optimization.cpp"


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

    /*-----------------------------------------------------------------------*/

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

    /*-----------------------------------------------------------------------*/

    int n2 = 5;
    globalData.N = n1;
    pair < double, double > cityGraph2[] = { make_pair(1, 0),
                                        make_pair(4, 1),
                                        make_pair(4, 2),
                                        make_pair(1, 3),
                                        make_pair(0, 1.5) };

    setCities(cityGraph2, n2);
    globalData.init();

    //sc_start();
    ACO colony2(alpha, beta, evaporation);
    colony2.run();


    return 0;
}