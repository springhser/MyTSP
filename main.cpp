/*
 * @Descripttion: 
 * @version: 
 * @Author: springhser
 * @Date: 2021-04-18 23:32:34
 * @LastEditors: springhser
 * @LastEditTime: 2021-06-13 19:59:19
 */
#include "DataReader.hpp"
#include "TSP.hpp"

// void testDataReader()
// {
//     TESTMODULE("read data")
//     TESTCASE("read data")
//     std::ifstream infile("Data/att48.tsp");
//     Points point_list;
//     DataReader::addContents(infile, point_list);

//     helptool::PrintContainer<Point2D> pprint;
//     pprint.printVecMultiLines(point_list);
// }

// void testInitTSP()
// {
//     TESTMODULE("Initialise TSP")
//     std::ifstream infile("Data/test.tsp");
//     Points point_list;
//     DataReader::addContents(infile, point_list);
//     TSP tsp(point_list);
//     TESTCASE("print node distance matrix")
//     tsp.printTour();
// }

void testDoOpt()
{
    TESTMODULE("Execute TSP")
    // std::ifstream infile("Data/test.tsp");
    std::ifstream infile("Data/att48.tsp");
    Points point_list;
    DataReader::addContents(infile, point_list);
    TSP tsp(point_list);
    tsp.printCost();
    TESTCASE("Optimise the Tour")
    tsp.optTour();
    tsp.printCost();
    tsp.printRes();
    
}
void testAll()
{
    // testDataReader();
    // testInitTSP();
    // testDoOpt();
#ifdef NO_TEST
    EdgeTest e_test;
    e_test.test();
    TourTest tour_test;
    tour_test.test();
    TSPTest tsp_test;
    tsp_test.test();
#endif
    helptool::Timer<std::chrono::milliseconds> timer;
    testDoOpt();


    
}

int main()
{
    testAll();
    return 0;
}
