/*
 * @Descripttion: 
 * @version: 
 * @Author: springhser
 * @Date: 2021-04-18 23:32:34
 * @LastEditors: springhser
 * @LastEditTime: 2021-05-26 22:06:32
 */
#include "DataReader.hpp"
#include "TSP.hpp"

void testDataReader()
{
    TESTMODULE("read data")
    TESTCASE("read data")
    std::ifstream infile("Data/att48.tsp");
    Points point_list;
    DataReader::addContents(infile, point_list);

    helptool::PrintContainer<Point2D> pprint;
    pprint.printVecMultiLines(point_list);
}

void testInitTSP()
{
    TESTMODULE("Initialise TSP")
    std::ifstream infile("Data/att48.tsp");
    Points point_list;
    DataReader::addContents(infile, point_list);
    TSP tsp(point_list);
    TESTCASE("print node distance matrix")
    tsp.printTour();
}
void testAll()
{
    testDataReader();
    testInitTSP();
}

int main()
{
    testAll();
    return 0;
}
