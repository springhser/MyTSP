/*
 * @Descripttion: 
 * @version: 
 * @Author: springhser
 * @Date: 2021-05-25 20:48:06
 * @LastEditors: springhser
 * @LastEditTime: 2021-05-26 21:36:29
 */
#ifndef DATA_READER_HPP
#define DATA_READER_HPP
#include <iostream>
#include <iomanip>  // for setw() and ws
#include <string>
#include <fstream>
#include <cstdlib>
#include<vector>

#include "Point.hpp"
#include "Utils/HelpTool.hpp"

/**
 * @brief read tsp data from file
 * @param ff file handle
 * @param point_list the result 
 * @return true:success, false:failed
 */
class DataReader
{
public:
    static bool addContents(std::ifstream& ff, Points& point_list)
    {
        try
        {
            int x_coord;
            int y_coord;
            int citi_no;
            int num_of_cities;
            ff >> num_of_cities;
            for(int i = 0; i<num_of_cities; i++)
            {
                ff >> citi_no; // useless
                ff >> x_coord;
                ff >> y_coord;
                Point2D p(x_coord, y_coord);
                point_list.push_back(p);
            }
        }
        catch(...)
        {
            PRINTN("read file error");
            return false;
        }
        
        return true;
    }

};

#endif