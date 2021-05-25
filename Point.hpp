/*
 * @Descripttion: 
 * @version: 
 * @Author: springhser
 * @Date: 2021-05-25 20:50:48
 * @LastEditors: springhser
 * @LastEditTime: 2021-05-25 21:28:50
 */
#ifndef POINT_HPP
#define POINT_HPP

#include <vector>

struct Point2D
{
    Point2D(double px, double py):x(px), y(py){}
    double x;
    double y;
};

using Points = std::vector<Point2D>;


#endif