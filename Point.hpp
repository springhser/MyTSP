/*
 * @Descripttion: 
 * @version: 
 * @Author: springhser
 * @Date: 2021-05-25 20:50:48
 * @LastEditors: springhser
 * @LastEditTime: 2021-06-05 21:45:37
 */
#ifndef POINT_HPP
#define POINT_HPP

#include <vector>
#include <sstream>
#include <iostream>

struct Point2D
{
    Point2D(double px = 0.0, double py = 0.0):x(px), y(py){}
    double x;
    double y;
    friend std::ostream& operator<<(std::ostream& os, const Point2D& p)
    {
        os << p.x << " " << p.y ;
        return os;
    }
};

using Points = std::vector<Point2D>;


#endif