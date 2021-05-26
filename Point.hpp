/*
 * @Descripttion: 
 * @version: 
 * @Author: springhser
 * @Date: 2021-05-25 20:50:48
 * @LastEditors: springhser
 * @LastEditTime: 2021-05-26 21:48:39
 */
#ifndef POINT_HPP
#define POINT_HPP

#include <vector>
#include <sstream>
#include <iostream>

struct Point2D
{
    Point2D(double px, double py):x(px), y(py){}
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