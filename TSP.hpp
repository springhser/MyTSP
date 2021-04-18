/*
 * @Descripttion: 
 * @version: 
 * @Author: springhser
 * @Date: 2020-12-21 22:35:55
 * @LastEditors: springhser
 * @LastEditTime: 2021-04-18 16:59:24
 */
#include <bits/stdc++.h>


struct Node
{
    int idx_;           // this index is an unique flag of a node.
    double x_;
    double y_;
    std::string name_; // the name of the node in real world, like city's name 
    Node(const double x = 0.0, const double y = 0.0):x_(x),y_(y){}
    
    friend bool operator==(const Node& lhs, const Node& rhs)
    {
        return lhs.x_ == rhs.x_ && lhs.y_ == rhs.y_;
    }

};
struct Edge
{
    Node ns_;
    Node ne_;
    double length_;
    explicit Edge():ns_(Node()),ne_(Node()){}
    Edge(const Node& ns, const Node& ne):ns_(ns),ne_(ne){}    
    double getLength() const
    {
        return length_;
    }
    
    friend bool operator==(const Edge& lhs, const Edge& rhs) 
    {
        return lhs.ns_ == rhs.ns_ && lhs.ne_ == rhs.ne_;
    }
};
struct Tour
{
    Tour(){}
    bool isPath() const
    {
        return false;
    }
    void calTourCost() 
    {
        cost_ = 0.0;
    }
    double getTourCost() const
    {
        return cost_;
    }

    // 
    bool isNodeConnected(const Node& n1, const Node& n2) const
    {
        return false;
    }

    bool isEdgeInTour(const Edge& e) const
    {
        return edges_.find(e) != edges_.end();
    }
    Node getFirstNode() const
    {
        return nodes_list_[0]; 
    }

    bool isTourEmpty()
    {
        return edges_.empty();
    }
    Node getNodeSucc(const Node& n) const 
    {
        return n;
    }

    Node getNodePrev(const Node& n) const 
    {
        return n;
    }

    std::vector<Node> getAdjacent(const Node& n)
    {
        return {getNodePrev(n), getNodeSucc(n)};
    }
    bool isEqual(const Tour& tour) const
    {
        return false;
    }
    std::vector<Node> nodes_list_;
    std::set<Edge> edges_; // to be removed
    double cost_;
};

/*
思路：
* 维护两个边的集合：remove集（待删除边集合，简称RSet）和add集（带加入的集合, 简称A_Set）。
算法简要流程：
1. 选择起点n1.
2. 选n1的前驱或者后继节点为n2，组成第一条待删除的边(n1,n2)进入R_Set.
3. 选n2的邻近节点n3
    (
    要求(n2,n3):
        （1）不属于原路径上的边。
        （2）不在R_Set中。
        （3）长度小于(n1,n2)
    )。
    (n2,n3)进入A_Set。
4. 从n3的前驱或者后继中选择节点n4，(n3,n4)进入R_Set.
5. 若n4和n1连接，即(n1,n2) (n3,n4)-> (n1,n4)(n2,n3) ，
   能形成一条路径，且使得dist(n1,n2)+ dist(n3,n4)> dist(n1,n4)+dist(n2,n3) ，
   则得到一条新的路径。
6. 否则，(n1,n4)进入A_Set。从n4出发，按照2，3，4的步骤重新找待删除的边和待添加边。
   这里，由于寻找越多，计算复杂度越大，通常当R_Set的大小超过5，退出搜索，
   表明从n1出发找不到更合适的路径。
*/

/*
Key function:
    1. get the previous and successor node of a specify node
    2. generate an edge
    3. calculate a distiance of an edge
    4. generate a path
    5. compare two edge
    6. compare two path
    7. get path length
    8. get edge length

*/

class TSP
{
public:
    TSP(const Tour& init_tour): tour_(init_tour)
    {
    
    }
    
    Tour getOptTour()
    {
    // select the first node n1;
    Node n1 = tour_.getFirstNode();
    
    // get succ or prev of n1;
    std::vector<Node> prv_succ= tour_.getAdjacent(n1);
    for(auto& n2: prv_succ)
    {
        Edge e1(n1, n2);
        if(set_R_.find(e1) == set_R_.end())
        {
            set_R_.insert(e1);
        }
        else
        {
            continue;
        }
        

        // get the neighbor node of n2.
        std::vector<Node> neighs = getNeighbor(n2);

        for(auto& neigh:neighs)
        {
            
        }
    }

    return lk_tour_;
    }

    std::vector<Node> getNeighbor(const Node& n2)
    {
        return std::vector<Node>();
    }

private:
    std::set<Edge> set_R_; // the set of edge that will be removed
    std::set<Edge> set_A_; // the set of edge that will be added
    
    int recur_deep_;  // the recurve deep
    Tour tour_;
    Tour lk_tour_;

};
