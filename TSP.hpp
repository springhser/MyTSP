/*
 * @Descripttion: 
 * @version: 
 * @Author: springhser
 * @Date: 2020-12-21 22:35:55
 * @LastEditors: springhser
 * @LastEditTime: 2021-05-23 12:11:04
 */
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <string>
#include <list>
#include <array>
#include <iostream>
#include <algorithm>
#include <math.h>



#define RECURSION_DEPTH 5

struct Node
{
    Node(const double x = 0.0, const double y = 0.0):x_(x),y_(y){}
    
    friend bool operator==(const Node& lhs, const Node& rhs)
    {
        return lhs.unq_idx_ == rhs.unq_idx_;
    }

    int unq_idx_;           // this index is an unique flag of a node.
    int order_;            // this is the order in a specify tour
    double x_;
    double y_;
    std::string name_; // the name of the node in real world, like city's name 

    Node* prev;
    Node* next;

};

static std::vector<Node> g_node_list;

using Matrix = std::vector<std::vector<int>>;

static Matrix g_node_dist_mat; 
struct Edge
{
    explicit Edge():ns_(Node()),ne_(Node())
    {
        calLength_();
    }
    Edge(const Node& ns, const Node& ne):ns_(ns),ne_(ne)
    {
        calLength_();
    }
        
    double getLength() const
    {
        return length_;
    }
    
    friend bool operator==(const Edge& lhs, const Edge& rhs) 
    {
        return (lhs.ns_ == rhs.ns_ && lhs.ne_ == rhs.ne_) ||
               (lhs.ne_ == rhs.ns_ && lhs.ns_ == rhs.ne_); // ignore the direction
    }

private:
    void calLength_()
    {
        length_ = sqrt((ns_.x_-ne_.x_)*(ns_.x_-ne_.x_)+(ns_.y_-ne_.y_)*(ns_.y_-ne_.y_));
    }
private:
    Node ns_;
    Node ne_;
    double length_;
};
struct Tour
{
    Tour(int node_size = 0):node_size_(node_size){}
    /**
     * @brief Determine whether it is a valid path
     * @return yes: true; no: false
     */
    bool isPath() const
    {
        Node* node_tra = head_node_;
        std::unordered_set<int> node_idx;
        while(node_tra)
        {
            if(node_idx.find(node_tra->unq_idx_)!=node_idx.end())
            {
                return node_idx.size() == node_size_;
            }
            node_idx.insert(node_tra->unq_idx_);
            node_tra = node_tra->next;
        }
        return node_idx.size() == node_size_;
    }

    /**
     * @brief calculate the cost of the tour.This function
     *        maybe will not be used very frequently 
     * @return the cost of the tour
     * 
     * @warning It should esure the tour is a vaild path, 
     *          or the cost will be set to -1
     */
    void calTourCost() // 
    {
        if(!isPath())
        {
            cost_ = -1;
            return;
        }
        double t_cost = 0.0;
        Node* node_tra = head_node_;
        while(node_tra && node_tra->next )
        {
            if(node_tra->next->unq_idx_ == head_node_->unq_idx_)
            {
                t_cost = t_cost + g_node_dist_mat[node_tra->next->unq_idx_][head_node_->unq_idx_];
                break;
            }
            t_cost = t_cost + g_node_dist_mat[node_tra->next->unq_idx_][node_tra->unq_idx_];
            node_tra = node_tra->next;
        }
        cost_ = t_cost;
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
        return false;
        // return edges_.find(e) != edges_.end();
    }                     
    Node getFirstNode() const
    {
        return nodes_list_[0]; 
    }

    bool isTourEmpty()
    {
        return false;
        // return edges_.empty();
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
    bool removeEdge(const Edge& e)
    {
        return true;
    }
    bool addEdge(const Edge& e)
    {
        return true;
    }
    std::vector<Node> nodes_list_;

    std::unordered_map<int, Node*> tour_map_; 
    Node* head_node_;
    
    double cost_;

    int node_size_;
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
//     TSP(const Tour& init_tour): tour_(init_tour),recur_depth_(0)
//     {
//     }
    
//     Tour getOptTour()
//     {
//     // select the first node n1;
//     Node n1 = tour_.getFirstNode();
    
//     doOpt(n1);
//     lk_tour_ = tour_;
//     return lk_tour_;
//     }

//     bool doOpt(const Node& n1)
//     {
//     // get succ or prev of n1;
//         std::vector<Node> prv_succ= tour_.getAdjacent(n1);

//         bool get_new_tour_flag = false;
//         for(auto& n2: prv_succ)
//         {
//             Edge e1(n1, n2);
//             if(set_R_.find(e1) == set_R_.end())
//             {
//                 set_R_.insert(e1);
//             }
//             else
//             {
//                 continue;
//             }
            
//             // get the neighbor node of n2.
//             std::vector<Node> neighs = getNeighbor(n2);
//             Node n3;

//             for(auto& neigh:neighs)
//             {
//                 Edge e(n2,neigh);
//                 if(tour_.isEdgeInTour(e))
//                 {
//                     continue;
//                 }
                
//                 if(isEdgeInRSet(e))
//                 {
//                     continue;
//                 }

//                 if(e.getLength() >= e1.getLength())
//                 {
//                     continue;
//                 }

//                 // add the new edge to set_A_
//                 set_A_.insert(e);
//                 n3 = neigh;
//                 if(get_new_tour_flag = doSelection(n3, n1))
//                 {
//                     break;
//                 }
//             }
//         }
//         return true;
//     }
    
//     bool doSelection(const Node& n3, const Node& n1 )
//     {
//         std::vector<Node> prv_succ= tour_.getAdjacent(n3);
//         for(auto& n4: prv_succ)
//         {
//             Edge e(n1, n4);
//             if(tour_.isEdgeInTour(e))
//             {
//                 continue;
//             }
            
//             if(isEdgeInRSet(e))
//             {
//                 continue;
//             }

//             Tour temp_tour = tour_;
//             for(auto& e_r: set_R_)
//             {
//                 temp_tour.removeEdge(e_r);
//             }

//             for(auto& e_a: set_A_)
//             {
//                 temp_tour.addEdge(e_a);
//             }
//             temp_tour.addEdge(e);
            
//             if(temp_tour.isPath() && temp_tour.getTourCost() < tour_.getTourCost())
//             {
//                 tour_ = temp_tour;
//                 recur_depth_ = 0;
//                 return true;
//             }
            
//             if(recur_depth_ > RECURSION_DEPTH)
//             {
//                 recur_depth_ = 0;
//                 return false;
//             }

//             set_A_.insert(e);
//             recur_depth_++;
//             doOpt(n4);
            
//         }
//     }
//     std::vector<Node> getNeighbor(const Node& n2)
//     {
//         return std::vector<Node>();
//     }

//     bool isEdgeInRSet(const Edge& e)
//     {
//         return set_R_.find(e) != set_R_.end();
//     }

//     bool isEdgeInASet(const Edge& e)
//     {
//         return set_A_.find(e) != set_A_.end();
//     }

// private:
//     std::set<Edge> set_R_; // the set of edge that to be removed
//     std::set<Edge> set_A_; // the set of edge that to be added
    
//     int recur_depth_;  // the recurve deep
//     Tour tour_;
//     Tour lk_tour_;

};


struct NodeMap
{

};