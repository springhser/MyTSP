/*
 * @Descripttion: 
 * @version: 
 * @Author: springhser
 * @Date: 2020-12-21 22:35:55
 * @LastEditors: springhser
 * @LastEditTime: 2021-05-25 08:12:48
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
#include <utility>
#include <functional>


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
using Edge = std::pair<int, int>;

struct hash_pair
{
    std::size_t operator()(const std::pair<int, int>& pair) const
    {
        std::size_t h1 = std::hash<int>{}(pair.first);
        std::size_t h2 = std::hash<int>{}(pair.second);
        return h1 ^ h2;
    }
};

static Matrix g_node_dist_mat; 

using EdgeSet = std::unordered_set<Edge, hash_pair>;
static double getEdgeLength(int idx1, int idx2)
{
    return g_node_dist_mat[idx1][idx2];
}
// struct Edge
// {
//     explicit Edge():ns_(Node()),ne_(Node())
//     {
//         calLength_();
//     }
//     Edge(const Node& ns, const Node& ne):ns_(ns),ne_(ne)
//     {
//         calLength_();
//     }
        
//     double getLength() const
//     {
//         return length_;
//     }
    
//     friend bool operator==(const Edge& lhs, const Edge& rhs) 
//     {
//         return (lhs.ns_ == rhs.ns_ && lhs.ne_ == rhs.ne_) ||
//                (lhs.ne_ == rhs.ns_ && lhs.ns_ == rhs.ne_); // ignore the direction
//     }

// private:
//     void calLength_()
//     {
//         length_ = sqrt((ns_.x_-ne_.x_)*(ns_.x_-ne_.x_)+(ns_.y_-ne_.y_)*(ns_.y_-ne_.y_));
//     }
// private:
//     Node ns_;
//     Node ne_;
//     double length_;
// };
struct Tour
{
    Tour(int node_size = 0):node_size_(node_size){}

    void deepCopy(const Tour& tour)
    {
        nodes_list_ = tour.nodes_list_;
        for(auto& n:nodes_list_)
        {
            tour_map_[n.unq_idx_] = &n;
        }
        for(auto& n : tour.nodes_list_ )
        {
            tour_map_[n.unq_idx_]->next = tour_map_[n.next->unq_idx_];
            tour_map_[n.unq_idx_]->prev = tour_map_[n.prev->unq_idx_];
        }
        head_node_ = tour_map_[tour.head_node_->unq_idx_];
        node_size_ = tour.node_size_;
        cost_ = tour.cost_;
    }
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
    

    /**
     * @brief get the cost of the tour.
     * @return the cost of the tour
     */
    double getTourCost() const
    {
        return cost_;
    }

    bool isNodeConnected(const Node& n1, const Node& n2) const
    {
        return (n1.prev->unq_idx_ == n2.unq_idx_ ||
                n1.next->unq_idx_ == n2.unq_idx_);
    }

    bool isEdgeInTour(const Edge& e) const
    {
        return false;
    }                     
    Node getFirstNode() const
    {
        return *head_node_; 
    }

    Node getNodeByIdx(int uique_idx)
    {
        return *tour_map_[uique_idx];
    }

    bool isTourEmpty()
    {
        return false;
    }
    Node getNodeSucc(const Node& n) const 
    {
        return *n.next;
    }

    Node getNodePrev(const Node& n) const 
    {
        return *n.prev;
    }

    std::vector<Node> getAdjacent(const Node& n)
    {
        return {getNodePrev(n), getNodeSucc(n)};
    }


    /**
     * @brief get the neighbor node of the cur_node.
     *        requirement:
     *        1. 
     *        2.
     *        3.
     * @return the cost of the tour
     */
    bool getNeighborNode(const Node& cur_node, const EdgeSet& R_set, double length, Node& nnode)
    {
        for(auto& n: nodes_list_)
        {
            // not the adjacent
            if(n.unq_idx_ == cur_node.unq_idx_ ||
               cur_node.next->unq_idx_ == n.unq_idx_ ||
               cur_node.prev->unq_idx_ == n.unq_idx_)
            {
                continue;
            }
            
            Edge t_e = n.unq_idx_ > cur_node.unq_idx_? 
                                std::pair<int, int>(cur_node.unq_idx_, n.unq_idx_):
                                std::pair<int, int>(n.unq_idx_, cur_node.unq_idx_);
            if(R_set.find(t_e) != R_set.end())
            {
                continue;
            }
            if(getEdgeLength(n.unq_idx_, cur_node.unq_idx_) >= length)
            {
                continue;
            }

            nnode = *tour_map_[n.unq_idx_];
            return true;

        }
        return false;
    }

    bool isEqual(const Tour& tour) const
    {
        return false;
    }

    bool relinkTour(const EdgeSet& R_set,
                    const EdgeSet& A_set)
    {
        Node* n1 = nullptr;
        Node* n2 = nullptr;
        // remove edge
        for(auto& e : R_set)
        {
            if(tour_map_.find(e.first) == tour_map_.end() || 
               tour_map_.find(e.second) == tour_map_.end())
            {
                return false;
            }
            n1 = tour_map_[e.first];
            n2 = tour_map_[e.second];
            
            if(n1->next->unq_idx_ == n2->prev->unq_idx_)
            {
                n1->next = nullptr;
                n2->prev = nullptr;
            }
            else
            {
                n1->prev = nullptr;
                n2->next = nullptr;
            }
            cost_ = cost_ - getEdgeLength(n1->unq_idx_, n2->unq_idx_);
        }

        // add edge
        for(auto& e : A_set)
        {
            if(tour_map_.find(e.first) == tour_map_.end() || 
               tour_map_.find(e.second) == tour_map_.end())
            {
                return false;
            }
            n1 = tour_map_[e.first];
            n2 = tour_map_[e.second];
            
            if(n1->next->unq_idx_ == n2->prev->unq_idx_)
            {
                n1->next = n2;
                n2->prev = n1;
            }
            else
            {
                n1->prev = n2;
                n2->next = n1;
            }
            cost_ = cost_ + getEdgeLength(n1->unq_idx_, n2->unq_idx_);
        }
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

class TSP
{
public:
    TSP(const Tour& init_tour): tour_(init_tour),recur_depth_(0)
    {
    }
    
    Tour getOptTour()
    {
        bool is_success = false;
        for(auto& n: tour_.nodes_list_)
        {
            Tour temp_tour;
            temp_tour.deepCopy(tour_);
            // select the first node n1;
            Node n1 = temp_tour.getNodeByIdx(n.unq_idx_);
            if(is_success = doOpt(n1, n1, temp_tour))
            {
                lk_tour_.deepCopy(temp_tour);
                break;
            }
        }
        
        return lk_tour_;
    }

    bool doOpt(const Node& n1, const Node& origin_node, Tour& temp_tour)
    {
        // get succ or prev of n1;
        std::vector<Node> prv_succ= tour_.getAdjacent(n1);

        bool get_new_tour_flag = false;
        for(auto& n2: prv_succ)
        {
            Edge e1(n1.unq_idx_, n2.unq_idx_);
            if(set_R_.find(e1) == set_R_.end())
            {
                set_R_.insert(e1);
            }
            else
            {
                continue;
            }
            
            // get the neighbor node of n2.
            Node n3;
            if(!tour_.getNeighborNode(n2, set_R_, getEdgeLength(n1.unq_idx_, n2.unq_idx_), n3))
            {
                return false;
            }

            Edge e(n2.unq_idx_,n3.unq_idx_);
            // add the new edge to set_A_
            set_A_.insert(e);

            if(get_new_tour_flag = doSelection(n3, n1, temp_tour))
            {
                break;
            }
            
        }
        return true;
    }
    
    bool doSelection(const Node& n3, const Node& n1, Tour& temp_tour)
    {
        std::vector<Node> prv_succ= tour_.getAdjacent(n3);
        bool get_new_tour_flag = false;
        for(auto& n4: prv_succ)
        {
            Edge e(n1.unq_idx_, n4.unq_idx_);
            if(tour_.isEdgeInTour(e))
            {
                continue;
            }
            
            if(isEdgeInRSet(e))
            {
                continue;
            }


            temp_tour.relinkTour(set_R_, set_A_);
            
            if(temp_tour.isPath() && temp_tour.getTourCost() < tour_.getTourCost())
            {
                tour_ = temp_tour;
                recur_depth_ = 0;
                return true;
            }
            
            if(recur_depth_ > RECURSION_DEPTH)
            {
                recur_depth_ = 0;
                return false;
            }

            set_A_.insert(e);
            recur_depth_++;
            if(get_new_tour_flag = doOpt(n4, n1, temp_tour))
            {
                break;
            }
        }

        return false;
    }


    bool isEdgeInRSet(const Edge& e)
    {
        return set_R_.find(e) != set_R_.end();
    }

    bool isEdgeInASet(const Edge& e)
    {
        return set_A_.find(e) != set_A_.end();
    }

private:
    EdgeSet set_R_; // the set of edge that to be removed
    EdgeSet set_A_; // the set of edge that to be added
    
    int recur_depth_;  // the recurve deep
    Tour tour_;
    Tour lk_tour_;

};
