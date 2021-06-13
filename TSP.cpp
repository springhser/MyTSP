/*
 * @Descripttion: 
 * @version: 
 * @Author: springhser
 * @Date: 2021-05-26 22:10:39
 * @LastEditors: springhser
 * @LastEditTime: 2021-06-13 21:00:19
 */
#include "TSP.hpp"
Matrix Tour::Node_Dist_Mat; 
NodeList Tour::Node_List;


/**
 * @brief Determine whether it is a valid path
 * @return yes: true; no: false
 */
bool Tour::isValidTour() const
{
    if(tour_edges_.size() != Node_List.size())
    {
        return false;
    }
    int cur_node = head_node_; 
    std::unordered_map<int, int> next_nodes;
    EdgeSet tour_edges_temp = tour_edges_;
    while(!tour_edges_temp.empty())
    {
        bool get_flag = false;
        for(auto& e: tour_edges_temp)
        {
            if(e.first == cur_node)
            {
                next_nodes[cur_node] = e.second;
                cur_node = e.second;
                tour_edges_temp.erase(e);
                get_flag = true;
                break;
            }
            if(e.second == cur_node)
            {
                next_nodes[cur_node] = e.first;
                cur_node = e.first;
                tour_edges_temp.erase(e);
                get_flag = true;
                break;
            }
        }
        if(!get_flag)
        {
            return false;
        }
        if(cur_node == head_node_)
        {
            break;
        }
    }
    
    if(next_nodes.size()!=Node_List.size())
    {
        return false;
    }
    
    return true;
}
bool Tour::relinkTour(const EdgeSet& R_set,
                const EdgeSet& A_set)
{
    if(R_set.size() != A_set.size())
    {
        return false;
    }
    EdgeSet tour_edges_temp = tour_edges_;
    for(const auto& e: R_set)
    {
        if(tour_edges_.find(e) != tour_edges_.end())
        {
            tour_edges_.erase(e);
        }
        else
        {
            return false;
        }
    }
    for(const auto& e: A_set)
    {
        if(tour_edges_.find(e) != tour_edges_.end())
        {
            return false;
        }
        else
        {
            tour_edges_.insert(e);
        }
    }
    if(isValidTour())
    {
        reGenTourMap();
    }
    else
    {
        tour_edges_ = tour_edges_temp;
        return false;
    }
    return true;
}
void Tour::reGenTourMap()
{
    int cur_node = head_node_; 
    EdgeSet tour_edges_temp = tour_edges_;
    while(!tour_edges_temp.empty())
    {
        bool get_flag = false;
        for(auto& e: tour_edges_temp)
        {
            if(e.first == cur_node)
            {
                int pre_node = cur_node;
                tour_map_[cur_node].next_idx_ = e.second;
                cur_node = e.second;
                tour_map_[cur_node].prev_idx_ = pre_node;
                
                tour_edges_temp.erase(e);
                break;
            }
            if(e.second == cur_node)
            {
                int pre_node = cur_node;
                tour_map_[cur_node].next_idx_ = e.first;
                cur_node = e.first;
                tour_map_[cur_node].prev_idx_ = pre_node;
                tour_edges_temp.erase(e);
                break;
            }
        }
        if(cur_node == head_node_)
        {
            break;
        }
    }
}




TSP::TSP(const Points& point_list):recur_depth_(0)
{
    // initialise global variable
    Tour::clearTour();
    Tour::initDistMat(point_list);
    Tour::initNodeList(point_list);
    // initialise member variable
    initTour(point_list);
}
void TSP::initTour(const Points& point_list)
{
    point_list_ = point_list;
    tour_.initTour();
}

void TSP::optTour()
{
    bool is_success = false;
    for(auto& n: Tour::Node_List)
    {
        Tour temp_tour;
        temp_tour = tour_;
        set_A_.clear();
        set_R_.clear();
        // select the first node n1;
        if(is_success = doOpt(n.unq_idx_, temp_tour))
        {
            tour_ = temp_tour;
        }
    }
}
bool TSP::doOpt(int n1_idx, Tour& temp_tour)
{
    // get succ or prev of n1;
    std::vector<int> prv_succ= tour_.getAdjacentIdxByIdx(n1_idx);
    bool get_new_tour_flag = false;
    for(auto& n2_idx: prv_succ)
    {
        // at the first the add set and remove set should be cleared.
        set_A_.clear();
        set_R_.clear();
        recur_depth_ = 0;
        Edge e_r(n1_idx, n2_idx);
        if(set_R_.find(e_r) == set_R_.end())
        {
            set_R_.insert(e_r);
        }
        else
        {
            continue;
        }
        
        if(doSelection(n2_idx, n1_idx, temp_tour))
        {
            return true;
        }
    }
    return false;
}

bool TSP::doSelection(int n2_idx, int origin_node_idx, Tour& temp_tour)
{
    if(set_R_.size() > RECURSION_DEPTH || set_A_.size() > RECURSION_DEPTH)
    {
        return false;
    }
    // get the neighbor node of n2.
    std::unordered_set<int> n3_idxes;
    if(!tour_.getNeighborNode(n2_idx, set_R_, set_A_, n3_idxes))
    {
        return false;
    } 
    for(auto n3_idx: n3_idxes)
    {
        if(origin_node_idx == n3_idx)
        {
            continue;
        }
        Edge e(n2_idx, n3_idx);
        
        if(tour_.isEdgeInTour(e))
        {
            continue;
        }
        
        if(isEdgeInRSet(set_R_, e))
        {
            continue;
        }
        if(isEdgeInASet(set_A_, e))
        {
            continue;
        }
        // add the new edge to set_A_
        set_A_.insert(e);
        if(doSelection2(n3_idx, origin_node_idx, temp_tour))
        {
            return true;
        }
        set_A_.erase(e);
    }
    return false;
}

bool TSP::doSelection2(int n3_idx, int origin_node_idx, Tour& temp_tour)
{
    if(set_R_.size() > RECURSION_DEPTH || set_A_.size() > RECURSION_DEPTH)
    {
        return false;
    }
    std::vector<int> prv_succ = tour_.getAdjacentIdxByIdx(n3_idx);
    for(auto& n4_idx : prv_succ)
    {
        if(origin_node_idx == n4_idx)
        {
            continue;
        }
        Edge e_final(origin_node_idx, n4_idx);
        if(tour_.isEdgeInTour(e_final))
        {
            continue;
        }
        
        if(isEdgeInRSet(set_R_, e_final))
        {
            continue;
        }
        if(isEdgeInASet(set_A_, e_final))
        {
            continue;
        }
        Edge e_r(n3_idx,n4_idx);
        set_R_.insert(e_r);
        set_A_.insert(e_final);
        // PRINTN("*******************************")
        // PRINTN("print middle result:"<<recur_depth_)
        // printRSet();
        // printASet();
        // PRINTN("*******************************")
        // PRINTN("")
        if(temp_tour.getEdgeSetLength(set_R_) > temp_tour.getEdgeSetLength(set_A_))
        {
            if(temp_tour.relinkTour(set_R_, set_A_))
            {
                recur_depth_ = 0;
                // PRINTN("Yes get it")
                return true;
            }
        }
        set_A_.erase(e_final);
        
        recur_depth_++;
        if(recur_depth_ > RECURSION_DEPTH)
        {
            // recur_depth_ = 0;
            return false;;
        }
        if(doSelection(n4_idx, origin_node_idx, temp_tour))
        {
            return true;
        }
        set_R_.erase(e_r);
    }
    return false;
}

