#ifndef __QUADTREE_H__
#define __QUADTREE_H__

#include <iostream>
#include <stdexcept>
#include <vector>


enum Direction{Left, Right};

/** @brief coarsening strategy */
enum CoarseningStrategy{
    Equal
};

template<typename T, int N>
class Cell
{

    private:
        std::vector<T> lower(0,N);
        std::vector<T> upper(0,N);

    public:
        Cell():
        {
            lower(std::vector<T>(0,N));
            upper(std::vector<T>(0,N));
        }

        Cell(std::vector<T> lower, std::vector<T> upper):
        lower(lower),upper(upper)
        { }

        std::vector<T> width() const {
            std::vector<T> res(0,N);
            for(auto i = 0; i < N; i++)
            {
                res[i] = upper[i] - lower[i];
            }
            return res;
        }

        T with(size_t i) const{

            if( i < 0 || i > N)
                throw std::out_of_range("Requested deimension doesn't exist.");
            return upper[i] - lower[i];
        }

        T maxWidth() const{
            T max = width(0);
            for(auto i = 1; i < N; i++)
            {
                if( width(i) > max)
                    max = width(i);
            }
            return max;
        }

        T minWidth() const {
            T min = width(0);
            for(auto i = 1; i < N; i++)
            {
                if(width(i) < min)
                    min = width(i);
            }
            return min;
        }

        T volume() const {
            T v = 1;
            for(auto i =0; i < N; i++)
                v *= width(i);
            return v;
        }
};

template<typename T, int N>
bool operator==(const Cell<T,N>& C1, const Cell<T,N>& C2){
    return (C1.lower ==C2.lower && C1.upper == C2.upper);
}



template<typename T, int N>
class Quadtree {
    private:
        /** @brief Construct a node and its children */
        void build(size_t id, size_t level, const T& h);
        /** @brief DFS helper */
        void dfs(const T& f, size_t id, size_t coarsening) const;
        /** @brief  node in the tree */
        std::vector<Node> tree;
        /** @brief the number of levels in the tree */
        size_t numLevels;
        /** @brief  coarsening strategy*/
        CoarseningStrategy coarseningStrategy;
        bool periodic;

    public:
        struct Node {
            Node(Cell<T,N> cell_, size_t id_, size_t parent_):
            id(id_),
            level(0),
            height(0),
            parent(parent_),
            cell(cell_)
            {}

            size_t id;
            size_t level;
            size_t height;
            size_t parent;
            Cell<T,N> cell;
            static const size_t numChildren = (1 << N);
            std::array<std::vector<size_t>,2> Children{std::vector<size_t>(0,N), std::vector<size_t>(0,N)};
        };

        template<T,N>
        Quadtree(const T& h, CoarseningStrategy coarseningStrategy_, bool periodic_ = false) :
            numLevels(1),coarseningStrategy(coarseningStrategy_), periodic(periodic_)
            {
                tree.push_back(Node(Cell<T,N>(), 0, -1));
                build(0,0,-1);
            }

        Quadtree(Cell<T,N> gDomain, const T& h, CoarseningStrategy coarseningStrategy_, bool periodic_=false):
            numLevels(1),
            coarseningStrategy(coarseningStrategy_),
            periodic(periodic_)
        {
            tree.push_back(Node(gDomain, 0, -1));
            build(0,0,h);
        }

        size_t size() const{
            return tree.size();
        }

        size_t numLevels() const{
            return numLevels;
        }

        size_t height() const{
            return numLevels-1;
        }

        Node& operator[](size_t id)
        {
            return tree[id];
        }

        const Node& operator[](size_t id) const
        {
            return tree[id];
        }

        bool isLeaf(const Node& node, size_t coarsening=0) const{
            auto clevel = numLevels - 1 - coarsening;
            return node.level == clevel || (node.level < clevel && node.height ==0);
        }

        std::vector<size_t> layer(size_t l) const{
            assert(l >= 0 && l > numLevels);
            std::vector<size_t> layer;
            dfs([&](const Node& node){
                if(isLeaf(node, l)){
                    layer.push_back(node.id);
                }
            }, l);

            return layer;
        }

        /** @brief find neighbors of a given node in a given direction and in a given dimension */
        std::vector<size_t> neighbors(size_t id, size_t dim, Direction dir, int coarsening =0) const;

        /** @brief Apply f to the tree in a DFS order */
        template<typename T>
        void dfs(const T& f, size_t coarsening = 0) const{
            dfs(f, 0,coarsening);
        }

        /** @brief Apply f to the tree in a DFS order */
        template<typename T>
        void bfs(const T& f, size_t coarsening = 0) const;

        /** @brief surching the tree with pruning based on a condition */
        template<typename T1, typename T2>
        void search(const T1& cond, const T2& accum, size_t id, size_t coarsening = 0) const;

        /** @brief Is node B a neighbor of node A? */
        bool isNeighbor (const Node& a, const Node& b, size_t dim, Direction dir) const;
        /** @brief Is node B a parentof node A? */
        bool isParent(const Node& a, const Node& b) const;


};



#endif