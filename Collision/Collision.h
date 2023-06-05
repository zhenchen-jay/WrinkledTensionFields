#ifndef COLLISION_H
#define COLLISION_H

#include <Eigen/Dense>
#include <vector>

struct BoundingBox
{
    double mins[3];
    double maxs[3];

    bool intersects(const BoundingBox &other) const
    {
        for (int axis = 0; axis < 3; axis++)
        {
            if (mins[axis] > other.maxs[axis] || maxs[axis] < other.mins[axis])
                return false;
        }
        return true;
    }
};

struct AABB
{
    BoundingBox bbox;
    virtual ~AABB() {}

    virtual void gatherIntersections(const BoundingBox &isectBox, std::vector<int> &isects) const = 0;
    void intersect(const BoundingBox &isectBox, std::vector<int> &isects) const
    {
        if (bbox.intersects(isectBox))
        {
            gatherIntersections(isectBox, isects);
        }
    }
};

struct AABBLeaf : public AABB
{    
    int face;

    virtual void gatherIntersections(const BoundingBox &isectBox, std::vector<int> &isects) const
    {
        isects.push_back(face);
    }
};

struct AABBNode : public AABB
{
    AABBNode() : left(nullptr), right(nullptr) {};
    ~AABBNode() { delete left; delete right; }    
    virtual void gatherIntersections(const BoundingBox &isectBox, std::vector<int> &isects) const
    {
        left->intersect(isectBox, isects);
        right->intersect(isectBox, isects);
    }
    
    AABB *left, *right;
};

bool insidebox (const Eigen::Vector3d q, const Eigen::MatrixXd &V1, const Eigen::MatrixXd &V2);
AABB *buildAABB(std::vector<std::pair<int, BoundingBox> > &boxes);
AABB *buildAABB(const Eigen::MatrixXd &Vstart, const Eigen::MatrixXd &Vend, const Eigen::MatrixXi &F, double padding);

BoundingBox wrapTriangle(const Eigen::MatrixXd &Vstart, const Eigen::MatrixXd &Vend, const Eigen::MatrixXi &F, int face, double padding);


#endif
