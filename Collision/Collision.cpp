#include "CTCD.h"
#include "Collision.h"
#include <vector>


bool insidebox (const Eigen::Vector3d q, const Eigen::MatrixXd &V1, const Eigen::MatrixXd &V2)
{
     Eigen::MatrixXd V(V1.rows()+ V2.rows(), 3);
     V << V1, V2;
     Eigen::Vector3d min = V.colwise().minCoeff();
     Eigen::Vector3d max = V.colwise().maxCoeff();
     
     bool inside = true;
     for (int i = 0; i < 3; i++)
     {
         if (q(i) < min(i) || q(i) > max(i))
             inside = false;
     }
 
     return inside;
}

AABB *buildAABB(std::vector<std::pair<int, BoundingBox> > &boxes)
{
    if (boxes.size() == 0)
        return nullptr;
    if (boxes.size() == 1)
    {
        AABBLeaf *result = new AABBLeaf;
        result->face = boxes[0].first;
        result->bbox = boxes[0].second;
        return result;
    }
    else
    {
        double mincentroids[3];
        double maxcentroids[3];
        BoundingBox enclosing;
        for (int j = 0; j < 3; j++)
        {
            mincentroids[j] = std::numeric_limits<double>::infinity();
            maxcentroids[j] = -std::numeric_limits<double>::infinity();
            enclosing.mins[j] = std::numeric_limits<double>::infinity();
            enclosing.maxs[j] = -std::numeric_limits<double>::infinity();
        }
        for (int i = 0; i < boxes.size(); i++)
        {
            for (int j = 0; j < 3; j++)
            {
                double centroid = 0.5*(boxes[i].second.mins[j] + boxes[i].second.maxs[j]);
                mincentroids[j] = std::min(mincentroids[j], centroid);
                maxcentroids[j] = std::max(maxcentroids[j], centroid);
                enclosing.mins[j] = std::min(enclosing.mins[j], boxes[i].second.mins[j]);
                enclosing.maxs[j] = std::max(enclosing.maxs[j], boxes[i].second.maxs[j]);
            }
        }

        int splitaxis = 0;
        double bestdist = 0;
        for (int j = 0; j < 3; j++)
        {
            double dist = maxcentroids[j] - mincentroids[j];
            if (dist > bestdist)
            {
                bestdist = dist;
                splitaxis = j;
            }
        }

        std::sort(boxes.begin(), boxes.end(),
            [splitaxis](const std::pair<int, BoundingBox> &b1, const std::pair<int, BoundingBox> &b2) -> bool
        {
            double c1 = 0.5*(b1.second.mins[splitaxis] + b1.second.maxs[splitaxis]);
            double c2 = 0.5*(b2.second.mins[splitaxis] + b2.second.maxs[splitaxis]);
            return c1 < c2;
        }
        );

        std::vector<std::pair<int, BoundingBox> > leftboxes;
        std::vector<std::pair<int, BoundingBox> > rightboxes;
        for (int i = 0; i < boxes.size() / 2; i++)
            leftboxes.push_back(boxes[i]);
        for (int i = boxes.size() / 2; i < boxes.size(); i++)
            rightboxes.push_back(boxes[i]);
        AABB *left = buildAABB(leftboxes);
        AABB *right = buildAABB(rightboxes);
        AABBNode *result = new AABBNode;
        result->bbox = enclosing;
        result->left = left;
        result->right = right;
        return result;
    }
}

BoundingBox wrapTriangle(const Eigen::MatrixXd &Vstart, const Eigen::MatrixXd &Vend, const Eigen::MatrixXi &F, int face, double padding)
{
    BoundingBox bb;
    for (int j = 0; j < 3; j++)
    {
        bb.mins[j] = std::numeric_limits<double>::infinity();
        bb.maxs[j] = -std::numeric_limits<double>::infinity();
    }
    for (int j = 0; j < 3; j++)
    {
        Eigen::Vector3d vert1 = Vstart.row(F(face, j)).transpose();
        Eigen::Vector3d vert2 = Vend.row(F(face, j)).transpose();
        for (int k = 0; k < 3; k++)
        {
            bb.mins[k] = std::min(bb.mins[k], vert1[k]-padding);
            bb.maxs[k] = std::max(bb.maxs[k], vert1[k]+padding);
            bb.mins[k] = std::min(bb.mins[k], vert2[k]-padding);
            bb.maxs[k] = std::max(bb.maxs[k], vert2[k]+padding);
        }
    }
    return bb;
}

AABB *buildAABB(const Eigen::MatrixXd &Vstart, const Eigen::MatrixXd &Vend, const Eigen::MatrixXi &F, double padding)
{
    std::vector<std::pair<int, BoundingBox> > boxes;
    int nfaces = F.rows();
    for (int i = 0; i < nfaces; i++)
    {
        BoundingBox bb = wrapTriangle(Vstart, Vend, F, i, padding);
        boxes.push_back(std::pair<int, BoundingBox>(i, bb));
    }
    return buildAABB(boxes);
}

