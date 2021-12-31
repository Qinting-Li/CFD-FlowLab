#include "GeometryUnit.h"

namespace PHSPACE
{
DYNode::DYNode()
{
    node2node  = nullptr;
    node2nodeK = nullptr;
}

DYNode::~DYNode()
{
    delete []node2node;
    delete []node2nodeK;
}

DYFace::DYFace()
{
    face2node = nullptr;
}
DYFace::~DYFace()
{
    delete []face2node;
}

DYCell::DYCell()
{
    cell2face = nullptr;
    cell2node = nullptr;
    cell2cell = nullptr;
}
DYCell::~DYCell()
{
    delete []cell2face;
    delete []cell2node;
    delete []cell2cell;
}

}