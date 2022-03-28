#include "GRHash.h"
#include "Geo_Face.h"

namespace GRIDER_SPACE
{
CAulong gethashkey(const Point3D& p)
{
    return p.hash_key();
}
CAulong gethashkey(const Point3D* p)
{
    return p->hash_key();
}
CAulong gethashkey(const Geo_Face& face)
{
    return face.hash_key();
}
CAulong gethashkey(const Geo_Face* face)
{
    return face->hash_key();
}
bool equal_hash(const Point3D& p1, const Point3D& p2)
{
    return (p1 == p2);
}
bool equal_hash(const Point3D* p1, const Point3D* p2)
{
    return (*p1==*p2);
}
bool equal_hash(const Geo_Face& face1, const Geo_Face& face2)
{
    return (face1 == face2);
}
bool equal_hash(const Geo_Face* face1, const Geo_Face* face2)
{
    return (*face1 == *face2);
}
}
