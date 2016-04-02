#ifndef RING
#define RING
#include <scene/geometry/geometry.h>

// ring is defined as a big disc with radius of 1 and a empty smaller disc with radius of 0.5. Center is (0,0,0)
// x - y plane defined shape

class Ring : public Geometry
{
public:
    Intersection GetIntersection(Ray r);
    virtual glm::vec2 GetUVCoordinates(const glm::vec3 &point);
    virtual glm::vec3 ComputeNormal(const glm::vec3 &P);
    void create();

    virtual void ComputeArea();
    virtual Intersection SampleOnGeometrySurface(const float &u, const float &v, const glm::vec3 &point);

    virtual void SetNormalTangentBitangent(const glm::vec3& point_local, Intersection& isx);
};

#endif // RING

