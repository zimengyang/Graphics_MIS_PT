#include <raytracing/intersection.h>
#include <raytracing/intersectionengine.h>

Intersection::Intersection():
    point(glm::vec3(0)),
    normal(glm::vec3(0)),
    tangent(glm::vec3(0)),
    bitangent(glm::vec3(0)),
    t(-1),
    texture_color(glm::vec3(1.0f))
{
    object_hit = NULL;
}

IntersectionEngine::IntersectionEngine()
{
    scene = NULL;
}

Intersection IntersectionEngine::GetIntersection(Ray r)
{
    Intersection nearest;
    for(Geometry* g : scene->objects)
    {
        Intersection isx = g->GetIntersection(r);
        if((isx.t < nearest.t && isx.object_hit != NULL) || nearest.t < 0)
        {
            nearest = isx;
        }
    }
    return nearest;
}

bool IntersectionComp(const Intersection &lhs, const Intersection &rhs)
{
    return lhs.t < rhs.t;
}

QList<Intersection> IntersectionEngine::GetAllIntersections(Ray r)
{
    QList<Intersection> result;
    for(Geometry* g : scene->objects)
    {
        Intersection isx = g->GetIntersection(r);
        if(isx.t > 0)
        {
            result.append(isx);
        }
    }
    std::sort(result.begin(), result.end(), IntersectionComp);
    return result;
}

glm::vec3 Intersection::ToLocalNormalCoordinate(const glm::vec3& w_world)
{
    glm::vec3 w_local;

    glm::mat4 M;
    M[0] = glm::vec4(tangent,0.0f);
    M[1] = glm::vec4(bitangent,0.0f);
    M[2] = glm::vec4(normal,0.0f);
    M[3] = glm::vec4(0,0,0,1);

    M = glm::transpose(M);

    w_local = glm::vec3(M * glm::vec4(w_world, 0.0f));

    return glm::normalize(w_local);
}

glm::vec3 Intersection::ToWorldNormalCoordinate(const glm::vec3& w_local)
{
    glm::vec3 w_world;

    glm::mat4 M;
    M[0] = glm::vec4(tangent,0.0f);
    M[1] = glm::vec4(bitangent,0.0f);
    M[2] = glm::vec4(normal,0.0f);
    M[3] = glm::vec4(0,0,0,1);

    //M = glm::transpose(M);

    w_world = glm::vec3(M * glm::vec4(w_local, 0.0f));

    return glm::normalize(w_world);
}
