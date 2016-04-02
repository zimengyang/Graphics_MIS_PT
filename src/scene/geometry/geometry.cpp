#include <scene/geometry/geometry.h>

float Geometry::RayPDF(const Intersection &isx, const Ray &ray)
{
    //The isx passed in was tested ONLY against us (no other scene objects), so we test if NULL
    //rather than if != this.
    if(isx.t <0 || isx.object_hit == NULL)
    {
        return 0;
    }

    //Add more here
    Intersection isxOnGeometry = this->GetIntersection(ray);
    float r = glm::distance(isx.point, isxOnGeometry.point);
    float cosTheta = glm::dot(-ray.direction, isxOnGeometry.normal);

    float pdf = r * r / cosTheta / area;
    return pdf;

//    float r2 = glm::length2(isx.point - ray.origin);
//    float costheta = glm::dot(isx.normal, -ray.direction);
//    float pdf = r2 / costheta / area;
//    return pdf;
}


glm::vec3 Geometry::toLocalDirection(const glm::vec3 &a)
{
    glm::vec4 direct(a, 0.0f);
    return glm::normalize(glm::vec3(this->transform.invT() * direct));
}

glm::vec3 Geometry::toLocalPoint(const glm::vec3 &a)
{
    glm::vec4 direct(a, 1.0f);
    return glm::vec3(this->transform.invT() * direct);
}


glm::vec3 Geometry::toWorldDirection(const glm::vec3 &a)
{
    glm::vec4 direct(a, 0.0f);
    return glm::normalize(glm::vec3(this->transform.T() * direct));
}

glm::vec3 Geometry::toWorldPoint(const glm::vec3 &a)
{
    glm::vec4 direct(a, 1.0f);
    return glm::vec3(this->transform.T() * direct);
}
