#include <scene/geometry/square.h>

void SquarePlane::ComputeArea()
{
    // use linear transformation matrix to calculate the area of shape
    float detT = glm::determinant(transform.T());
    float r = 1.0f;
    float originArea = r*r;

    area = detT*originArea;
    if(transform.getScale().z != 1.0f)
        area = area/transform.getScale().z;
}

Intersection SquarePlane::SampleOnGeometrySurface(const float &u, const float &v, const glm::vec3 &point)
{
    glm::vec3 hitPoint(u - 0.5f, v - 0.5f, 0.0f);
    Intersection result;
    result.point = glm::vec3(transform.T() * glm::vec4(hitPoint,1.0f));
    result.object_hit = this;
    result.texture_color = Material::GetImageColorInterp(GetUVCoordinates(hitPoint), material->texture);

    //tangent and bitangent
    SetNormalTangentBitangent(hitPoint,result);

    // what's the purpose of t of a sample intersection on light surface?
    result.t = 1.0f;

    return result;
}

void SquarePlane::SetNormalTangentBitangent(const glm::vec3 &point_local, Intersection &isx)
{
    glm::vec3 up(0,1,0);

    glm::vec3 normal_local = glm::vec3(0,0,1);
    glm::vec3 tangent_local = glm::normalize(glm::cross(up,normal_local));
    glm::vec3 bitangent_local = glm::normalize(glm::cross(normal_local,tangent_local));

    // transform to world coordinates
    isx.normal = glm::normalize(glm::vec3(transform.invTransT() * glm::vec4(normal_local,0.0f)));
    isx.tangent = glm::normalize(glm::vec3(transform.T() * glm::vec4(tangent_local,0.0f)));
    isx.bitangent = glm::normalize(glm::vec3(transform.T() * glm::vec4(bitangent_local,0.0f)));
}

Intersection SquarePlane::GetIntersection(Ray r)
{
    //Transform the ray
    Ray r_loc = r.GetTransformedCopy(transform.invT());
    Intersection result;

    //Ray-plane intersection
    float t = glm::dot(glm::vec3(0,0,1), (glm::vec3(0.5f, 0.5f, 0) - r_loc.origin)) / glm::dot(glm::vec3(0,0,1), r_loc.direction);
    glm::vec4 P = glm::vec4(t * r_loc.direction + r_loc.origin, 1);
    //Check that P is within the bounds of the square
    if(t > 0 && P.x >= -0.5f && P.x <= 0.5f && P.y >= -0.5f && P.y <= 0.5f)
    {
        result.point = glm::vec3(transform.T() * P);
        result.object_hit = this;
        result.t = glm::distance(result.point, r.origin);
        result.texture_color = Material::GetImageColorInterp(GetUVCoordinates(glm::vec3(P)), material->texture);

        //Store the tangent and bitangent
        SetNormalTangentBitangent(glm::vec3(P), result);
        return result;
    }
    return result;
}


glm::vec2 SquarePlane::GetUVCoordinates(const glm::vec3 &point)
{
    return glm::vec2(point.x + 0.5f, point.y + 0.5f);
}

glm::vec3 SquarePlane::ComputeNormal(const glm::vec3 &P)
{
        return glm::vec3(0,0,1);
}

void SquarePlane::create()
{
    GLuint cub_idx[6];
    glm::vec3 cub_vert_pos[4];
    glm::vec3 cub_vert_nor[4];
    glm::vec3 cub_vert_col[4];

    cub_vert_pos[0] = glm::vec3(-0.5f, 0.5f, 0);  cub_vert_nor[0] = glm::vec3(0, 0, 1); cub_vert_col[0] = material->base_color;
    cub_vert_pos[1] = glm::vec3(-0.5f, -0.5f, 0); cub_vert_nor[1] = glm::vec3(0, 0, 1); cub_vert_col[1] = material->base_color;
    cub_vert_pos[2] = glm::vec3(0.5f, -0.5f, 0);  cub_vert_nor[2] = glm::vec3(0, 0, 1); cub_vert_col[2] = material->base_color;
    cub_vert_pos[3] = glm::vec3(0.5f, 0.5f, 0);   cub_vert_nor[3] = glm::vec3(0, 0, 1); cub_vert_col[3] = material->base_color;

    cub_idx[0] = 0; cub_idx[1] = 1; cub_idx[2] = 2;
    cub_idx[3] = 0; cub_idx[4] = 2; cub_idx[5] = 3;

    count = 6;

    bufIdx.create();
    bufIdx.bind();
    bufIdx.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufIdx.allocate(cub_idx, 6 * sizeof(GLuint));

    bufPos.create();
    bufPos.bind();
    bufPos.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufPos.allocate(cub_vert_pos, 4 * sizeof(glm::vec3));

    bufNor.create();
    bufNor.bind();
    bufNor.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufNor.allocate(cub_vert_nor, 4 * sizeof(glm::vec3));

    bufCol.create();
    bufCol.bind();
    bufCol.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufCol.allocate(cub_vert_col, 4 * sizeof(glm::vec3));
}
