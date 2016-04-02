#include <scene/geometry/ring.h>

Intersection Ring::SampleOnGeometrySurface(const float &u, const float &v, const glm::vec3 &point)
{
    float r =  0.5f + u * 0.5f;
    float theta = v * 2.0f * PI;

    glm::vec3 hitPoint(r * cos(theta),
                       r * sin(theta),
                       0.0f);

    Intersection result;
    result.point = glm::vec3(transform.T() * glm::vec4(hitPoint, 1.0f));
    result.object_hit = this;
    result.texture_color = Material::GetImageColorInterp(GetUVCoordinates(hitPoint), material->texture);

    //Store the tangent and bitangent
    SetNormalTangentBitangent(hitPoint,result);

    result.t = 1.0f;

    return result;
}

Intersection Ring::GetIntersection(Ray r)
{
    //Transform the ray
    Ray r_loc = r.GetTransformedCopy(transform.invT());
    Intersection result;

    //Ray-plane intersection
    float t = glm::dot(glm::vec3(0,0,1), (glm::vec3(0.5f, 0.5f, 0) - r_loc.origin)) / glm::dot(glm::vec3(0,0,1), r_loc.direction);
    glm::vec4 P = glm::vec4(t * r_loc.direction + r_loc.origin, 1);

    float dist2 = (P.x * P.x + P.y * P.y);
    if(t > 0 && dist2 <= 1.0f && dist2 >= 0.25f)
    {
        result.point = glm::vec3(transform.T() * P);
        result.object_hit = this;
        result.t = glm::distance(result.point, r.origin);
        result.texture_color = Material::GetImageColorInterp(GetUVCoordinates(glm::vec3(P)), material->texture);

        //Store the tangent and bitangent
        SetNormalTangentBitangent(glm::vec3(P),result);

        return result;
    }
    return result;
}

void Ring::SetNormalTangentBitangent(const glm::vec3 &point_local, Intersection &isx)
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
void Ring::ComputeArea()
{
    // use linear transformation matrix to calculate the area of shape
    float detT = glm::determinant(transform.T());
    float r1 = 0.5f, r2=1.0f;
    float originArea = PI*(r2*r2 - r1*r1);

    area = detT*originArea;
    if(transform.getScale().z != 1.0f)
        area = area/transform.getScale().z;
}

glm::vec2 Ring::GetUVCoordinates(const glm::vec3 &point)
{
    return glm::vec2(point.x/2.0f + 0.5f, point.y/2.0f + 0.5f);
}

glm::vec3 Ring::ComputeNormal(const glm::vec3 &P)
{
    return glm::vec3(0,0,1);
}

void Ring::create()
{
    GLuint idx[120];

    glm::vec3 vert_pos[40];
    glm::vec3 vert_nor[40];
    glm::vec3 vert_col[40];


    glm::vec4 pt1(0.5f, 0, 0, 1);
    glm::vec4 pt2(1.0f, 0, 0, 1);

    float angle = 18.0f * DEG2RAD;
    glm::vec3 axis(0,0,1);
    for(int i = 0; i < 20; i++)
    {
        //Position
        glm::vec3 new_pt1 = glm::vec3(glm::rotate(glm::mat4(1.0f), angle * i, axis) * pt1);
        vert_pos[i*2] = new_pt1;

        glm::vec3 new_pt2 = glm::vec3(glm::rotate(glm::mat4(1.0f), angle * i, axis) * pt2);
        vert_pos[i*2+1] = new_pt2;

        //Normal
        vert_nor[i*2] = glm::vec3(0,0,1);vert_nor[i*2+1] = glm::vec3(0,0,1);
        //Color
        vert_col[i*2] = material->base_color; vert_col[i*2+1] = material->base_color;
    }

    //Fill the indices.
    int index = 0;
    for(int i = 0; i <= 18; i++)
    {
        idx[index++] = 2 * i;
        idx[index++] = 2 * i + 1;
        idx[index++] = 2 * i + 2;

        idx[index++] = 2 * i + 1;
        idx[index++] = 2 * i + 2;
        idx[index++] = 2 * i + 3;
    }


    idx[index++] = 38;
    idx[index++] = 39;
    idx[index++] = 0;

    idx[index++] = 0;
    idx[index++] = 39;
    idx[index++] = 1;


    count = 120;

    bufIdx.create();
    bufIdx.bind();
    bufIdx.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufIdx.allocate(idx, 120 * sizeof(GLuint));

    bufPos.create();
    bufPos.bind();
    bufPos.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufPos.allocate(vert_pos, 40 * sizeof(glm::vec3));

    bufNor.create();
    bufNor.bind();
    bufNor.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufNor.allocate(vert_nor, 40 * sizeof(glm::vec3));

    bufCol.create();
    bufCol.bind();
    bufCol.setUsagePattern(QOpenGLBuffer::StaticDraw);
    bufCol.allocate(vert_col, 40 * sizeof(glm::vec3));
}
