#include <scene/materials/bxdfs/phongbxdf.h>

glm::vec3 PhongBxDF::EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi, float& pdf) const
{
    pdf = PDF(wo, wi);

    if(wo.z < 0 || wi.z < 0 )
        return glm::vec3(0);

    glm::vec3 N(0,0,1);
    glm::vec3 wo_reflect = 2.0f*glm::dot(N,wo) * N - wo;

    float p = glm::abs(glm::dot(wo_reflect,wi));

    return (diffuse_color / PI + glm::pow(p,specular_power) * specular_color);
}

glm::vec3 PhongBxDF::EvaluateHemisphereScatteredEnergy(const glm::vec3 &wo, int num_samples, const glm::vec2 *samples) const
{
    return glm::vec3(0);
}

glm::vec3 PhongBxDF::SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const
{
    glm::vec3 N(0,0,1);
    wi_ret = 2.0f * glm::dot(N,wo)*N - wo;

    float costheta = glm::pow(rand1, 0.5f);
    float sintheta = glm::sqrt(glm::max(0.0f, 1.0f - costheta*costheta));
    float phi = rand2 * 2* PI;

    glm::vec3 delta = SphericalDirection(sintheta,costheta,phi);

    wi_ret = wi_ret + 0.2f * delta;
    wi_ret = glm::normalize(wi_ret);

    return EvaluateScatteredEnergy(wo,wi_ret,pdf_ret);
}

float PhongBxDF::PDF(const glm::vec3 &wo, const glm::vec3 &wi) const
{
    if(isPerfectReflective(wo,wi))
        return 1.0f;
    else
        return 0.2f/PI;
}
