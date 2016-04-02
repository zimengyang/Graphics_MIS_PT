#include <scene/materials/bxdfs/specularreflectionBxDF.h>

glm::vec3 SpecularReflectionBxDF::EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi, float& pdf) const
{
    pdf = PDF(wo, wi);
    if(isPerfectReflective(wo,wi))
        return reflection_color;
    else
        return glm::vec3(0);
}
glm::vec3 SpecularReflectionBxDF::EvaluateHemisphereScatteredEnergy(const glm::vec3 &wo, int num_samples, const glm::vec2 *samples) const
{
    //TODO
    return glm::vec3(0);
}

glm::vec3 SpecularReflectionBxDF::SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const
{
    // specular reflection wi = delta(w-wo) where w is perfectly reflective direction
    glm::vec3 N(0,0,1);
    wi_ret = 2.0f*glm::dot(N,wo)*N - wo;
    wi_ret = glm::normalize(wi_ret);
    //pdf_ret = PDF(wo, wi_ret);

    return EvaluateScatteredEnergy(wo, wi_ret, pdf_ret);
}

float SpecularReflectionBxDF::PDF(const glm::vec3 &wo, const glm::vec3 &wi) const
{
    if(isPerfectReflective(wo,wi))
        return INFINITY;
    else
        return 0.0f;
}
