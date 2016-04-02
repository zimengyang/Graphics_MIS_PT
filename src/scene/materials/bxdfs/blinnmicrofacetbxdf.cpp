#include <scene/materials/bxdfs/blinnmicrofacetbxdf.h>

glm::vec3 BlinnMicrofacetBxDF::EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi, float &pdf) const
{
    pdf = PDF(wo,wi);

    if(wo.z < 0 || wo.z*wi.z<0)
        return glm::vec3(0);

    glm::vec3 N(0,0,1);
    glm::vec3 wh = glm::normalize(wi+wo);

    float theta_i = glm::abs(glm::dot(N,wi));
    float theta_o = glm::abs(glm::dot(N,wo));
    if(theta_i == 0.0f || theta_o == 0.0f)
        return glm::vec3(0.0f);

    // for D(wh) term
    float costhetaH = glm::abs(wh.z);
    float D = (exponent + 2.0f) / (2.0f*PI) * glm::pow(costhetaH, exponent);

    // for fresnel term
    float F = fresnel->Evaluate(glm::abs(wh.z));//1.0f;

    // for geometric term
    float M = glm::abs(2.0f * glm::dot(N,wh)*glm::dot(N,wo) / glm::dot(wo,wh));
    float S = glm::abs(2.0f * glm::dot(N,wh)*glm::dot(N,wi) / glm::dot(wo,wh));
    float G = glm::min(1.0f,glm::min(M,S));

    float BrDF = (D * F * G) / (4.0f * theta_i * theta_o);

    return reflection_color * BrDF;
}
glm::vec3 BlinnMicrofacetBxDF::EvaluateHemisphereScatteredEnergy(const glm::vec3 &wo, int num_samples, const glm::vec2 *samples) const
{
    //TODO
    return glm::vec3(0);
}

glm::vec3 BlinnMicrofacetBxDF::SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const
{
    // sample half-angle vector wi_ret
    float costheta = glm::pow(rand1, 1.0f / (exponent + 1));
    float sintheta = glm::sqrt(glm::max(0.0f, 1.0f - costheta*costheta));
    float phi = rand2 * 2* PI;

    glm::vec3 wh = glm::normalize(SphericalDirection(sintheta,costheta,phi));

    if(wo.z * wh.z < 0)
        wh = -wh;

    wi_ret = 2.0f * glm::dot(wo,wh) * wh - wo;

    return EvaluateScatteredEnergy(wo, wi_ret, pdf_ret);
}

float BlinnMicrofacetBxDF::PDF(const glm::vec3 &wo, const glm::vec3 &wi) const
{
    glm::vec3 wh = glm::normalize(wo+wi);

    if(glm::dot(wh,wo) <= 0.0f || wo.z < 0)
        return 0.0f;

    float pdf;
    float costheta = glm::abs(wh.z);

    pdf = (exponent + 1.0f) / (2.0f*PI*4.0f*glm::dot(wo,wh)) * glm::pow(costheta, exponent);

    return pdf;
}
