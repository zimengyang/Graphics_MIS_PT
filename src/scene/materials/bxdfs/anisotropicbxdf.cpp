#include <scene/materials/bxdfs/anisotropicbxdf.h>

// for anisotropic bxdf
glm::vec3 AnisotropicBxDF::EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi, float &pdf) const
{
    glm::vec3 N(0,0,1);
    glm::vec3 wh = glm::normalize(wi+wo);

    float theta_i = glm::abs(glm::dot(N,wi));
    float theta_o = glm::abs(glm::dot(N,wo));
    if(theta_i == 0.0f || theta_o == 0.0f)
        return glm::vec3(0.0f);

    // for D(wh) term
    float costhetaH = glm::abs(wh.z);
    float d=1.f-costhetaH*costhetaH;
    if(d==0) return glm::vec3(0);
    float e = e1 * wh.x * wh.x + e2 * wh.y * wh.y;
    e = e / d;
    float D = glm::sqrt((e1 + 2.0f) * (e2 + 2.0f) )/ 2.0f / PI * glm::pow(costhetaH,e);

    // for fresnel term
    float F = fresnel->Evaluate(glm::abs(wh.z));//1.0f;

    // for geometric term
    float M = glm::abs(2.0f * glm::dot(N,wh)*glm::dot(N,wo) / glm::dot(wo,wh));
    float S = glm::abs(2.0f * glm::dot(N,wh)*glm::dot(N,wi) / glm::dot(wi,wh));
    float G = glm::min(1.0f,glm::min(M,S));

    float BrDF = (D * F * G) / (4.0f * theta_i * theta_o);

    return reflection_color * BrDF;
}
glm::vec3 AnisotropicBxDF::SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const
{
    float u1 = rand1;
    float u2 = rand2;
    float phi,costheta;
    if(u1<0.25)
    {
        SampleFirstQuadrant(4.0f*u1,u2,phi,costheta);
    }
    else if(u1<0.5f)
    {
        SampleFirstQuadrant(4.0f*(0.5f-u1),u2,phi,costheta);
        phi = PI - phi;
    }
    else if(u1<0.75f)
    {
        SampleFirstQuadrant(4.0f*(u1-0.5f),u2,phi,costheta);
        phi += PI;
    }
    else
    {
        SampleFirstQuadrant(4.0f*(1.0f-u1),u2,phi,costheta);
        phi = 2.0f * PI - phi;
    }
    float sintheta = glm::sqrt(glm::max(0.f, 1.f-costheta*costheta));
    glm::vec3 wh = SphericalDirection(sintheta,costheta,phi);
    if(wh.z * wo.z <0)
        wh = -wh;

    wi_ret = 2.0f * glm::dot(wo,wh) * wh - wo;

    pdf_ret = PDF(wo,wi_ret);

    return EvaluateScatteredEnergy(wo,wi_ret,pdf_ret);
}

float AnisotropicBxDF::PDF(const glm::vec3 &wo, const glm::vec3 &wi) const
{
    glm::vec3 wh = glm::normalize(wo+wi);

    float costheta = glm::abs(wh.z);
    float ds = 1- costheta*costheta;

    float pdf = 0.0f;
    if(ds > 0.0f && glm::dot(wo,wh) > 0)
    {
        float e = e1 * wh.x * wh.x + e2 * wh.y * wh.y;
        e = e / ds;
        float d = glm::sqrt((e1 + 1.0f) * (e2 + 1.0f) / 2.0f / PI * glm::pow(costheta,e));
        pdf = d / (4.0f * glm::dot(wo,wh));
    }
    return pdf;
}

void AnisotropicBxDF::SampleFirstQuadrant(const float u1, const float u2, float &phi, float &costheta) const
{
    if(e1 == e2)
    {
        phi = PI * u1 * 0.5f;
    }
    else
        phi = glm::atan(glm::sqrt(e1 +1) / (e2+1))*glm::tan(PI*u1*0.5f);

    float costhi = glm::cos(phi),sinthi=glm::sin(phi);
    costheta = glm::pow(u2,1.0f/(e1*costhi*costhi + e2*sinthi*sinthi));
}
