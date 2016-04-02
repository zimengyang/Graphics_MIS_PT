#include <scene/materials/bxdfs/bxdf.h>

glm::vec3 BxDF::SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const
{
    //TODO
    wi_ret = glm::vec3(0);
    pdf_ret = 0.0f;
    return glm::vec3(0);
}

glm::vec3 BxDF::EvaluateHemisphereScatteredEnergy(const glm::vec3 &wo, int num_samples, const glm::vec2* samples) const
{
    //TODO
    return glm::vec3(0);
}

glm::vec3 BxDF::EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi, float &pdf) const
{
    return glm::vec3(0);
}

float BxDF::PDF(const glm::vec3 &wo, const glm::vec3 &wi) const
{
    //TODO
    return 0.0f;
}

bool BxDF::isPerfectReflective(const glm::vec3 &wo, const glm::vec3 &wi) const
{
    glm::vec3 N(0,0,1);
    glm::vec3 ref = 2.0f * glm::dot(wo,N) * N - wo;
    ref = glm::normalize(ref);

    float len = glm::length(ref-wi);
    if(len < 0.01f)
        return true;
    else
        return false;
}


glm::vec3 BxDF::SphericalDirection(float sintheta, float costheta, float phi) const
{
    glm::vec3 direction(sintheta * glm::cos(phi), sintheta * glm::sin(phi), costheta);
    return glm::normalize(direction);
}

float FresnelDielectric::Evaluate(const float& costheta)
{
    float et=etat,ei=etai;
    float cosi = glm::clamp(costheta,-1.0f,1.0f);
    if( cosi < 0)
    {
        float temp = et;
        et = ei;
        ei =temp;
    }

    float sint = ei/et*glm::sqrt(glm::max(0.0f,1.0f-cosi*cosi));
    if(sint >= 1.0f)
        return 1.0f;
    else
    {
        float cost = glm::sqrt(glm::max(0.0f,1.0f - sint*sint));

        return FreDiel(ei,glm::abs(cosi),et,glm::abs(cost));
    }
}

float FresnelDielectric::FreDiel(const float&ei, const float &cosi, const float& et, const float &cost)
{
    float Rpar = (et*cosi - ei * cost) / ( et*cosi + ei * cost);
    float Rper = (ei*cosi - et * cost) / ( ei*cosi + et * cost);
    return (Rpar * Rpar + Rper * Rper)/2.0f;
}
