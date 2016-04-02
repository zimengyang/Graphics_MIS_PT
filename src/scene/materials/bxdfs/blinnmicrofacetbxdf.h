#pragma once
#include <scene/materials/bxdfs/bxdf.h>

class BlinnMicrofacetBxDF : public BxDF
{
public:
//Constructors/Destructors
    BlinnMicrofacetBxDF() : BlinnMicrofacetBxDF(glm::vec3(0.5f), 1.0f)
    {

    }
    BlinnMicrofacetBxDF(const glm::vec3 &color) : BlinnMicrofacetBxDF(color, 1.0f)
    {

    }
    BlinnMicrofacetBxDF(const glm::vec3 &color, float exp) : BxDF(BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)), reflection_color(color), exponent(exp)
    {
        fresnel = new FresnelNo();
    }
//Functions
    virtual glm::vec3 EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi,float& pdf) const;
    virtual glm::vec3 EvaluateHemisphereScatteredEnergy(const glm::vec3 &wo, int num_samples, const glm::vec2 *samples) const;
    virtual glm::vec3 SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const;
    virtual float PDF(const glm::vec3 &wo, const glm::vec3 &wi) const;


//Member variables
    glm::vec3 reflection_color;
    float exponent;
};
