#ifndef ANISOTROPICBXDF
#define ANISOTROPICBXDF
#include <scene/materials/bxdfs/bxdf.h>

class AnisotropicBxDF : public BxDF
{
public:

    AnisotropicBxDF() : AnisotropicBxDF(glm::vec3(0.5f), 4.0f, 20.0f)
    {

    }
    AnisotropicBxDF(const glm::vec3 &color) : AnisotropicBxDF(color, 4.0f, 20.0f)
    {

    }
    AnisotropicBxDF(const glm::vec3 &color, float exp1,float exp2) :
        BxDF(BxDFType(BSDF_REFLECTION | BSDF_SPECULAR)),
        reflection_color(color), e1(exp1), e2(exp2)
    {
        fresnel = new FresnelNo();
    }
//Functions
    virtual glm::vec3 EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi,float& pdf) const;

    virtual glm::vec3 SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const;
    virtual float PDF(const glm::vec3 &wo, const glm::vec3 &wi) const;

    void SampleFirstQuadrant(const float u1,const float u2,float& phi,float& costheta) const;

//Member variables
    glm::vec3 reflection_color;
    float e1,e2;
};

#endif // ANISOTROPICBXDF

