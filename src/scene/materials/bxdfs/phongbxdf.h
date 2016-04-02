#ifndef PHONGBXDF
#define PHONGBXDF

#include <scene/materials/bxdfs/bxdf.h>

class PhongBxDF : public BxDF
{
public:
//Constructors/Destructors
    PhongBxDF(const glm::vec3 &diffuseColor, const glm::vec3 &specularColor, const float& specularPower) :
        BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)),
        diffuse_color(diffuseColor),
        specular_color(specularColor),
        specular_power(specularPower)
    {}
//Functions
    virtual glm::vec3 EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi, float& pdf) const;
    virtual glm::vec3 EvaluateHemisphereScatteredEnergy(const glm::vec3 &wo, int num_samples, const glm::vec2 *samples) const;
    virtual glm::vec3 SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const;
    virtual float PDF(const glm::vec3 &wo, const glm::vec3 &wi) const;

//Member variables
    glm::vec3 diffuse_color;
    glm::vec3 specular_color;
    float specular_power;
};

#endif // PHONGBXDF

