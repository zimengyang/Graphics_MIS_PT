#ifndef TRANSMISSIONBXDF
#define TRANSMISSIONBXDF
#include <scene/materials/bxdfs/bxdf.h>

class TransmissionBxDF : public BxDF
{
public:
    float etai,etat;
    TransmissionBxDF():
        TransmissionBxDF(1.0f, 1.0f, glm::vec3(1.0f))
    {
    }

    TransmissionBxDF(const float ei,const float et,const glm::vec3 &color) :
        BxDF(BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR)),
        trans_color(color),
        etai(ei),etat(et)
    {
        fresnel = new FresnelDielectric(etai,etat);
    }
//Functions
    virtual glm::vec3 EvaluateScatteredEnergy(const glm::vec3 &wo, const glm::vec3 &wi,float& pdf) const;

    virtual glm::vec3 SampleAndEvaluateScatteredEnergy(const glm::vec3 &wo, glm::vec3 &wi_ret, float rand1, float rand2, float &pdf_ret) const;
    virtual float PDF(const glm::vec3 &wo, const glm::vec3 &wi) const;

    glm::vec3 trans_color;
};

#endif // TRANSMISSIONBXDF

