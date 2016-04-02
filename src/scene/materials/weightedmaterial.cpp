#include <scene/materials/weightedmaterial.h>
#include <random>

WeightedMaterial::WeightedMaterial() : Material(){}
WeightedMaterial::WeightedMaterial(const glm::vec3 &color) : Material(color){}

glm::vec3 WeightedMaterial::EvaluateScatteredEnergy(const Intersection &isx, const glm::vec3 &woW, const glm::vec3 &wiW, float &pdf, BxDFType flags)
{
    std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<float> uniform_distribution(0.0f,1.0f);

    float x = uniform_distribution(generator);

    float allMatchedBxDFWeights=0.0f;
    for(int i=0;i < bxdfs.size(); i++)
        if(bxdfs[i]->MatchesFlags(flags))
        {
           allMatchedBxDFWeights +=  bxdf_weights[i];
        }

    if(fequal(allMatchedBxDFWeights,0.0f)) // no bxdf mathes the flag
    {
        return glm::vec3(0);
    }
    else  // make the range of x between [0,allMatchedBxDFWeights]
    {
        x = x * allMatchedBxDFWeights;
    }

    int index = 0;
    float sum = 0.0f;
    for(int i=0; i < this->bxdfs.size(); i++)
        if(bxdfs[i]->MatchesFlags(flags))
        {
            if(sum <= x && x < (sum+this->bxdf_weights[i]))
            {
                index = i;
                break;
            }
            else
            {
                sum += this->bxdf_weights[i];
            }
        }

    return base_color * isx.texture_color * bxdfs[index]->EvaluateScatteredEnergy(woW, wiW, pdf);
}

glm::vec3 WeightedMaterial::SampleAndEvaluateScatteredEnergy(const Intersection &isx, const glm::vec3 &woW, glm::vec3 &wiW_ret, float &pdf_ret, BxDFType flags)
{
    std::mt19937 generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<float> uniform_distribution(0.0f,1.0f);

    float x = uniform_distribution(generator);

    float allMatchedBxDFWeights=0.0f;
    for(int i=0;i < bxdfs.size(); i++)
        if(bxdfs[i]->MatchesFlags(flags))
        {
           allMatchedBxDFWeights +=  bxdf_weights[i];
        }

    if(fequal(allMatchedBxDFWeights,0.0f)) // no bxdf mathes the flag
    {
        return glm::vec3(0);
    }
    else  // make the range of x between [0,allMatchedBxDFWeights]
    {
        x = x * allMatchedBxDFWeights;
    }

    int index = 0;
    float sum = 0.0f;
    for(int i=0; i < this->bxdfs.size(); i++)
        if(bxdfs[i]->MatchesFlags(flags))
        {
            if(sum <= x && x < (sum+this->bxdf_weights[i]))
            {
                index = i;
                break;
            }
            else
            {
                sum += this->bxdf_weights[i];
            }
        }

    float rand1 = uniform_distribution(generator);
    float rand2 = uniform_distribution(generator);

    return base_color * isx.texture_color * bxdfs[index]->SampleAndEvaluateScatteredEnergy(woW,wiW_ret,rand1,rand2,pdf_ret);
}
