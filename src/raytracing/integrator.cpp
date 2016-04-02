#include <raytracing/integrator.h>
#include <QStack>

Integrator::Integrator():
    max_depth(5),
    Number_Light(10),
    Number_BRDF(10)
{
    //random number generator and uniform distribution
    generator = std::mt19937(std::chrono::system_clock::now().time_since_epoch().count());
    uniform_distribution = std::uniform_real_distribution<float>(0.0f,1.0f);

    scene = NULL;
    intersection_engine = NULL;

    //std::cout<<"Integrator constructor here!\n";
}

glm::vec3 ComponentMult(const glm::vec3 &a, const glm::vec3 &b)
{
    return glm::vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}

void Integrator::SetDepth(unsigned int depth)
{
    max_depth = depth;
}

float Integrator::PowerHeuristic(const float &pdf_s, const float &n_s, const float &pdf_f, const float &n_f)
{
    if(isinf(pdf_s))
        return 1.0f;
    if(isinf(pdf_f))
        return 0.0f;
    else if(fequal(pdf_f,0.0f) && fequal(pdf_s,0.0f))
        return 0.f;
    else
        return glm::pow(pdf_s * n_s,2.0f) / (glm::pow(pdf_s * n_s,2.0f) + glm::pow(pdf_f * n_f,2.0f));
}

//Russian roulette
bool Integrator::RussianRoulette(const glm::vec3 &color,const int& depth)
{
    if(depth <= 2)
        return false;

    float maxElemet = glm::max(glm::max(color.x,color.y),color.z);

    throughput *= maxElemet;

    if(uniform_distribution(generator) > throughput)
        return true;

    return false;

}

glm::vec3 Integrator::EstimateLight(Geometry *& light, Ray r,unsigned int depth)
{
    QStack<glm::vec3> FS;
    QStack<glm::vec3> directLight;
    FS.clear();directLight.clear();


    throughput = 1.0f;
    while(depth < max_depth)
    {

        Intersection intersection = intersection_engine->GetIntersection(r);
        if(intersection.t > 0 && intersection.object_hit->material->is_light_source)
        {   
            //F_absdot.push(glm::vec3(1));
            //directLight.push(intersection.object_hit->material->base_color * intersection.texture_color);
            break;
        }

        if(intersection.t > 0)
        {
            glm::vec3 wo_world = -r.direction;
            glm::vec3 wo_local = intersection.ToLocalNormalCoordinate(wo_world);
            glm::vec3 wj_world;
            glm::vec3 wj_local;

            glm::vec3 F;
            float absDot;
            float pdf;

            glm::vec3 Ld = EstimateDirectLight(intersection,r,light,wj_world);

            //wj_local = intersection.ToLocalNormalCoordinate(wj_world);

            F = intersection.object_hit->material->SampleAndEvaluateScatteredEnergy(intersection,
                                                                           wo_local,
                                                                           wj_local,pdf);
            wj_world = intersection.ToWorldNormalCoordinate(wj_local);

            absDot  = glm::abs(glm::dot(intersection.normal,wj_world));

            directLight.push(Ld);
            if(isinf(pdf))
                FS.push(F);
            else
                FS.push(F*absDot);

            r = Ray(intersection.point + glm::sign(glm::dot(wj_world,intersection.normal)) * 1e-3f*intersection.normal,wj_world);

            //Russian Roulette
            if(RussianRoulette(F*absDot, depth))
                break;
        }
        else
        {
            break;
        }

        depth++;
    }

    glm::vec3 indirectLight(0);
    while(!directLight.empty() && !FS.empty())
    {
        glm::vec3 F = FS.pop();
        glm::vec3 direct = directLight.pop();

        indirectLight = F * (indirectLight) + direct;
    }

    return indirectLight;
}
//Basic ray trace
glm::vec3 Integrator::TraceRay(Ray r, unsigned int depth)
{

    Intersection intersection = intersection_engine->GetIntersection(r);
    if(intersection.t > 0 && intersection.object_hit->material->is_light_source)
    {
        return intersection.texture_color * intersection.object_hit->material->base_color;
    }

    glm::vec3 color(0);
    for(Geometry* light : scene->lights)
       color = color + EstimateLight(light,r,depth);

    return color;

}

DirectLightingIntegrator::DirectLightingIntegrator()
{
}

// MIS: sampling light source
glm::vec3 Integrator::MIS_SampleLight(Intersection &intersection, Ray &r, Geometry* &light)
{
    if(Number_Light == 0)
        return glm::vec3(0);

    // Direct light estimator: sample Light source
    glm::vec3 sum_light_sample(0);
    for(int i = 0; i < Number_Light; i++)
    {
        float u = uniform_distribution(generator);
        float v = uniform_distribution(generator);

        Intersection lightSample = light->SampleOnGeometrySurface(u, v, intersection.point + 1e-3f * intersection.normal);
        glm::vec3 wj = glm::normalize(lightSample.point - intersection.point);
        glm::vec3 wo = - r.direction;
        glm::vec3 P = intersection.point;
        glm::vec3 N = intersection.normal;

        float pdf_light = light->RayPDF(intersection, Ray(P + float(1e-3)*N, wj));

        Intersection lightIntersection = intersection_engine->GetIntersection(Ray(P + float(1e-3)*N, wj));
        float temp, pdf_brdf;

        glm::vec3 wo_local = intersection.ToLocalNormalCoordinate(wo);
        glm::vec3 wj_local = intersection.ToLocalNormalCoordinate(wj);

        glm::vec3 F = intersection.object_hit->material->EvaluateScatteredEnergy(intersection, wo_local, wj_local, pdf_brdf);
        // reach light directly && pdf(wj) > 0
        if(lightIntersection.t > 0 && lightIntersection.object_hit == light && pdf_light > 0 && pdf_brdf > 0)
        {

            glm::vec3 Ld = light->material->EvaluateScatteredEnergy(lightSample, wo, -wj, temp);
            float W = PowerHeuristic(pdf_light, float(Number_Light), pdf_brdf, float(Number_BRDF)); // cause shadow in center

            sum_light_sample = sum_light_sample +
                                  W * F * Ld * float(fabs(glm::dot(wj, N))) / pdf_light;

        }
    }
    return sum_light_sample / float(Number_Light);
}

// MIS: sampling BRDF
glm::vec3 Integrator::MIS_SampleBRDF(Intersection &intersection, Ray &r, Geometry* &light)
{
    if(Number_BRDF == 0)
        return glm::vec3(0);

    // Direct light estimator: sample BRDF
    glm::vec3 sum_brdf_sample(0.0f);
    for(int i = 0; i < Number_BRDF; i++)
    {
        glm::vec3 wo_local = intersection.ToLocalNormalCoordinate(-r.direction);
        glm::vec3 wj_local;
        float pdf_brdf;
        glm::vec3 F = intersection.object_hit->material->SampleAndEvaluateScatteredEnergy(intersection,wo_local,wj_local,pdf_brdf);

        glm::vec3 wj_world = intersection.ToWorldNormalCoordinate(wj_local);
        glm::vec3 wo_world = - r.direction;

        Intersection isxOnLight = intersection_engine->GetIntersection(Ray(intersection.point+float(1e-3)*intersection.normal, wj_world));

        if(isxOnLight.t > 0 && isxOnLight.object_hit == light && pdf_brdf > 0)
        {
            float temp,pdf_light = light->RayPDF(intersection, Ray(intersection.point, wj_world));
            float W = PowerHeuristic(pdf_brdf,float(Number_BRDF),pdf_light,float(Number_Light));
            glm::vec3 Ld = light->material->EvaluateScatteredEnergy(isxOnLight,wo_world,-wj_world,temp);

            if(pdf_light > 0 )
            {
                if(isinf(pdf_brdf)) // delta specular surface
                {
                    sum_brdf_sample = sum_brdf_sample +
                            F * Ld * float(fabs(glm::dot(wj_world, intersection.normal))) / pdf_light;
                }
                else
                {
                    sum_brdf_sample = sum_brdf_sample +
                            W * F * Ld * float(fabs(glm::dot(wj_world,intersection.normal))) / pdf_brdf;
                }
            }
        }

    }
    return sum_brdf_sample / float(Number_BRDF);
}

glm::vec3 Integrator::MIS_SampleLight_Ld(Intersection &intersection, Ray &r, Geometry* &light)
{
    if(Number_Light == 0)
        return glm::vec3(0);

    // Direct light estimator: sample Light source
    glm::vec3 sum_light_sample(0);
    for(int i = 0; i < Number_Light; i++)
    {
        float u = uniform_distribution(generator);
        float v = uniform_distribution(generator);

        Intersection lightSample = light->SampleOnGeometrySurface(u, v, intersection.point + 1e-3f * intersection.normal);
        glm::vec3 wj = glm::normalize(lightSample.point - intersection.point);
        glm::vec3 wo = - r.direction;
        glm::vec3 P = intersection.point;
        glm::vec3 N = intersection.normal;

        float pdf_light = light->RayPDF(intersection, Ray(P + float(1e-3)*N, wj));

        Intersection lightIntersection = intersection_engine->GetIntersection(Ray(P + float(1e-3)*N, wj));
        float temp, pdf_brdf;
        glm::vec3 wo_local = intersection.ToLocalNormalCoordinate(wo);
        glm::vec3 wi_local = intersection.ToLocalNormalCoordinate(wj);

        glm::vec3 F = intersection.object_hit->material->EvaluateScatteredEnergy(intersection, wo_local, wi_local, pdf_brdf);
        // reach light directly && pdf(wj) > 0
        if(lightIntersection.t > 0 && lightIntersection.object_hit == light && pdf_light > 0 && pdf_brdf > 0)
        {

            glm::vec3 Ld = light->material->EvaluateScatteredEnergy(lightSample, wo, -wj, temp);
            float W = PowerHeuristic(pdf_light, float(Number_Light), pdf_brdf, float(Number_BRDF)); // cause shadow in center

            sum_light_sample = sum_light_sample +
                                  F * W * Ld * glm::abs(glm::dot(N,wj))/ pdf_light ;

        }
    }
    return sum_light_sample / float(Number_Light);
}

// MIS: sampling brdf and assigned the sampled direction wj
glm::vec3 Integrator::MIS_SampleBRDF_Ld(Intersection &intersection, Ray &r, Geometry* &light, glm::vec3& wj)
{
    if(Number_BRDF == 0)
        return glm::vec3(0);

    // Direct light estimator: sample BRDF
    glm::vec3 sum_brdf_sample(0.0f);

    glm::vec3 wo_local = intersection.ToLocalNormalCoordinate(-r.direction);
    glm::vec3 wj_local;
    float pdf_brdf;
    glm::vec3 F = intersection.object_hit->material->SampleAndEvaluateScatteredEnergy(intersection,wo_local,wj_local,pdf_brdf);

    glm::vec3 wj_world = intersection.ToWorldNormalCoordinate(wj_local);
    glm::vec3 wo_world = - r.direction;

    Intersection isxOnLight = intersection_engine->GetIntersection(Ray(intersection.point+float(1e-3)*intersection.normal, wj_world));
    if(isxOnLight.t > 0 && isxOnLight.object_hit == light && pdf_brdf > 0)
    {
        float temp,pdf_light = light->RayPDF(intersection, Ray(intersection.point, wj_world));
        float W = PowerHeuristic(pdf_brdf,float(Number_BRDF),pdf_light,float(Number_Light));
        glm::vec3 Ld = light->material->EvaluateScatteredEnergy(isxOnLight,wo_world,-wj_world,temp);

        if(pdf_light > 0)
        {
            if(isinf(pdf_brdf) ) // delta specular surface
            {
                sum_brdf_sample = sum_brdf_sample +
                        F * Ld * float(fabs(glm::dot(wj_world, intersection.normal)))/ pdf_light;
            }
            else
            {
                sum_brdf_sample = sum_brdf_sample +
                        W * F* Ld * float(fabs(glm::dot(wj_world, intersection.normal))) / pdf_brdf;
            }
        }
    }

    wj = wj_world;
    return sum_brdf_sample / float(Number_BRDF);
}

glm::vec3 Integrator::EstimateDirectLight(Intersection &isx, Ray &ray, Geometry* &light, glm::vec3 &wj)
{
    glm::vec3 Ld(0);

    Ld = MIS_SampleBRDF_Ld(isx,ray,light,wj) + MIS_SampleLight(isx,ray,light);

//    Ld = MIS_SampleBRDF_Ld(isx,ray,light,wj);

//    glm::vec3 wo_local = isx.ToLocalNormalCoordinate(-ray.direction);
//    glm::vec3 wj_local = isx.ToLocalNormalCoordinate(wj);
//    float pdf;
//    glm::vec3 F_brdf = isx.object_hit->material->EvaluateScatteredEnergy(isx,wo_local,wj_local,pdf);

//    float absDot = glm::abs(glm::dot(isx.normal,wj));

//    Ld = Ld + MIS_SampleLight(isx,ray,light) / F_brdf / absDot;

    return Ld;
}

//direct lighting integrator
glm::vec3 DirectLightingIntegrator::TraceRay(Ray r, unsigned int depth)
{
    if(depth >= max_depth)
        return glm::vec3(0);

    Intersection intersection = intersection_engine->GetIntersection(r);
    if(intersection.t > 0)
    {
        if(intersection.object_hit->material->is_light_source)
            return intersection.object_hit->material->base_color * intersection.texture_color;

        glm::vec3 resultColor(0);

        // light sampling
        glm::vec3 sample_light(0);
        glm::vec3 sample_brdf(0);

        for(Geometry* light : scene->lights)
        {
            glm::vec3 sample_light_i = MIS_SampleLight(intersection, r, light);
            sample_light = sample_light + sample_light_i;

            glm::vec3 sample_brdf_i = MIS_SampleBRDF(intersection, r, light);
            sample_brdf = sample_brdf + sample_brdf_i;
        }

        // combined result
        resultColor =  sample_brdf + sample_light;
        return resultColor;
    }
    else
        return glm::vec3(0);

}
