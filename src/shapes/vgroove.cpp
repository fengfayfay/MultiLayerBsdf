
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// shapes/heightfield.cpp*

#include <stdlib.h>
#include "microfacet.h"
#include "shapes/vgroove.h"
#include "shapes/triangle.h"
#include "paramset.h"

namespace pbrt {


void 
makeGroove(std::unique_ptr<Point3f[]>& P, std::unique_ptr<Point2f[]>& uvs, 
           std::unique_ptr<int[]>& indices, int grooveIndex, float grooveAngle, 
           int grooveRes) {


    float grooveSize = 1.0;
    float res = grooveRes;
    float xinterval = grooveSize / res;
    float ymin = -grooveSize * 0.5;
    float yrange = grooveSize;

    //each groove has 4 triangles
    int pos = 6 * grooveIndex;
    float xmin = xinterval * grooveIndex -grooveSize * 0.5;
    P[pos].x = xmin; 
    P[pos].y = ymin;
    P[pos].z = 0;
    P[pos+1].x = P[pos].x; 
    P[pos+1].y = P[pos].y + yrange;
    P[pos+1].z = 0;
    P[pos+2].x = P[pos].x + xinterval; 
    P[pos+2].y = P[pos+1].y;
    P[pos+2].z = 0;
    P[pos+3].x = P[pos].x + xinterval; 
    P[pos+3].y = P[pos].y;
    P[pos+3].z = 0;
    P[pos+4].x = (P[pos].x + P[pos+3].x ) *0.5;
    P[pos+4].y = P[pos].y;
    P[pos+4].z = -xinterval * 0.5 / tan(Radians(grooveAngle));
    P[pos+5].x = P[pos+4].x;
    P[pos+5].y = P[pos+1].y;
    P[pos+5].z = P[pos+4].z;

    for (int i = 0; i < 6; i++) {
        uvs[pos+i].x = P[pos+i].x;
        uvs[pos+i].y = P[pos+i].y;
    }
   
    int *vp = indices.get() + 3 * 4 * grooveIndex;
    *vp++ = pos;
    *vp++ = pos+1;
    *vp++ = pos+5;
    *vp++ = pos+5;
    *vp++ = pos+4;
    *vp++ = pos;
    
    *vp++ = pos+4;
    *vp++ = pos+5;
    *vp++ = pos+2;
    *vp++ = pos+2;
    *vp++ = pos+3;
    *vp++ = pos+4;
}


// Heightfield Definitions
std::vector<std::shared_ptr<Shape>> CreateVGroove(
    const Transform *ObjectToWorld, const Transform *WorldToObject,
    bool reverseOrientation, const ParamSet &params) {
    int nx = params.FindOneInt("nu", 100);
    int ny = params.FindOneInt("nv", 2);

    float roughu = params.FindOneFloat("roughu", 0);
    std::cout << "roughu " << roughu << "\n";
    float roughv = params.FindOneFloat("roughv", roughu);

    bool useBeckmann = params.FindOneBool("useBeckmann", true);
    float grooveAngle = params.FindOneFloat("grooveAngle", 45.0);

    MicrofacetDistribution *dist;
    if (useBeckmann) {
        //dist = new BeckmannDistribution(roughu, roughv, false, true); 
        dist = new BeckmannDistribution(roughu, roughv, false); 
    } else {
        //dist = new TrowbridgeReitzDistribution(roughu, roughv, false, true); 
        dist = new TrowbridgeReitzDistribution(roughu, roughv, false); 
    }

    srand(time(NULL)); 

    int ntris = 4 * (nx);
    int nverts = 6 * (nx);
    std::unique_ptr<int[]> indices(new int[3 * ntris]);
    std::unique_ptr<Point3f[]> P(new Point3f[nverts]);
    std::unique_ptr<Point2f[]> uvs(new Point2f[nverts]);
    // Compute heightfield vertex positions
    int pos = 0;
    Vector3f wo(0, 0, 1);
    for (int x = 0; x < nx; ++x) {
        if (roughu > 0) {
            float u1 = ((float) rand()/(RAND_MAX));
            float u2 = ((float) rand()/(RAND_MAX));
            Point2f u (u1, u2);
            Vector3f wh;
            wh = dist->Sample_wh(wo, u);
            float sampleAngle = 90.0 - Degrees(acos(wh.z));
            makeGroove(P, uvs, indices, x, sampleAngle, nx);
        } else {
            makeGroove(P, uvs, indices, x, grooveAngle, nx);
        }
    }

    delete dist;

    return CreateTriangleMesh(ObjectToWorld, WorldToObject, reverseOrientation,
                              ntris, indices.get(), nverts, P.get(), nullptr,
                              nullptr, uvs.get(), nullptr, nullptr);
}

}  // namespace pbrt
