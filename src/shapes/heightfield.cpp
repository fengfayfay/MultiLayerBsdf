
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
#include "shapes/heightfield.h"
#include "shapes/triangle.h"
#include "paramset.h"

#include <fstream>

namespace pbrt {

// Heightfield Definitions
std::vector<std::shared_ptr<Shape>> CreateHeightfield(
    const Transform *ObjectToWorld, const Transform *WorldToObject,
    bool reverseOrientation, const ParamSet &params) {

    // read Pz values from file 
    const std::string filename = params.FindOneFilename("filename", "");
    int nx = params.FindOneInt("nu", -1);
    int ny = params.FindOneInt("nv", -1);
    int nitems = nx * ny;
    float* z = new float[nitems];
    std::fstream myfile(filename, std::ios_base::in);
    float a;
    int co = 0;
    while (myfile >> a){
      //printf("%f ", a);
      z[co] = a;
      co++;
    }  
    
    
    // int nx = params.FindOneInt("nu", -1);
    // int ny = params.FindOneInt("nv", -1);
    // int nitems;
    // const Float *z = params.FindFloat("Pz", &nitems);
    // CHECK_EQ(nitems, nx * ny);
    CHECK(nx != -1 && ny != -1 && z != nullptr);

    int ntris = 2 * (nx - 1) * (ny - 1);
    std::unique_ptr<int[]> indices(new int[3 * ntris]);
    std::unique_ptr<Point3f[]> P(new Point3f[nx * ny]);
    std::unique_ptr<Point2f[]> uvs(new Point2f[nx * ny]);
//     for vertex normal
    std::unique_ptr<Normal3f[]> N(new Normal3f[nx*ny]);
    int nverts = nx * ny;
    // Compute heightfield vertex positions
    int pos = 0;
    float domain = 80;
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            P[pos].x = uvs[pos].x = (float)x / (float)(nx - 1)*domain-domain/2;
            P[pos].y = uvs[pos].y = (float)y / (float)(ny - 1)*domain-domain/2;
            P[pos].z = z[pos];
            ++pos;
        }
    }

//     Fill in heightfield vertex offset array
    int *vp = indices.get();
    int* count = new int[nitems];
    Point3f p0, p1, p2;
    Vector3f dp02, dp12;
    Normal3f normal1,normal2;
    for (int y = 0; y < ny - 1; ++y) {
        for (int x = 0; x < nx - 1; ++x) {
#define VERT(x, y) ((x) + (y)*nx)
            *vp++ = VERT(x, y);
            *vp++ = VERT(x + 1, y);
            *vp++ = VERT(x + 1, y + 1);

            *vp++ = VERT(x, y);
            *vp++ = VERT(x + 1, y + 1);
            *vp++ = VERT(x, y + 1);

	    // calculate vertex normals for normal interpolation
	    // surface normal
	    p0 = P[VERT(x, y)];
	    p1 = P[VERT(x + 1, y)];
	    p2 = P[VERT(x + 1, y + 1)];
	    dp02 = p0 - p2;
	    dp12 = p1 - p2;
	    normal1 = Normal3f(Normalize(Cross(dp02, dp12)));

	    p0 = P[VERT(x, y)];
	    p1 = P[VERT(x + 1, y+1)];
	    p2 = P[VERT(x, y + 1)];
	    dp02 = p0 - p2;
	    dp12 = p1 - p2;
	    normal2 = Normal3f(Normalize(Cross(dp02, dp12)));
	    N[VERT(x, y)] += normal1 + normal2; 
	    N[VERT(x + 1, y)] += normal1;
	    N[VERT(x + 1, y + 1)] += normal1 + normal2;
	    N[VERT(x, y + 1)] += normal2;
	    // increment # of contribution
	    count[VERT(x, y)] += 2;
	    count[VERT(x+1, y)]++;
	    count[VERT(x+1, y+1)] += 2;
	    count[VERT(x, y+1)]++;
	    
        }
#undef VERT
    }

//     // average vertex normal
    pos = 0;
    for (int y = 0; y<ny; ++y){
      for (int x = 0; x<nx; ++x){
    	N[pos] /= count[pos];
    	++pos;
      }	
    }

    delete[] z;
    delete[] count;
    
    // return CreateTriangleMesh(ObjectToWorld, WorldToObject, reverseOrientation,
    //                           ntris, indices.get(), nverts, P.get(), nullptr,
    //                           nullptr, uvs.get(), nullptr, nullptr);
    std::cout<<"finish heightfield -> triangle"<<std::endl;
    return CreateTriangleMesh(ObjectToWorld, WorldToObject, reverseOrientation,
                              ntris, indices.get(), nverts, P.get(), nullptr,
                              N.get(), uvs.get(), nullptr, nullptr);
}

}  // namespace pbrt
