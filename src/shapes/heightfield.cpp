
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

    float alpha = params.FindOneFloat("alpha", 0.5);

    float rL = 20;            // heightfield edge length
    float cl = 0.1;           // auto correlation length
    float h = alpha * cl;     // height

    std::unique_ptr<Point3f[]> P(new Point3f[nx * ny]);
    randomsurface(nx,rL,h,cl,P);

    CHECK(nx != -1 && ny != -1);

    int ntris = 2 * (nx - 1) * (ny - 1);
    std::unique_ptr<int[]> indices(new int[3 * ntris]);
    std::unique_ptr<Point2f[]> uvs(new Point2f[nx * ny]);
    //     for vertex normal
    std::unique_ptr<Normal3f[]> N(new Normal3f[nx*ny]);
    int nverts = nx * ny;

    int pos = 0;
    for (int y = 0; y < ny; ++y) {
      for (int x = 0; x < nx; ++x) {
        // move set x y z in randomsurface function
        //P[pos].x =  (float)x / (float)(nx - 1)*2*rL-rL;
        //P[pos].y =  (float)y / (float)(ny - 1)*2*rL-rL;
        //P[pos].z = z[pos];

        uvs[pos].x = (float)x / (float) (nx - 1);
        uvs[pos].y = (float)y / (float) (ny - 1);
        ++pos;
      }
    }

    // Fill in heightfield vertex offset array
    int *vp = indices.get();
    // for vertex normal computation
    int* count = new int[nitems];
    //std::unique_ptr<int[]> count(new int[nitems]);
    Point3f p0, p1, p2;
    Vector3f dp02, dp12;
    Normal3f normal1,normal2;

    float slopesum = 0;

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

        float costheta = Dot(normal1,Normal3f(0.f,0.f,1.f));
        if (costheta<0){
        std::cout<<"downside normal!"<<std::endl;
        }
        float sintheta = sqrt(1 - pow(costheta,2.f));
        slopesum += sintheta/costheta;

        p0 = P[VERT(x, y)];
        p1 = P[VERT(x + 1, y+1)];
        p2 = P[VERT(x, y + 1)];
        dp02 = p0 - p2;
        dp12 = p1 - p2;
        normal2 = Normal3f(Normalize(Cross(dp02, dp12)));

        costheta = Dot(normal2,Normal3f(0.f,0.f,1.f));
        if (costheta<0){
        std::cout<<"downside normal!"<<std::endl;
        }
        sintheta = sqrt(1 - pow(costheta,2.f));
        slopesum += sintheta/costheta;

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

    slopesum /= ntris;
    std::cout<<"average slope: "<<slopesum<<std::endl;

    // average vertex normal
    pos = 0;
    for (int y = 0; y<ny; ++y){
      for (int x = 0; x<nx; ++x){
        N[pos] /= count[pos];
        ++pos;
      }
    }

    delete[] count;

    // return CreateTriangleMesh(ObjectToWorld, WorldToObject, reverseOrientation,
    //                           ntris, indices.get(), nverts, P.get(), nullptr,
    //                           nullptr, uvs.get(), nullptr, nullptr);
    std::cout<<"finish heightfield -> triangle mesh"<<std::endl;
    return CreateTriangleMesh(ObjectToWorld, WorldToObject, reverseOrientation,
                              ntris, indices.get(), nverts, P.get(), nullptr,
                              N.get(), uvs.get(), nullptr, nullptr);
  }

  void randomsurface(int N, float rL, float h, float cl, std::unique_ptr<Point3f[]>& P){
    // isotropic surface
    // calculate x, y coordinates
    int k =0;
    for (int i = 0; i<N; ++i){
      for (int j = 0; j<N; ++j){
      P[k].x =  -rL/2.0f + float (j * rL)/float(N-1);
      P[k].y =  -rL/2.0f + float (i * rL)/float(N-1); // y is same as x
      k++;
      }
    }

    /////////////////////////////////////////////////////////////////
    // FFT part
    fftw_complex *in, *out1, *out2, *out3;
    fftw_plan p1,p2,p3;

    // calculate the transform fft2(F)
    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N);
    out1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N);
    p1 = fftw_plan_dft_2d(N, N, in, out1, FFTW_FORWARD, FFTW_MEASURE);

    // set real part of input matrix F ( Gaussian filter )
    k=0;
    while(k<N*N){
        in[k][0] = exp( - (P[k].x * P[k].x + P[k].y * P[k].y) / (cl * cl / 2));
        in[k][1] = 0;
        k++;
    }

      // calculate the transform
      fftw_execute(p1);

      // calculate the transform fft2(Z)
      out2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N);
      p2 = fftw_plan_dft_2d(N, N, in, out2, FFTW_FORWARD, FFTW_MEASURE);

      // set real part of input matrix Z
      // Z is uncorrelated Gaussian random rough surface distribution
      //  with mean 0 and standard deviation h
      srand(time(NULL)); // set random seed
      k=0;
      while (k < N*N){
        in[k][0] =  h * ((float) rand() / (RAND_MAX));
        in[k++][1] = 0;
      }

      // calculate the transform
      fftw_execute(p2);

      // calculate (ifft2(fft2(Z).*fft2(F))
      out3 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * N * N);
      p3 = fftw_plan_dft_2d(N, N, in, out3, FFTW_BACKWARD, FFTW_MEASURE);

      k=0;
      while (k < N*N){
        in[k][0] = out1[k][0] * out2[k][0] - out1[k][1] * out2[k][1];
        in[k][1] = out1[k][0] * out2[k][1] + out1[k][1] * out2[k][0];
        k++;
      }

      // calculate the transform
      fftw_execute(p3);

      // normalizing prefactors
      float factor = 1.0/sqrt(M_PI) * rL/(float)N /cl;
      k=0;
      while (k < N*N){
        P[k].z = factor * out3[k][0] / ( N * N );
        k++;
      }

      // deallocation
      fftw_destroy_plan(p1);
      fftw_destroy_plan(p2);
      fftw_destroy_plan(p3);
      fftw_free(in);
      fftw_free(out1); fftw_free(out2); fftw_free(out3);
  }

}  // namespace pbrt
