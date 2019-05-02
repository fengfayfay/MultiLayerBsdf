#ifndef _MLBRDF_H
#define _MLBRDF_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tensorflow/c/c_api.h>
#include "geometry.h"

namespace pbrt{

class MLFresnel {
public:
    MLFresnel(const std::string& modelPathPrefix=".", const std::string& fresnelOpName="dense_3_1/Elu");
    ~MLFresnel();
    bool loadAndRestore();
    bool setupEvalTensors();
    
    Spectrum eval(float alpha, float muI, float muH);

private:
    std::string modelPathPrefix;
    std::string fresnelOpName;
    TF_Session* sess;
    TF_Graph* graph;
    TF_Status* status;
    TF_Buffer* graph_def;
    
    TF_Tensor *eval_input_tensors[1];
    TF_Tensor *eval_output_tensors[1];
    TF_Output eval_run_inputs[1];
    TF_Output eval_run_outputs[1];
};

class RealNVPScatter {
public:
    RealNVPScatter();
    ~RealNVPScatter();
    void init(const std::string& modelPathPrefix=".");
    bool loadAndRestore();
    bool setupSampleTensors();
    bool setupEvalTensors();
    
    float eval(float thetaI, float alpha, const Vector2f& sampleN);
    pbrt::Vector2f  sample(float thetaI, float alpha);      

private:
    std::string modelPathPrefix;
    TF_Session* sess;
    TF_Graph* graph;
    TF_Status* status;
    TF_Buffer* graph_def;

    TF_Tensor *eval_input_tensors[3];
    TF_Tensor *sample_input_tensors[3];
    TF_Tensor *eval_output_tensors[1];
    TF_Tensor *sample_output_tensors[1];
    TF_Output sample_run_inputs[3];
    TF_Output eval_run_inputs[3];
    TF_Output sample_run_outputs[1];
    TF_Output eval_run_outputs[1];
};

class RealNVPScatterSpectrum {
public:
    RealNVPScatterSpectrum(Vector3f& energyRatio, const std::string& modelPathPrefix=".", 
        const std::string& fresnelPrefix="None", int numChannels = 3);
    ~RealNVPScatterSpectrum();

    Spectrum eval(Float alpha, Float muI, Float muH, Float thetaI, const Vector2f &sampleN);
    Float eval(float thetaI, float alpha, const Vector2f &sampleN);
    Spectrum evalFresnel(Float alpha, Float muI, Float muH);
    pbrt::Vector2f sample(float thetaI, float alpha);
private:
    const std::string modelPathPrefix;
    const std::string fresnelPrefix;
    int numChannels;
    MLFresnel *mlFresnel;
    RealNVPScatter *nvpScatter;
    Vector3f energyRatio;
};

}
#endif
