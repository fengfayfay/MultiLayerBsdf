#ifndef _MLBRDF_H
#define _MLBRDF_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tensorflow/c/c_api.h>
#include "geometry.h"

namespace pbrt{
class RealNVPScatter {
public:
    RealNVPScatter(const std::string& modelPathPrefix=".");
    ~RealNVPScatter();
    bool loadAndRestore();
    bool setupSampleTensors();
    bool setupEvalTensors();
    
    float eval(float thetaI, float alpha, const Vector2f& sampleN);
    pbrt::Vector2f  sample(float thetaI, float alpha);      

private:
    std::string modelPath;
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
}
#endif
