#include "mlbrdf.h"
#include "spectrum.h"
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <tensorflow/c/c_api.h>

namespace pbrt{

TF_Buffer* read_file(const char* file);

void free_buffer(void* data, size_t length) { free(data); }

void deallocator(void* ptr, size_t len, void* arg) { free((void*)ptr); }

RealNVPScatterSpectrum::RealNVPScatterSpectrum(Vector3f& energyRatio, 
            const std::string& modelPathPrefix, const std::string& fresnelPrefix,  int numChannels):
            energyRatio(energyRatio),
            modelPathPrefix(modelPathPrefix), mlFresnel(NULL), numChannels(numChannels) {

    if (fresnelPrefix.compare("None") != 0) mlFresnel = new MLFresnel(fresnelPrefix);
    nvpScatter = new RealNVPScatter[1];
    for (int i = 0; i < 1; i++) {
        nvpScatter[i].init(modelPathPrefix);
    }
}

RealNVPScatterSpectrum::~RealNVPScatterSpectrum() {

    if (mlFresnel) delete mlFresnel;
    delete [] nvpScatter;
}

//Spectrum 
Float
RealNVPScatterSpectrum::eval(float thetaI, float alpha, const Vector2f &sampleN){
    
    Float value[3];
    value[0] = nvpScatter[0].eval(thetaI, alpha, sampleN);
    if (isnan(value[0])) value[0] = 0;
    return value[0] * energyRatio[0];

    /*
    //add fresnel evaluation
    for (int i = 0; i < numChannels; i++) {
        value[i] = value[0];
        value[i] *= energyRatio[i];
    }
    //convert rgb value to spectrum
    
    return Spectrum::FromRGB(value);
    */
}

Spectrum
RealNVPScatterSpectrum::eval(Float alpha, Float muI, Float muH, Float thetaI, const Vector2f &sampleN){

    Float dist = eval(thetaI, alpha, sampleN);
    if (dist < 1e-6) return 0;
    Spectrum F (1.0);
    if (mlFresnel) F = mlFresnel->eval(alpha, muI, muH);
    return F * dist;  

}

pbrt::Vector2f RealNVPScatterSpectrum::sample(float thetaI, float alpha) {
    return nvpScatter[0].sample(thetaI, alpha);
}

RealNVPScatter::RealNVPScatter(){
}

void
RealNVPScatter::init(const std::string& pathPrefix) {
    modelPathPrefix = pathPrefix;
    loadAndRestore();
    setupSampleTensors();
    setupEvalTensors();
}

bool 
RealNVPScatter::loadAndRestore() {
    // load graph
    // ================================================================================
    std::string filePath = modelPathPrefix + "/graph.pb";
    graph_def = read_file(filePath.c_str());
    graph = TF_NewGraph();
    
    status = TF_NewStatus();
    TF_ImportGraphDefOptions* opts = TF_NewImportGraphDefOptions();
    TF_GraphImportGraphDef(graph, graph_def, opts, status);
    TF_DeleteImportGraphDefOptions(opts);
    if (TF_GetCode(status) != TF_OK) {
        fprintf(stderr, "ERROR: Unable to import graph %s\n", TF_Message(status));
        return 1;
    }
    fprintf(stdout, "Successfully imported graph\n");

    // create session
    // ================================================================================
    TF_SessionOptions* opt = TF_NewSessionOptions();
    sess = TF_NewSession(graph, opt, status);
    TF_DeleteSessionOptions(opt);
    if (TF_GetCode(status) != TF_OK) {
        fprintf(stderr, "ERROR: Unable to create session %s\n", TF_Message(status));
        return 1;
    }
    fprintf(stdout, "Successfully created session\n");

    // run init operation
    // ================================================================================
    const TF_Operation* init_op = TF_GraphOperationByName(graph, "init");
    const TF_Operation* const* targets_ptr = &init_op;

      TF_SessionRun(sess,
                  /* RunOptions */ NULL,
                  /* Input tensors */ NULL, NULL, 0,
                  /* Output tensors */ NULL, NULL, 0,
                  /* Target operations */ targets_ptr, 1,
                  /* RunMetadata */ NULL,
                  /* Output status */ status);
    if (TF_GetCode(status) != TF_OK) {
      fprintf(stderr, "ERROR: Unable to run init_op: %s\n", TF_Message(status));
      return 1;
    }
    fprintf(stdout, "Successfully init session\n");

    // run restore
    // ================================================================================
    TF_Operation* checkpoint_op = TF_GraphOperationByName(graph, "save/Const");
    TF_Operation* restore_op = TF_GraphOperationByName(graph, "save/restore_all");


    filePath = modelPathPrefix + "/model";
    const char* checkpoint_path_str = filePath.c_str();
    size_t checkpoint_path_str_len = strlen(checkpoint_path_str);
    size_t encoded_size = TF_StringEncodedSize(checkpoint_path_str_len);

    // The format for TF_STRING tensors is:
    //   start_offset: array[uint64]
    //   data:         byte[...]
    size_t total_size = sizeof(int64_t) + encoded_size;
    char* input_encoded = (char*)malloc(total_size);
    memset(input_encoded, 0, total_size);
    TF_StringEncode(checkpoint_path_str, checkpoint_path_str_len,
                    input_encoded + sizeof(int64_t), encoded_size, status);
    if (TF_GetCode(status) != TF_OK) {
      fprintf(stderr, "ERROR: something wrong with encoding: %s",
              TF_Message(status));
      return 1;
    }

    TF_Tensor* path_tensor = TF_NewTensor(TF_STRING, NULL, 0, input_encoded,
                                          total_size, &deallocator, 0);

    TF_Output* run_path = (TF_Output*)malloc(1 * sizeof(TF_Output));
    run_path[0].oper = checkpoint_op;
    run_path[0].index = 0;

    TF_Tensor** run_path_tensors = (TF_Tensor**)malloc(1 * sizeof(TF_Tensor*));
    run_path_tensors[0] = path_tensor;

    TF_SessionRun(sess,
                  /* RunOptions */ NULL,
                  /* Input tensors */ run_path, run_path_tensors, 1,
                  /* Output tensors */ NULL, NULL, 0,
                  /* Target operations */ &restore_op, 1,
                  /* RunMetadata */ NULL,
                  /* Output status */ status);
    if (TF_GetCode(status) != TF_OK) {
      fprintf(stderr, "ERROR: Unable to run restore_op: %s\n",
              TF_Message(status));
      return 1;
    }
    fprintf(stdout, "Successfully restore session\n");
    
    free((void*)input_encoded);
    return true;
}

bool 
RealNVPScatter::setupSampleTensors() {

    // gerenate input
    // ================================================================================
    TF_Status* status = TF_NewStatus();
    TF_Operation* input_op_x = TF_GraphOperationByName(graph, "x_test_placeholder");
    if (input_op_x == NULL) {
      printf("input_op_x not found\n");
      exit(0);
    } else {
      printf("input_op_x has %i inputs\n", TF_OperationNumOutputs(input_op_x));
    }
    
    TF_Operation* input_op_y = TF_GraphOperationByName(graph, "y_test_placeholder");
    if (input_op_y == NULL) {
      printf("input_op_y not found\n");
      exit(0);
    } else {
      printf("input_op_y has %i inputs\n", TF_OperationNumOutputs(input_op_y));
    }
    
    TF_Operation* input_op_a = TF_GraphOperationByName(graph, "alpha_placeholder");
    if (input_op_a == NULL) {
      printf("input_op_a not found\n");
      exit(0);
    } else {
      printf("input_op_a has %i inputs\n", TF_OperationNumOutputs(input_op_a));
    }
    

    float* raw_input_y = (float*)malloc(sizeof(float));
    raw_input_y[0] = 30.0/180.0 * 3.14159265;
    int64_t* raw_input_dims = (int64_t*)malloc(2 * sizeof(int64_t));
    raw_input_dims[0] = 1;
    raw_input_dims[1] = 1;
    
    float* raw_input_a = (float*)malloc(sizeof(float));
    raw_input_a[0] = 0.5;
 
    int* raw_input_x = (int*)malloc(sizeof(int));
    raw_input_x[0] = 1;


    // prepare inputs

    TF_Tensor* input_tensor_x =
        TF_NewTensor(TF_INT32, NULL, 0, raw_input_x,
                     sizeof(int), &deallocator, NULL);

    TF_Tensor* input_tensor_y =
        TF_NewTensor(TF_FLOAT, raw_input_dims, 2, raw_input_y,
                     sizeof(float), &deallocator, NULL);

    TF_Tensor* input_tensor_a =
        TF_NewTensor(TF_FLOAT, NULL, 0, raw_input_a,
                     sizeof(float), &deallocator, NULL);

    // void* input_data = TF_TensorData(input_tensor);
    // printf("input_data[0] = %f\n", ((float*)input_data)[0]);
    // printf("input_data[1] = %f\n", ((float*)input_data)[1]);

    sample_run_inputs[0].oper = input_op_x;
    sample_run_inputs[0].index = 0;
    sample_run_inputs[1].oper = input_op_y;
    sample_run_inputs[1].index = 0;
    sample_run_inputs[2].oper = input_op_a;
    sample_run_inputs[2].index = 0;

    sample_input_tensors[0] = input_tensor_x;
    sample_input_tensors[1] = input_tensor_y;
    sample_input_tensors[2] = input_tensor_a;

    // prepare outputs
    // ================================================================================
    TF_Operation* output_op = TF_GraphOperationByName(graph, "feng_test_samples");
    printf("output_op has %i outputs\n", TF_OperationNumOutputs(output_op));


    sample_run_outputs[0].oper = output_op;
    sample_run_outputs[0].index = 0;

    float* raw_output_data = (float*)malloc(2 * sizeof(float));
    raw_output_data[0] = 0.f;
    raw_output_data[1] = 0.f;
    int64_t* raw_output_dims = (int64_t*)malloc(1 * sizeof(int64_t));
    raw_output_dims[0] = 2;

    TF_Tensor* output_tensor =
        TF_NewTensor(TF_FLOAT, raw_output_dims, 1, raw_output_data,
                     2 * sizeof(float), &deallocator, NULL);
    sample_output_tensors[0] = output_tensor;

    // run network
    // ================================================================================
    TF_SessionRun(sess,
                  /* RunOptions */ NULL,
                  /* Input tensors */ sample_run_inputs, sample_input_tensors, 3,
                  /* Output tensors */ sample_run_outputs, sample_output_tensors, 1,
                  /* Target operations */ NULL, 0,
                  /* RunMetadata */ NULL,
                  /* Output status */ status);
    if (TF_GetCode(status) != TF_OK) {
      fprintf(stderr, "ERROR: Unable to run output_op: %s\n", TF_Message(status));
      return 1;
    }

    // printf("output-tensor has %i dims\n", TF_NumDims(run_output_tensors[0]));

    float* output_data = (float*) TF_TensorData(sample_output_tensors[0]);

    printf("output %f %f\n", output_data[0], output_data[1] );
   
    /* 
    free((void*)raw_input_x);
    free((void*)raw_input_y);
    free((void*)raw_input_a);
    free((void*)raw_input_dims);
    */

    return true;

}

bool RealNVPScatter::setupEvalTensors() {
    
    int64_t* raw_input_dims = (int64_t*)malloc(2 * sizeof(int64_t));
    
    TF_Operation* input_op_prob_x = TF_GraphOperationByName(graph, "prob_placeholder");
    if (input_op_prob_x == NULL) {
      printf("input_op_prob_x not found\n");
      exit(0);
    } else {
      printf("input_op_prob_x has %i inputs\n", TF_OperationNumOutputs(input_op_prob_x));
    }
 
    TF_Operation* output_prob_op = TF_GraphOperationByName(graph, "dist_op/feng_prob/add");
    if (output_prob_op == NULL) {
      printf("output_prob_op not found\n");
      exit(0);
    } else {
      printf("output_prob_op has %i outputs\n", TF_OperationNumOutputs(output_prob_op));
    }

    float *prob_output = (float*)malloc( sizeof(float)); 
    TF_Tensor* output_prob_tensor =
        TF_NewTensor(TF_FLOAT, NULL, 0, prob_output,
                     sizeof(float), &deallocator, NULL);


    eval_output_tensors[0] = output_prob_tensor;
    float* raw_input_prob_x = (float*)malloc(2*sizeof(float));
    raw_input_prob_x[0] = 0.5;
    raw_input_prob_x[1] = 0.5;
    raw_input_dims[0] = 1;
    raw_input_dims[1] = 2;
    TF_Tensor* input_tensor_prob_x =
        TF_NewTensor(TF_FLOAT, raw_input_dims, 2, raw_input_prob_x,
                     sizeof(float)*2, &deallocator, NULL);

 
    eval_run_outputs[0].oper = output_prob_op;
    eval_run_outputs[0].index = 0;

    eval_run_inputs[0].oper = input_op_prob_x;
    eval_run_inputs[0].index = 0;
    eval_run_inputs[1].oper = sample_run_inputs[1].oper;
    eval_run_inputs[1].index = 0;
    eval_run_inputs[2].oper = sample_run_inputs[2].oper;
    eval_run_inputs[2].index = 0;

    eval_input_tensors[0] = input_tensor_prob_x;
    eval_input_tensors[1] = sample_input_tensors[1];
    eval_input_tensors[2] = sample_input_tensors[2];

    printf("start running prob session\n");
    TF_SessionRun(sess,
                  /* RunOptions */ NULL,
                  /* Input tensors */ eval_run_inputs, eval_input_tensors, 3,
                  /* Output tensors */ eval_run_outputs, eval_output_tensors, 1,
                  /* Target operations */ NULL, 0,
                  /* RunMetadata */ NULL,
                  /* Output status */ status);
    if (TF_GetCode(status) != TF_OK) {
      fprintf(stderr, "ERROR: Unable to run output_op: %s\n", TF_Message(status));
      return 1;
    }
    printf("finished running prob\n");
    
    float* output_data = (float*) TF_TensorData(eval_output_tensors[0]);

    printf("output %f\n", output_data[0]);
    printf("finished running prob\n");


    // you do not want see me creating all the other tensors; Enough lines for
    // this simple example!

    // free up stuff
    // ================================================================================
    // I probably missed something here

    //free((void*)raw_input_x);
    //free((void*)raw_input_y);
    //free((void*)raw_input_a);
    //free((void*)raw_input_dims);
    return true;
}

float 
RealNVPScatter::eval(float thetaI, float alpha, const pbrt::Vector2f& sampleN) {
    float* alphaPtr = (float*) TF_TensorData(eval_input_tensors[2]);
    float* thetaIPtr = (float*) TF_TensorData(eval_input_tensors[1]);
    float* samplePtr = (float*) TF_TensorData(eval_input_tensors[0]);
    samplePtr[0] = sampleN[0];
    samplePtr[1] = sampleN[1];
    alphaPtr[0] = alpha;
    thetaIPtr[0] = thetaI;
    TF_SessionRun(sess,
                  /* RunOptions */ NULL,
                  /* Input tensors */ eval_run_inputs, eval_input_tensors, 3,
                  /* Output tensors */ eval_run_outputs, eval_output_tensors, 1,
                  /* Target operations */ NULL, 0,
                  /* RunMetadata */ NULL,
                  /* Output status */ status);
    if (TF_GetCode(status) != TF_OK) {
      fprintf(stderr, "ERROR: Unable to run output_op: %s\n", TF_Message(status));
      return 0;
    }
    float* output_data = (float*) TF_TensorData(eval_output_tensors[0]);
    if (isNaN(output_data[0])) return 0;
    float prob = expf(output_data[0]);
    return prob;
}

pbrt::Vector2f 
RealNVPScatter::sample(float thetaI, float alpha) {
    float* alphaPtr = (float*) TF_TensorData(sample_input_tensors[2]);
    float* thetaIPtr = (float*) TF_TensorData(sample_input_tensors[1]);
    alphaPtr[0] = alpha;
    thetaIPtr[0] = thetaI;
    TF_SessionRun(sess,
                  /* RunOptions */ NULL,
                  /* Input tensors */ eval_run_inputs, eval_input_tensors, 3,
                  /* Output tensors */ eval_run_outputs, eval_output_tensors, 1,
                  /* Target operations */ NULL, 0,
                  /* RunMetadata */ NULL,
                  /* Output status */ status);
    if (TF_GetCode(status) != TF_OK) {
      fprintf(stderr, "ERROR: Unable to run output_op: %s\n", TF_Message(status));
      return pbrt::Vector2f(0, 0);
    }
    float* output_data = (float*) TF_TensorData(eval_output_tensors[0]);
    return pbrt::Vector2f(output_data[0], output_data[1]);
}



RealNVPScatter::~RealNVPScatter() {
    TF_CloseSession(sess, status);
    TF_DeleteSession(sess, status);

    TF_DeleteStatus(status);
    TF_DeleteBuffer(graph_def);

    TF_DeleteGraph(graph);
}


//MLFresnel Object

MLFresnel::MLFresnel(const std::string& pathPrefix, const std::string& fresnelOpName):
            modelPathPrefix(pathPrefix), fresnelOpName(fresnelOpName){
    std::cout<<"path to fresnel: "<< pathPrefix;
    loadAndRestore();
    setupEvalTensors();
}

bool 
MLFresnel::loadAndRestore() {
    // load graph
    // ================================================================================
    std::string filePath = modelPathPrefix + "/graph.pb";
    std::cout << "path to fresnel: " << filePath;
    graph_def = read_file(filePath.c_str());
    graph = TF_NewGraph();
    
    status = TF_NewStatus();
    TF_ImportGraphDefOptions* opts = TF_NewImportGraphDefOptions();
    TF_GraphImportGraphDef(graph, graph_def, opts, status);
    TF_DeleteImportGraphDefOptions(opts);
    if (TF_GetCode(status) != TF_OK) {
        fprintf(stderr, "ERROR: Unable to import graph %s\n", TF_Message(status));
        return 1;
    }
    fprintf(stdout, "Successfully imported fresnel graph\n");

    // create session
    // ================================================================================
    TF_SessionOptions* opt = TF_NewSessionOptions();
    sess = TF_NewSession(graph, opt, status);
    TF_DeleteSessionOptions(opt);
    if (TF_GetCode(status) != TF_OK) {
        fprintf(stderr, "ERROR: Unable to create session %s\n", TF_Message(status));
        return 1;
    }
    fprintf(stdout, "Successfully created session\n");

    // run init operation
    // ================================================================================
    const TF_Operation* init_op = TF_GraphOperationByName(graph, "init");
    const TF_Operation* const* targets_ptr = &init_op;

    TF_SessionRun(sess,
                  /* RunOptions */ NULL,
                  /* Input tensors */ NULL, NULL, 0,
                  /* Output tensors */ NULL, NULL, 0,
                  /* Target operations */ targets_ptr, 1,
                  /* RunMetadata */ NULL,
                  /* Output status */ status);
    if (TF_GetCode(status) != TF_OK) {
      fprintf(stderr, "ERROR: Unable to run init_op: %s\n", TF_Message(status));
      return 1;
    }
    fprintf(stdout, "Successfully init session\n");

    // run restore
    // ================================================================================
    TF_Operation* checkpoint_op = TF_GraphOperationByName(graph, "save/Const");
    TF_Operation* restore_op = TF_GraphOperationByName(graph, "save/restore_all");

    filePath = modelPathPrefix + "/model";
    const char* checkpoint_path_str = filePath.c_str();
    size_t checkpoint_path_str_len = strlen(checkpoint_path_str);
    size_t encoded_size = TF_StringEncodedSize(checkpoint_path_str_len);

    // The format for TF_STRING tensors is:
    //   start_offset: array[uint64]
    //   data:         byte[...]
    size_t total_size = sizeof(int64_t) + encoded_size;
    char* input_encoded = (char*)malloc(total_size);
    memset(input_encoded, 0, total_size);
    TF_StringEncode(checkpoint_path_str, checkpoint_path_str_len,
                    input_encoded + sizeof(int64_t), encoded_size, status);
    if (TF_GetCode(status) != TF_OK) {
      fprintf(stderr, "ERROR: something wrong with encoding: %s",
              TF_Message(status));
      return 1;
    }

    TF_Tensor* path_tensor = TF_NewTensor(TF_STRING, NULL, 0, input_encoded,
                                          total_size, &deallocator, 0);

    TF_Output* run_path = (TF_Output*)malloc(1 * sizeof(TF_Output));
    run_path[0].oper = checkpoint_op;
    run_path[0].index = 0;

    TF_Tensor** run_path_tensors = (TF_Tensor**)malloc(1 * sizeof(TF_Tensor*));
    run_path_tensors[0] = path_tensor;

    TF_SessionRun(sess,
                  /* RunOptions */ NULL,
                  /* Input tensors */ run_path, run_path_tensors, 1,
                  /* Output tensors */ NULL, NULL, 0,
                  /* Target operations */ &restore_op, 1,
                  /* RunMetadata */ NULL,
                  /* Output status */ status);
    if (TF_GetCode(status) != TF_OK) {
      fprintf(stderr, "ERROR: Unable to run restore_op: %s\n",
              TF_Message(status));
      return 1;
    }
    fprintf(stdout, "Successfully restore session\n");
    
    free((void*)input_encoded);
    return true;
}

bool MLFresnel::setupEvalTensors() {
    
    int64_t* raw_input_dims = (int64_t*)malloc(2 * sizeof(int64_t));
    
    TF_Operation* input_op_prob_x = TF_GraphOperationByName(graph, "fresnel_eval_placeholder");
    if (input_op_prob_x == NULL) {
      printf("input_op_prob_x not found\n");
      exit(0);
    } else {
      printf("input_op_prob_x has %i inputs\n", TF_OperationNumOutputs(input_op_prob_x));
    }
 
    TF_Operation* output_prob_op = TF_GraphOperationByName(graph, fresnelOpName.c_str());
    if (output_prob_op == NULL) {
      printf("output_prob_op not found\n");
      exit(0);
    } else {
      printf("output_prob_op has %i outputs\n", TF_OperationNumOutputs(output_prob_op));
    }

    float *prob_output = (float*)malloc(3 * sizeof(float)); 
    TF_Tensor* output_prob_tensor =
        TF_NewTensor(TF_FLOAT, NULL, 0, prob_output,
                     sizeof(float)*3, &deallocator, NULL);


    eval_output_tensors[0] = output_prob_tensor;
    float* raw_input_prob_x = (float*)malloc(2 * sizeof(float));
    raw_input_prob_x[0] = 0.5;
    raw_input_prob_x[1] = 0.5;
    raw_input_dims[0] = 1;
    raw_input_dims[1] = 2;
    TF_Tensor* input_tensor_prob_x =
        TF_NewTensor(TF_FLOAT, raw_input_dims, 2, raw_input_prob_x,
                     sizeof(float)*2, &deallocator, NULL);

 
    eval_run_outputs[0].oper = output_prob_op;
    eval_run_outputs[0].index = 0;

    eval_run_inputs[0].oper = input_op_prob_x;
    eval_run_inputs[0].index = 0;

    eval_input_tensors[0] = input_tensor_prob_x;

    printf("start running prob session\n");
    TF_SessionRun(sess,
                  /* RunOptions */ NULL,
                  /* Input tensors */ eval_run_inputs, eval_input_tensors, 1,
                  /* Output tensors */ eval_run_outputs, eval_output_tensors, 1,
                  /* Target operations */ NULL, 0,
                  /* RunMetadata */ NULL,
                  /* Output status */ status);
    if (TF_GetCode(status) != TF_OK) {
      fprintf(stderr, "ERROR: Unable to run output_op: %s\n", TF_Message(status));
      return 1;
    }
    printf("finished running prob\n");
    
    float* output_data = (float*) TF_TensorData(eval_output_tensors[0]);

    printf("fresnel output %f %f %f\n", output_data[0], output_data[1], output_data[2]);
    printf("finished running prob\n");

    return true;
}

Spectrum 
MLFresnel::eval(Float alpha, Float muI, Float muH) {
    //float* alphaPtr = (float*) TF_TensorData(eval_input_tensors[2]);
    //float* thetaIPtr = (float*) TF_TensorData(eval_input_tensors[1]);
    float* samplePtr = (float*) TF_TensorData(eval_input_tensors[0]);
    samplePtr[0] = muI;
    samplePtr[1] = muH;
    //alphaPtr[0] = alpha;
    //thetaIPtr[0] = thetaI;
    TF_SessionRun(sess,
                  /* RunOptions */ NULL,
                  /* Input tensors */ eval_run_inputs, eval_input_tensors, 1,
                  /* Output tensors */ eval_run_outputs, eval_output_tensors, 1,
                  /* Target operations */ NULL, 0,
                  /* RunMetadata */ NULL,
                  /* Output status */ status);
    if (TF_GetCode(status) != TF_OK) {
      fprintf(stderr, "ERROR: Unable to run output_op: %s\n", TF_Message(status));
      return 0;
    }
    float* output_data = (float*) TF_TensorData(eval_output_tensors[0]);
    if (isNaN(output_data[0])) return 0;
    //std::cout<<"ML Fresnel: "<< output_data[0] << " " << output_data[1] << " " << output_data[2] << "\n";
    return Spectrum::FromRGB(output_data);
}


MLFresnel::~MLFresnel() {
    TF_CloseSession(sess, status);
    TF_DeleteSession(sess, status);

    TF_DeleteStatus(status);
    TF_DeleteBuffer(graph_def);
    TF_DeleteGraph(graph);
}

///////////////end MLFresnel///////////////////////////////////
//Required MISC functions

TF_Buffer* read_file(const char* file) {

    FILE* f = fopen(file, "rb");
    if (!f) {
        std::cout << file << " missing\n";
        return NULL;
    }
    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);  // same as rewind(f);

    void* data = malloc(fsize);
    fread(data, fsize, 1, f);
    fclose(f);

    TF_Buffer* buf = TF_NewBuffer();
    buf->data = data;
    buf->length = fsize;
    buf->data_deallocator = free_buffer;
    return buf;
}
}
