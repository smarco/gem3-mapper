/*
 *  GEM-Cutter "Highly optimized genomic resources for GPUs"
 *  Copyright (c) 2011-2018 by Alejandro Chacon    <alejandro.chacond@gmail.com>
 *
 *  Licensed under GNU General Public License 3.0 or later.
 *  Some rights reserved. See LICENSE, AUTHORS.
 *  @license GPL-3.0+ <http://www.gnu.org/licenses/gpl-3.0.en.html>
 */

#ifndef GPU_FMI_PRIMITIVES_DECODE_C_
#define GPU_FMI_PRIMITIVES_DECODE_C_

#include "../include/gpu_fmi_primitives.h"
#include "../include/gpu_sa_primitives.h"

/************************************************************
Functions to get the GPU FMI buffers
************************************************************/

gpu_fmi_decode_init_pos_t* gpu_fmi_decode_buffer_get_init_pos_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.decode.initPositions.h_initBWTPos);
}

gpu_fmi_decode_end_pos_t* gpu_fmi_decode_buffer_get_end_pos_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.decode.endPositions.h_endBWTPos);
}

gpu_fmi_decode_text_pos_t* gpu_fmi_decode_buffer_get_ref_pos_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.decode.textPositions.h_textPos);
}

/************************************************************
Functions to get the maximum elements of the buffers
************************************************************/

uint32_t gpu_fmi_decode_buffer_get_max_positions_(const void* const fmiBuffer){
  const gpu_buffer_t* const mBuff = (gpu_buffer_t *) fmiBuffer;
  return(mBuff->data.decode.numMaxInitPositions);
}

/************************************************************
Functions to initialize the buffers (DECODE)
************************************************************/

size_t gpu_fmi_decode_input_size()
{
  return(sizeof(gpu_fmi_decode_init_pos_t) + sizeof(gpu_fmi_decode_end_pos_t) + sizeof(gpu_fmi_decode_text_pos_t));
}

void gpu_fmi_decode_reallocate_host_buffer_layout(gpu_buffer_t* const mBuff)
{
  const void* rawAlloc = mBuff->h_rawData;
  //Adjust the host buffer layout (input)
  mBuff->data.decode.initPositions.h_initBWTPos = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.decode.initPositions.h_initBWTPos + mBuff->data.decode.numMaxInitPositions);
  //Adjust the host buffer layout (intermediate-results)
  mBuff->data.decode.endPositions.h_endBWTPos   = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.decode.endPositions.h_endBWTPos + mBuff->data.decode.numMaxEndPositions);
  //Adjust the host buffer layout (output)
  mBuff->data.decode.textPositions.h_textPos    = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.decode.textPositions.h_textPos + mBuff->data.decode.numMaxTextPositions);
}

void gpu_fmi_decode_reallocate_device_buffer_layout(gpu_buffer_t* const mBuff)
{
  const void* rawAlloc = mBuff->d_rawData;
  //Adjust the host buffer layout (input)
  mBuff->data.decode.initPositions.d_initBWTPos = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.decode.initPositions.d_initBWTPos + mBuff->data.decode.numMaxInitPositions);
  //Adjust the host buffer layout (output)
  mBuff->data.decode.endPositions.d_endBWTPos   = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.decode.endPositions.d_endBWTPos + mBuff->data.decode.numMaxEndPositions);
  //Adjust the host buffer layout (output)
  mBuff->data.decode.textPositions.d_textPos    = GPU_ALIGN_TO(rawAlloc,16);
  rawAlloc = (void *) (mBuff->data.decode.textPositions.d_textPos + mBuff->data.decode.numMaxTextPositions);
}

void gpu_fmi_decode_init_buffer_(void* const fmiBuffer)
{
  gpu_buffer_t* const mBuff            = (gpu_buffer_t *) fmiBuffer;
  const size_t        sizeBuff         = mBuff->sizeBuffer * 0.95;
  const uint32_t      numMaxPositions  = sizeBuff / gpu_fmi_decode_input_size();

  //set the type of the buffer
  mBuff->typeBuffer = GPU_FMI_DECODE_POS | GPU_SA_DECODE_POS;

  //Set real size of the input
  mBuff->data.decode.numMaxInitPositions = numMaxPositions;
  mBuff->data.decode.numMaxEndPositions  = numMaxPositions;
  mBuff->data.decode.numMaxTextPositions = numMaxPositions;
  gpu_fmi_decode_reallocate_host_buffer_layout(mBuff);
  gpu_fmi_decode_reallocate_device_buffer_layout(mBuff);
}

void gpu_fmi_decode_init_and_realloc_buffer_(void* const fmiBuffer, const uint32_t numDecodes)
{
  gpu_buffer_t* const mBuff                   = (gpu_buffer_t *) fmiBuffer;
  const uint32_t      idSupDevice             = mBuff->idSupportedDevice;
  const float         resizeFactor            = 2.0;
  const size_t        bytesPerDecodeBuffer    = numDecodes * gpu_fmi_decode_input_size();

  //printf("RESIZE[FMI_DECODE] %d %d \n",  mBuff->sizeBuffer, bytesPerDecodeBuffer * resizeFactor);
  //Recalculate the minimum buffer size
  mBuff->sizeBuffer = bytesPerDecodeBuffer * resizeFactor;

  //FREE HOST AND DEVICE BUFFER
  GPU_ERROR(gpu_buffer_free(mBuff));

  //Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));

  //ALLOCATE HOST AND DEVICE BUFFER
  CUDA_ERROR(cudaHostAlloc((void**) &mBuff->h_rawData, mBuff->sizeBuffer, cudaHostAllocMapped));
  CUDA_ERROR(cudaMalloc((void**) &mBuff->d_rawData, mBuff->sizeBuffer));

  gpu_fmi_decode_init_buffer_(fmiBuffer);
}


/************************************************************
Functions to transfer data HOST <-> DEVICE (Decode)
************************************************************/

gpu_error_t gpu_fmi_decode_transfer_CPU_to_GPU(gpu_buffer_t* const mBuff)
{
  const gpu_fmi_decode_init_pos_buffer_t* initPosBuff = &mBuff->data.decode.initPositions;
  const cudaStream_t                      idStream    =  mBuff->listStreams[mBuff->idStream];
  const size_t                            cpySize     =  initPosBuff->numDecodings * sizeof(gpu_fmi_decode_init_pos_t);

  //Transfer seeds from CPU to the GPU
  CUDA_ERROR(cudaMemcpyAsync(initPosBuff->d_initBWTPos, initPosBuff->h_initBWTPos, cpySize, cudaMemcpyHostToDevice, idStream));

  return (SUCCESS);
}

gpu_error_t gpu_fmi_decode_transfer_GPU_to_CPU(gpu_buffer_t* const mBuff)
{
  const gpu_fmi_decode_end_pos_buffer_t   *endPosBuff  = &mBuff->data.decode.endPositions;
  const gpu_fmi_decode_text_pos_buffer_t  *textPosBuff = &mBuff->data.decode.textPositions;
  const cudaStream_t                      idStream     =  mBuff->listStreams[mBuff->idStream];
        size_t                            cpySize      =  0;

  //Transfer SA intervals (occurrence results) from CPU to the GPU
  if(mBuff->index->activeModules & GPU_SA_DECODE_POS){
    cpySize = endPosBuff->numDecodings * sizeof(gpu_fmi_decode_text_pos_t);
    CUDA_ERROR(cudaMemcpyAsync(textPosBuff->h_textPos, textPosBuff->d_textPos, cpySize, cudaMemcpyDeviceToHost, idStream));
  }else{
    cpySize = endPosBuff->numDecodings * sizeof(gpu_fmi_decode_end_pos_t);
    CUDA_ERROR(cudaMemcpyAsync(endPosBuff->h_endBWTPos, endPosBuff->d_endBWTPos, cpySize, cudaMemcpyDeviceToHost, idStream));
  }

  return (SUCCESS);
}

void gpu_fmi_decode_send_buffer_(void* const fmiBuffer, const uint32_t numDecodings, const uint32_t samplingRate)
{
  gpu_buffer_t* const mBuff  = (gpu_buffer_t *) fmiBuffer;
  const uint32_t idSupDevice = mBuff->idSupportedDevice;

  //Set real size of the input
  mBuff->data.decode.initPositions.numDecodings = numDecodings;
  mBuff->data.decode.endPositions.numDecodings  = numDecodings;
  mBuff->data.decode.textPositions.numDecodings = numDecodings;
  mBuff->data.decode.samplingRate               = samplingRate;

  //Select the device of the Multi-GPU platform
  CUDA_ERROR(cudaSetDevice(mBuff->device[idSupDevice]->idDevice));

  GPU_ERROR(gpu_fmi_decode_transfer_CPU_to_GPU(mBuff));
  GPU_ERROR(gpu_fmi_decode_process_buffer(mBuff));
  if(mBuff->index->activeModules & GPU_SA_DECODE_POS)
    GPU_ERROR(gpu_sa_decode_process_buffer(mBuff));
  GPU_ERROR(gpu_fmi_decode_transfer_GPU_to_CPU(mBuff));
}

void gpu_fmi_decode_receive_buffer_(const void* const fmiBuffer)
{
  const gpu_buffer_t* const mBuff    = (gpu_buffer_t *) fmiBuffer;
  const cudaStream_t        idStream =  mBuff->listStreams[mBuff->idStream];
  //Synchronize Stream (the thread wait for the commands done in the stream)
  CUDA_ERROR(cudaStreamSynchronize(idStream));
}

#endif /* GPU_FMI_PRIMITIVES_DECODE_C_ */

