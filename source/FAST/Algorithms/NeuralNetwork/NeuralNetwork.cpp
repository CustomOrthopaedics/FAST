#include "NeuralNetwork.hpp"
#include "FAST/Data/Image.hpp"
#include "FAST/Algorithms/ImageResizer/ImageResizer.hpp"

#include <tensorflow/core/framework/step_stats.pb.h>
#include <tensorflow/core/framework/tensor.h>
#include <tensorflow/core/framework/types.pb.h>
#include <tensorflow/core/lib/strings/stringprintf.h>
#include <tensorflow/core/platform/env.h>
#include <tensorflow/core/platform/logging.h>
#include <tensorflow/core/platform/mutex.h>
#include <tensorflow/core/platform/types.h>
#include <tensorflow/core/public/session.h>

namespace fast {

// See here for reference: https://github.com/tensorflow/tensorflow/blob/86f5ab7474825da756838b34e1b4eac93f5fc68a/tensorflow/contrib/android/jni/tensorflow_inference_jni.cc

void NeuralNetwork::load(std::string networkFilename) {

	tensorflow::SessionOptions options;
	tensorflow::ConfigProto &config = options.config;
	mSession.reset(tensorflow::NewSession(options));
	tensorflow::GraphDef tensorflow_graph;

    ReadBinaryProto(tensorflow::Env::Default(), networkFilename, &tensorflow_graph);

	reportInfo() << "Creating session." << reportEnd();
	tensorflow::Status s = mSession->Create(tensorflow_graph);
	if (!s.ok()) {
		throw Exception("Could not create TensorFlow Graph");
	}

	// Clear the proto to save memory space.
	tensorflow_graph.Clear();
	reportInfo() << "TensorFlow graph loaded from: " << networkFilename << reportEnd();

	mModelLoaded = true;
}


NeuralNetwork::NeuralNetwork() {
	createInputPort<Image>(0, true, INPUT_STATIC_OR_DYNAMIC, true);
	mModelLoaded = false;
	mMeanImageLoaded = false;
	createOpenCLProgram(std::string(FAST_SOURCE_DIR) + "Algorithms/NeuralNetwork/NeuralNetwork.cl");
}

void NeuralNetwork::execute() {
    if(!mModelLoaded)
		throw Exception("Network and weights must be loaded in NeuralNetwork before execution.");

	std::vector<Image::pointer> images = getMultipleStaticInputData<Image>();

	/*
	// Assuming only one input layer here
	caffe::Blob<float>* input_layer = mNetwork->input_blobs()[0];

    if(input_layer->num() != images.size()) {
		// Only reshape if necessary
		// batch size x channels x width x height
		input_layer->Reshape(images.size(), input_layer->channels(), input_layer->height(), input_layer->width());
		mNetwork->Reshape();
		reportInfo() << "Net reshaped" << reportEnd();
	}

	// Pre processing, subtract mean, resize, crop, convert to normalized float etc.

    // TODO Resize images if necessary
    images = resizeImages(images);

	// Subtract mean
    if(mMeanImageLoaded) {
		// Mean image has been loaded, subtract it from each image
		images = subtractMeanImage(images);
	}

	// TODO normalization?
	 */

	executeNetwork(images);
}

std::vector<float> NeuralNetwork::getNetworkOutput() {
	/*
	caffe::Blob<float>* layerData = mNetwork->output_blobs()[0];
	std::vector<float> result(
			layerData->cpu_data(),
            layerData->cpu_data() + layerData->num()*layerData->channels()*layerData->height()*layerData->width()
	);
    return result;
    */
}

std::vector<float> NeuralNetwork::getNetworkOutput(std::string layerName) {
	/*
    boost::shared_ptr<caffe::Blob<float> > layerData = mNetwork->blob_by_name(layerName);
    std::vector<float> result(
			layerData->cpu_data(),
			layerData->cpu_data() + layerData->num()*layerData->channels()*layerData->height()*layerData->width()
	);
	return result;
	 */
}

void NeuralNetwork::executeNetwork(const std::vector<Image::pointer>& images) {
	// Create input tensor
	tensorflow::Tensor input_tensor(
			tensorflow::DT_FLOAT,
			tensorflow::TensorShape({1, mHeight, mWidth, 1})
	);

	auto input_tensor_mapped = input_tensor.tensor<float, 4>();

    // TODO support multiple images
	// TODO need to reshape network to support this
	reportInfo() << "TensorFlow: Copying Data." << reportEnd();
    Image::pointer image = images[0];
	ImageAccess::pointer access = image->getImageAccess(ACCESS_READ);
	for (int i = 0; i < mHeight; ++i) { // y
		for (int j = 0; j < mWidth; ++j) { // x
			input_tensor_mapped(0, i, j, 0) = access->getScalar(Vector2i(j, i)) / 255;
		}
	}

    // TODO Need to know names of inputs and outputs in advance
	// Input: Only single for now
	// Output: Can be multiple

	std::vector <std::pair<std::string, tensorflow::Tensor>> input_tensors(
			{{mInputName, input_tensor}});

	std::vector <tensorflow::Tensor> output_tensors;
	std::vector <std::string> output_names({mOutputName});

	tensorflow::Status s;
	s = mSession->Run(input_tensors, output_names, {}, &output_tensors);

	if (!s.ok()) {
		reportError() << "Error during inference: " << s << reportEnd();
	}

	tensorflow::Tensor *output = &output_tensors[0];

	const auto outputData = output->flat<float>(); // This is some sort of Eigen tensor type
	std::vector<float> resultData;
	for (int i = 0; i < outputData.size(); ++i) {
		resultData.push_back(outputData(i));
	}


	/*
 	// TODO, gives images directly to GPU
	// Set images to input layer
	reportInfo() << "Adding input data to input layer.." << reportEnd();
	int counter = 0;
	caffe::Blob<float>* input_layer = mNetwork->input_blobs()[0];
	float* input_data = input_layer->mutable_cpu_data(); // This is the input data layer
	for(Image::pointer image : images) {
        // Some sanity checks
		Vector3ui size = image->getSize();
		if(size.x() != input_layer->width() || size.y() != input_layer->height())
			throw Exception("Mismatch between size of input layer and an image given to NeuralNetwork::executeNetwork()");
		if(image->getDataType() != TYPE_FLOAT)
			throw Exception("Data type of images feed to executeNetwork() must be float");
		ImageAccess::pointer access = image->getImageAccess(ACCESS_READ);
		float* pixels = (float*)access->get();
		for(int i = 0; i < size.x()*size.y(); ++i) {
			input_data[counter + i] = pixels[i];
		}
		counter += size.x()*size.y();
	}

	// Do a forward pass
	reportInfo() << "Neural network foward pass executing..." << reportEnd();
	mNetwork->Forward();
	reportInfo() << "Neural network foward pass finished" << reportEnd();
	 */
}

std::vector<SharedPointer<Image>> NeuralNetwork::resizeImages(const std::vector<SharedPointer<Image>> &images) {
	/*
	reportInfo() << "Resizing images.." << reportEnd();
	caffe::Blob<float>* input_layer = mNetwork->input_blobs()[0];
    std::vector<Image::pointer> resizedImages;
	for(Image::pointer image : images) {
		// Resize image to fit input layer
		if(input_layer->width() != image->getWidth() || input_layer->height() != image->getHeight()) {
			// Only resize if needed
            ImageResizer::pointer resizer = ImageResizer::New();
            resizer->setWidth(input_layer->width());
            resizer->setHeight(input_layer->height());
            resizer->setInputData(image);
            resizer->update();
            Image::pointer resizedImage = resizer->getOutputData<Image>();
            resizedImages.push_back(resizedImage);
		} else {
			resizedImages.push_back(image);
		}
	}

	return resizedImages;
	 */
}

void NeuralNetwork::loadBinaryMeanImage(std::string filename) {
	/*
	// Network must be loaded first
	if(!mNetwork.isValid()) {
		throw Exception("Must load network definition before loading weights");
	}

	reportInfo() << "Loading mean image file.." << reportEnd();
	caffe::BlobProto blob_proto;
	caffe::ReadProtoFromBinaryFileOrDie(filename.c_str(), &blob_proto);

	caffe::Blob<float> mMeanBlob;
	mMeanBlob.FromProto(blob_proto);
	OpenCLDevice::pointer device = getMainDevice();
	caffe::Blob<float>* input_layer = mNetwork->input_blobs()[0];
    // Assuming float image here
	mMeanImage = cl::Image2D(
			device->getContext(),
			CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			getOpenCLImageFormat(device, CL_MEM_OBJECT_IMAGE2D, TYPE_FLOAT, input_layer->channels()),
			input_layer->width(), input_layer->height(),
			0,
			mMeanBlob.mutable_cpu_data()
	);
	reportInfo() << "Finished loading mean image file." << reportEnd();

	mMeanImageLoaded = true;
	 */
}

std::vector<SharedPointer<Image> > NeuralNetwork::subtractMeanImage(const std::vector<SharedPointer<Image> >& images) {
	if(!mMeanImageLoaded) {
		throw Exception("Mean image not loaded, cannot subtract mean image from images");
	}

	reportInfo() << "Subtracting mean image.." << reportEnd();
	OpenCLDevice::pointer device = getMainDevice();
	cl::Program program = getOpenCLProgram(device);
	cl::Kernel normalizationKernel(program, "subtractMeanImage");
	normalizationKernel.setArg(1, mMeanImage);

	std::vector<Image::pointer> preProcessedImages;
	for(auto &&image : images) {
        Image::pointer preProcessedImage = Image::New();
		preProcessedImage->create(image->getSize(), TYPE_FLOAT, 1);
		OpenCLImageAccess::pointer access = image->getOpenCLImageAccess(ACCESS_READ, device);
		OpenCLImageAccess::pointer access2 = preProcessedImage->getOpenCLImageAccess(ACCESS_READ_WRITE, device);
		normalizationKernel.setArg(0, *(access->get2DImage()));
		normalizationKernel.setArg(2, *(access2->get2DImage()));

		device->getCommandQueue().enqueueNDRangeKernel(
				normalizationKernel,
				cl::NullRange,
				cl::NDRange(image->getWidth(), image->getHeight()),
				cl::NullRange
		);

		preProcessedImages.push_back(preProcessedImage);

	}

	return preProcessedImages;
}


};
