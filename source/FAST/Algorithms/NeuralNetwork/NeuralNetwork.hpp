#ifndef NEURAL_NETWORK_HPP_
#define NEURAL_NETWORK_HPP_

#include "FAST/ProcessObject.hpp"
#include <tensorflow/core/public/session.h>

namespace fast {

class Image;

class NeuralNetwork : public ProcessObject {
public:
    void load(std::string networkFilename);
    void setInputParameters(std::string inputNodeNames, int width, int height);
    void setOutputParameters(std::vector<std::string> > outputNodeNames);

    // Use this if only 1 network output layer
    std::vector<float> getNetworkOutput();

    // Get output by layer name
    std::vector<float> getNetworkOutput(std::string layerName);

protected:
    NeuralNetwork();
    UniquePointer<tensorflow::Session> mSession;
    bool mModelLoaded;
    int mWidth;
    int mHeight;
    std::string mInputName;
    std::vector<std::string> > mOutputNames;

    void execute();

    void executeNetwork(const std::vector<SharedPointer<Image> >& images);
    std::vector<SharedPointer<Image> > resizeImages(const std::vector<SharedPointer<Image> >& images);
    std::vector<SharedPointer<Image> > subtractMeanImage(const std::vector<SharedPointer<Image> >& images);
};

}

#endif