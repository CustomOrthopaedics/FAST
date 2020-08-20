#ifndef CENTERLINE_EXTRACTION_HPP_
#define CENTERLINE_EXTRACTION_HPP_

#include "FAST/ProcessObject.hpp"

namespace fast {

class Image;

class FAST_EXPORT  CenterlineExtraction : public ProcessObject {
	FAST_OBJECT(CenterlineExtraction)
    public:
    void setOffset(Vector3f off);
    private:
		CenterlineExtraction();
    Vector3f offset;
		void execute();
        SharedPointer<Image> calculateDistanceTransform(SharedPointer<Image> input);
};

}

#endif
