#include "AirwaySegmentation.hpp"
#include "FAST/Algorithms/GaussianSmoothingFilter/GaussianSmoothingFilter.hpp"
#include "FAST/Data/Segmentation.hpp"
#include <unordered_set>
#include <stack>
#include <queue>
#include <vector>
#include <algorithm>
#include <math.h>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>

using namespace boost::accumulators;

namespace fast {

Vector3f icohalfMm[21] = {
    Vector3f(0.276388, 0.447220, 0.850649),
    Vector3f(-0.723607, 0.447220, 0.525725),
    Vector3f(-0.723607, 0.447220, -0.525725),
    Vector3f(0.276388, 0.447220, -0.850649),
    Vector3f(0.894426, 0.447216, 0.000000),
    Vector3f(0.000000, 1.000000, 0.000000),
    Vector3f(0.951058, 0.000000, -0.309013),
    Vector3f(-0.587786, 0.000000, -0.809017),
    Vector3f(-0.951058, 0.000000, -0.309013),
    Vector3f(0.587786, 0.000000, -0.809017),
    Vector3f(0.000000, 0.000000, -1.000000),
    Vector3f(0.688189, 0.525736, 0.499997),
    Vector3f(-0.262869, 0.525738, 0.809012),
    Vector3f(-0.850648, 0.525736, 0.000000),
    Vector3f(-0.262869, 0.525738, -0.809012),
    Vector3f(0.688189, 0.525736, -0.499997),
    Vector3f(0.162456, 0.850654, 0.499995),
    Vector3f(0.525730, 0.850652, 0.000000),
    Vector3f(-0.425323, 0.850654, 0.309011),
    Vector3f(-0.425323, 0.850654, -0.309011),
    Vector3f(0.162456, 0.850654, -0.499995)
};

float deltaW = 200.0;
float dr = 0.5;
float rMax = 20.0;
float maxRadiusIncrease = 2.0;
float maxAirwayDensity = -550.0;
int maxVoxelVal = 1000;

// https://stackoverflow.com/questions/19271568/trilinear-interpolation
float interpolate1D(float v1, float v2, float x){
    return v1*(1-x) + v2*x;
}

float interpolate2D(float v1, float v2, float v3, float v4, float x, float y) {
    float s = interpolate1D(v1, v2, x);
    float t = interpolate1D(v3, v4, x);
    return interpolate1D(s, t, y);
}

float interpolate3D(float *verts, Vector3f point) {
    float s = interpolate2D(verts[0], verts[4], verts[2], verts[6], point.x(), point.y());
    float t = interpolate2D(verts[1], verts[5], verts[3], verts[7], point.x(), point.y());
    return interpolate1D(s, t, point.z());
}

int AirwaySegmentation::getIndex(Vector3i point) {
	return point.x() + point.y()*width + point.z()*width*height;
}

int AirwaySegmentation::getIndex(int x, int y, int z) {
	return x + y*width + z*width*height;
}

static float calcDistance(Vector3f start, Vector3f end) {
	return sqrt(pow(end.x() - start.x(), 2) + pow(end.y() - start.y(), 2) + pow(end.z() - start.z(), 2));
}

int AirwaySegmentation::interp3D(short *vol, Vector3f point) {
	Vector3i pointIdx(floor(point.x()), floor(point.y()), floor(point.z()));

	float xPoint = point.x() - floor(point.x());
	float yPoint = point.y() - floor(point.y());
	float zPoint = point.z() - floor(point.z());
	Vector3f localPoint(xPoint, yPoint, zPoint);

	float vals[8];
	float *verts = vals;
	int valIdx = 0;

	for (int x = 0; x <= 1; ++x) {
		for (int y = 0; y <= 1; ++y) {
			for (int z = 0; z <= 1; ++z) {
				Vector3i neighborIdx = Vector3i(x, y, z) + pointIdx;

				// reached edge of volume, return max voxel value to stop ray from growing any further
				if(neighborIdx.x() >= width || neighborIdx.y() >= height || neighborIdx.z() >= depth) {
					return maxVoxelVal;
				}

				int voxValue = vol[getIndex(neighborIdx)];
				vals[valIdx++] = voxValue;
			}
		}
	}

	return interpolate3D(verts, localPoint);
}

float AirwaySegmentation::findDistanceToWall(short *vol, Vector3f dir, Vector3i startPoint) {
	int startVal = vol[getIndex(startPoint)];
	int endVal = 0;
	float distance = 0.0;

	do {
		distance += dr;
		Vector3f endPoint = (dir * distance) + Vector3f((float)startPoint.x(), (float)startPoint.y(), (float)startPoint.z());
		endVal = interp3D(vol, endPoint);
	} while(endVal - startVal < deltaW && distance < rMax);

	return distance * mmPerVx;
}

Voxel AirwaySegmentation::getVoxelData(short *vol, Vector3i point) {
	std::vector<float> diameters;
	for (int i = 0; i < 21; i++) {
		float dist1 = findDistanceToWall(vol, icohalfVx[i], point);
		float dist2 = findDistanceToWall(vol, -1.0 * icohalfVx[i], point);

		diameters.push_back(dist1 + dist2);
	}

	std::sort(diameters.begin(), diameters.end());

	// find lower half of diameters
	accumulator_set<float, stats<tag::variance> > acc;
	for (int i = 0; i < 10; i++) {
		acc(diameters[i]);
	}

	float radiiMean = mean(acc) / 2.0;
	float radiStdDev = sqrt(variance(acc)) / 2.0;

	float centricity = 1 - (radiStdDev / radiiMean);
	float smallestRadius = diameters[0] / 2.0;

	return Voxel(point, centricity, smallestRadius, radiiMean);
}

AirwaySegmentation::AirwaySegmentation() {
	createInputPort<Image>(0);
	createOutputPort<Segmentation>(0);

	createOpenCLProgram(Config::getKernelSourcePath() + "Algorithms/AirwaySegmentation/AirwaySegmentation.cl");
}

Vector3i AirwaySegmentation::findSeedVoxel(Image::pointer volume) {

	ImageAccess::pointer access = volume->getImageAccess(ACCESS_READ);
	short* data = (short*)access->get();

    int slice = volume->getDepth()*0.6;

    int threshold = -700;
    int Tseed = -950;
    float minArea = 7.5f*7.5f*3.14f; // min radius of 7.5
    float maxArea = 25.0f*25.0f*3.14f; // max radius of 25.0
    Vector3i currentSeed(0,0,0);
    float currentCentricity = 99999.0f; // Distance from center
    float spacing = 1.0f;//volume->getSpacing().x();

    std::unordered_set<int> visited;
    int width = volume->getWidth();
    int height = volume->getHeight();
    int depth = volume->getDepth();

    for(int x = width*0.25; x < width*0.75; ++x) {
        for(int y = height*0.25; y < height*0.75; ++y) {
            Vector3i testSeed(x,y,slice);
            if(data[testSeed.x() + testSeed.y()*width + testSeed.z()*width*height] > Tseed)
                continue;
            std::stack<Vector3i> stack;
            stack.push(testSeed);
            int perimenter = 0;
            bool invalid = false;
            visited.clear();
            visited.insert(testSeed.x()+testSeed.y()*width);

            while(!stack.empty() && !invalid) {
                Vector3i v = stack.top();
                stack.pop();

                for(int a = -1; a < 2 && !invalid; ++a) {
                for(int b = -1; b < 2; ++b) {
                    Vector3i c(v.x()+a, v.y()+b, v.z());
                    if(c.x() < 0 || c.y() < 0 ||
                    		c.x() >= width || c.y() >= height) {
                        invalid = true;
                        break;
                    }

                    if(data[c.x() + c.y()*width + c.z()*width*height] <= threshold && visited.find(c.x()+c.y()*volume->getWidth()) == visited.end()) {
                        visited.insert(c.x()+c.y()*volume->getWidth());
                        stack.push(c);
                        if(visited.size()*spacing*spacing > maxArea) {
                            invalid = true;
                            break;
                        }
                    }
                }}
            }

            //float compaction = (4.0f*3.14*area)/(perimenter*perimenter);
            if(!invalid && visited.size()*spacing*spacing > minArea) {
                float centricity = sqrt(pow(testSeed.x()-volume->getWidth()*0.5f,2.0f)+pow(testSeed.y()-volume->getHeight()*0.5f,2.0f));
                if(centricity < currentCentricity) {
                    // Accept as new seed
                    currentSeed = testSeed;
                    currentCentricity = centricity;
                }
            }
        }
    }

    return currentSeed;
}

int AirwaySegmentation::grow(Vector3i seed, uchar* mask, std::vector<Vector3i> neighbors, short* data, float threshold) {
    std::priority_queue<Voxel> queue;

	float pathMinRadius = 9999.0;
	queue.push(Voxel(seed, 1.0, pathMinRadius, 1.0));

	int maskSize = 0;
	int volSize = height * width * depth;

	// stop at half vol size in case of major leaks
    while (!queue.empty() && maskSize < volSize / 2) {
        Voxel currVox = queue.top();
        queue.pop();

		pathMinRadius = currVox.minRadius;

		Vector3i x = currVox.point;
        mask[getIndex(x)] = 1;

        // Add 26 neighbors
        for (int i = 0; i < 25; ++i) {
        	Vector3i neighbor = neighbors[i];
            Vector3i y(x.x()+neighbor.x(), x.y()+neighbor.y(), x.z()+neighbor.z());
			int volIdx = getIndex(y);

			// outside of volume
			if(y.x() < 0 || y.y() < 0 || y.z() < 0 || y.x() >= width || y.y() >= height || y.z() >= depth) {
                continue;
            }

			// vox already in mask
			if (mask[volIdx] != 0) {
				continue;
			}

			// above threshold
			if (data[volIdx] > threshold) {
				continue;
			}

			// add vox to mask
			mask[volIdx] = 1;
			maskSize++;

			Voxel vox = getVoxelData(data, y);

			// radius is too large, voxel may be leaking
			if (vox.meanRadii > pathMinRadius * maxRadiusIncrease) {
				continue;
			}

			if (pathMinRadius < vox.minRadius) {
				vox.minRadius = pathMinRadius;
			}

			queue.push(vox);
        }
    }

    return maskSize;
}

void AirwaySegmentation::regionGrowing(Image::pointer volume, Segmentation::pointer segmentation, const std::vector<Vector3i> seeds) {
	segmentation->createFromImage(volume);
	ImageAccess::pointer volAccess = volume->getImageAccess(ACCESS_READ);
	short* volData = (short*)volAccess->get();

	ImageAccess::pointer segAccess = segmentation->getImageAccess(ACCESS_READ_WRITE);
	uchar* segData = (uchar*)segAccess->get();
	memset(segData, 0, width*height*depth);

	uchar* seedMask = (uchar*) malloc(width*height*depth);

	// Create neighbor list
	std::vector<Vector3i> neighborList;
	for(int a = -1; a < 2; a++) {
	for(int b = -1; b < 2; b++) {
	for(int c = -1; c < 2; c++) {
		if(a == 0 && b == 0 && c == 0)
			continue;
		neighborList.push_back(Vector3i(a,b,c));
	}}}

	for (Vector3i seed : seeds) {
		// reset mask
		memset(seedMask, 0, width*height*depth);

		Reporter::info() << "Segmenting Seed: " << seed.transpose() << Reporter::end();

		// perform region growing, store results in seedMask
		grow(seed, seedMask, neighborList, volData, maxAirwayDensity);

		// add mask to final seg
		for (int x = 0; x < width; ++x) {
			for (int y = 0; y < height; ++y) {
				for (int z = 0; z < depth; ++z) {
					int index = getIndex(x, y, z);
					if (seedMask[index] == 1) {
						segData[index] = 1;
					}
				}
			}
		}
	}
	free(seedMask);
}

Image::pointer AirwaySegmentation::convertToHU(Image::pointer image) {
	OpenCLDevice::pointer device = std::dynamic_pointer_cast<OpenCLDevice>(getMainDevice());
	cl::Program program = getOpenCLProgram(device);

	OpenCLImageAccess::pointer input = image->getOpenCLImageAccess(ACCESS_READ, device);
	Image::pointer newImage = Image::New();
	newImage->create(image->getSize(), TYPE_INT16, 1);
	newImage->setSpacing(image->getSpacing());
	SceneGraph::setParentNode(newImage, image);

	cl::Kernel kernel(program, "convertToHU");

	kernel.setArg(0, *input->get3DImage());
	if(device->isWritingTo3DTexturesSupported()) {
		OpenCLImageAccess::pointer output = newImage->getOpenCLImageAccess(ACCESS_READ_WRITE, device);
		kernel.setArg(1, *output->get3DImage());
	} else {
		OpenCLBufferAccess::pointer output = newImage->getOpenCLBufferAccess(ACCESS_READ_WRITE, device);
		kernel.setArg(1, *output->get());
	}

	device->getCommandQueue().enqueueNDRangeKernel(
			kernel,
			cl::NullRange,
			cl::NDRange(image->getWidth(), image->getHeight(), image->getDepth()),
			cl::NullRange
    );

	return newImage;
}

void AirwaySegmentation::morphologicalClosing(Segmentation::pointer segmentation) {
	width = segmentation->getWidth();
	height = segmentation->getHeight();
	depth = segmentation->getDepth();

	// TODO need support for no 3d write
	OpenCLDevice::pointer device = std::dynamic_pointer_cast<OpenCLDevice>(getMainDevice());
	cl::Program program = getOpenCLProgram(device);

	Segmentation::pointer segmentation2 = Segmentation::New();
	segmentation2->create(segmentation->getSize(), TYPE_UINT8, 1);
	ImageAccess::pointer access = segmentation2->getImageAccess(ACCESS_READ_WRITE);
	uchar* data = (uchar*)access->get();
	memset(data, 0, width*height*depth);
	access->release();

	cl::Kernel dilateKernel(program, "dilate");
	cl::Kernel erodeKernel(program, "erode");

	if(device->isWritingTo3DTexturesSupported()) {
		OpenCLImageAccess::pointer input = segmentation->getOpenCLImageAccess(ACCESS_READ_WRITE, device);
		OpenCLImageAccess::pointer input2 = segmentation2->getOpenCLImageAccess(ACCESS_READ_WRITE, device);
		dilateKernel.setArg(0, *input->get3DImage());
		dilateKernel.setArg(1, *input2->get3DImage());
	} else {
		OpenCLImageAccess::pointer input = segmentation->getOpenCLImageAccess(ACCESS_READ_WRITE, device);
		OpenCLBufferAccess::pointer input2 = segmentation2->getOpenCLBufferAccess(ACCESS_READ_WRITE, device);
		dilateKernel.setArg(0, *input->get3DImage());
		dilateKernel.setArg(1, *input2->get());
	}

	device->getCommandQueue().enqueueNDRangeKernel(
			dilateKernel,
			cl::NullRange,
			cl::NDRange(width, height, depth),
			cl::NullRange
	);

	if(device->isWritingTo3DTexturesSupported()) {
		OpenCLImageAccess::pointer input = segmentation->getOpenCLImageAccess(ACCESS_READ_WRITE, device);
		OpenCLImageAccess::pointer input2 = segmentation2->getOpenCLImageAccess(ACCESS_READ_WRITE, device);
		erodeKernel.setArg(0, *input2->get3DImage());
		erodeKernel.setArg(1, *input->get3DImage());
	} else {
		OpenCLBufferAccess::pointer input = segmentation->getOpenCLBufferAccess(ACCESS_READ_WRITE, device);
		OpenCLImageAccess::pointer input2 = segmentation2->getOpenCLImageAccess(ACCESS_READ_WRITE, device);
		erodeKernel.setArg(0, *input2->get3DImage());
		erodeKernel.setArg(1, *input->get());
	}

	device->getCommandQueue().enqueueNDRangeKernel(
			erodeKernel,
			cl::NullRange,
			cl::NDRange(width, height, depth),
			cl::NullRange
	);

}

void AirwaySegmentation::execute() {
	Image::pointer image = getInputData<Image>();

	width = image->getWidth();
    height = image->getHeight();
    depth = image->getDepth();

	// Convert to signed HU if unsigned
	if(image->getDataType() == TYPE_UINT16) {
		image = convertToHU(image);
	}
	if(image->getDataType() != TYPE_INT16) {
		throw Exception("Input image to airway segmentation must be of data type INT16");
	}

	// Smooth image
	if(mSmoothingSigma > 0) {
		GaussianSmoothingFilter::pointer filter = GaussianSmoothingFilter::New();
        filter->setInputData(image);
        filter->setStandardDeviation(mSmoothingSigma);
        DataChannel::pointer port = filter->getOutputPort();
        filter->update();
        image = port->getNextFrame<Image>();
    }

	Segmentation::pointer segmentation = getOutputData<Segmentation>();
	
	if (mSeedPoints.size() > 0) {
		for (int i = 0; i < mSeedPoints.size(); ++i) {
			Vector3i seed = mSeedPoints[i];

			// Validate seed point
			if(seed.x() < 0 || seed.y() < 0 || seed.z() < 0 ||
					seed.x() >= image->getWidth() || seed.y() >= image->getHeight() || seed.z() >= image->getDepth()) {
				throw Exception("Seed point was not inside image in AirwaySegmentation");
			}
		}

		regionGrowing(image, segmentation, mSeedPoints);
	} else {
		Vector3i seed = findSeedVoxel(image);

		if(seed == Vector3i::Zero()) {
			throw Exception("No seed found.");
		}

		autoSeed = seed;

		reportInfo() << "Using seed point: " << seed.transpose() << reportEnd();

		regionGrowing(image, segmentation, std::vector<Vector3i>{seed});
	}

	// Do morphological closing to remove holes in segmentation
	morphologicalClosing(segmentation);
}

void AirwaySegmentation::addSeedPoint(int x, int y, int z) {
    addSeedPoint(Vector3i(x, y, z));
}

void AirwaySegmentation::addSeedPoint(Vector3i seed) {
	mSeedPoints.push_back(seed);
}

void AirwaySegmentation::setSmoothing(float sigma) {
	mSmoothingSigma = sigma;
}

void AirwaySegmentation::setVoxSpacing(Vector3f spacing) {
	voxSpacing = spacing;
	
	// distort icosphere to match voxel aspect ratio
	for (int i = 0; i < 21; i++) {
		icohalfVx[i] = icohalfMm[i];
		icohalfVx[i] = icohalfVx[i].cwiseQuotient(spacing);
	}

	mmPerVx = calcDistance(Vector3f(0.0, 0.0, 0.0), spacing);
}

}
