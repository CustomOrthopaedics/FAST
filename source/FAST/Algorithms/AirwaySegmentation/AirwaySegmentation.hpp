#ifndef AIRWAY_SEGMENTATION_HPP_
#define AIRWAY_SEGMENTATION_HPP_

#include "FAST/Algorithms/SegmentationAlgorithm.hpp"

namespace fast {

struct Voxel {
	Vector3i point;
	float centricity;
	float minRadius;
	float meanRadii;

	Voxel(Vector3i p, float cent, float minRad, float meanRad) {
		point = p;
		centricity = cent;
		minRadius = minRad;
		meanRadii = meanRad;
	}

	// https://stackoverflow.com/a/5712235
	// overload operator so Voxel can be used in priority queue
	bool operator <(const struct Voxel& other) const {
		return centricity < other.centricity;
	}
};

class Image;
class Segmentation;

class FAST_EXPORT  AirwaySegmentation : public SegmentationAlgorithm {
	FAST_OBJECT(AirwaySegmentation)
	public:
	    void addSeedPoint(int x, int y, int z);
		void addSeedPoint(Vector3i seed);
		void setVoxSpacing(Vector3f);
		/**
		 * Set the sigma value of the gaussian smoothing performed before segmentation.
		 * Default is 0.5. A higher value can be used for low dose CT.
		 * @param sigma
		 */
		void setSmoothing(float sigma);
		Vector3i autoSeed;
	private:
		AirwaySegmentation();
		void execute();
		static Vector3i findSeedVoxel(SharedPointer<Image> volume);
		SharedPointer<Image> convertToHU(SharedPointer<Image> image);
		void morphologicalClosing(SharedPointer<Segmentation> segmentation);
		int interp3D(short *vol, Vector3f point);
		int getIndex(Vector3i);
		int getIndex(int x, int y, int z);
		float findDistanceToWall(short *vol, Vector3f dir, Vector3i startPoint);
		Voxel getVoxelData(short *vol, Vector3i point);
		int grow(Vector3i seed, uchar* mask, std::vector<Vector3i> neighbors, short* data, float threshold);
		void regionGrowing(Image::pointer volume, Segmentation::pointer segmentation, const std::vector<Vector3i> seeds);

		std::vector<Vector3i> mSeedPoints;
		Vector3f voxSpacing = Vector3f(1, 1, 1);
		float mSmoothingSigma = 0.5;
        bool mUseManualSeedPoint = false;
		int height;
		int width;
		int depth;
		Vector3f icohalfVx[21];
		float mmPerVx = 1.0;
};

}


#endif
