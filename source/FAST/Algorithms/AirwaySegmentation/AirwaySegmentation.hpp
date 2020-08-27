#ifndef AIRWAY_SEGMENTATION_HPP_
#define AIRWAY_SEGMENTATION_HPP_

#include "FAST/Algorithms/SegmentationAlgorithm.hpp"

namespace fast {

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
		float findDistanceToWall(short *vol, Vector3f dir, Vector3i startPoint);
		std::vector<float> getVoxelData(short *vol, Vector3i point);
		int grow(uchar* segmentation, std::vector<Vector3i> neighbors, std::vector<Vector3i>& voxels, short* data, float threshold, int width, int height, int depth, float previousVolume, float volumeIncreaseLimit, int volumeMinimum);
		void regionGrowing(Image::pointer volume, Segmentation::pointer segmentation, const std::vector<Vector3i> seeds);

		std::vector<Vector3i> mSeedPoints;
		Vector3f voxSpacing = Vector3f(1, 1, 1);
		float mSmoothingSigma = 0.5;
        bool mUseManualSeedPoint = false;
		int height;
		int width;
		int depth;
};

}


#endif
