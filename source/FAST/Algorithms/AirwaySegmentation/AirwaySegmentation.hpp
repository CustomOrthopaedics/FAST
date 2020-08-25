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
		// static int getIndex(Vector3i point, int height, int width);
		// static float findDistanceToWall(short *vol, Vector3f dir, Vector3i startPoint, int height, int width);
		// static bool canSpawnNewVoxel(short *vol, Vector3i point, int height, int width, float maxRadius);

		std::vector<Vector3i> mSeedPoints;
		float mSmoothingSigma = 0.5;
        bool mUseManualSeedPoint = false;
};

}


#endif
