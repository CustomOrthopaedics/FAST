#ifndef AIRWAY_SEGMENTATION_HPP_
#define AIRWAY_SEGMENTATION_HPP_

#include "FAST/Algorithms/SegmentationAlgorithm.hpp"

namespace fast {
struct VoxelRay {
	Vector3f direction;
	float length;
	int startHU;
	int endHU;

	VoxelRay(Vector3f dir, float len, int sHU, int eHU) {
		direction = dir;
		length = len;
		startHU = sHU;
		endHU = eHU;
	}
};

struct Voxel {
	Vector3i point;
	float centricity;
	float minRadius;
	float meanRadii;
	int pathVoxelIdx;

	// debugging code
	std::vector<VoxelRay> rays;
	int maskIdx;

	Voxel(Vector3i p, float cent, float minRad, float meanRad) {
		point = p;
		centricity = cent;
		minRadius = minRad;
		meanRadii = meanRad;
	}

	void addRay(VoxelRay ray) {
		rays.push_back(ray);
	}

	// https://stackoverflow.com/a/5712235
	// overload operator so Voxel can be used in priority queue
	bool operator <(const struct Voxel& other) const {
		return centricity < other.centricity;
	}
};

struct SegAlgConfig {
	float maxRadiusIncrease = 2.3;
	float maxAirwayDensity = -600.0;
	int minNeighbors = 9;
	float branchEndMaxRadius = 1.3;
	float mSmoothingSigma = 0.5;
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
		std::vector<Voxel> maskVoxels;
		void setSensitivity(int sensitivity);
		void setDebug(bool);
		void setBB(Vector3i, Vector3i);
		SegAlgConfig getAlgConfig();        
	private:
		AirwaySegmentation();
		void execute();
		static Vector3i findSeedVoxel(SharedPointer<Image> volume);
		SharedPointer<Image> convertToHU(SharedPointer<Image> image);
		void morphologicalClosing(SharedPointer<Segmentation> segmentation);
		int interp3D(short *vol, Vector3f point);
		int getIndex(Vector3i);
		int getIndex(int x, int y, int z);
		VoxelRay findDistanceToWall(short *vol, Vector3f dir, Vector3i startPoint);
		Voxel getVoxelData(short *vol, Vector3i point);
		int grow(Vector3i seed, uchar* mask, std::vector<Vector3i> neighbors, short* data, float threshold, int& maskIdx);
		void regionGrowing(Image::pointer volume, Segmentation::pointer segmentation, const std::vector<Vector3i> seeds);
		int numNeighbors(uchar *, Vector3i);

		std::vector<Vector3i> mSeedPoints;
		Vector3f voxSpacing = Vector3f(1, 1, 1);
    bool mUseManualSeedPoint = false;
		int height;
		int width;
		int depth;
		Vector3f icohalfVx[21];
		float mmPerVx = 1.0;

		bool debug = false;

		// alg constants
		float deltaW = 200.0;
		float dr = 0.5;
		float rMax = 20.0;
		float minCent = 0.0;
		int maxVoxelVal = 1000;
		float branchEndMinCentricity = 0.5;

		// Voxels identified to be leakage by "maxRadiusIncrease" param must have a centricity above this value
		// in order to be kept in the mask. This is to reduce leakage at the ends of branches and make the mask "cleaner".
		float minLeakageCentricity = 0.5;

		// path length leakage detection
		float maxPathRadiusIncrease = 0.8;
		int pathLengthMinVoxels = 20;

		SegAlgConfig algConfig;

		// bounding box in voxel coordinates
		Vector3i bbMin = Vector3i(0, 0, 0);
		Vector3i bbMax = Vector3i(100000, 100000, 100000);
};

}


#endif
