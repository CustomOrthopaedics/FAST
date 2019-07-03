#pragma once

#include "SpatialDataObject.hpp"

// Forward declare
typedef struct _openslide openslide_t;

namespace fast {

class Image;


class FAST_EXPORT ImagePyramid : public SpatialDataObject {
    FAST_OBJECT(ImagePyramid)
    public:

        typedef struct Patch {
            std::shared_ptr<uchar[]> data;
            int width;
            int height;
            int offsetX;
            int offsetY;
        } Patch;

        typedef struct Level {
            int width;
            int height;
            int patches;
        } Level;

        void create(openslide_t* fileHandle, std::vector<Level> levels);
        Patch getPatch(std::string tile);
        Patch getPatch(int level, int patchX, int patchY);
        int getNrOfLevels();
        int getLevelWidth(int level);
        int getLevelHeight(int level);
        int getLevelPatches(int level);
        int getFullWidth();
        int getFullHeight();
        SharedPointer<Image> getLevelAsImage(int level);
        SharedPointer<Image> getPatchAsImage(int level, int offsetX, int offsetY, int width, int height);
        void free(ExecutionDevice::pointer device) override;
        void freeAll() override;
        ~ImagePyramid();
    private:
        std::vector<Level> m_levels;

        std::unordered_map<std::string, Patch> m_patchCache;

        uint64_t m_patchCacheMemoryUsage = 0;

        openslide_t* m_fileHandle;
};

}