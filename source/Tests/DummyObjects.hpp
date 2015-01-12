#ifndef DUMMYOBJECTS_HPP_
#define DUMMYOBJECTS_HPP_

#include "ProcessObject.hpp"
#include "DataObject.hpp"
#include <boost/unordered_map.hpp>

namespace fast {

// Create a dummy class that extends the DataObject which is an abstract class
class DummyDataObject : public DataObject {
    FAST_OBJECT(DummyDataObject)
    public:
        bool hasBeenFreed(ExecutionDevice::pointer device) {
            if(mFreed.count(device) > 0) {
                return mFreed[device];
            }
            return false;
        };
    private:
        DummyDataObject() {};
        void free(ExecutionDevice::pointer device) {
            mFreed[device] = true;
        };
        void freeAll() {};
        boost::unordered_map<ExecutionDevice::pointer, bool> mFreed;
};


// Create a dummy class that extends the ProcessObject which is an abstract class
class DummyProcessObject : public ProcessObject {
    FAST_OBJECT(DummyProcessObject)
    public:
        void setIsModified() { mIsModified = true; };
        void setInput(DataObject::pointer data) {  };
        void setInput(DataObject::pointer data, bool required) { ProcessObject::setInputRequired(0, required); };
        void setInputRequired(uint number) { ProcessObject::setInputRequired(number, true); };
        bool hasExecuted() { return mHasExecuted; };
        void setHasExecuted(bool value) { mHasExecuted = value; };
        void updateDataTimestamp() { getOutputData<DummyDataObject>(0)->updateModifiedTimestamp(); };
    private:
        DummyProcessObject() : mHasExecuted(false) {};
        void execute() {
            std::cout << "executing" << std::endl;
            mHasExecuted = true;
            getOutputData<DummyDataObject>(0);
        };
        bool mHasExecuted;
};

};




#endif /* DUMMYOBJECTS_HPP_ */
