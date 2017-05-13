#include "acq/impl/cloudPrimitive.hpp"
#include <iostream>

namespace acq {

// ******** This class is inspired from "cloudManager" 

void CloudPrimitive::addPrimitive(Primitive* primitive) {
    _primitives.push_back(primitive);
} 

void CloudPrimitive::setPrimitive(Primitive* primitive, int index) {
    if (index >= _primitives.size()) {
        if (index != _primitives.size())
            std::cerr << "[CloudPrimitive::setPrimitive] "
                      << "Warning, creating " << index - _primitives.size()
                      << " empty primitive when inserting to index " << index
                      << ", current size is " << _primitives.size()
                      << "...why not use addPrimitive?\n";
        _primitives.resize(index + 1);
    }

    _primitives.at(index) = primitive;
} 

Primitive* CloudPrimitive::getPrimitive(int index) {
    if (index < _primitives.size())
        return _primitives.at(index);
    else {
        std::cerr << "Cannot return primitive with id " << index
                  << ", only have " << _primitives.size()
                  << " primitives ...returning empty primitive\n";
        throw new std::runtime_error("No such primitive");
    }
}

int CloudPrimitive::findBestScore() {
    int numberPrim = _primitives.size() ;
    double bestScore = 0 ;
    double thisScore ;
    int bestPrimIdx ;

    // go over all the primitives to find the best score
    for (int i=0; i<numberPrim ; i++) {
        // get back information for this primitive 
        Primitive* thisPrim = this->getPrimitive(i) ;
        thisScore = thisPrim->getScore() ;

        // compare with previous result 
        if (thisScore > bestScore) {
            bestPrimIdx = i ;
            bestScore = thisScore ;
        }
    }

    return bestPrimIdx ;
}

// delete the primitive at the position index
void CloudPrimitive::deletePrimitive(int index) {
        delete this->getPrimitive(index);
        _primitives.erase(_primitives.begin() + index);
}

// delete all the primitives and the vector 
void CloudPrimitive::clearAllPrimitives() {
    std::vector<Primitive*>::iterator thisPrimIt;

    for ( thisPrimIt = _primitives.begin(); thisPrimIt != _primitives.end(); ) {
        delete * thisPrimIt;  
        thisPrimIt = _primitives.erase(thisPrimIt);
    }
}



} //...ns acq
