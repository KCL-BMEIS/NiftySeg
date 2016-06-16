/**
 * @file _seg_PatchMatchResult.cpp
 * @author Ferran Prados
 * @date 16/06/2016
 *
 * Copyright (c) 2016, University College London. All rights reserved.
 * Centre for Medical Image Computing (CMIC)
 * See the LICENSE.txt file in the nifty_seg root folder
 *
 */
 
#include "_seg_PatchMatchResult.h"

PatchMatchResult::PatchMatchResult() {
    ANN_0=std::numeric_limits<float>::max();
    ANN_1=std::numeric_limits<float>::max();
    patchmatch_0=-1;
    patchmatch_1=-1;
    image_0=-1;
    image_1=-1;
}

PatchMatchResult::~PatchMatchResult(){
}

PatchMatchResult::PatchMatchResult(PatchMatchResult &copy) {
    setANN(copy.getANN());
    setPatchMatch(copy.getPatchMatch());
    setImage(copy.getImage());
    buffer();
}

PatchMatchResult& PatchMatchResult::operator = (PatchMatchResult &copy) {
    setANN(copy.getANN());
    setPatchMatch(copy.getPatchMatch());
    setImage(copy.getImage());
    buffer();

    return *this;
}

bool PatchMatchResult::operator < (PatchMatchResult &rhs) {
    bool smallerThan=getANN()<=rhs.getANN();
    return smallerThan;
}

bool PatchMatchResult::operator > (PatchMatchResult &rhs) {
    bool biggerThan=getANN()>rhs.getANN();
    return biggerThan;
}

bool PatchMatchResult::operator >= (PatchMatchResult &rhs) {
    bool biggerEqualThan=getANN()>=rhs.getANN();
    return biggerEqualThan;
}

bool PatchMatchResult::operator <= (PatchMatchResult &rhs) {
    bool smallerEqualThan=getANN()<=rhs.getANN();
    return smallerEqualThan;
}

bool PatchMatchResult::operator == (PatchMatchResult &rhs) {
    bool equal=getANN()==rhs.getANN();
    return equal;
}

void PatchMatchResult::setANN(float ann) {
    ANN_0=ann;
}

void PatchMatchResult::setImage(int img) {
    image_0=img;
}

void PatchMatchResult::setPatchMatch(long voxel) {
    patchmatch_0=voxel;
}

void PatchMatchResult::buffer(){
    ANN_0=ANN_1;
    image_0=image_1;
    patchmatch_0=patchmatch_1;
}

float PatchMatchResult::getANN() {
    return ANN_0;
}

int PatchMatchResult::getImage() {
    return image_0;
}

long PatchMatchResult::getPatchMatch() {
    return patchmatch_0;
}



