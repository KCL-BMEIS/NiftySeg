/**
 * @file _seg_PatchMatchResult.h
 * @author Ferran Prados
 * @date 16/06/2016
 *
 * Copyright (c) 2016, University College London. All rights reserved.
 * Centre for Medical Image Computing (CMIC)
 * See the LICENSE.txt file in the nifty_seg root folder
 *
 */
 
#ifndef _SEG_PATCHMATCHRESULT_H
#define _SEG_PATCHMATCHRESULT_H

#include "_seg_tools.h"

class PatchMatchResult
{
    private:
        float ANN_0;
        float ANN_1;
        long patchmatch_0;
        long patchmatch_1;
        int image_0;
        int image_1;

    public:
        PatchMatchResult();
        ~PatchMatchResult();
        PatchMatchResult(PatchMatchResult &);
        PatchMatchResult& operator = (PatchMatchResult &);
        bool operator < (PatchMatchResult &);
        bool operator > (PatchMatchResult &);
        bool operator >= (PatchMatchResult &);
        bool operator <= (PatchMatchResult &);
        bool operator == (PatchMatchResult &);
        void setANN(float);
        float getANN();
        void setImage(int);
        void setPatchMatch(long);
        int getImage();
        long getPatchMatch();
        void buffer();
};
#endif // _SEG_PATCHMATCHRESULT_H
