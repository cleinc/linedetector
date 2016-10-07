/*
 * system.h
 *
 *  Created on: 2011. 6. 17.
 *  Modified on: 2013. 10. 1.
 *      Author: Jinhan Lee
*/

#ifndef SYSTEM_H_
#define SYSTEM_H_

#include <iostream>
#include <fstream>
#include <cstdio>
#include <stdio.h>
#include <string>
#include <limits>
#include <memory>

#include <math.h>
#include <vector>

#include <opencv2/opencv.hpp>
// #include <opencv2/features2d/features2d.hpp>
// #include <opencv2/nonfree/features2d.hpp>
// #include <opencv2/imgproc/imgproc.hpp>
// #include <opencv2/core/core.hpp>
// #include <opencv2/highgui/highgui.hpp>
// #include <opencv2/calib3d/calib3d.hpp>
// #include <opencv2/video/video.hpp>
// #include <opencv2/gpu/gpu.hpp>
// #include <opencv2/core/internal.hpp>

#include <gflags/gflags.h>

using namespace std;
using namespace cv;

struct SEGMENT {
    float x1, y1, x2, y2, angle;
    int label;
};

struct Point4f {
	float x1,y1,x2,y2;
};

#endif
