## This implementation has been imported as an OpenCV extra module (ximgproc)
https://docs.opencv.org/4.x/df/ded/group__ximgproc__fast__line__detector.html

# Straight line segment extractor

Simple but efficient/effective line segement detector.
This detector, used in works listed below, extracts line segments from images more effectively than classical Hough transform or LSD.
This detector also has a function that merges noisy-broken short line segments into one segment for more reliable detection.

Please cite one of these papers if you use this code in your research:

Lee, Jin Han, et al. "Place recognition using straight lines for vision-based SLAM." 2013 IEEE International Conference on Robotics and Automation (ICRA). IEEE, 2013.

Lee, Jin Han, et al. "Outdoor place recognition in urban environments using straight lines." 2014 IEEE International Conference on Robotics and Automation (ICRA). IEEE, 2014.

Zhang, Guoxuan, et al. "Building a 3-D Line-Based Map Using Stereo SLAM." IEEE Transactions on Robotics 31.6 (2015): 1364-1377.

# Dependency
Opencv 2.4.x and upper versions

Google Flags 2.1.0 ($ sudo apt-get install libgflags-dev)

# Usage

$ mkdir build

$ cd build

$ cmake ../

$ make

$ ./linedetection -i ../img/squre.jpg

You can see the simple usage of this detector in linedetection.cpp.

If you have any question, feel free to contact me jhlee@cle.vision
