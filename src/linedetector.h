/*
 * LineDetector.h
 *
 *  Created on: 2012. 5. 23.
 *  Modified on: 2013. 8. 17.
 *
 *      Author: Jinhan Lee
 */

#ifndef LINEDETECTOR_H_
#define LINEDETECTOR_H_

#include "system.h"

DECLARE_int32(lt);
DECLARE_double(dt);

class LineDetector
{
public:
  LineDetector() { threshold_dist=(float)FLAGS_dt; threshold_length=(int)FLAGS_lt;
    init_label=0; imagewidth=720; imageheight=405;
  }
  ~LineDetector(){};
	template<class tType>
    void incidentPoint( tType * pt, Mat & l );
	void mergeLines(SEGMENT * Seg1, SEGMENT * Seg2, SEGMENT * SegMerged);
  bool getPointChain( const Mat & img, const Point pt, Point * chained_pt, int
      & direction, int step );
  double dist_point_line( const Mat & p, Mat & l );
	bool mergeSegments( SEGMENT * seg1, SEGMENT * seg2, SEGMENT * seg_merged );
	void extractSegments( vector<Point2i> * points, vector<SEGMENT> * segments );
	void lineDetection( Mat & src, vector<SEGMENT> * segments_all );
	void pointInboardTest(Mat & src, Point2i * pt);
  void getAngle(SEGMENT *seg);
	void additionalOperationsOnSegments(Mat & src, SEGMENT * seg);
	void drawArrow( Mat & mat, const SEGMENT * seg, Scalar bgr=Scalar(0,255,0));
	void set( Size sz, float th_dist=(float)FLAGS_dt, int th_length=FLAGS_lt);

  int init_label;

private:
	int imagewidth, imageheight, threshold_length;
  float threshold_dist;

};

#endif
