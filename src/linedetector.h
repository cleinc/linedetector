#ifndef LINEDETECTOR_H_
#define LINEDETECTOR_H_

#include "system.h"

DECLARE_double(dt);
DECLARE_int32(lt);

class LineDetector
{
public:
  LineDetector() {
    threshold_dist=FLAGS_dt;
    threshold_length=FLAGS_lt;
    init_label=0;
  }
  ~LineDetector(){};
	template<class tType>
    void incidentPoint( tType * pt, Mat & l );
	void mergeLines(SEGMENT * Seg1, SEGMENT * Seg2, SEGMENT * SegMerged);
  bool getPointChain( const Mat & img, const Point pt, Point * chained_pt, int
      & direction, int step );
  // bool getPointChain( const Mat & img, const Mat & scales, const Point pt, Point * chained_pt, int
  // & direction, int step );
  double distPointLine( const Mat & p, Mat & l );
	bool mergeSegments( SEGMENT * seg1, SEGMENT * seg2, SEGMENT * seg_merged );
	void extractSegments( vector<Point2i> * points, vector<SEGMENT> * segments );
	void lineDetection( Mat & src, vector<SEGMENT> & segments_all, bool merge = true );
	void pointInboardTest(Mat & src, Point2i * pt);
  void getAngle(SEGMENT *seg);
	void additionalOperationsOnSegments(Mat & src, SEGMENT * seg);
  void drawArrow( Mat& mat, const SEGMENT* seg, Scalar bgr=Scalar(0,255,0),
      int thickness=1, bool directed=true);

private:
	int init_label, imagewidth, imageheight, threshold_length;
  float threshold_dist;
};

#endif
