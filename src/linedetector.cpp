#include "linedetector.h"

void LineDetector::mergeLines(SEGMENT * seg1, SEGMENT * seg2, SEGMENT * seg_merged)
{
	double xg = 0.0, yg = 0.0;
	double delta1x = 0.0, delta1y = 0.0, delta2x = 0.0, delta2y = 0.0;
	float ax = 0, bx = 0, cx = 0, dx = 0;
	float ay = 0, by = 0, cy = 0, dy = 0;
	double li = 0.0, lj = 0.0;
	double thi = 0.0, thj = 0.0, thr = 0.0;
	double axg = 0.0, bxg = 0.0, cxg = 0.0, dxg = 0.0, delta1xg = 0.0, delta2xg = 0.0;

	ax = seg1->x1;
	ay = seg1->y1;

	bx = seg1->x2;
	by = seg1->y2;
	cx = seg2->x1;
	cy = seg2->y1;

	dx = seg2->x2;
	dy = seg2->y2;

	float dlix = (bx - ax);
	float dliy = (by - ay);
	float dljx = (dx - cx);
	float dljy = (dy - cy);

	li = sqrt((double) (dlix * dlix) + (double) (dliy * dliy));
	lj = sqrt((double) (dljx * dljx) + (double) (dljy * dljy));

	xg = (li * (double) (ax + bx) + lj * (double) (cx + dx))
					/ (double) (2.0 * (li + lj));
	yg = (li * (double) (ay + by) + lj * (double) (cy + dy))
					/ (double) (2.0 * (li + lj));

	if(dlix == 0.0f) thi = CV_PI / 2.0;
	else thi = atan(dliy / dlix);

	if(dljx == 0.0f) thj = CV_PI / 2.0;
	else thj = atan(dljy / dljx);

	if (fabs(thi - thj) <= CV_PI / 2.0)
	{
		thr = (li * thi + lj * thj) / (li + lj);
	}
	else
	{
		double tmp = thj - CV_PI * (thj / fabs(thj));
		thr = li * thi + lj * tmp;
		thr /= (li + lj);
	}

	axg = ((double) ay - yg) * sin(thr) + ((double) ax - xg) * cos(thr);

	bxg = ((double) by - yg) * sin(thr) + ((double) bx - xg) * cos(thr);

	cxg = ((double) cy - yg) * sin(thr) + ((double) cx - xg) * cos(thr);

	dxg = ((double) dy - yg) * sin(thr) + ((double) dx - xg) * cos(thr);

	delta1xg = min(axg,min(bxg,min(cxg,dxg)));
	delta2xg = max(axg,max(bxg,max(cxg,dxg)));

	delta1x = delta1xg * cos(thr) + xg;
	delta1y = delta1xg * sin(thr) + yg;
	delta2x = delta2xg * cos(thr) + xg;
	delta2y = delta2xg * sin(thr) + yg;

	seg_merged->x1 = (float) delta1x;
	seg_merged->y1 = (float) delta1y;
	seg_merged->x2 = (float) delta2x;
	seg_merged->y2 = (float) delta2y;
}

double LineDetector::distPointLine( const Mat & p, Mat & l )
{
	double x, y, w;

	x = l.at<double>(0,0);
	y = l.at<double>(1,0);

	w = sqrt(x*x+y*y);

	l.at<double>(0,0) = x  / w;
	l.at<double>(1,0) = y  / w;
	l.at<double>(2,0) = l.at<double>(2,0)  / w;

	return l.dot(p);
}

bool LineDetector::mergeSegments( SEGMENT * seg1, SEGMENT * seg2, SEGMENT * seg_merged )
{
	double o[] = { 0.0, 0.0, 1.0 };
	double a[] = { 0.0, 0.0, 1.0 };
	double b[] = { 0.0, 0.0, 1.0 };
	double c[3];

	o[0] = ( seg2->x1 + seg2->x2 ) / 2.0;
	o[1] = ( seg2->y1 + seg2->y2 ) / 2.0;

	a[0] = seg1->x1;
	a[1] = seg1->y1;
	b[0] = seg1->x2;
	b[1] = seg1->y2;

	Mat ori = Mat(3, 1, CV_64FC1, o).clone();
	Mat p1 = Mat(3, 1, CV_64FC1, a).clone();
	Mat p2 = Mat(3, 1, CV_64FC1, b).clone();
	Mat l1 = Mat(3, 1, CV_64FC1, c).clone();

	l1 = p1.cross(p2);

	Point2f seg1mid, seg2mid;
	seg1mid.x = (seg1->x1 + seg1->x2) /2.0f;
	seg1mid.y = (seg1->y1 + seg1->y2) /2.0f;
	seg2mid.x = (seg2->x1 + seg2->x2) /2.0f;
	seg2mid.y = (seg2->y1 + seg2->y2) /2.0f;

	double seg1len, seg2len;
	seg1len = sqrt((seg1->x1 - seg1->x2)*(seg1->x1 - seg1->x2)+(seg1->y1 - seg1->y2)*(seg1->y1 - seg1->y2));
	seg2len = sqrt((seg2->x1 - seg2->x2)*(seg2->x1 - seg2->x2)+(seg2->y1 - seg2->y2)*(seg2->y1 - seg2->y2));

	double middist = sqrt((seg1mid.x - seg2mid.x)*(seg1mid.x - seg2mid.x) + (seg1mid.y - seg2mid.y)*(seg1mid.y - seg2mid.y));

	float angdiff = seg1->angle - seg2->angle;
	angdiff = fabs(angdiff);

	double dist = distPointLine( ori, l1 );

	if ( fabs( dist ) <= threshold_dist * 2.0
			&& middist <= seg1len / 2.0 + seg2len / 2.0 + 20.0
			&& angdiff <= CV_PI / 180.0f * 5.0f) {
		mergeLines(seg1, seg2, seg2);
		return true;
	} else {
		return false;
	}
}

template<class tType>
void LineDetector::incidentPoint( tType * pt, Mat & l )
{
	double a[] = { (double)pt->x, (double)pt->y, 1.0 };
	double b[] = { l.at<double>(0,0), l.at<double>(1,0), 0.0 };
	double c[3];

	Mat xk = Mat(3, 1, CV_64FC1, a).clone();
	Mat lh = Mat(3, 1, CV_64FC1, b).clone();
	Mat lk = Mat(3, 1, CV_64FC1, c).clone();

	lk = xk.cross(lh);
	xk = lk.cross(l);

	double s = 1.0 / xk.at<double>(2,0);
	xk.convertTo(xk, -1, s);

	pt->x = (float)xk.at<double>(0,0) < 0.0f ? 0.0f : (float)xk.at<double>(0,0)
			>= (imagewidth - 1.0f) ? (imagewidth - 1.0f) : (float)xk.at<double>(0,0);
	pt->y = (float)xk.at<double>(1,0) < 0.0f ? 0.0f : (float)xk.at<double>(1,0)
			>= (imageheight - 1.0f) ? (imageheight - 1.0f) : (float)xk.at<double>(1,0);

}

void LineDetector::extractSegments( vector<Point2i> * points, vector<SEGMENT> * segments )
{
	bool is_line;

	int i, j;
	SEGMENT seg;
	Point2i ps, pe, pt;

	vector<Point2i> l_points;

	int total = points->size();

	for ( i = 0; i + threshold_length < total; i++ ) {
		ps = points->at(i);
		pe = points->at(i + threshold_length);

		double a[] = { (double)ps.x, (double)ps.y, 1 };
		double b[] = { (double)pe.x, (double)pe.y, 1 };
		double c[3], d[3];

		Mat p1 = Mat(3, 1, CV_64FC1, a).clone();
		Mat p2 = Mat(3, 1, CV_64FC1, b).clone();
		Mat p = Mat(3, 1, CV_64FC1, c).clone();
		Mat l = Mat(3, 1, CV_64FC1, d).clone();
		l = p1.cross(p2);

		is_line = true;

		l_points.clear();
		l_points.push_back(ps);

		for ( j = 1; j < threshold_length; j++ ) {
			pt.x = points->at(i+j).x;
			pt.y = points->at(i+j).y;

			p.at<double>(0,0) = (double)pt.x;
			p.at<double>(1,0) = (double)pt.y;
			p.at<double>(2,0) = 1.0;

			double dist = distPointLine( p, l );

			if ( fabs( dist ) > threshold_dist ) {
				is_line = false;
				break;
			}
			l_points.push_back(pt);
		}

		// Line check fail, test next point
		if ( is_line == false )
			continue;

		l_points.push_back(pe);

		Vec4f line;
		fitLine( Mat(l_points), line, CV_DIST_L2, 0, 0.01, 0.01);
		a[0] = line[2];
		a[1] = line[3];
		b[0] = line[2] + line[0];
		b[1] = line[3] + line[1];

		p1 = Mat(3, 1, CV_64FC1, a).clone();
		p2 = Mat(3, 1, CV_64FC1, b).clone();

		l = p1.cross(p2);

		incidentPoint( &ps, l );

		// Extending line
		for ( j = threshold_length + 1; i + j < total; j++ ) {
			pt.x = points->at(i+j).x;
			pt.y = points->at(i+j).y;

			p.at<double>(0,0) = (double)pt.x;
			p.at<double>(1,0) = (double)pt.y;
			p.at<double>(2,0) = 1.0;

			double dist = distPointLine( p, l );

			if ( fabs( dist ) > threshold_dist ) {
				j--;
				break;
			}

			pe = pt;
			l_points.push_back(pt);
		}
		fitLine( Mat(l_points), line, CV_DIST_L2, 0, 0.01, 0.01);
		a[0] = line[2];
		a[1] = line[3];
		b[0] = line[2] + line[0];
		b[1] = line[3] + line[1];

		p1 = Mat(3, 1, CV_64FC1, a).clone();
		p2 = Mat(3, 1, CV_64FC1, b).clone();

		l = p1.cross(p2);

		Point2f e1, e2;
		e1.x = ps.x;
		e1.y = ps.y;
		e2.x = pe.x;
		e2.y = pe.y;

		incidentPoint( &e1, l );
		incidentPoint( &e2, l );
		seg.x1 = e1.x;
		seg.y1 = e1.y;
		seg.x2 = e2.x;
		seg.y2 = e2.y;

		segments->push_back(seg);
		i = i + j;
	}
}

void LineDetector::pointInboardTest(Mat & src, Point2i * pt)
{
	pt->x = pt->x <= 5.0f ? 5.0f : pt->x >= src.cols - 5.0f ? src.cols - 5.0f : pt->x;
	pt->y = pt->y <= 5.0f ? 5.0f : pt->y >= src.rows - 5.0f ? src.rows - 5.0f : pt->y;
}

bool LineDetector::getPointChain( const Mat & img, const Point pt, Point * chained_pt,
    int & direction, int step )
{
	int ri, ci;
	int indices[8][2]={ {1,1}, {1,0}, {1,-1}, {0,-1}, {-1,-1},{-1,0}, {-1,1}, {0,1} };

	for ( int i = 0; i < 8; i++ ) {
		ci = pt.x + indices[i][1];
		ri = pt.y + indices[i][0];

		if ( ri < 0 || ri == img.rows || ci < 0 || ci == img.cols )
			continue;

		if ( img.at<unsigned char>(ri, ci) == 0 )
			continue;

		if(step == 0) {
			chained_pt->x = ci;
			chained_pt->y = ri;
			direction = i;
			return true;
		} else {
      if(abs(i-direction) <= 2 || abs(i-direction) >= 6)
      {
				chained_pt->x = ci;
				chained_pt->y = ri;
				direction = i;
				return true;
			} else
				continue;
		}
	}
	return false;
}

void LineDetector::lineDetection( Mat & src, vector<SEGMENT> & segments_all, bool merge )
{
	int r, c;
	imageheight=src.rows; imagewidth=src.cols;

	vector<Point2i> points;
	vector<SEGMENT> segments, segments_tmp, segments_tmp2;
	Mat canny = src.clone();
	Canny(src, canny, 50, 50, 3);

	for(int i=0; i<src.rows;i++) {
		for(int j=0; j<src.cols;j++) {
			if( i < 5 || i > src.rows-5 || j < 5 || j > src.cols - 5)
				canny.at<unsigned char>(i,j) = 0;
		}
	}

	SEGMENT seg, seg1, seg2;

	for ( r = 0; r < imageheight; r++ ) {
		for ( c = 0; c < imagewidth; c++ ) {
			// Find seeds - skip for non-seeds
			if ( canny.at<unsigned char>(r,c) == 0 )
				continue;

			// Found seeds
			Point2i pt;
			pt.x = c;
			pt.y = r;

			points.push_back(pt);
			canny.at<unsigned char>(pt.y, pt.x) = 0;

			int direction = 0;
			int step = 0;
			while (getPointChain( canny, pt, &pt, direction, step)) {
				points.push_back(pt);
				step++;
				canny.at<unsigned char>(pt.y, pt.x) = 0;
			}

			if ( points.size() < (unsigned int)threshold_length + 1 ) {
				points.clear();
				continue;
			}

			extractSegments( &points, &segments );

			if ( segments.size() == 0 ) {
				points.clear();
				continue;
			}
			for ( int i = 0; i < (int)segments.size(); i++ ) {
				seg = segments.at(i);
				float length = sqrt((seg.x1 - seg.x2)*(seg.x1 - seg.x2) + (seg.y1 - seg.y2)*(seg.y1 - seg.y2));
				if(length < threshold_length) continue;
				if( (seg.x1 <= 5.0f && seg.x2 <= 5.0f)
						|| (seg.y1 <= 5.0f && seg.y2 <= 5.0f)
						|| (seg.x1 >= imagewidth - 5.0f && seg.x2 >= imagewidth - 5.0f)
						|| (seg.y1 >= imageheight - 5.0f && seg.y2 >= imageheight - 5.0f) )
					continue;

				additionalOperationsOnSegments(src, &seg);
        if(!merge) {
          segments_all.push_back(seg);
        }
        segments_tmp.push_back(seg);
			}
			points.clear();
			segments.clear();
		}
	}
  if(!merge)
    return;

  bool is_merged = false;
  int ith = segments_tmp.size() - 1;
  int jth = ith - 1;
  while(true)
  {
    seg1 = segments_tmp[ith];
    seg2 = segments_tmp[jth];
		is_merged = mergeSegments(&seg1, &seg2, &seg2);
    if(is_merged == true)
    {
      additionalOperationsOnSegments(src, &seg2);
      vector<SEGMENT>::iterator it = segments_tmp.begin() + ith;
      *it = seg2;
      segments_tmp.erase(segments_tmp.begin()+jth);
      ith--;
      jth = ith - 1;
    }
    else
    {
      jth--;
    }
    if(jth < 0) {
      ith--;
      jth = ith - 1;
    }
    if(ith == 1 && jth == 0)
      break;
  }
  segments_all = segments_tmp;
}

void LineDetector::getAngle(SEGMENT *seg)
{
	float dx = (float)(seg->x2 - seg->x1);
	float dy = (float)(seg->y2 - seg->y1);
	double ang=0.0;

	if(dx == 0.0f) {
		if(dy > 0)
			ang = CV_PI / 2.0;
		else
			ang = -1.0 * CV_PI / 2.0;
	}
	else if(dy == 0.0f) {
		if(dx > 0)
			ang = 0.0;
		else
			ang = CV_PI;
	}
	else if(dx < 0.0f && dy > 0.0f)
		ang = CV_PI + atan( dy/dx );
	else if(dx > 0.0f && dy < 0.0f)
		ang = 2*CV_PI + atan( dy/dx );
	else if(dx < 0.0f && dy < 0.0f)
		ang = CV_PI + atan( dy/dx );
	else
		ang = atan( dy/dx );

	if(ang > 2.0 * CV_PI)
		ang -= 2.0 * CV_PI;
	seg->angle = (float)ang;
}

void LineDetector::additionalOperationsOnSegments(Mat & src, SEGMENT * seg)
{
	if(seg->x1 == 0.0f && seg->x2 == 0.0f && seg->y1 == 0.0f && seg->y2 == 0.0f)
		return;

	getAngle(seg);
	double ang = (double)seg->angle;

	Point2f start = Point2f(seg->x1, seg->y1);
	Point2f end = Point2f(seg->x2, seg->y2);

	double dx = 0.0, dy = 0.0;
	dx = (double) end.x - (double) start.x;
	dy = (double) end.y - (double) start.y;

	int num_points = 10;
	Point2f *points = new Point2f[num_points];

	points[0] = start;
	points[num_points - 1] = end;
	for (int i = 0; i < num_points; i++) {
		if (i == 0 || i == num_points - 1)
			continue;
		points[i].x = points[0].x + (dx / double(num_points - 1) * (double) i);
		points[i].y = points[0].y + (dy / double(num_points - 1) * (double) i);
	}

	Point2i *points_right = new Point2i[num_points];
	Point2i *points_left = new Point2i[num_points];
	double gap = 1.0;

	for(int i = 0; i < num_points; i++) {
		points_right[i].x = cvRound(points[i].x + gap*cos(90.0 * CV_PI / 180.0 + ang));
		points_right[i].y = cvRound(points[i].y + gap*sin(90.0 * CV_PI / 180.0 + ang));
		points_left[i].x = cvRound(points[i].x - gap*cos(90.0 * CV_PI / 180.0 + ang));
		points_left[i].y = cvRound(points[i].y - gap*sin(90.0 * CV_PI / 180.0 + ang));
		pointInboardTest(src, &points_right[i]);
		pointInboardTest(src, &points_left[i]);
	}

	int iR = 0, iL = 0;
	for(int i = 0; i < num_points; i++) { 
		iR += src.at<unsigned char>(points_right[i].y, points_right[i].x);
		iL += src.at<unsigned char>(points_left[i].y, points_left[i].x);
	}

	if(iR > iL)
	{
    std::swap(seg->x1, seg->x2);
    std::swap(seg->y1, seg->y2);
		ang = ang + CV_PI;
		if(ang >= 2.0*CV_PI)
			ang = ang - 2.0 * CV_PI;
		seg->angle = (float)ang;
	}

	delete[] points; delete[] points_right; delete[] points_left;

	seg->label = init_label++;
	return;
}

void LineDetector::drawArrow( Mat & mat, const SEGMENT * seg, Scalar bgr, int thickness, bool directed)
{
	Point2i p1;

	double gap = 10.0;
	double ang = (double)seg->angle;
	double arrow_angle = 30.0;

	p1.x = round(seg->x2 - gap*cos(arrow_angle * CV_PI / 180.0 + ang));
	p1.y = round(seg->y2 - gap*sin(arrow_angle * CV_PI / 180.0 + ang));
	pointInboardTest(mat, &p1);

  line(mat, Point(round(seg->x1), round(seg->y1)), Point(round(seg->x2),
        round(seg->y2)), bgr, thickness, 1);
  if(directed)
    line(mat, Point(round(seg->x2), round(seg->y2)), p1, bgr, thickness, 1);
}
