//
//  raytracer.cpp
//  raytrace2
//
//  Created by Björn Dagerman on 2012-05-22.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
#ifndef _TRACER_
#define _TRACER_


#include "Vec3.h"
#include "Image.h"
#include "Vector"
#include "xmmintrin.h"
using namespace std;

#define NR_OF_REGIONS 200
#define NR_OF_THREADS 8
#define MAX_DEPTH 10

#define SHUFFLE_REGIONS 0

class Scene;
class Tracer;

struct Region{
    Region(int xs, int xe,int ys,int ye) : xstart(xs), ystart(ys), xend(xe), yend(ye) {}
    int xstart, xend, ystart, yend;  
};

class Ray{
public:
	Ray() : o(Vec3(0,0,0)), d(Vec3(0,0,0)) {};
	Ray(Vec3& o, Vec3& d) :o(o), d(d) {};
	Vec3 o;
	Vec3 d;
};

class Tracer{
    static Tracer *instance;
public:
    static Tracer *getInstance(){ if (!instance) instance = new Tracer(); return instance; }
    
	Tracer();
	void setTarget(int *width, int *height,int *w_width, int*w_heightfloat, void*rgb_b, float *sr);

	void trace(Ray &ray, int depth, float &dist, int&nrOfIntersectionTests, float &r, float&g, float&b, bool useSSE, float coef);
    void trace_sse(__m128 &rox, __m128 &roy,__m128 &roz,__m128 &rdx,__m128 &rdy,__m128 &rdz, int depth, float& aDist, int&nrOfIntersectionTests, __m128 &r, __m128&g, __m128&b, bool useSSE, float coef);
	void init();
    void render(Vec3 &cam,int&nrOfIntersectionTests, int xstart, int xend, int ystart, int yend);
    void render(Vec3&cam);
    void render_sse(int ystart,int yend, int xstart, int xend, int maxdepth);
    
    bool useSuperSampling;
    void shuffle();
    void resetRegions();
    
    void*rgb_buffer;
	float xs, ys, xe, ye, dx, dy, xStart, yStart;
    float *samplerate;
	Scene* scene;
	int* width;
    int* height;
    int *window_height;
    int* window_width;
    bool renderChecker;
    int maxDepth;
    
    bool *showBenchmarking;
    bool shouldRender;

    Vec3 campos;
    
    pthread_t threads[NR_OF_THREADS];
};


#endif