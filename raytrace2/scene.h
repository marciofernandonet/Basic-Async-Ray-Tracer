//
//  scene.h
//  raytrace2
//
//  Created by Björn Dagerman on 2012-05-22.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.

#ifndef scn
#define scn

#include "raytracer.h"
#include "xmmintrin.h"

#define HIT		 1		
#define MISS	 0		
#define INSIDE	-1		

class Material{
public:
	Material() : color(Color(0,0,0)), refl(0), diff(0), spec(0) {}
    
	void setColor(Color& clr) { color = clr; }
	void setDiffuse(float d) { diff = d; }
	void setReflection(float r) { refl = r; }
    void setSpecular(float s) { spec = s;}

	float getSpecular() { return spec; }
	float getDiffuse() { return diff; }
	float getReflection() { return refl; }
    Color getColor() { return color; }
private:
	Color color;
	float refl, diff, spec;
};

class Primitive{
public:
	enum{
		SPHERE = 1,
		PLANE,
        CHECKERBOARD
	};
	Primitive() : isLight(false) {};
	Material* getMaterial() { return &material; }
	virtual int getType() = 0;

	virtual int intersect(Ray &ray, float &t) = 0;
    virtual void intersect_sse(__m128 &rox, __m128 &roy, __m128 &roz,__m128 &rdx, __m128 &rdy, __m128 &rdz, __m128 &t, __m128 &hit_res) = 0;

	virtual Vec3 getNormal(Vec3 &p) = 0;
    bool isLight;
protected:
	Material material;
};


struct Sphere : public Primitive{
	int getType() { return SPHERE; }
	Sphere(Vec3 &c, float r) : c(c), r(r) {};
	Vec3& getCenter() { return c; }
    bool hitSphere(Ray &r, float &t);
    
    int intersect(Ray &ray, float &t);
    void intersect_sse(__m128 &rox, __m128 &roy, __m128 &roz,__m128 &rdx, __m128 &rdy, __m128 &rdz, __m128 &t, __m128 &hit_res);
	Vec3 getNormal(Vec3 &p) { 
        Vec3 n = p - c;
        float dot = DOT(n,n);
        dot = InvSqrt(dot);
        n *= dot;
        return n;
    }
private:
	Vec3 c;
    float r;
};

struct PlanePrim : public Primitive{
	int getType() { return PLANE; }
	PlanePrim(Vec3 &n, float d) : plane(Plane(n, d)) {};
	Vec3& getNormal() { return plane.N; }
	int intersect(Ray& ray, float &t);
    void intersect_sse(__m128 &rox, __m128 &roy, __m128 &roz,__m128 &rdx, __m128 &rdy, __m128 &rdz, __m128 &t, __m128 &hit_res);

	Vec3 getNormal( Vec3& a_Pos );
	Plane plane;
};

struct CheckerBoard : public PlanePrim{
    int getType() { return CHECKERBOARD; }
    CheckerBoard(Vec3 &n, float d) : PlanePrim(n, d) {};
    Vec3 getNormal(Vec3 &p);
};

struct Scene{
	Scene() : nrOfPrimitives(0), primitives(0) { init(); };
	~Scene();
	void init();
	Primitive* getPrimitive(int i) { return primitives[i]; }
    int nrOfPrimitives;
private:
    Primitive** primitives;
};

#endif