#include "Vec3.h"
#include "scene.h"
#include "raytracer.h"

inline void PlanePrim::intersect_sse(__m128 &rox, __m128 &roy, __m128 &roz,__m128 &rdx, __m128 &rdy, __m128 &rdz, __m128 &t, __m128 &hit_res){
    __m128 nx = _mm_set1_ps(plane.N.x);
    __m128 ny = _mm_set1_ps(plane.N.y);
    __m128 nz = _mm_set1_ps(plane.N.z);
    
    __m128 D = dot_sse(nx, ny, nz, rdx, rdy, rdz);
    
    __m128 dist = -_mm_div_ps(_mm_add_ps(dot_sse(nx, ny, nz, rox, roy, roz),  _mm_set1_ps(plane.D)) , D);
    
    hit_res = _mm_set1_epi32(MISS);
    
    for (int i=0;i<4;++i){
        if (t[i] < dist[i]){
            t[i] = dist[i];
            hit_res[i] = HIT;
        }
    }
}

inline void Sphere::intersect_sse(__m128 &rox, __m128 &roy, __m128 &roz,__m128 &rdx, __m128 &rdy, __m128 &rdz, __m128 &t, __m128 &hit_res){
    __m128 dist_x = _mm_sub_ps(_mm_set1_ps(c.x), rox);
    __m128 dist_y = _mm_sub_ps(_mm_set1_ps(c.y), roy);
    __m128 dist_z = _mm_sub_ps(_mm_set1_ps(c.z), roz);
    
    __m128 B = dot_sse(rdx,rdy,rdz, dist_x, dist_y, dist_z);
    
	__m128 D = _mm_add_ps(_mm_sub_ps(_mm_mul_ps(B,B), dot_sse(dist_x, dist_y, dist_z, dist_x, dist_y, dist_z)), _mm_set_ps1(r * r));  
    
    __m128 sq = _mm_sqrt_ps(D);
    
    __m128 t0 = _mm_sub_ps(B, sq);
    __m128 t1 = _mm_add_ps(B, sq);
    
    hit_res = _mm_set1_epi32(MISS);
    
    for (int i=0;i<4;++i){
        if ((t0[i] > 0.1f) && (t0[i] < t[i])){
            t[i] = t0[i];
            hit_res[i] = HIT;
        }
    }
    
    for (int i=0;i<4;++i){
        if ((t1[i] > 0.1f) && (t1[i] < t[i])){
            t[i] = t1[i];
            hit_res[i] = HIT;
        }
    }
}

inline bool Sphere::hitSphere(Ray &r, float &t){
	Vec3 dist = c - r.o;
	float B = (r.d.x * dist.x + r.d.y * dist.y + r.d.z * dist.z);
	float D = B*B - DOT( dist, dist ) + this->r*this->r;
	if (D < 0.0f) return false;
    
    float sq = sqrtf(D);
    
	float t0 = B - sq;
	float t1 = B + sq;
	bool retvalue = false;
	if ((t0 > 0.1f ) && (t0 < t)){
		t = t0;
		retvalue = true;
	}
	if ((t1 > 0.1f ) && (t1 < t)){
		t = t1;
		retvalue = true;
	}
	return retvalue;
}

Vec3 CheckerBoard::getNormal(Vec3 &pos){
    return PlanePrim::getNormal(pos);
}


inline int Sphere::intersect( Ray &ray, float &t){
    return hitSphere(ray, t);
}

inline int PlanePrim::intersect(Ray &ray, float &t){
	float d = DOT( plane.N, ray.d );
	if (d != 0){
		float dist = -(DOT( plane.N, ray.o) + plane.D) / d;
		if (dist > 0){
			if (dist < t) {
				t = dist;
				return HIT;
			}
		}
	}
	return MISS;
}



Vec3 PlanePrim::getNormal(Vec3 &pos){
	return plane.N;
}

void Scene::init(){
	primitives = new Primitive*[25];
	nrOfPrimitives = 0;
    
    //wall at0,0,1
    Vec3 v121( 0,0, 1 );
	primitives[nrOfPrimitives] = new PlanePrim( v121, -13.0f );
	primitives[nrOfPrimitives]->getMaterial()->setReflection( 1 );
	primitives[nrOfPrimitives]->getMaterial()->setDiffuse( 0.0f );
    primitives[nrOfPrimitives]->getMaterial()->setSpecular(1.0f);
    Color c121( 1,1,1 );
	primitives[nrOfPrimitives++]->getMaterial()->setColor(c121);
    
    
    //wall at 0 0 -1
    Vec3 v12( 0,0, -1 );
	primitives[nrOfPrimitives] = new PlanePrim( v12, 13.0f );
	primitives[nrOfPrimitives]->getMaterial()->setReflection( 1.0f );
	primitives[nrOfPrimitives]->getMaterial()->setDiffuse( 0.0f );
    primitives[nrOfPrimitives]->getMaterial()->setSpecular(1.0f);
    Color c12( 0.3f,0.3f,0.7f );
	primitives[nrOfPrimitives++]->getMaterial()->setColor(c12);
    //-
    
    //the ground
    Vec3 v1( 0,1, 0 );
	primitives[nrOfPrimitives] = new CheckerBoard( v1, 3.4f );
	primitives[nrOfPrimitives]->getMaterial()->setReflection( 0.7f );
	primitives[nrOfPrimitives]->getMaterial()->setDiffuse( 1.0f );
    primitives[nrOfPrimitives]->getMaterial()->setSpecular(0.0f);
    Color c1( 0.3f, 0.3f, 0.7f );
	primitives[nrOfPrimitives++]->getMaterial()->setColor(c1);
    //
    
    
    // 0 -1 0
    Vec3 v13( 0,-1, 0 );
	primitives[nrOfPrimitives] = new PlanePrim( v13, -3.4f );
	primitives[nrOfPrimitives]->getMaterial()->setReflection( 1 );
	primitives[nrOfPrimitives]->getMaterial()->setDiffuse( 1.0f );
    primitives[nrOfPrimitives]->getMaterial()->setSpecular(1.0f - primitives[nrOfPrimitives]->getMaterial()->getDiffuse());
    Color c13( 1,1,1 );
	primitives[nrOfPrimitives++]->getMaterial()->setColor(c13);
    //
    
	//middle sphere
    Vec3 v2( 0, 1, 8 );
	primitives[nrOfPrimitives] = new Sphere( v2, 1.0f );
	primitives[nrOfPrimitives]->getMaterial()->setReflection( 1 );
    primitives[nrOfPrimitives]->getMaterial()->setDiffuse( 0.2 );
    primitives[nrOfPrimitives]->getMaterial()->setSpecular(1.0f - primitives[nrOfPrimitives]->getMaterial()->getDiffuse());
    Color c2( 1,1,1);
	primitives[nrOfPrimitives++]->getMaterial()->setColor( c2 );
	//-
    
    //left sphere
    Vec3 v3( -5.5f, -0.5, 7 );
    Color c3( 0.7f, 0.7f, 0 );
	primitives[nrOfPrimitives] = new Sphere( v3, 2 );
	primitives[nrOfPrimitives]->getMaterial()->setReflection( 0.6f );
	primitives[nrOfPrimitives]->getMaterial()->setDiffuse( 0.9f );
    primitives[nrOfPrimitives]->getMaterial()->setSpecular(1.0f - primitives[nrOfPrimitives]->getMaterial()->getDiffuse());
	primitives[nrOfPrimitives++]->getMaterial()->setColor( c3 );
    //-
    
    //right sphere
    Vec3 v31( 5.5f, -0.5, 7 );
    Color c31( 0.7f, 0.0f, 0 );
	primitives[nrOfPrimitives] = new Sphere( v31, 2 );
	primitives[nrOfPrimitives]->getMaterial()->setReflection( 0.6f );
	primitives[nrOfPrimitives]->getMaterial()->setDiffuse( 0.9f );
    primitives[nrOfPrimitives]->getMaterial()->setSpecular(1.0f);
	primitives[nrOfPrimitives++]->getMaterial()->setColor( c31 );
    //-
    
    //a light
    Vec3 v4( 3, 2.5f, 3 );
    Color c4(0.7f,0.7f,0.7f);
	primitives[nrOfPrimitives] = new Sphere( v4, 0.5f );
	primitives[nrOfPrimitives]->isLight = true;
	primitives[nrOfPrimitives++]->getMaterial()->setColor( c4 );
	//--
    
    //another light
    Vec3 v5( 0, 5, 10 );
    Color c5( 0.6f, 0.6f, 0.8f );
	primitives[nrOfPrimitives] = new Sphere( v5, 0.5f );
	primitives[nrOfPrimitives]->isLight = true;
	primitives[nrOfPrimitives++]->getMaterial()->setColor( c5 );
	//--    
}
