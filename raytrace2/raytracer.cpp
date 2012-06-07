
//
//  raytracer.cpp
//  raytrace2
//
//  Created by Björn Dagerman on 2012-05-22.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.

#include </usr/include/GL/glew.h>
#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/GLUT.h>

#else
#include <GL/glut.h>
#endif

#include "raytracer.h"
#include "scene.h"
#include "stdint.h"
#include "stdio.h"
#include <iostream>
#include <pthread.h>

void *threaded_render(void *cam);

Tracer*Tracer::instance = 0;

Tracer::Tracer(){
	scene = new Scene();
    campos = Vec3(0,0,0);
    renderChecker = true;
    maxDepth = MAX_DEPTH;
}

void Tracer::setTarget(int *width, int *height,int *ww, int*wh, void*rgb_b, float*sr){
	this->width = width;
	this->height = height;
    this->window_width = ww;
    this->window_height = wh;
    samplerate = sr;
    campos = Vec3(0,0,-5);
    useRandom = true;
    rgb_buffer = rgb_b;
    
    resetRegions();
    shuffle();
}

//Color ambient(0.01f,0.01f,0.01f);
Color ambient(0.0f,0.0f,0.0f);
Color bg(0.01f,0.01f,0.01f);

void Tracer::trace_sse(__m128 &rox, __m128 &roy,__m128 &roz,__m128 &rdx,__m128 &rdy,__m128 &rdz, int depth, float &aDist, int&nrOfIntersectionTests, __m128 &r, __m128&g, __m128&b, bool useSSE, float coef) {
	if (depth > maxDepth || coef <= 0.001f) return;
    
    __m128 bgcolor_r, bgcolor_g, bgcolor_b;
    bgcolor_r = bgcolor_g = bgcolor_b = _mm_set1_ps(0.2f);
	aDist = 1000000.0f;
    
    __m128 dist = _mm_set1_ps(aDist);
    
	Vec3 pi;
	Primitive* prim = 0;
    
    int hit = -1;
    
    __m128 hit_res;
    __m128 pix, piy, piz;
    
	int result;
	for ( int s = 0; s < scene->nrOfPrimitives; s++ )
	{
		Primitive* pr = scene->getPrimitive(s);
        pr->intersect_sse(rox, roy, roz, rdx, rdy, rdz, dist, hit_res);
        
        
        nrOfIntersectionTests++;
        for (int i=0;i<4;++i){
            int res = hit_res[i];
            if (res) {
                prim = pr;
                result = res;
                hit = s;
            }
        }
	}
	if (!prim || prim->isLight) {
        for (int i =0;i<4;++i){
            if (r[i] == g[i] == b[i] == 0){
                r = bgcolor_r;
                g = bgcolor_g;
                b = bgcolor_b;
            }
        }
        
        return;
    }
	else{
        pix = _mm_add_ps( rox,  _mm_mul_ps(rdx, dist) );
        piy = _mm_add_ps( roy,  _mm_mul_ps(rdy, dist) );
        piz = _mm_add_ps( roz,  _mm_mul_ps(rdz, dist) );
        
        
		for ( int l = 0; l < scene->nrOfPrimitives; l++ )
		{
			Primitive* p = scene->getPrimitive(l);
            
            __m128 shade = _mm_set1_ps(1.0f);

			if (p->isLight) {
				Primitive* light = p;
                Vec3 c = ((Sphere*)light)->getCenter();
                __m128 Lx = _mm_sub_ps(_mm_set1_ps(c.x), pix);
                __m128 Ly = _mm_sub_ps(_mm_set1_ps(c.y), piy);
                __m128 Lz = _mm_sub_ps(_mm_set1_ps(c.z), piz);
                
				if (light->getType() == Primitive::SPHERE){
                    __m128 tdist = length_sse(Lx, Ly,Lz);
                    
                    __m128 rsq = _mm_div_ps(_mm_set1_ps(1.0f), tdist);
                    
                    Lx = _mm_mul_ps(Lx, rsq);
                    Ly = _mm_mul_ps(Ly, rsq);
                    Lz = _mm_mul_ps(Lz, rsq);
                    
                    __m128 rvx = _mm_add_ps(pix, _mm_mul_ps(Lx, _mm_set1_ps(EPSILON)));
                    __m128 rvy = _mm_add_ps(piy, _mm_mul_ps(Ly, _mm_set1_ps(EPSILON)));
                    __m128 rvz = _mm_add_ps(piz, _mm_mul_ps(Lz, _mm_set1_ps(EPSILON)));
                    
                    for ( int s = 0; s < scene->nrOfPrimitives; s++ ){
						Primitive* pr = scene->getPrimitive(s);
						
                        nrOfIntersectionTests++;
                        if (!pr->isLight){
                            __m128 hit_res2;
                            pr->intersect_sse(rvx, rvy, rvz, Lx, Ly, Lz, tdist, hit_res2);
                            
                            for (int i = 0; i<4;++i){
                                if (hit_res2[i] == HIT){
                                    shade[i] = 0.0f;
                                }
                            }                            
						}
					}
				}
				normalize_sse(Lx, Ly, Lz);
				Vec3 n = prim->getNormal(pi);
                
                __m128 Nx = _mm_set1_ps(n.x);
                __m128 Ny = _mm_set1_ps(n.y);
                __m128 Nz = _mm_set1_ps(n.z);
                
                
				if (prim->getMaterial()->getDiffuse() > 0){
                    __m128 dot = dot_sse(Nx, Ny, Nz, Lx, Ly, Lz);
                    
                    __m128 diff = _mm_mul_ps(_mm_mul_ps(_mm_set1_ps(prim->getMaterial()->getDiffuse()), dot) , shade);
                    
                    
                    Vec3 pmc = prim->getMaterial()->getColor();
                    
                    __m128 pmc_r = _mm_set1_ps(pmc.r), pmc_g = _mm_set1_ps(pmc.g), pmc_b = _mm_set1_ps(pmc.b);
                    
                    if (prim->getType() == Primitive::CHECKERBOARD && renderChecker){
                        for (int i =0; i < 4;++i){
                            bool horz = ((int)(pix[i]*100.0f/80.0f) % 2 == 0); 
                            bool vert = ((int)(piz[i]*100.0f/80.0f) % 2 == 0); 
                            
                            if( (horz && !vert) || (!horz && vert)) {
                                pmc_r[i] *= 1.3f;
                                pmc_g[i] *= 1.3f;
                                pmc_b[i] *= 1.3f;
                            } 
                        }
                    }
                    
                    Vec3 lmc = light->getMaterial()->getColor();
                    
                    __m128 lmc_r = _mm_set1_ps(lmc.x), lmc_g = _mm_set1_ps(lmc.y), lmc_b = _mm_set1_ps(lmc.z);
                    
                    for (int i=0;i<4;++i){
                        if (diff[i] < 0) diff[i] = 0;
                    }
                    
                    r = _mm_add_ps(r, _mm_mul_ps(lmc_r, diff));
                    g = _mm_add_ps(g, _mm_mul_ps(lmc_g, diff));
                    b = _mm_add_ps(b, _mm_mul_ps(lmc_b, diff));
                }
                //todo: add specular
			}
		}
        
        //todo: add reflections
    }
    
    for (int i =0;i<4;++i){
        if (r[i] > 1.0f) r[i] = 1.0f;
        if (g[i] > 1.0f) g[i] = 1.0f;
        if (b[i] > 1.0f) b[i] = 1.0f;
    }
    
	return;
}

void Tracer::trace(Ray &ray, int depth, float &dist, int&nrOfIntersectionTests, float &r, float&g, float&b, bool useSSE, float coef) {
	if (depth > maxDepth || coef <= 0.001f) return;
    
    Color bgcolor = bg;
	dist = 1000000.0f;
	Vec3 pi;
	Primitive *prim = 0;

    if (r == g == b == 0){
        r = ambient.r;
        g = ambient.g;
        b = ambient.b;   
    }
    
    int hit = -1;
	int result;
	for (int s=0;s<scene->nrOfPrimitives;++s) {
		Primitive* pr = scene->getPrimitive(s);
		int res = pr->intersect(ray, dist);
        //nrOfIntersectionTests++;
		if (res) {
			prim = pr;
			result = res;
            hit = s;
		}
	}
	if (!prim || prim->isLight) {
        if (r == ambient.r && g == ambient.g && b == ambient.b){
            r += bgcolor.r;
            g += bgcolor.g;
            b += bgcolor.b;   
        }
        return;
    }
	else{
		pi = ray.o + ray.d * dist;

		for ( int l = 0; l < scene->nrOfPrimitives; l++ )
		{
			Primitive* p = scene->getPrimitive(l);
            float shade = 1.0f;
			if (p->isLight) {
				Primitive* light = p;
                Vec3 L = ((Sphere*)light)->getCenter() - pi;
                float ldist = LENGTH(L);
                L*=(1.0f/ldist);
                
                Vec3 rv = pi + L * EPSILON;
                Ray r2 = Ray( rv, L );
                
                for ( int s = 0; s < scene->nrOfPrimitives; s++ ){
                    Primitive* pr = scene->getPrimitive(s);
                    
                    //nrOfIntersectionTests++;
                    if ((!pr->isLight) && (pr->intersect( r2, ldist ))){
                        if (hit != s){
                            shade = 0;
                            break;
                        }
                        
                    }
                }
                
				Vec3 N = prim->getNormal(pi);
				if (prim->getMaterial()->getDiffuse() > 0){
					float dot = DOT(N, L);
					if (dot > 0){                        
						float diff = dot * prim->getMaterial()->getDiffuse() * shade;

                        Vec3 pmc = prim->getMaterial()->getColor();
                        if (prim->getType() == Primitive::CHECKERBOARD && renderChecker){
                            bool h = ((int)(pi.x*100.0f/80.0f) % 2 == 0); 
                            bool v = ((int)(pi.z*100.0f/80.0f) % 2 == 0); 
                            if((h && !v) || (!h && v))
                              pmc *= 1.3f;
                        }
                        
                        
                        Vec3 lmc = light->getMaterial()->getColor();
                        
                        r += diff *pmc.x*lmc.x;
                        g += diff *pmc.y*lmc.y;
                        b += diff *pmc.z*lmc.z;
					}
				}

				if (prim->getMaterial()->getSpecular() > 0 && depth < (maxDepth-2)){
					Vec3 R = L - 2.0f * DOT( L, N ) * N;
                    float phongTerm = MAX(DOT(R, ray.d), 0.0f) ;
                    phongTerm = powf(phongTerm, 4) * prim->getMaterial()->getSpecular() * shade;
                    Vec3 lmc = light->getMaterial()->getColor();
                    r += phongTerm * lmc.r;
                    g += phongTerm * lmc.g;
                    b += phongTerm * lmc.b;
				}
			}
		}
        float refl = prim->getMaterial()->getReflection();
        if (refl > 0.0f){
            Vec3 N = prim->getNormal( pi );
            
            Vec3 R = ray.d - 2.0f * N * DOT(ray.d, N);
            float r2 = 0.0f;
            float g2 = 0.0f;
            float b2 = 0.0f;
            
            float dist;
            
            Vec3 rv = pi + R *EPSILON;
            Ray ray(rv, R);
            
            trace(ray, depth + 1, dist, nrOfIntersectionTests, r2, g2, b2, useSSE, coef * refl);
            Vec3 pmc = prim->getMaterial()->getColor();
            Color rcol(r2,g2,b2);
            Vec3 reflection(refl*rcol.x*pmc.x, refl*rcol.y*pmc.y, refl*rcol.z*pmc.z);
            
            r= r*(1.0f - refl) + reflection.r;
            g= g*(1.0f - refl) + reflection.g;
            b= b*(1.0f - refl) + reflection.b;
        }
	}
    if (r > 1.0f) r = 1.0f;
    if (g > 1.0f) g = 1.0f;
    if (b > 1.0f) b = 1.0f;
    
	return;
}

const float MULT = 3.0f;

void Tracer::init(){ 
    shouldRender = true;    
    float w = ((MULT*(float)*window_width)/((float)*window_height));
	xs = -w, xe = w, ys = yStart = -MULT, ye = MULT;

	dx = (xe - xs) / *window_width;
	dy = (ye - ys) / *window_height;
}

vector<Region> vr_regions;

void Tracer::resetRegions(){
    vr_regions.clear();
    
    float dy = *window_height/(float)NR_OF_REGIONS;
    float dx = *window_width/(float)NR_OF_REGIONS;
    
    for (int y=0;y<NR_OF_REGIONS;++y){
        for (int x = 0; x<NR_OF_REGIONS;++x){
            int xss = x*dx;
            int yss = y*dy;
            int xee = (x+1)*dx*1.01f;// + (int)(dx+10.0f);
            int yee = (y+1)*dy*1.01f;// + (int)(dy+10.0f);
            Region r(xss, xee, yss, yee);
            
            vr_regions.push_back(r);
        }
    }
}
int lastregion = -1;

void Tracer::shuffle(){
    lastregion = -1;
    
    for (int i=0; i<vr_regions.size(); i++) {
        unsigned long idx = i + (rand() % (vr_regions.size()-i));
        Region temp = vr_regions.at(i);
        vr_regions.at(i) = vr_regions.at(idx);
        vr_regions.at(idx) = temp;
    }
}

int complete = 0;
time_t start = 0;

void *threaded_render(void *cam){
    int *id = (int*)cam;
    int i = *id;
    delete id;
    
    Tracer *tracer = Tracer::getInstance();
    
    Region r = vr_regions.at((++lastregion) % vr_regions.size());
    if (complete >= vr_regions.size() && complete <= vr_regions.size()+(NR_OF_THREADS/2)  && *tracer->showBenchmarking){
        clock_t end = clock();
        float total = end - start;
        
        cout << "complete rendering took: " << total/CLOCKS_PER_SEC << " seconds.\n";
    }
    int depth = complete >= NR_OF_THREADS ? 1 : tracer->maxDepth-2;
    tracer->render_sse(r.ystart, r.yend, r.xstart, r.xend, depth);
    ++complete;

    pthread_join(tracer->threads[i], 0);
    return 0;
}




void Tracer::render(Vec3 &cam){
    if (start == 0) start = clock();
    
    if (campos.x != cam.x ||campos.y != cam.y ||campos.z != cam.z){
        campos = cam;
        if (useRandom)
            shuffle();    
        complete = 0;
        start = clock();
    }

    shouldRender = true;
    
    for (int i = 0; i < NR_OF_THREADS;++i){
        int msg = pthread_kill(threads[i], 0);
        if (msg != 0){
            int *id = new int();
            *id = i;
            pthread_create(threads+i, NULL, threaded_render, id);
        }
    }
}

void Tracer::render_sse(int ystart,int yend, int xstart, int xend, int maxdepth){
    //float ww = *window_width/(*width);
    //float hh = *window_height/(*height);
    
//    const float MULT = 3.0f;
    float w = ((MULT*(float)*window_width)/((float)*window_height));
	float m_WX1 = -w,m_SY = -MULT;
    m_SY += dy * ystart;
    xStart += dx * xstart;
    
    Vec3 cam = campos;
    __m128 defx = _mm_set_ps(m_WX1 + 3*dx, m_WX1 + 2*dx, m_WX1 + dx, m_WX1);
    defx = _mm_add_ps(defx, _mm_set_ps1(xstart*dx));
    __m128 dx(_mm_set1_ps(4*this->dx));
    
    __m128 cx;
    __m128 cy;
    
    __m128 dir_z = _mm_set1_ps(-cam.z);
    __m128 dir_z_sq =  _mm_mul_ps(dir_z, dir_z);

    

    float *clrs = new float[12];
    /*
    int clrs_size = 3*4**samplerate**samplerate;
    float tmp_clrs[clrs_size];
    */
        for (int yi=ystart;(yi < yend && yi < *window_height); ++yi){
            int c = 3*(yi* (*window_width) + xstart);
            cx = defx;
            cy = _mm_set1_ps(m_SY);
            __m128 dir_y = _mm_sub_ps(cy, _mm_set1_ps(cam.y));
            __m128 dir_y_sq =  _mm_mul_ps(dir_y, dir_y);
            //int xlaps = 0;
            
            __m128 dir_zy_sq = _mm_add_ps(dir_y_sq, dir_z_sq);
            for (int xi=xstart; (xi < (xend-3) &&(xi < (*window_width-3))); xi+=4){
                
                if (!shouldRender){delete clrs; return;}
                
                __m128 dir_x = _mm_sub_ps(cx, _mm_set1_ps(cam.x));
                __m128 dir_x_sq =  _mm_mul_ps(dir_x, dir_x);          
                __m128 dir_xyz_sq = _mm_add_ps(dir_x_sq, dir_zy_sq);
                
                //using rsqrt&mult
                
                /*__m128 mult = _mm_rsqrt_ps(dir_xyz_sq);
                
                 __m128  dir_x_n = _mm_mul_ps(dir_x, mult); 
                 __m128  dir_y_n = _mm_mul_ps(dir_y, mult); 
                 __m128  dir_z_n = _mm_mul_ps(dir_z, mult); */
                
                //-
                //using sqrt&div
                
                
                __m128 mult1 = _mm_sqrt_ps(dir_xyz_sq);
                
                __m128 dir_x_n = _mm_div_ps(dir_x, mult1);
                __m128 dir_y_n = _mm_div_ps(dir_y, mult1);
                __m128 dir_z_n = _mm_div_ps(dir_z, mult1);
                
                //-
                //wow, using rsqrt is way slower, why?
                
                Vec3 r1_d(dir_x_n[0], dir_y_n[0],dir_z_n[0]);
                Vec3 r2_d(dir_x_n[1], dir_y_n[1],dir_z_n[1]);
                Vec3 r3_d(dir_x_n[2], dir_y_n[2],dir_z_n[2]);
                Vec3 r4_d(dir_x_n[3], dir_y_n[3],dir_z_n[3]);
                
                float d1, d2, d3 ,d4; 
                memset(clrs, 0, sizeof(float)*12);
                
                int intersectiontests = 0; //why do you exist?
               
                
                if (true){
                    Ray ray1(cam, r1_d);
                    Ray ray2(cam, r2_d);
                    Ray ray3(cam, r3_d);
                    Ray ray4(cam, r4_d);
                    
                    trace(ray1, maxdepth, d1, intersectiontests, *(clrs), *(clrs+1),*(clrs+2), true, 1);
                    trace(ray2, maxdepth, d2, intersectiontests, *(clrs+3), *(clrs+4),*(clrs+5), true, 1);
                    trace(ray3, maxdepth, d3, intersectiontests, *(clrs+6), *(clrs+7),*(clrs+8), true, 1);
                    trace(ray4, maxdepth, d4, intersectiontests, *(clrs+9), *(clrs+10),*(clrs+11), true, 1);
                    
                    memcpy(((float*)rgb_buffer + c), clrs, sizeof(float)*12);
                    
                }
                else{
                    __m128 ox = _mm_set1_ps(cam.x);
                    __m128 oy = _mm_set1_ps(cam.y);
                    __m128 oz = _mm_set1_ps(cam.z);
                    
                    
                    __m128 rr, gg, bb;
                    
                    trace_sse(ox, oy, oz, dir_x_n, dir_x_n, dir_z_n, maxdepth, d1, intersectiontests, rr, gg, bb, true, 1);
                    
                    for (int i=0;i<4;++i){
                        int idx = c+i*3;
                        *((float*)rgb_buffer +idx) = rr[i];
                        *((float*)rgb_buffer +idx+1) = gg[i];
                        *((float*)rgb_buffer +idx+2) = bb[i];
                    }
                }
                c+=12;
                
                cx = _mm_add_ps(cx, dx);
            }
            m_SY += dy;
        }
        delete clrs;
}

void Tracer::render(Vec3 &cam,int&nrOfIntersectionTests, int xstart, int xend, int ystart, int yend){
	Vec3 o = cam;
    
    //float ww = *window_width/(*width);
    //float hh = *window_height/(*height);
    
    const float MULT = 3.0f;
    float w = ((MULT*(float)*window_width)/((float)*window_height));
    
	float m_WX1 = -w, m_SY = -MULT;
    m_SY += dy * ystart;
    xStart += dx * xstart;
    
    unsigned int index = 0;
	for ( int y = ystart; y < (yend); y++ ) {
		xStart = m_WX1;

		for ( int x = xstart; x < xend; x++ ){
            float r = 0;
            float g = 0;
            float b = 0;
            
            for(float fragmentx = xStart; fragmentx < xStart + dx; fragmentx += dx*0.5f)
                for(float fragmenty = yStart; fragmenty < yStart + dy; fragmenty += dy*0.5f)
                {
                    if (!shouldRender) return;
                    
                    Vec3 dir = Vec3(fragmentx, fragmenty,0) -o;
                    NORMALIZE(dir);
                    Ray ray(o, dir);
                    float dist = 1000000.0f;
                    
                    float coef = 0.25f;
                    float red = 0, green = 0, blue = 0; 
                    trace(ray, 1, dist, nrOfIntersectionTests, red, green, blue, false, 1);
                    r += red *coef;
                    g += green *coef;
                    b += blue*coef;
            }
            
            *(((float*)rgb_buffer + index++)) = r;
            *(((float*)rgb_buffer + index++)) = g;
            *(((float*)rgb_buffer + index++)) = b;
            
			xStart += dx;
		}
		m_SY += dy;
	}
}