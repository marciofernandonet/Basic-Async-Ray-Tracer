//
//  Vec3.cpp
//  raytrace2
//
//  Created by Björn Dagerman on 2012-05-30.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include "Vec3.h"

float MAX(float A, float B){
    return (A > B) ? A : B;
}

__m128 dot_sse(__m128 &x1, __m128 &y1, __m128 &z1,__m128 &x2, __m128 &y2, __m128 &z2){
    return _mm_add_ps(_mm_add_ps(_mm_mul_ps(x1, x2), _mm_mul_ps(y1, y2)), _mm_mul_ps(z1, z2)); 
}

__m128 length_sse(__m128 &x, __m128 &y, __m128 &z){
    return _mm_sqrt_ps(dot_sse(x, y, z, x, y, z));
}

void normalize_sse(__m128 &x, __m128 &y,__m128 &z){
    
    __m128 l = _mm_div_ps( _mm_set1_ps( 1.0f ) , _mm_sqrt_ps( length_sse(x, y, z) ));
    
    x = _mm_mul_ps(x, l);
    y = _mm_mul_ps(y, l);
    z = _mm_mul_ps(z, l);
}

float InvSqrt(float x){
    float xhalf = 0.5f*x;
    int i = *(int*)&x; // get bits for floating value
    i = 0x5f3759df - (i>>1); // gives initial guess y0
    x = *(float*)&i; // convert bits back to float
    x = x*(1.5f-xhalf*x*x); // Newton step, repeating increases accuracy
    return x;
}


float Q_rsqrt(float number){//Quake 3 fast sqrt
    long i;
    float x2, y;
    const float threehalfs = 1.5F;
    
    x2 = number * 0.5F;
    y  = number;
    i  = * ( long * ) &y;                       // evil floating point bit level hacking
    i  = 0x5f3759df - ( i >> 1 );               // what the fuck?
    y  = * ( float * ) &i;
    y  = y * ( threehalfs - ( x2 * y * y ) );   // 1st iteration
    //      y  = y * ( threehalfs - ( x2 * y * y ) );   // 2nd iteration, this can be removed
    
    return y;
}