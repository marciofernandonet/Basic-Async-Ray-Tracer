//Supplied class (from demo project)

#ifndef _IMAGE_H_
#define _IMAGE_H_

#include "Vec3.h"
#include "iostream"

class Image { 
private:
	int width, height;
	Vec3 *pixels;

public:
	Image(const int w, const int h) : width(w), height(h) {
		pixels = new Vec3[width * height];
	}
	~Image(void) {
		delete [] pixels;
	}
	int getWidth(void) const { return width; }
	int getHeight(void) const { return height; }
	void setPixel(const int x, const int y, const Vec3 & c) { pixels[y * width + x] = c; }
	Vec3 getPixel(const int x, const int y) { return pixels[y * width + x]; }
	Vec3 * getPixelBufferPtr(void) { return pixels; }
    
    //method taken from 
    //http://stackoverflow.com/questions/3322705/difficulty-writing-a-bmp-file-in-c-manually-using-xcode 
    //2012-05-22
    void save(const char*fileName){ 
        FILE * out = fopen("ss.ppm", "wb");
        fprintf(out, "P6 %d %d 255\n", width, height);
        
        for(int y=0; y<height; y++)
            for(int x=0; x<width; x++)
            {
                Vec3 c = getPixel(x,height - y);
                putc(c.x*255, out);
                putc(c.y*255, out);
                putc(c.z*255, out);
            }
        
        fclose(out);
        
        std::cout << "Screenshot saved.\n";
    }
};

#endif