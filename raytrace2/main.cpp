//
//  main.cpp
//  raytrace2
//
//  Created by Björn Dagerman on 2012-05-22.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include </usr/include/GL/glew.h>
#ifdef __APPLE__
#include <OpenGL/OpenGL.h>
#include <GLUT/GLUT.h>

#else
#include <GL/glut.h>
#endif
#include <iostream>
#include "Vec3.h"
#include "raytracer.h"
#include "scene.h"
#include "Image.h"
#include "scene.h"

#include <stdio.h>
#include <time.h>
float samplerate = 4.0f;
int width = 1280, height = 720;
int window_width = 1280, window_height = 720;
//int width = 640, height = 480;
//int window_width = 640, window_height = 480;

Tracer *tracer = 0;
void *rgb_buffer = 0;

bool useSSE = true;
bool frameLess = true;

Image *image = new Image(window_width, window_height);	
Vec3 camPos(0,0,-5);
int nrOfIntersectionTests = 0;

bool showBenchmarking = false;
unsigned int text;


void init();
void keypress(unsigned char key, int x, int y);
void changeSize(int w, int h);
void display(void);
void update(int i);
void saveImg();


 
int main (int argc, char * argv[])
{
    // insert code here...
    std::cout << "Hello, World!\n";
    
    glutInit(&argc, argv);
    
    
    init();
    update(0);

    Tracer::getInstance()->init();

    glutMainLoop();

    return 0;
}

void update(int i){
    time_t start = clock();
    
    Tracer::getInstance()->render(camPos);

    glutPostRedisplay();
    glutTimerFunc( MAX(16 - ((clock() - start)/1000.0f), 0) , update, 0);
}

void init(void){
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutInitWindowSize(window_width, window_height);
	glutCreateWindow("");
	glutDisplayFunc(display);
	glutReshapeFunc(changeSize);
	glutKeyboardFunc(keypress);
    
    glEnable(GL_TEXTURE_2D);
    glGenTextures(1,&text);
    
    glBindTexture(GL_TEXTURE_2D, text); 
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); 
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    

	glClearColor(0.0f, 1.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    
    
    rgb_buffer = malloc(sizeof(float)*window_width*window_height*3);
    
    tracer = Tracer::getInstance();
    tracer->showBenchmarking = &showBenchmarking;
    tracer->setTarget(&width, &height, &window_width, &window_height, rgb_buffer, &samplerate);
}

int lastdrawn = 0;
int lastc = 0;
time_t totaltime = 0;
int drawNr = 0;

void display(void){
    time_t start = clock();
    
    glClear(GL_COLOR_BUFFER_BIT);
    
    int c = 0;
    
    glBegin(GL_POINTS);
    for (int y=0;y<window_height;++y){
        for (int x=0;x<window_width;++x){
            glColor3fv((float*)rgb_buffer + c);
            glVertex2i(x, y);
            c+=3;
        }
    }
    glEnd();
     
    glutSwapBuffers();
    
    totaltime += (clock() - start);
    drawNr++;
    if (drawNr >= 100 && showBenchmarking){
        float total = totaltime*1.0f/CLOCKS_PER_SEC;
        cout << "Last " << drawNr << "draws took: " << total  << " sekunder. " << total/(float)drawNr <<  " seconds per draw (avg).\n";
        drawNr=0;
        totaltime = 0;
    }
    
}

void saveImg(){
    int c = 0;
    for (int y=0;y<window_height;++y){
        for (int x=0;x<window_width;++x){
            float *frgbb = (float*)rgb_buffer;
            Vec3 clr(frgbb[c],frgbb[c+1],frgbb[c+2]);
            image->setPixel(x, y, clr);
            c+=3;
        }
    }
    
    image->save("img.bmp");
}

float movementspeed = 0.05f;
void keypress(unsigned char key, int x, int y) {
    switch (key) {
        case 'd':
            camPos.x-=movementspeed;
            break;
        case 'a':
            camPos.x+=movementspeed;
            break;
        case 'w':
            camPos.z-=movementspeed;
            break;
        case 's':
            camPos.z+=movementspeed;
            break;
        case '<':
            camPos.y+=movementspeed;
            break;
        case '>':
            camPos.y-=movementspeed;
            break;
        case 'p':
            saveImg();
            break;
        case 'i':
            tracer->useRandom = !tracer->useRandom;
            std::cout << "Random region" << (tracer->useRandom ? " enabled\n" : " disabled\n");
            tracer->resetRegions();
            break;
        case 'o':
            tracer->renderChecker = !tracer->renderChecker;
            break;
        case 'f':
            frameLess = !frameLess;
            std::cout << "frameless " << (frameLess ? "enabled\n" : "disabled\n");
            break;
        case 'u':
            showBenchmarking = !showBenchmarking;
            std::cout << (showBenchmarking ? "Showing " : "Hiding ") << "benchmarking results.\n";
            break;
        case 'm':
            useSSE = !useSSE;
            break;
        case '+':
            samplerate++;
            width = window_width*samplerate;
            //tracer->resetRegions();
           // tracer->SetTarget(&width, &height, &window_width, &window_height, buffer, x_buffer, y_buffer, z_buffer, rgb_buffer, c_rgb_buffer);
            break;
            
        case '-':
            samplerate--;
            height = window_height*samplerate;
            //tracer->resetRegions();
            //tracer->SetTarget(&width, &height, &window_width, &window_height, buffer, x_buffer, y_buffer, z_buffer, rgb_buffer, c_rgb_buffer);
            break;
        default:
            break;
    }
    tracer->shouldRender = false;
}

void changeSize(int w, int h){
    window_width = w;
    window_height = h;
    delete image;
    image = new Image(w,h);    
    rgb_buffer = realloc(rgb_buffer, sizeof(float)*w*h*3);
    
    glViewport(0,0,w,h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0,(GLfloat)w,0,(GLfloat)h);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

