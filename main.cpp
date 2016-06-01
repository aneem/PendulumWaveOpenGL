#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <GL/glut.h>
#include <cstring>
#include "math.h"
#include <conio.h>

using namespace std;

const int ms_per_frame = 20;  // min ms per frame
int time_step_counter = -1, pause_flag = 1;
double **times, **thetas, **vs;

// pendulum parameters
const int ln=15;
double length[]={0.2115,0.2182,0.2252,0.2325,0.2402,0.2482,0.2567,0.2657,0.2751,0.285,0.2954,0.3065,0.3181,0.3305,0.3436};
double initial_length=5;
const double two_pi=2*3.1416;

//set lighting intensity and color
GLfloat qaAmbientLight[]    = {0.2, 0.2, 0.2, 2.0};
GLfloat qaDiffuseLight[]    = {0.8, 0.8, 0.8, 2.0};
GLfloat qaSpecularLight[]    = {1.0, 1.0, 1.0, 2.0};
GLfloat qaLightPosition[]    = {5, 10, 1, .1};

//colors
GLfloat black[] = {0.0, 0.0, 0.0, 1.0}; //Black Color
GLfloat black30[] = {0.0, 0.0, 0.0, .30}; //30% black Color
GLfloat white[] = {1.0, 1.0, 1.0, 1.0}; //White Color
GLfloat red[] = {1.0, 0.0, 0.0, 1.0};
GLfloat green[] = {0.0, 1.0, 0.0, 1.0};
GLfloat blue[] = {.0, 0.0, 1.0, 1.0};

//for material
GLfloat ambientMaterial[] = {.6, 0.65, .65, 1};
GLfloat diffuseMaterial[] = {1, 1, 1, 1};
GLfloat specularMaterial[] = {0.5, .50,.50, 1};
GLfloat specularMaterialshiness =0.6 ;

GLdouble x_camera=0,y_camera=0,z_camera=0,theta_camera=4.71,phi_camera=0.14,theta_camera_inc=0.01,radius_camera_inc=0.1,phi_camera_inc=0.01,radius_camera=11.7;

/* g= gravitational constant
theta0= initial angle
v0= intital angular velocity
ti= initial time
tf= final time duration
alpha= damping factor
m= mass of bob */
int N = 3000;
double g = 9.8, theta0 = 3.1416/4, v0 = 0, ti = 0, tf = 250, alpha = initial_length*1.6, m = 10;

 // opengl options 
void initgLOptions() {
   glEnable(GL_BLEND);
   glEnable(GL_LINE_SMOOTH);
   glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
}

void initLighting(){
    // Enable lighting
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    // Set lighting intensity and color
    glLightfv(GL_LIGHT0, GL_AMBIENT, qaAmbientLight);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, qaDiffuseLight);
    glLightfv(GL_LIGHT0, GL_SPECULAR, qaSpecularLight);
    // Set the light position
     glLightfv(GL_LIGHT0, GL_POSITION, qaLightPosition);
}

void drawaxes(void){
    glBegin(GL_LINES);
        glColor3fv(red);                //x axis
        glVertex3d(0.0,  0.0,  0.0);
        glVertex3d(10,  0.0,  0.0);
        glColor3fv(green);              //y axis
        glVertex3d(0.0,  0.0,  0.0);
        glVertex3d(00.0,  10.0,  0.0);
        glColor3fv(blue);               //z axis
        glVertex3d(0.0,  0.0,  0.0);
        glVertex3d(0.0,  0.0,  10.0);
    glEnd();
}

void camera(void){
    //drawaxes();
    glLoadIdentity();
    x_camera=radius_camera*cos(theta_camera);
    y_camera=radius_camera*sin(phi_camera);
    z_camera=radius_camera*sin(theta_camera);
    gluLookAt(x_camera,y_camera,z_camera,0,0,0,0,1,0);
}

void timer(int id) {
    if(pause_flag == 0){
        time_step_counter++;
        if(time_step_counter >= N)
            time_step_counter = 0;
    }
    glutPostRedisplay();
}

void rk4_sys_integrator(int pendno,double g, double l, double alpha, double m, double theta0, double v0, int N, double ti, double tf, double **times, double **thetas, double **vs){
    double h = (tf - ti)/N;
    double t = ti;
    double theta, v, k11, k12, k21, k22, k31, k32, k41, k42;
    int i;
    alpha=alpha/l;
    theta = theta0;
    v = v0;

    times[pendno] = (double *) malloc(sizeof(double)*N);
    thetas[pendno] = (double *) malloc(sizeof(double)*N);
    vs[pendno] = (double *) malloc(sizeof(double)*N);

    thetas[pendno][0] = theta;
    vs[pendno][0] = v;
    times[pendno][0] = t;

    for(i = 1; i<N; i++){
        k11 = h*v;
        k12 = -h*(alpha/m)*v - h*(g/l)*sin(theta);

        k21 = h*(v + k11/2);
        k22 = -h*(alpha/m)*(v + k12/2) - h*(g/l)*sin(theta + k12/2);

        k31 = h*(v + k21/2);
        k32 = -h*(alpha/m)*(v + k22/2) -h*(g/l)*sin(theta + k22/2);

        k41 = h*(v + k31);
        k42 =  -h*(alpha/m)*(v + k32) - h*(g/l)*sin(theta + k32);

        theta = theta + (k11 + 2*k21 + 2*k31 + k41)/6;
        v = v + (k12 + 2*k22 + 2*k32 + k42)/6;
        t = t + h;

        thetas[pendno][i] = theta;
        vs[pendno][i] = v;
        times[pendno][i] = t;
   }
}

void drawPendulum(int pendno,double l,double time, double theta, double v) {
    double pend_dist_multiplier=0.4;
    camera();
    glTranslatef(0,2,0);    //for centering the overall system of pendulum
    glBegin(GL_QUADS);
    //1st face
    glNormal3f( 0.0,  0.0,  1.0);
    glVertex3d(-0.5,  0.0,  0.0);
    glVertex3d( 0.5,  0.0,  0.0);
    glVertex3d( 0.5,  0.5,  0.0);
    glVertex3d(-0.5,  0.5,  0.0);
    //2nd face
    glVertex3d(-0.5,  0.0,  -6.0);
    glVertex3d( 0.5,  0.0,  -6.0);
    glVertex3d( 0.5,  0.5,  -6.0);
    glVertex3d(-0.5,  0.5,  -6.0);
    //3rd face
    glVertex3d(-0.5,  0.5,  -6.0);
    glVertex3d( -0.5,  0.5,  0.0);
    glVertex3d( -0.5,  0.0,  0.0);
    glVertex3d(-0.5,  0.0,  -6.0);
    //4th face
    glVertex3d(0.5,  0.5,  -6.0);
    glVertex3d( 0.5,  0.5,  0.0);
    glVertex3d( 0.5,  0.0,  0.0);
    glVertex3d(0.5,  0.0,  -6.0);
    //5th face
    glVertex3d(0.5,  0.5,  -6.0);
    glVertex3d( 0.5,  0.5,  0.0);
    glVertex3d( -0.5,  0.5,  0.0);
    glVertex3d(-0.5,  0.5,  -6.0);
    //6th face
    glVertex3d(0.5,  0.0,  -6.0);
    glVertex3d( 0.5,  0.0,  0.0);
    glVertex3d( -0.5,  0.0,  0.0);
    glVertex3d(-0.5,  0.0,  -6.0);
    glEnd();

    // The pendulum bob
    double pi = 3.1416,x1 = 0, y1 = 0,x2, y2;

    if(0 <= theta < pi/2){
            x2 = l*cos(pi/2-theta);
            y2 = -l*sin(pi/2-theta);
    }
    else if(pi/2 <= theta < pi){
            x2 = l*cos(theta - pi/2);
            y2 = l*sin(theta - pi/2);
    }
    else if(-pi/2 <= theta < 0){
            x2 = -l*cos(pi/2-abs(theta));
            y2 = -l*sin(pi/2-abs(theta));
    }
    else if(-pi <= theta < -pi/2){
            x2 = -l*cos(abs(theta) - pi/2);
            y2 = l*sin(abs(theta) - pi/2);
    }
    else{
        printf("invalid theta\n");
    }
    // string connecting the bob with the wooden plank
    glLineWidth(1);
    glColor3f(0.3, 0.3, 0.3);
    glBegin(GL_LINES);
    glVertex3f(x1,y1,- pendno*pend_dist_multiplier);
    glVertex3f(x2,y2,- pendno*pend_dist_multiplier);
    glEnd();

    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, ambientMaterial);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuseMaterial);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, specularMaterial);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, specularMaterialshiness);

    glPushMatrix();
    glTranslated(x2, y2,- pendno*pend_dist_multiplier);
    glutSolidSphere(0.2, 50, 50);
    glPopMatrix();
}

void printString(const char *str, double x, double y, double size) {
   glPushMatrix();
   glTranslatef(x,y,0);
   glScalef(size/153.0,size/153.0,1.0);
   int lineCt = 0;
   int len = strlen(str);
   for (int i = 0; i < len; i++) {
      if (str[i] == '\n') {
         lineCt++;
         glPopMatrix();
         glPushMatrix();
         glTranslatef(x,y-size*1.15*lineCt,0);
         glScalef(size/153.0,size/153.0,1.0);
      }
      else {
         glutStrokeCharacter(GLUT_STROKE_ROMAN,str[i]);
      }
   }
   glPopMatrix();
}

void display(void) {
    glutTimerFunc(ms_per_frame,timer,1);
    glClearColor (0.0,0.0,0.0,1.0);
    glClear (GL_COLOR_BUFFER_BIT);
    
    for(int i=0; i<ln;i++){
        drawPendulum(i,length[i],times[i][time_step_counter],thetas[i][time_step_counter],vs[i][time_step_counter]);
    }
    printString("Press w,a,s,d to pan screen.", 0, 0, 0.5);
    printString("Press +,- to z00m in/out.", 0, -.5, .5);
    printString("Press Spacebar to play/pause Animation and Esc to exit.", 0, -1, .5);
        glFlush();
        glutSwapBuffers();
}

void reshape (int w, int h) {
    glViewport (0, 0, (GLsizei)w, (GLsizei)h);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    gluPerspective (60, (GLfloat)w / (GLfloat)h, 0.1, 100.0);
    glMatrixMode (GL_MODELVIEW);
}

/* quit on escape key press and toggle start/stop of simulation with space key */
void keyboard(unsigned char key, int x, int y){
    switch (key) {
        case 27:       //esc button
            exit(0);
            break;
        case 32:     //space button
            if(pause_flag == 0) 
                pause_flag = 1;
            else
                pause_flag = 0;
            break;
        case 119:       //w
        case 87:
            phi_camera=phi_camera+phi_camera_inc;
            break;
        case 97:        //a
        case 65:
            theta_camera=theta_camera+theta_camera_inc;
            break;
        case 115:       //s
        case 83:
            phi_camera=phi_camera-phi_camera_inc;
            break;
        case 100:       //d
        case 68:
            theta_camera=theta_camera-theta_camera_inc;
            break;
        case 43:
            radius_camera=radius_camera-radius_camera_inc;
            break;
        case 45:
            radius_camera=radius_camera+radius_camera_inc;
            break;
    }
}

void init(void ){
    times = (double **) malloc(sizeof(double *)*ln);
    thetas = (double **) malloc(sizeof(double *)*ln);
    vs = (double **) malloc(sizeof(double *)*ln);

    for (int i=0;i<ln;i++){
        length[i]=length[i]*10;
    }
}

int main (int argc, char **argv) {
    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_DOUBLE);
    glutInitWindowSize (600, 600);
    glutInitWindowPosition (0, 0);
    glutCreateWindow ("Pendulum Wave Simulation");
    init();
    for(int i=0;i<ln;i++){
        rk4_sys_integrator(i,g, length[i], alpha, m, theta0, v0, N, ti, tf, &times[0], &thetas[0], &vs[0]);
    }
    initLighting();
    glutDisplayFunc (display);
    glutReshapeFunc (reshape);
    glutKeyboardFunc (keyboard);
    initgLOptions();
    glutMainLoop ();
    return 0;
}
