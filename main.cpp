#include <iostream>
#include "freeglut/glut.h"
#include "gjk_epa.hpp"

using namespace std;

using namespace std;


#define GLUT_WHEEL_UP   3
#define GLUT_WHEEL_DOWN 4



struct Shape
{
    Vector3 mPosition;

    Shape()
    {
    }

    virtual const Vector3 support_world( const Vector3& direction ) const {}
    virtual void draw() {}
};


struct ShpereShape : public Shape
{
    float   mRadius;

    ShpereShape( float radius )
     : Shape() , mRadius(radius)
    {
    }


    const Vector3 support_world( const Vector3& direction ) const
    {
        Vector3 res =  direction.getUnit() * mRadius;
        return res + mPosition;
    }

    void draw()
    {
        glPushMatrix();
        glTranslatef( mPosition.x , mPosition.y , mPosition.z );
        glutWireSphere(mRadius , 10 , 10);
        glPopMatrix();
    }

};




struct BoxShape : public Shape
{
    Vector3 mExtent;

    BoxShape( const Vector3& extent )
    : mExtent(extent)
    {

    }


    const Vector3 support_world( const Vector3& direction ) const
    {
        Vector3 res = Vector3(direction.x < 0.0 ? -mExtent.x : mExtent.x,
                              direction.y < 0.0 ? -mExtent.y : mExtent.y,
                              direction.z < 0.0 ? -mExtent.z : mExtent.z);
        return res + mPosition;
    }

    void draw()
    {
         glPushMatrix();
         glTranslatef( mPosition.x , mPosition.y , mPosition.z );
         glScalef(mExtent.x , mExtent.y , mExtent.z);
         glutWireCube(2.0);
         glPopMatrix();

    }
};



Shape *shape1 = new BoxShape(Vector3(2,1,2));
Shape *shape0 = new ShpereShape(2);


//=================================  Unit Test ==================================//

float Width = 600;
float Height = 400;

void print_text(int x, int y, string String)
{
   //(x,y) is from the bottom left of the window
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, Width , 0, Height, -1.0f, 1.0f);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glPushAttrib(GL_DEPTH_TEST);
    glDisable(GL_DEPTH_TEST);
    glRasterPos2i(x,y);
    for(unsigned int i=0; i<String.size(); i++)
    {
        glutBitmapCharacter(GLUT_BITMAP_9_BY_15, String[i]);
    }
    glPopAttrib();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}


/// Camera value
Vector3 mEye;
Vector3 mCenter;
Vector3 mUp;


float TimeStep = 1.0 / 60.0;
bool Pause = true;

// Initialization
void initCamera()
{

    float aspect = Width / Height;
    float zNear  = 1.0;
    float zFar   = 1024;
    float fov    = 45.0;

    mEye    =  Vector3(0,0,0);
    mCenter =  Vector3(0,0,0);
    mUp     =  Vector3(0,1,0);


    glLoadIdentity();

    glMatrixMode(GL_PROJECTION);
    gluPerspective( fov , aspect , zNear , zFar );


    glMatrixMode(GL_MODELVIEW);
    //	  gluLookAt(0.0, 0.0, 5.0,  /* eye is at (0,0,5) */
    //	    0.0, 0.0, 0.0,      /* center is at (0,0,0) */
    //	    0.0, 1.0, 0.);      /* up is in positive Y direction */

}




// Display the scene
void display()
{

    glClearColor(0.2f, 0.2f, 0.2f, 0.0f);

    glViewport(0, 0, Width , Height );




    glLoadIdentity();

    //	float aspect = 1.0;//Width / Height;
    //	float zNear  = 1.0;
    //	float zFar   = 1024;
    //	float fov    = 45.0;
    //
    //	glMatrixMode(GL_PROJECTION);
    //	gluPerspective( fov , aspect , zNear , zFar );



    glMatrixMode(GL_MODELVIEW);
    gluLookAt( -mEye.x , -mEye.y , -mEye.z ,  mCenter.x , mCenter.y , mCenter.z , mUp.x , mUp.y , mUp.z );



    glEnable(GL_DEPTH_TEST);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    ContactBasicData Date;
    if( GJKEPAGenerator::CollisionSolveGJK_EPA(shape0 , shape1 , &Date))
    {

         Date.normal.normalize();

         Vector3 a = Date.point;
         Vector3 b = Date.point + Date.normal * Date.penetration;

         glPushMatrix();
         glColor3f(0,1,1);
         glLineWidth(8);
         glBegin(GL_LINES);
          glVertex3f(a.x,a.y,a.z);
          glVertex3f(b.x,b.y,b.z);
         glEnd();
         glLineWidth(1);
         glColor3f(0,1,0);
         glPopMatrix();

         print_text( 20 , Height - 60 , " COLLISION : Yes");
    }
    else
    {
        glColor3f(1,1,1);
        print_text( 20 , Height - 60 , " COLLISION : No");

    }


    shape0->draw();
    shape1->draw();



    print_text( 20 , Height - 20 , " KEY   : (UP,DOWN,RIGHT,LEFT) : Move-Object");
    print_text( 20 , Height - 40 , " MOUSE : Move-Camera");


    glutSwapBuffers();


}






float zoom_distance = 25;

float oldX = 0;
float oldY = 0;

float angle_X = 0;
float angle_Y = 0;

int   mouse_button = -1;



void UpdateTime(void)
{

    //**********  Camera Position *******//
    Vector3 DircetionLeft(1,0,0);
    Vector3 DircetionLook(0,0,zoom_distance);
    DircetionLook.rotateAroundAxis(Vector3(0,1,0) , angle_X);
    DircetionLeft.rotateAroundAxis(Vector3(0,1,0) , angle_X);
    DircetionLook.rotateAroundAxis(DircetionLeft  , angle_Y);
    mEye = DircetionLook;//+ Vector3(-10,-15,0);


    //************ Render Scene *********//
    display();


};


// Reshape function
void reshape(int width, int height)
{
    Width  = width;
    Height = height;
}


// Called when a mouse button event occurs
void mouseButton(int button, int state, int x, int y)
{
    mouse_button = button;

    float aspect = Width / Height;
    float m_x = ((x / Width ) - 0.5f) * aspect * 0.834;
    float m_y = ((y / Height) - 0.5f) * 0.834;

    oldX = m_x;
    oldY = m_y;

    if (state == GLUT_UP )
    {
        if ( button == GLUT_WHEEL_UP )
        {
            zoom_distance -= 0.5;
        }
        else if( button == GLUT_WHEEL_DOWN )
        {
            zoom_distance += 0.5;
        }
    }


    if (button == GLUT_RIGHT_BUTTON )
    {

    }
}

// Called when a mouse motion event occurs
void mouseMotion(int x, int y)
{
   float aspect = Width / Height;
   float m_x = ((x / Width ) - 0.5f) * aspect * 0.834;
   float m_y = ((y / Height) - 0.5f) * 0.834;

   float speedX = (m_x - oldX);
   float speedY = (m_y - oldY);

   if( mouse_button == GLUT_LEFT_BUTTON )
   {
      float coff = 5.5f;
      angle_X += speedX * coff;
      angle_Y += speedY * coff;
   }

   oldX = m_x;
   oldY = m_y;

}

void processNormalKeys(unsigned char key, int x, int y)
{

    switch(key)
    {
      case 27 : exit(0); break;
      case ' ' : Pause = !Pause; break;
    }

}


void processSpecialKeys(int key, int x, int y) {

    switch(key)
    {
    case GLUT_KEY_UP:
        shape0->mPosition = shape0->mPosition + Vector3(0,0,0.1);
    break;
    case GLUT_KEY_DOWN:
        shape0->mPosition = shape0->mPosition - Vector3(0,0,0.1);
    break;
    case GLUT_KEY_LEFT:
        shape0->mPosition = shape0->mPosition + Vector3(0.1,0,0);
    break;
    case GLUT_KEY_RIGHT:
        shape0->mPosition = shape0->mPosition - Vector3(0.1,0,0);
    break;
    }
}




// Main function
int main(int argc, char** argv)
{
    // Initialize GLUT
    glutInit(&argc, argv);

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    // Initialize the size of the GLUT windows
    glutInitWindowSize( Width ,  Height );
    // Initialize the position of the GLUT windows
    glutInitWindowPosition( 0 , 0 );
    // Initialize the create window
    glutCreateWindow("Demo: GJK-EPA");

    initCamera();

    glutDisplayFunc(display);

    // Glut Idle function that is continuously called
    glutReshapeFunc(reshape);
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMotion);

    // here are the new entries
    glutKeyboardFunc(processNormalKeys);
    glutSpecialFunc(processSpecialKeys);


    glutIdleFunc(UpdateTime);

    // Glut main looop
    glutMainLoop();


    /**/
    return 0;
}
