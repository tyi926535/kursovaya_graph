#include "Render.h"

#include <sstream>
#include <iostream>

#include <windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>
#include <glut.h>
#include "MyOGL.h"

#include "Camera.h"
#include "Light.h"
#include "Primitives.h"

#include "GUItextRectangle.h"
#include <math.h>          // Заголовочный файл для математической библиотеки
#include <stdio.h>         // Заголовочный файл для стандартного ввода/вывода
#include <stdlib.h>        // Заголовочный файл для стандартной библиотеки
#include <gl\glu.h>        // Заголовочный файл для библиотеки GLu32

//#include <glm/glm.hpp>
//#include <glm/gtc/matrix_transform.hpp>
//#include <glm/gtc/type_ptr.hpp>



bool textureMode = true;
bool lightMode = true;
GLfloat rtri;           // Угол для треугольник
GLfloat rquad;
GLfloat efp;
GLfloat tyu;
GLfloat h1[] = {0,7,14,7,0 ,-7,-14,-7,0};
GLfloat h2[] = { 0,7,14,7,0 ,-7,-14,-7,0 };
GLfloat h3[] = { 0,-4, -8, -4, 0,4,8,4,0};

Vector3 points[4][4];
Vector3* dragablePoint = 0;
void anim();
//класс для настройки камеры
class CustomCamera : public Camera
{
public:
	//дистанция камеры
	double camDist;
	//углы поворота камеры
	double fi1, fi2;

	
	//значния масеры по умолчанию
	CustomCamera()
	{
		camDist = 15;
		fi1 = 1;
		fi2 = 1;
	}

	
	//считает позицию камеры, исходя из углов поворота, вызывается движком
	void SetUpCamera()
	{
		//отвечает за поворот камеры мышкой
		lookPoint.setCoords(0, 0, 0);

		pos.setCoords(camDist*cos(fi2)*cos(fi1),
			camDist*cos(fi2)*sin(fi1),
			camDist*sin(fi2));

		if (cos(fi2) <= 0)
			normal.setCoords(0, 0, -1);
		else
			normal.setCoords(0, 0, 1);

		LookAt();
	}

	void CustomCamera::LookAt()
	{
		//функция настройки камеры
		gluLookAt(pos.X(), pos.Y(), pos.Z(), lookPoint.X(), lookPoint.Y(), lookPoint.Z(), normal.X(), normal.Y(), normal.Z());
	}



}  camera;   //создаем объект камеры


//Класс для настройки света
class CustomLight : public Light
{
public:
	CustomLight()
	{
		//начальная позиция света
		pos = Vector3(1, 1, 3);
	}

	
	//рисует сферу и линии под источником света, вызывается движком
	void  DrawLightGhismo()
	{
		glDisable(GL_LIGHTING);

		
		glColor3d(0.9, 0.8, 0);
		Sphere s;
		s.pos = pos;
		s.scale = s.scale*0.08;
		s.Show();
		
		if (OpenGL::isKeyPressed('G'))
		{
			glColor3d(0, 0, 0);
			//линия от источника света до окружности
			glBegin(GL_LINES);
			glVertex3d(pos.X(), pos.Y(), pos.Z());
			glVertex3d(pos.X(), pos.Y(), 0);
			glEnd();

			//рисуем окруность
			Circle c;
			c.pos.setCoords(pos.X(), pos.Y(), 0);
			c.scale = c.scale*1.5;
			c.Show();
		}

	}

	void SetUpLight()
	{
		GLfloat amb[] = { 0.2, 0.2, 0.2, 0 };
		GLfloat dif[] = { 1.0, 1.0, 1.0, 0 };
		GLfloat spec[] = { .7, .7, .7, 0 };
		GLfloat position[] = { pos.X(), pos.Y(), pos.Z(), 1. };

		// параметры источника света
		glLightfv(GL_LIGHT0, GL_POSITION, position);
		// характеристики излучаемого света
		// фоновое освещение (рассеянный свет)
		glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
		// диффузная составляющая света
		glLightfv(GL_LIGHT0, GL_DIFFUSE, dif);
		// зеркально отражаемая составляющая света
		glLightfv(GL_LIGHT0, GL_SPECULAR, spec);

		glEnable(GL_LIGHT0);
	}


} light;  //создаем источник света




//старые координаты мыши
int mouseX = 0, mouseY = 0;

void mouseEvent(OpenGL *ogl, int mX, int mY)
{
	int dx = mouseX - mX;
	int dy = mouseY - mY;
	mouseX = mX;
	mouseY = mY;

	//меняем углы камеры при нажатой левой кнопке мыши
	if (OpenGL::isKeyPressed(VK_RBUTTON))
	{
		camera.fi1 += 0.01*dx;
		camera.fi2 += -0.01*dy;
	}

	
	//двигаем свет по плоскости, в точку где мышь
	if (OpenGL::isKeyPressed('G') && !OpenGL::isKeyPressed(VK_LBUTTON))
	{
		LPPOINT POINT = new tagPOINT();
		GetCursorPos(POINT);
		ScreenToClient(ogl->getHwnd(), POINT);
		POINT->y = ogl->getHeight() - POINT->y;

		Ray r = camera.getLookRay(POINT->x, POINT->y);

		double z = light.pos.Z();

		double k = 0, x = 0, y = 0;
		if (r.direction.Z() == 0)
			k = 0;
		else
			k = (z - r.origin.Z()) / r.direction.Z();

		x = k*r.direction.X() + r.origin.X();
		y = k*r.direction.Y() + r.origin.Y();

		light.pos = Vector3(x, y, z);
	}

	if (OpenGL::isKeyPressed('G') && OpenGL::isKeyPressed(VK_LBUTTON))
	{
		light.pos = light.pos + Vector3(0, 0, 0.02*dy);
	}

	if (dragablePoint && OpenGL::isKeyPressed(VK_LBUTTON)) {
		(*dragablePoint) = (*dragablePoint) + Vector3(0.02*dx, 0, 0.02 * dy);
	}

	
}

void mouseWheelEvent(OpenGL *ogl, int delta)
{

	if (delta < 0 && camera.camDist <= 1)
		return;
	if (delta > 0 && camera.camDist >= 100)
		return;

	camera.camDist += 0.01*delta;

}

void keyDownEvent(OpenGL *ogl, int key)
{
	if (key == 'L')
	{
		lightMode = !lightMode;
	}

	if (key == 'T')
	{
		textureMode = !textureMode;
	}

	if (key == 'R')
	{
		camera.fi1 = 1;
		camera.fi2 = 1;
		camera.camDist = 15;

		light.pos = Vector3(1, 1, 3);
	}

	if (key == 'F')
	{
		light.pos = camera.pos;
	}
	if (key == 'P') {
		//anim();
	}
	if (key == VK_LBUTTON) {
		double mvMatrix[16], prMatrix[16];
		int viewPort[4];
		double x, y, z;
		glGetDoublev(GL_MODELVIEW_MATRIX, mvMatrix);
		glGetDoublev(GL_PROJECTION_MATRIX, prMatrix);
		glGetIntegerv(GL_VIEWPORT, viewPort);

		LPPOINT POINT = new tagPOINT();
		GetCursorPos(POINT);
		ScreenToClient(ogl->getHwnd(), POINT);
		POINT->y = ogl->getHeight() - POINT->y;

		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				int res = gluProject(points[i][j].X(), points[i][j].Y(), points[i][j].Z(), mvMatrix, prMatrix, viewPort, &x, &y, &z);
				if ((x - POINT->x) * (x - POINT->x) + (y - POINT->y) * (y - POINT->y) < 49) {
					dragablePoint = &points[i][j];
					break;
				}
			}
		}
		delete POINT;

	}

}

void keyUpEvent(OpenGL *ogl, int key)
{
	if (key == VK_LBUTTON) {
		dragablePoint = 0;
	}
}



GLuint texId;
GLuint tex2;
//GLuint texId2;

//выполняется перед первым рендером
void initRender(OpenGL *ogl)
{
	//настройка текстур

	//4 байта на хранение пикселя
	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

	//настройка режима наложения текстур
	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	//включаем текстуры
	glEnable(GL_TEXTURE_2D);
	

	//массив трехбайтных элементов  (R G B)
	RGBTRIPLE *texarray;
	//RGBTRIPLE* texarray2;

	//массив символов, (высота*ширина*4      4, потомучто   выше, мы указали использовать по 4 байта на пиксель текстуры - R G B A)
	char *texCharArray;
	//char* texCharArray2;
	int texW, texH;
	//int texW2, texH2;
	OpenGL::LoadBMP("text1.bmp", &texW, &texH, &texarray);
	OpenGL::RGBtoChar(texarray, texW, texH, &texCharArray);

	//OpenGL::LoadBMP("text2.bmp", &texW2, &texH2, &texarray2);
	//OpenGL::RGBtoChar(texarray2, texW2, texH2, &texCharArray2);
	
	
	//генерируем ИД для текстуры
	glGenTextures(1, &texId);
	//glGenTextures(2, &texId[2]);

	//биндим айдишник, все что будет происходить с текстурой, будте происходить по этому ИД
	glBindTexture(GL_TEXTURE_2D, texId);
	//glBindTexture(GL_TEXTURE_2D, texId[2]);

	//загружаем текстуру в видеопямять, в оперативке нам больше  она не нужна
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texW, texH, 0, GL_RGBA, GL_UNSIGNED_BYTE, texCharArray);
	//glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, texW2, texH2, 0, GL_RGBA, GL_UNSIGNED_BYTE, texCharArray2);
	
	//отчистка памяти
	free(texCharArray);
	free(texarray);
	//free(texCharArray2);
	//free(texarray2);

	//наводим шмон
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST); 
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	/*glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);*/
	//камеру и свет привязываем к "движку"
	ogl->mainCamera = &camera;
	ogl->mainLight = &light;

	// нормализация нормалей : их длины будет равна 1
	glEnable(GL_NORMALIZE);

	// устранение ступенчатости для линий
	glEnable(GL_LINE_SMOOTH); 


	//   задать параметры освещения
	//  параметр GL_LIGHT_MODEL_TWO_SIDE - 
	//                0 -  лицевые и изнаночные рисуются одинаково(по умолчанию), 
	//                1 - лицевые и изнаночные обрабатываются разными режимами       
	//                соответственно лицевым и изнаночным свойствам материалов.    
	//  параметр GL_LIGHT_MODEL_AMBIENT - задать фоновое освещение, 
	//                не зависящее от сточников
	// по умолчанию (0.2, 0.2, 0.2, 1.0)

	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 0);

	camera.fi1 = -1.3;
	camera.fi2 = 0.8;


	srand(time(0));

	
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			points[i][j] = Vector3(10.0 * i / 4 - 1, 10.0 * j / 4 - 1, rand() % 5 - 2.5);
		}
	}

}

double fact(int num) {
	int r = 1;
	for (int i = 2; i <= num; i++)
		r *= i;
	return r;
}

double bern(int i, int n, double u) {
	return 1.0 * fact(n) * pow(u, i) * pow(1-u, n - i) / fact(i) / fact(n - i);
}
double* norm3(double a[])
{
	double x, y, z;
	x = a[1];
	y = a[2];
	z = a[3];
	double l = sqrt(x * x + y * y + z * z);
	double n[] = { x / l, y / l, z / l };
	return a;
}
double* norm2(double A[], double B[])
{
	double x, y, z;
	double a[] = { A[0] - B[0],A[1] - B[1],A[2] - B[2] };
	x = a[1];
	y = a[2];
	z = a[3];
	double l = sqrt(x * x + y * y + z * z);
	double n[] = { x / l, y / l, z / l };
	return n;
}
double* norm(double A[], double B[], double C[])
{
	double x, y, z;
	double a[] = { A[0] - B[0],A[1] - B[1],A[2] - B[2] };
	double c[] = { C[0] - B[0],C[1] - B[1],C[2] - B[2] };
	x = a[1] * c[2] - c[1] * a[2];
	y = c[0] * a[2] - a[0] * c[2];
	z = a[0] * c[1] - c[0] * a[1];
	double l = sqrt(x * x + y * y + z * z);
	double n[] = { x / l, y / l, z / l };
	return n;
}
double f(double a, double b, double c, double t)
{
	return a * (1 - t) * (1 - t) + 2 * b * t * (1 - t) + c * t * t;
}
void shester() {
	double xn, yn, xn1, yn1, xn2 = 0, yn2 = 0, x1 = 1, x2 = 1, x3 = 1, x4 = 1, y1 = 0, y2 = 0, y3 = 0, y4 = 0, y5 = 0, x5 = 1, z1 = -4, z2 = 0, z3, z4, fi2 = 0, fi1;
	double r = 1; x1 = 1; y1 = 0; fi2 = 0; int rt = -3, h = 0; z1 = 0; double r2 = 1.5, r3 = 2;
	while (rt < 1) {
		xn = r * cos(fi2);
		yn = r * sin(fi2);
		if (h == 0) {
			h = 1; x2 = r2 * cos(fi2); y2 = r2 * sin(fi2);
		}
		if (h == 1) {
			h = 0;
			x5 = r3 * cos(fi2);
			y5 = r3 * sin(fi2);
		}
		if (xn >= 0.9) { rt++; }
		{
			double A1[] = { x1,y1,z1 };
			double B1[] = { x1,y1,z1 - 0.5 };
			double C1[] = { xn,yn,z1 - 0.5 };
			double D1[] = { xn,yn,z1 };
			double* N = norm(A1, B1, C1);
			glBegin(GL_QUADS);
			glNormal3dv(N);
			glVertex3dv(B1);
			glVertex3dv(A1);
			glVertex3dv(D1);
			glVertex3dv(C1);
			glEnd();
		}
		{
			double A1[] = { x1,y1,z1 };
			double B1[] = { xn,yn,z1 };
			double C1[] = { x2,y2,z1 };
			double D1[] = { x4,y4,z1 };
			double* N = norm(A1, B1, C1);
			glBegin(GL_QUADS);
			glNormal3dv(N);
			glVertex3dv(B1);
			glVertex3dv(A1);
			glVertex3dv(D1);
			glVertex3dv(C1);
			glEnd();
		}
		{
			double A1[] = { x1,y1,z1 - 0.5 };
			double B1[] = { xn,yn,z1 - 0.5 };
			double C1[] = { x2,y2,z1 - 0.5 };
			double D1[] = { x4,y4,z1 - 0.5 };
			double* N = norm(A1, B1, C1);
			glBegin(GL_QUADS);
			glNormal3dv(N);
			glVertex3dv(B1);
			glVertex3dv(A1);
			glVertex3dv(D1);
			glVertex3dv(C1);
			glEnd();
		}
		{
			double A1[] = { x5,y5,z1 };
			double B1[] = { x4,y4,z1 };
			double C1[] = { x5,y5,z1 - 0.5 };
			double D1[] = { x4,y4,z1 - 0.5 };
			double* N = norm( A1,B1, C1);
			glBegin(GL_QUADS);
			glNormal3dv(N);
			glVertex3dv(B1);
			glVertex3dv(A1);
			glVertex3dv(D1);
			glVertex3dv(C1);
			glEnd();
		}
		{
			double A1[] = { x2,y2,z1 };
			double B1[] = { x5,y5,z1 };
			double C1[] = { x2,y2,z1 - 0.5 };
			double D1[] = { x5,y5,z1 - 0.5 };
			double* N = norm(A1, B1, C1);
			glBegin(GL_QUADS);
			glNormal3dv(N);
			glVertex3dv(B1);
			glVertex3dv(A1);
			glVertex3dv(D1);
			glVertex3dv(C1);
			glEnd();
		}
		{
			double A1[] = { x2,y2,z1 };
			double B1[] = { x5,y5,z1 };
			double C1[] = { x4,y4,z1 };
			double* N = norm(A1, B1, C1);
			glBegin(GL_TRIANGLES);
			glNormal3dv(N);
			glVertex3dv(B1);
			glVertex3dv(A1);
			glVertex3dv(C1);
			glEnd();
		}
		{
			double A1[] = { x2,y2,z1 - 0.5 };
			double B1[] = { x5,y5,z1 - 0.5 };
			double C1[] = { x4,y4,z1 - 0.5 };
			double* N = norm(A1, B1, C1);
			glBegin(GL_TRIANGLES);
			glNormal3dv(N);
			glVertex3dv(B1);
			glVertex3dv(A1);
			glVertex3dv(C1);
			glEnd();

		}

		fi2 = fi2 + 0.35; x1 = xn; y1 = yn; x4 = x2; y4 = y2;
	}

}
void vert(GLfloat tyu)
{
	glPushMatrix();
	glRotatef(tyu, 0.0f, 0.0f, 1.0f);
	shester();
	glPopMatrix();
}
void anim() {
	double LIN1[3]; efp = 0;
	glTranslatef(0, 0, -3);
	
	for (int i = 0; i < 9; i+2) {
		for (double t = 0; t <= 1.0001; t += 0.2)
		{glPushMatrix();
			LIN1[0] = f(h1[i], h2[i], h3[i], t);
			LIN1[1] = f(h1[i+1], h2[i+1], h3[i+1], t);
			LIN1[2] = f(h1[i+2], h2[i+2], h3[i+2], t);
			glTranslatef(LIN1[0], LIN1[1],LIN1[2]);
		glRotatef(efp, 0.0f, 0.0f, 1.0f);
		shester();
		efp = efp + 0.1;
		}
		glPopMatrix();

	}
	

}

void line(double A1[], double B1[], double C1[])
{
	double LIN1[3];
	glBegin(GL_LINE_STRIP);
	for (double t = 0; t <= 1.0001; t += 0.02)
	{
		
		LIN1[0] = f(A1[0], B1[0], C1[0], t);
		LIN1[1] = f(A1[1], B1[1], C1[1], t);
		LIN1[2] = f(A1[2], B1[2], C1[2], t);
		double* N = norm3(LIN1);
		glNormal3dv(N);
		glVertex3dv(LIN1);
		
	}
	glEnd();
}
	
void strelkiM() {
		glLineWidth(2);
		{double A1[] = { 0, 0, 0.3 };
		double B1[] = { 6, 0, 0.3 };
		double* N = norm2(A1, B1);
		glBegin(GL_LINES);
		glNormal3dv(N);
		glVertex3d(0, 0, 0.3);
		glVertex3d(6, 0, 0.3);
		glEnd(); }
	}
void strelkiCh() {
		glLineWidth(5);
		{double A1[] = { 0, 0, 0.3 };
		double B1[] = { 4, 0, 0.3 };
		double* N = norm2(A1, B1);
		glBegin(GL_LINES);
		glNormal3dv(N);
		glVertex3d(0, 0, 0.3);
		glVertex3d(4, 0, 0.3);
		glEnd(); }
	}



double increment_delta_time = 0.005;
double t_max = 0;
void Render(OpenGL *ogl)
{
		glDisable(GL_TEXTURE_2D);
		//glDisable(GL_TEXTURE_1D);
		glDisable(GL_LIGHTING);

		glEnable(GL_DEPTH_TEST);
		if (textureMode)
		{
			glEnable(GL_TEXTURE_2D); //glEnable(GL_TEXTURE_1D);
		}

		if (lightMode)
			glEnable(GL_LIGHTING);


		//альфаналожение
		//glEnable(GL_BLEND);
		//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


		//настройка материала
		GLfloat amb[] = { 0.2, 0.2, 0.1, 1. };
		GLfloat dif[] = { 0.4, 0.65, 0.5, 1. };
		GLfloat spec[] = { 0.9, 0.8, 0.3, 1. };
		GLfloat sh = 0.1f * 256;


		//фоновая
		glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
		//дифузная
		glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
		//зеркальная
		glMaterialfv(GL_FRONT, GL_SPECULAR, spec); \
			//размер блика
			glMaterialf(GL_FRONT, GL_SHININESS, sh);
	
	//чтоб было красиво, без квадратиков (сглаживание освещения)
	glShadeModel(GL_SMOOTH);
	//===================================
	//Прогать тут  

	
		glPushMatrix();
		glRotatef(90, 1, 0, 0);
		//glRotatef(180, 0, 1, 0);

		{
			//Циферблат
			
			{
				glColor3d(1, 1, 1);
				double x = 0, y = 0, z = 0;
				float f, l, u, n;
				double xn, yn, x1 = 8, x2 = 0, x3, x4, y1 = 0, y2 = 0, y3, y4, z1 = 0, z2 = 0, z3, z4, fi2 = 0, fi1;


				int th = -3;
				double r = 8;
				while (th < 1) {
					xn = r * cos(fi2);
					yn = r * sin(fi2);
					x2 = xn; y2 = yn;
					{
						if (xn >= 7.8) { th++; }
						double A1[] = { 0,0,0 };
						double B1[] = { x2,y2,z1 };
						double C1[] = { x1,y1,z1 };
						double* N = norm(A1, B1, C1);
						glBegin(GL_TRIANGLES);
						glNormal3dv(N);
						glColor3d(1, 1, 1);
						glVertex3dv(A1);
						glVertex3dv(B1);
						glVertex3dv(C1);
						glEnd();
					}
					fi2 = fi2 + 0.2; x1 = x2; y1 = y2;
				}

				z1 = 0.2;
				glLineWidth(1);
				glColor3d(0, 0, 0);
				{
					double A1[] = { 4,0,z1 };
					double B1[] = { 3.7,3.7,z1 };
					double C1[] = { 0,4,z1 };
					line(A1, B1, C1);
				}
				{
					double A1[] = { 4,0,z1 };
					double B1[] = { 3.7,-3.7,z1 };
					double C1[] = { 0,-4,z1 };
					line(A1, B1, C1);
				}
				{
					double A1[] = { -4,0,z1 };
					double B1[] = { -3.7,-3.7,z1 };
					double C1[] = { 0,-4,z1 };
					line(A1, B1, C1);
				}
				{
					double A1[] = { -4,0,z1 };
					double B1[] = { -3.7,3.7,z1 };
					double C1[] = { 0,4,z1 };
					line(A1, B1, C1);
				}
				int rty = 0;
				{
					x1 = 4; x2 = 0; y1 = 0; y2 = 0;  z2 = 0; fi2 = 0;
					th = -1;
					r = 4, xn = 4;
					while (th < 1) {
						xn = r * cos(fi2);
						yn = r * sin(fi2);
						if (xn > 3.8) { th++; }
						if (th < 1)
						{
							double A1[] = { xn,yn,z1 };
							glPointSize(6);
							double* N = norm3(A1);
							glBegin(GL_POINTS);
							glNormal3dv(N);
							glVertex3d(xn, yn, z1);
							glEnd();
						}
						fi2 = fi2 + 0.525;

					}
				}
				glLineWidth(2);
				glColor3d(0, 0, 0);
				{
					double A1[] = { 7,0,z1 };
					double B1[] = { 6.5,6.5,z1 };
					double C1[] = { 0,7,z1 };
					line(A1, B1, C1);
				}
				{
					double A1[] = { 7,0,z1 };
					double B1[] = { 6.5,-6.5,z1 };
					double C1[] = { 0,-7,z1 };
					line(A1, B1, C1);
				}
				{
					double A1[] = { -7,0,z1 };
					double B1[] = { -6.5,-6.5,z1 };
					double C1[] = { 0,-7,z1 };
					line(A1, B1, C1);
				}
				{
					double A1[] = { -7,0,z1 };
					double B1[] = { -6.5,6.5,z1 };
					double C1[] = { 0,7,z1 };
					line(A1, B1, C1);
				}
				glColor3d(0, 0, 0);
				{
					double A1[] = { 7.5,0,z1 };
					double B1[] = { 6.9,6.9,z1 };
					double C1[] = { 0,7.5,z1 };
					line(A1, B1, C1);
				}
				{
					double A1[] = { 7.5,0,z1 };
					double B1[] = { 6.9,-6.9,z1 };
					double C1[] = { 0,-7.5,z1 };
					line(A1, B1, C1);
				}
				{
					double A1[] = { -7.5,0,z1 };
					double B1[] = { -6.9,-6.9,z1 };
					double C1[] = { 0,-7.5,z1 };
					line(A1, B1, C1);
				}
				{
					double A1[] = { -7.5,0,z1 };
					double B1[] = { -6.9,6.9,z1 };
					double C1[] = { 0,7.5,z1 };
					line(A1, B1, C1);
				}

				{
					x1 = 7; x2 = 7.5; y1 = 0; y2 = 0;  z2 = 0; fi2 = 0;
					th = -3;
					r = 7; double r1 = 7.5, xn1, yn1;
					while (fi2 < 360) {
						xn = r * cos(fi2);
						yn = r * sin(fi2);
						xn1 = r1 * cos(fi2);
						yn1 = r1 * sin(fi2);
						if (xn > 6.9) { th++; }
						if (th < 1)
						{
							double A1[] = { xn, yn, z1 };
							double B1[] = { xn1, yn1, z1 };
							double* N = norm2(A1, B1);
							glPointSize(6);
							glBegin(GL_LINES);
							glNormal3dv(N);
							glVertex3d(xn, yn, z1);
							glVertex3d(xn1, yn1, z1);
							glEnd();
						}
						fi2 = fi2 + 0.10525;
					}
				}
				glLineWidth(4);
				
				{
					{
						double A1[] = { -0.4, 5.7, z1 };
						double B1[] = { -0.1, 5, z1 };
						double* N = norm2(A1, B1);
						glBegin(GL_LINES);
						glNormal3dv(N);
						glVertex3d(-0.4, 5.7, z1);
						glVertex3d(-0.1, 5, z1);
						glEnd(); }

					{
						double A1[] = { -0.4, 5, z1 };
						double B1[] = { -0.1, 5.7, z1 };
						double* N = norm2(A1, B1);
						glBegin(GL_LINES);
						glNormal3dv(N);
						glVertex3d(-0.4, 5, z1);
						glVertex3d(-0.1, 5.7, z1);
						glEnd(); }

					{double A1[] = { 0.1, 5, z1 };
					double B1[] = { 0.1, 5.7, z1 };
					double* N = norm2(A1, B1);
					glBegin(GL_LINES);
					glNormal3dv(N);
					glVertex3d(0.1, 5, z1);
					glVertex3d(0.1, 5.7, z1);
					glEnd(); }


					{double A1[] = { 0.3, 5, z1 };
					double B1[] = { 0.3, 5.7, z1 };
					double* N = norm2(A1, B1);
					glBegin(GL_LINES);
					glNormal3dv(N);
					glVertex3d(0.3, 5, z1);
					glVertex3d(0.3, 5.7, z1);
					glEnd(); }

				}
				{

					{double A1[] = { 5.5, 0.35, z1 };
					double B1[] = { 5.5, -0.35, z1 };
					double* N = norm2(A1, B1);
					glBegin(GL_LINES);
					glNormal3dv(N);
					glVertex3d(5.5, 0.35, z1);
					glVertex3d(5.5, -0.35, z1);
					glEnd(); }
					{double A1[] = { 5.7, 0.35, z1 };
					double B1[] = { 5.7, -0.35, z1 };
					double* N = norm2(A1, B1);
					glNormal3dv(N);
					glBegin(GL_LINES);
					glVertex3d(5.7, 0.35, z1);
					glVertex3d(5.7, -0.35, z1);
					glEnd(); }
					{double A1[] = { 5.3, 0.35, z1 };
					double B1[] = { 5.3, -0.35, z1 };
					double* N = norm2(A1, B1);
					glBegin(GL_LINES);
					glNormal3dv(N);
					glVertex3d(5.3, 0.35, z1);
					glVertex3d(5.3, -0.35, z1);
					glEnd(); }
				}
				{
					{double A1[] = { -0.2, -5.7, z1 };
					double B1[] = { -0.2, -5, z1 };
					double* N = norm2(A1, B1);
					glBegin(GL_LINES);
					glNormal3dv(N);
					glVertex3d(-0.2, -5.7, z1);
					glVertex3d(-0.2, -5, z1);
					glEnd(); }
					{double A1[] = { 0.2, -5.7, z1 };
					double B1[] = { 0, -5, z1 };
					double* N = norm2(A1, B1);
					glBegin(GL_LINES);
					glNormal3dv(N);
					glVertex3d(0.2, -5.7, z1);
					glVertex3d(0, -5, z1);
					glEnd(); }
					{double A1[] = { 0.2, -5.7, z1 };
					double B1[] = { 0.4, -5, z1 };
					double* N = norm2(A1, B1);
					glBegin(GL_LINES);
					glNormal3dv(N);
					glVertex3d(0.2, -5.7, z1);
					glVertex3d(0.4, -5, z1);
					glEnd(); }
				}
				{
					{double A1[] = { 0.2, -5.7, z1 };
					double B1[] = { 0.4, -5, z1 };
					double* N = norm2(A1, B1);
					glBegin(GL_LINES);
					glNormal3dv(N);
					glVertex3d(0.2, -5.7, z1);
					glVertex3d(0.4, -5, z1);
					glEnd(); }
					{double A1[] = { -5.4, -0.35, z1 };
					double B1[] = { -5.7, 0.35, z1 };
					double* N = norm2(A1, B1);
					glBegin(GL_LINES);
					glNormal3dv(N);
					glVertex3d(-5.4, -0.35, z1);
					glVertex3d(-5.7, 0.35, z1);
					glEnd(); }
					{double A1[] = { -5.7, -0.35, z1 };
					double B1[] = { -5.4, 0.35, z1 };
					double* N = norm2(A1, B1);
					glBegin(GL_LINES);
					glNormal3dv(N);
					glVertex3d(-5.7, -0.35, z1);
					glVertex3d(-5.4, 0.35, z1);
					glEnd(); }
				}
			}
			glColor3d(1, 0, 0);
			double xn, yn, xn1, yn1, xn2 = 0, yn2 = 0, x1 = 8, x2 = 0, x3 = 8, x4 = 9, y1 = 0, y2 = 0, y3 = 0, y4 = 0, y5 = 0, x5 = 8.7, z1 = 0, z2 = 0, z3, z4, fi2 = 0, fi1;
			int th = -5;
			double r = 8, r1 = 9, r2 = 8.7;
			while (th < 1) {
				xn = r * cos(fi2);
				yn = r * sin(fi2);
				xn2 = r2 * cos(fi2);
				yn2 = r2 * sin(fi2);
				xn1 = r1 * cos(fi2);
				yn1 = r1 * sin(fi2);
				x2 = xn1; y2 = yn1; x1 = xn; y1 = yn;
				{
					double A1[] = { x1,y1,z1 };
					double B1[] = { xn2,yn2,z1 + 3 };
					double C1[] = { x2,y2,z1 };
					double A2[] = { x3,y3,z1 };
					double B2[] = { x5,y5,z1 + 3 };
					double C2[] = { x4,y4,z1 };
					double LIN1[3], LIN2[3];
					if (xn > 7.8) { th++; }
					glBegin(GL_LINE_STRIP);
					for (double t = 0; t <= 1.0001; t += 0.001)
					{
						LIN1[0] = f(A1[0], B1[0], C1[0], t);
						LIN1[1] = f(A1[1], B1[1], C1[1], t);
						LIN1[2] = f(A1[2], B1[2], C1[2], t);
						glNormal3dv(LIN1);
						glVertex3dv(LIN1);
						LIN2[0] = f(A2[0], B2[0], C2[0], t);
						LIN2[1] = f(A2[1], B2[1], C2[1], t);
						LIN2[2] = f(A2[2], B2[2], C2[2], t);
						glNormal3dv(LIN2);
						glVertex3dv(LIN2);
					}
					glEnd();
				}
				fi2 = fi2 + 0.1; x4 = x2; y4 = y2; x3 = x1; y3 = y1; x5 = xn2; y5 = yn2;
			}
			th = 0;
			r = 9; x1 = 9; y1 = 0; fi2 = 0; int rt = -2;
			while (rt < 1) {
				xn = r * cos(fi2);
				yn = r * sin(fi2);
				{
					double A1[] = { x1,y1,z1 };
					double B1[] = { x1,y1,z1 - 6 };
					double C1[] = { xn,yn,z1 - 6 };
					double D1[] = { xn,yn,z1 };
					if (xn >= 8.8) { rt++; }
					double* N = norm(B1,A1,  C1);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3dv(B1);
					glVertex3dv(A1);
					glVertex3dv(D1);
					glVertex3dv(C1);
					glEnd();

				}
				fi2 = fi2 + 0.225; x1 = xn; y1 = yn;
			}

			Vector3 p1, p2, p3, p4;
			Vector3 prev1, prev2;

			glColor3d(1, 0, 0);
			p1 = { -4, -8, -0.5 };
			p2 = { -4, -8, -5.5 };
			p3 = { -5, -10.5, -5.5 };
			p4 = { -5, -10.5, -0.5 };

			{
				double A[] = { -4, -8, -0.5 };
				double B[] = { -4, -8, -5.5 };
				double C[] = { -5, -10.5, -5.5 };
				double* N = norm(A, B, C);
				glBegin(GL_QUADS);
				glNormal3dv(N);
				glVertex3d(-4, -8, -0.5);
				glVertex3d(-4, -8, -5.5);
				glVertex3d(-5, -10.5, -5.5);
				glVertex3d(-5, -10.5, -0.5);
				glEnd(); }
			{
				double A[] = { -4, -8, -0.5 };
				double B[] = { -5, -10.5, -0.5 };
				double C[] = { -7, -10.5, -0.5 };
				double* N = norm(A, B, C);
				glBegin(GL_QUADS);
				glNormal3dv(N);
				glVertex3d(-4, -8, -0.5);
				glVertex3d(-5, -10.5, -0.5);
				glVertex3d(-7, -10.5, -0.5);
				glVertex3d(-6, -7, -0.5);
				glEnd(); }
			{
				double A[] = { -4, -8, -5.5 };
				double B[] = { -5, -10.5, -5.5 };
				double C[] = { -7, -10.5, -5.5 };
				double* N = norm(B,A,  C);
				glBegin(GL_QUADS);
				glNormal3dv(N);
				glVertex3d(-4, -8, -5.5);
				glVertex3d(-5, -10.5, -5.5);
				glVertex3d(-7, -10.5, -5.5);
				glVertex3d(-6, -7, -5.5);
				glEnd(); }
			{
				double A[] = { -7, -10.5, -0.5 };
				double B[] = { -6, -7, -0.5 };
				double C[] = { -6, -7, -5.5 };
				double* N = norm(B,A,  C);
				glBegin(GL_QUADS);
				glNormal3dv(N);
				glVertex3d(-7, -10.5, -0.5);
				glVertex3d(-6, -7, -0.5);
				glVertex3d(-6, -7, -5.5);
				glVertex3d(-7, -10.5, -5.5);
				glEnd(); }
			{
				double A[] = { -7, -10.5, -0.5 };
				double B[] = { -7, -10.5, -5.5 };
				double C[] = { -5, -10.5, -5.5 };
				double* N = norm( B,A, C);
				glBegin(GL_QUADS);
				glNormal3dv(N);
				glVertex3d(-7, -10.5, -0.5);
				glVertex3d(-7, -10.5, -5.5);
				glVertex3d(-5, -10.5, -5.5);
				glVertex3d(-5, -10.5, -0.5);
				glEnd(); }
			{
				double A[] = { -4, -8, -0.5 };
				double B[] = { -4, -8, -5.5 };
				double C[] = {-5, -10.5, -5.5};
				double* N = norm(A, B, C);
				glBegin(GL_QUADS);
				glNormal3dv(N);
				glVertex3d(-4, -8, -0.5);
				glVertex3d(-4, -8, -5.5);
				glVertex3d(-5, -10.5, -5.5);
				glVertex3d(-5, -10.5, -0.5);
				glEnd(); }
			{
				double A[] = { -4, -8, -0.5 };
				double B[] = { -5, -10.5, -0.5 };
				double C[] = { -7, -10.5, -0.5 };
				double* N = norm(A, B, C);
				glBegin(GL_QUADS);
				glNormal3dv(N);
				
				glVertex3d(-4, -8, -0.5);
				glVertex3d(-5, -10.5, -0.5);
				glVertex3d(-7, -10.5, -0.5);
				glVertex3d(-6, -7, -0.5);
				glEnd(); }
			{
				double A[] = { -4, -8, -5.5 };
				double B[] = { -5, -10.5, -5.5 };
				double C[] = { -7, -10.5, -5.5 };
				double* N = norm(A, B, C);
				glBegin(GL_QUADS);
				glNormal3dv(N);
				
				glVertex3d(-4, -8, -5.5);
				glVertex3d(-5, -10.5, -5.5);
				glVertex3d(-7, -10.5, -5.5);
				glVertex3d(-6, -7, -5.5);
				glEnd(); }
			{
				double A[] = { -7, -10.5, -0.5 };
				double B[] = { -6, -7, -0.5 };
				double C[] = { -6, -7, -5.5 };
				double* N = norm(A, B, C);
				glBegin(GL_QUADS);
				glNormal3dv(N);
				glVertex3d(-7, -10.5, -0.5);
				glVertex3d(-6, -7, -0.5);
				glVertex3d(-6, -7, -5.5);
				glVertex3d(-7, -10.5, -5.5);
				glEnd(); }
			{
				double A[] = { -7, -10.5, -0.5 };
				double B[] = { -7, -10.5, -5.5 };
				double C[] = { -5, -10.5, -5.5 };
				double* N = norm(A, B, C);
				glBegin(GL_QUADS);
				glNormal3dv(N);
				glVertex3d(-7, -10.5, -0.5);
				glVertex3d(-7, -10.5, -5.5);
				glVertex3d(-5, -10.5, -5.5);
				glVertex3d(-5, -10.5, -0.5);
				glEnd(); }


			{
				double A[] = { 4, -8, -0.5 };
				double B[] = { 4, -8, -5.5 };
				double C[] = { 5, -10.5, -5.5 };
				double* N = norm(A, B, C);
				glBegin(GL_QUADS);
				glNormal3dv(N);
				glVertex3d(4, -8, -0.5);
				glVertex3d(4, -8, -5.5);
				glVertex3d(5, -10.5, -5.5);
				glVertex3d(5, -10.5, -0.5);
				glEnd(); }
			{
				double A[] = { -4, -8, -0.5 };
				double B[] = { -5, -10.5, -0.5 };
				double C[] = { -7, -10.5, -0.5 };
				double* N = norm(A, B, C);
				glBegin(GL_QUADS);
				glNormal3dv(N);
				glVertex3d(-4, -8, -0.5);
				glVertex3d(-5, -10.5, -0.5);
				glVertex3d(-7, -10.5, -0.5);
				glVertex3d(-6, -7, -0.5);
				glEnd(); }
			{
				double A[] = { 4, -8, -5.5 };
				double B[] = { 5, -10.5, -5.5 };
				double C[] = { 7, -10.5, -5.5 };
				double* N = norm(A, B, C);
				glBegin(GL_QUADS);
				glNormal3dv(N);
				glVertex3d(4, -8, -5.5);
				glVertex3d(5, -10.5, -5.5);
				glVertex3d(7, -10.5, -5.5);
				glVertex3d(6, -7, -5.5);
				glEnd(); }
			{
				double A[] = { 7, -10.5, -0.5 };
				double B[] = { 6, -7, -0.5 };
				double C[] = { 6, -7, -5.5 };
				double* N = norm(A, B, C);
				glBegin(GL_QUADS);
				glNormal3dv(N);
				glVertex3d(7, -10.5, -0.5);
				glVertex3d(6, -7, -0.5);
				glVertex3d(6, -7, -5.5);
				glVertex3d(7, -10.5, -5.5);
				glEnd(); }
			{
				double A[] = { 7, -10.5, -0.5 };
				double B[] = { 7, -10.5, -5.5 };
				double C[] = { 5, -10.5, -5.5 };
				double* N = norm(A, B, C);
				glBegin(GL_QUADS);
				glNormal3dv(N);
				glVertex3d(7, -10.5, -0.5);
				glVertex3d(7, -10.5, -5.5);
				glVertex3d(5, -10.5, -5.5);
				glVertex3d(5, -10.5, -0.5);
				glEnd(); }

		}
		
		{ {
				glColor3d(1, 1, 1);
				glBindTexture(GL_TEXTURE_2D, 1);
				glTexCoord2d(1.0f, 4.5f);

				GLdouble plane[4] = { -0.35f, 0.5f, 0.0f, 0.0f };
				GLUquadricObj* cylinder;
				cylinder = gluNewQuadric();
				GLUquadricObj* colp;
				colp = gluNewQuadric();


				glPushMatrix();       // сохраняем текущие координаты
				glEnableClientState(GL_TEXTURE_COORD_ARRAY);
				glTranslated(-6, 8, -2);  // сдвигаемся по оси Х на единицу 

				glEnable(GL_CLIP_PLANE0);
				glClipPlane(GL_CLIP_PLANE0, plane);
				glutSolidSphere(4.0f, 32, 32);
				glDisable(GL_CLIP_PLANE0);
				gluQuadricDrawStyle(cylinder, GLU_FILL);
				glRotatef(270, 1, 0, 0);
				glRotatef(320, 0, 1, 0);
				glTranslated(0.2, 0, -1.1);


				//////////////////////////////////////////////
				gluCylinder(cylinder, 0.4, 0.4, 5, 10, 10);

				glTranslated(0.2, 0, 5);
				gluCylinder(colp, 1.2, 0, 1, 10, 10);
				glPopMatrix();  // возвращаемся к старой системе координат

				glPushMatrix();       // сохраняем текущие координаты
				glTranslated(6, 8, -2);  // сдвигаемся по оси Х на единицу 
				glRotatef(180, 0, 1, 0);
				glEnable(GL_CLIP_PLANE0);
				glClipPlane(GL_CLIP_PLANE0, plane);
				glutSolidSphere(4.0f, 32, 32);
				glDisable(GL_CLIP_PLANE0);
				gluQuadricDrawStyle(cylinder, GLU_FILL);
				glRotatef(270, 1, 0, 0);
				glRotatef(320, 0, 1, 0);
				glTranslated(0.2, 0, -1.1);
				gluCylinder(cylinder, 0.4, 0.4, 5, 10, 10);
				glTranslated(0.2, 0, 5);
				gluCylinder(colp, 1.2, 0, 1, 10, 10);
				glPopMatrix();  // возвращаемся к старой системе координат


				glColor3d(1, 1, 1);

				{
					double A[] = { -8, 11.3, -1.2 };
					double B[] = { -8, 11.3, -2.8 };
					double C[] = { -6, 13, -2.8 };
					double* N = norm(A, B, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(-8, 11.3, -1.2);
					glVertex3d(-8, 11.3, -2.8);
					glVertex3d(-6, 13, -2.8);
					glVertex3d(-6, 13, -1.2);
					glEnd(); }
				{
					double A[] = { -8, 11.3, -1.2 };
					double B[] = { -8, 11.1, -1.2 };
					double C[] = { -6.3, 12.8, -1.2 };
					double* N = norm(B, A, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(-8, 11.3, -1.2);
					glVertex3d(-8, 11.1, -1.2);
					glVertex3d(-6.3, 12.8, -1.2);
					glVertex3d(-6, 13, -1.2);
					glEnd(); }
				{
					double A[] = { -8, 11.3, -2.8 };
					double B[] = { -8, 11.1, -2.8 };
					double C[] = { -6.3, 12.8, -2.8 };
					double* N = norm(A, B, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(-8, 11.3, -2.8);
					glVertex3d(-8, 11.1, -2.8);
					glVertex3d(-6.3, 12.8, -2.8);
					glVertex3d(-6, 13, -2.8);
					glEnd(); }
				{
					double A[] = { -8, 11.1, -1.2 };
					double B[] = { -8, 11.1, -2.8 };
					double C[] = { -6.3, 12.8, -2.8 };
					double* N = norm(A, B, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(-8, 11.1, -1.2);
					glVertex3d(-8, 11.1, -2.8);
					glVertex3d(-6.3, 12.8, -2.8);
					glVertex3d(-6.3, 12.8, -1.2);
					glEnd(); }
				{
					double A[] = { -6, 13, -2.8 };
					double B[] = { -6, 13, -1.2 };
					double C[] = { -3.5, 15.5, -1.2 };
					double* N = norm(B, A, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(-6, 13, -2.8);
					glVertex3d(-6, 13, -1.2);
					glVertex3d(-3.5, 15.5, -1.2);
					glVertex3d(-3.5, 15.5, -2.8);
					glEnd(); }
				{
					double A[] = { -6.3, 12.8, -1.2 };
					double B[] = { -3.2, 15.2, -1.2 };
					double C[] = { -3.5, 15.5, -1.2 };
					double* N = norm(B, A, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(-6.3, 12.8, -1.2);
					glVertex3d(-3.2, 15.2, -1.2);
					glVertex3d(-3.5, 15.5, -1.2);
					glVertex3d(-6, 13, -1.2);
					glEnd(); }
				{
					double A[] = { -6.3, 12.8, -2.8 };
					double B[] = { -3.2, 15.2, -2.8 };
					double C[] = { -3.5, 15.5, -2.8 };
					double* N = norm(A, B, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(-6.3, 12.8, -2.8);
					glVertex3d(-3.2, 15.2, -2.8);
					glVertex3d(-3.5, 15.5, -2.8);
					glVertex3d(-6, 13, -2.8);
					glEnd(); }
				{
					double A[] = { -6.3, 12.8, -2.8 };
					double B[] = { -6.3, 12.8, -1.2 };
					double C[] = { -3.2, 15.2, -1.2 };
					double* N = norm(B, A, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(-6.3, 12.8, -2.8);
					glVertex3d(-6.3, 12.8, -1.2);
					glVertex3d(-3.2, 15.2, -1.2);
					glVertex3d(-3.2, 15.2, -2.8);
					glEnd(); }
				{
					double A[] = { -3.5, 15.5, -1.2 };
					double B[] = { -3.5, 15.5, -2.8 };
					double C[] = { -2.5, 16, -2.8 };
					double* N = norm(A, B, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(-3.5, 15.5, -1.2);
					glVertex3d(-3.5, 15.5, -2.8);
					glVertex3d(-2.5, 16, -2.8);
					glVertex3d(-2.5, 16, -1.2);
					glEnd(); }
				{
					double A[] = { -3.5, 15.5, -1.2 };
					double B[] = { -2.5, 16, -1.2 };
					double C[] = { -2.5, 15.6, -1.2 };
					double* N = norm(A, B, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(-3.5, 15.5, -1.2);
					glVertex3d(-2.5, 16, -1.2);
					glVertex3d(-2.5, 15.6, -1.2);
					glVertex3d(-3.2, 15.2, -1.2);
					glEnd(); }
				{
					double A[] = { -3.5, 15.5, -2.8 };
					double B[] = { -2.5, 16, -2.8 };
					double C[] = { -2.5, 15.6, -2.8 };
					double* N = norm(B, A, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(-3.5, 15.5, -2.8);
					glVertex3d(-2.5, 16, -2.8);
					glVertex3d(-2.5, 15.6, -2.8);
					glVertex3d(-3.2, 15.2, -2.8);
					glEnd(); }
				{
					double A[] = { -3.2, 15.2, -1.2 };
					double B[] = { -3.2, 15.2, -2.8 };
					double C[] = { -2.5, 15.6, -2.8 };
					double* N = norm(A, B, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(-3.2, 15.2, -1.2);
					glVertex3d(-3.2, 15.2, -2.8);
					glVertex3d(-2.5, 15.6, -2.8);
					glVertex3d(-2.5, 15.6, -1.2);
					glEnd(); }

				//////
				{
					double A[] = { 8, 11.3, -1.2 };
					double B[] = { 8, 11.3, -2.8 };
					double C[] = { 6, 13, -2.8 };
					double* N = norm(B, A, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(8, 11.3, -1.2);
					glVertex3d(8, 11.3, -2.8);
					glVertex3d(6, 13, -2.8);
					glVertex3d(6, 13, -1.2);
					glEnd(); }
				{
					double A[] = { 8, 11.3, -1.2 };
					double B[] = { 8, 11.1, -1.2 };
					double C[] = { 6.3, 12.8, -1.2 };
					double* N = norm(C,A, B);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(8, 11.3, -1.2);
					glVertex3d(8, 11.1, -1.2);
					glVertex3d(6.3, 12.8, -1.2);
					glVertex3d(6, 13, -1.2);
					glEnd(); }
				{
					double A[] = { 8, 11.3, -2.8 };
					double B[] = { 8, 11.1, -2.8 };
					double C[] = { 6.3, 12.8, -2.8 };
					double* N = norm(B, A, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(8, 11.3, -2.8);
					glVertex3d(8, 11.1, -2.8);
					glVertex3d(6.3, 12.8, -2.8);
					glVertex3d(6.3, 13, -2.8);
					glEnd(); }
				{
					double A[] = { 8, 11.1, -1.2 };
					double B[] = { 8, 11.1, -2.8 };
					double C[] = { 6.3, 12.8, -2.8 };
					double* N = norm(B, A, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(8, 11.1, -1.2);
					glVertex3d(8, 11.1, -2.8);
					glVertex3d(6.3, 12.8, -2.8);
					glVertex3d(6.3, 12.8, -1.2);
					glEnd(); }
				{
					double A[] = { 6, 13, -2.8 };
					double B[] = { 6, 13, -1.2 };
					double C[] = { 3.5, 15.5, -1.2 };
					double* N = norm(A, B, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(6, 13, -2.8);
					glVertex3d(6, 13, -1.2);
					glVertex3d(3.5, 15.5, -1.2);
					glVertex3d(3.5, 15.5, -2.8);
					glEnd(); }
				{
					double A[] = { 6.3, 12.8, -1.2 };
					double B[] = { 3.2, 15.2, -1.2 };
					double C[] = { 3.5, 15.5, -1.2 };
					double* N = norm(A, B, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(6.3, 12.8, -1.2);
					glVertex3d(3.2, 15.2, -1.2);
					glVertex3d(3.5, 15.5, -1.2);
					glVertex3d(6, 13, -1.2);
					glEnd(); }
				{
					double A[] = { 6.3, 12.8, -2.8 };
					double B[] = { 3.2, 15.2, -2.8 };
					double C[] = { 3.5, 15.5, -2.8 };
					double* N = norm(B, A, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(6.3, 12.8, -2.8);
					glVertex3d(3.2, 15.2, -2.8);
					glVertex3d(3.5, 15.5, -2.8);
					glVertex3d(6, 13, -2.8);
					glEnd(); }
				{
					double A[] = { 6.3, 12.8, -2.8 };
					double B[] = { 6.3, 12.8, -1.2 };
					double C[] = { 3.2, 15.2, -1.2 };
					double* N = norm(A, B, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(6.3, 12.8, -2.8);
					glVertex3d(6.3, 12.8, -1.2);
					glVertex3d(3.2, 15.2, -1.2);
					glVertex3d(3.2, 15.2, -2.8);
					glEnd(); }

				{
					double A[] = { 3.5, 15.5, -1.2 };
					double B[] = { 3.5, 15.5, -2.8 };
					double C[] = { 2.5, 16, -2.8 };
					double* N = norm(B, A, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(3.5, 15.5, -1.2);
					glVertex3d(3.5, 15.5, -2.8);
					glVertex3d(2.5, 16, -2.8);
					glVertex3d(2.5, 16, -1.2);
					glEnd(); }
				{
					double A[] = { 3.5, 15.5, -1.2 };
					double B[] = { 2.5, 16, -1.2 };
					double C[] = { 2.5, 15.6, -1.2 };
					double* N = norm(B, A, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(3.5, 15.5, -1.2);
					glVertex3d(2.5, 16, -1.2);
					glVertex3d(2.5, 15.6, -1.2);
					glVertex3d(3.2, 15.2, -1.2);
					glEnd(); }
				{
					double A[] = { 3.5, 15.5, -2.8 };
					double B[] = { 2.5, 16, -2.8 };
					double C[] = { 2.5, 15.6, -2.8 };
					double* N = norm(A, B, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(3.5, 15.5, -2.8);
					glVertex3d(2.5, 16, -2.8);
					glVertex3d(2.5, 15.6, -2.8);
					glVertex3d(3.2, 15.2, -2.8);
					glEnd(); }
				{
					double A[] = { 3.2, 15.2, -1.2 };
					double B[] = { 3.2, 15.2, -2.8 };
					double C[] = { 2.5, 15.6, -2.8 };
					double* N = norm(B, A, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glBegin(GL_QUADS);
					glVertex3d(3.2, 15.2, -1.2);
					glVertex3d(3.2, 15.2, -2.8);
					glVertex3d(2.5, 15.6, -2.8);
					glVertex3d(2.5, 15.6, -1.2);
					glEnd(); }


				{
					double A[] = { -2.5, 15.6, -1.2 };
					double B[] = { -2.5, 15.6, -2.8 };
					double C[] = { 2.5, 15.6, -1.2 };
					double* N = norm(B,A, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(-2.5, 15.6, -1.2);
					glVertex3d(-2.5, 15.6, -2.8);
					glVertex3d(2.5, 15.6, -2.8);
					glVertex3d(2.5, 15.6, -1.2);
					glEnd(); }
				{
					double A[] = { 2.5, 16, -1.2 };
					double B[] = { 2.5, 15.6, -1.2 };
					double C[] = { -2.5, 15.6, -1.2 };
					double* N = norm( C, A,B);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(2.5, 16, -1.2);
					glVertex3d(2.5, 15.6, -1.2);
					glVertex3d(-2.5, 15.6, -1.2);
					glVertex3d(-2.5, 16, -1.2);

					glEnd(); }
				{
					double A[] = { 2.5, 16, -2.8 };
					double B[] = { 2.5, 15.6, -2.8 };
					double C[] = { -2.5, 15.6, -2.8 };
					double* N = norm(B, A, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glVertex3d(2.5, 16, -2.8);
					glVertex3d(2.5, 15.6, -2.8);
					glVertex3d(-2.5, 15.6, -2.8);
					glVertex3d(-2.5, 16, -2.8);
					glEnd(); }
				{
					double A[] = { -2.5, 16, -1.2 };
					double B[] = { -2.5, 16, -2.8 };
					double C[] = { 2.5, 16, -2.8 };
					double* N = norm(A, B, C);
					glBegin(GL_QUADS);
					glNormal3dv(N);
					glBegin(GL_QUADS);
					glVertex3d(-2.5, 16, -1.2);
					glVertex3d(-2.5, 16, -2.8);
					glVertex3d(2.5, 16, -2.8);
					glVertex3d(2.5, 16, -1.2);
					glEnd(); }
				glColor3d(-1, -1, -1);
			}}

	
			int r = 0;
		glPushMatrix();
		
		
		glTranslatef(0.0f, 0.0f, -3.0f);
		vert(tyu);
		glTranslatef(2.0f, -4.0f, 0.0f);
		vert(tyu);
		glTranslatef(-6.0f, 2.0f, 0.0f);
		vert(tyu);
		
		glTranslatef(8.0f, 4.0f, 0.0f);
		vert(tyu);
		
		glTranslatef(-9.0f, 2.0f, 0.0f);
		vert(tyu);
		glTranslatef(5.0f, 2.0f, 0.0f);
		vert(tyu);
		glTranslatef(2.0f, -10.0f, 0.0f);
		vert(tyu);
		glTranslatef(5.0f, 3.0f, 0.0f);
		vert(tyu);
		tyu = tyu+1.3f;
		glPopMatrix();

		//anim();

		glPushMatrix();
		glRotatef(rquad, 0.0f, 0.0f, 1.0f);
		strelkiCh();
		glRotatef(rtri, 0.0f, 0.0f, 1.0f);
		strelkiM();
		rtri -= 1.8f;                    
		rquad -= 0.3f;
		glPopMatrix();
		


		

		

					


glPopMatrix();

		
		


glDisable(GL_TEXTURE_2D);
	
	
		
	//Сообщение вверху экрана
	glMatrixMode(GL_PROJECTION);	//Делаем активной матрицу проекций. 
	                                //(всек матричные операции, будут ее видоизменять.)
	glPushMatrix();   //сохраняем текущую матрицу проецирования (которая описывает перспективную проекцию) в стек 				    
	glLoadIdentity();	  //Загружаем единичную матрицу
	glOrtho(0, ogl->getWidth(), 0, ogl->getHeight(), 0, 1);	 //врубаем режим ортогональной проекции

	glMatrixMode(GL_MODELVIEW);		//переключаемся на модел-вью матрицу
	glPushMatrix();			  //сохраняем текущую матрицу в стек (положение камеры, фактически)
	glLoadIdentity();		  //сбрасываем ее в дефолт

	glDisable(GL_LIGHTING);





	GuiTextRectangle rec;		   //классик моего авторства для удобной работы с рендером текста.
	rec.setSize(300, 150);
	rec.setPosition(10, ogl->getHeight() - 150 - 10);


	std::stringstream ss;
	ss << "T - вкл/выкл текстур" << std::endl;
	ss << "L - вкл/выкл освещение" << std::endl;
	ss << "F - Свет из камеры" << std::endl;
	ss << "G - двигать свет по горизонтали" << std::endl;
	ss << "G+ЛКМ двигать свет по вертекали" << std::endl;
	ss << "Коорд. света: (" << light.pos.X() << ", " << light.pos.Y() << ", " << light.pos.Z() << ")" << std::endl;
	ss << "Коорд. камеры: (" << camera.pos.X() << ", " << camera.pos.Y() << ", " << camera.pos.Z() << ")" << std::endl;
	ss << "Параметры камеры: R="  << camera.camDist << ", fi1=" << camera.fi1 << ", fi2=" << camera.fi2 << std::endl;
	
	rec.setText(ss.str().c_str());
	rec.Draw();

	glMatrixMode(GL_PROJECTION);	  //восстанавливаем матрицы проекции и модел-вью обратьно из стека.
	glPopMatrix();


	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

}