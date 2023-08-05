//--------------------------------------------------------
// JUST MATH:
// Bilinear Patch
//
// This implements ray intersection with a bilinear patch, 
// also known as an irregular quadrilateral, composed of four 3D points.
// Based on the code from:
//  "Cool Patches: A Geometric Approach to Ray/Bilinear Patch Intersections"
//  from Ray Tracing Gems, NVIDIA, 2019
//
//--------------------------------------------------------

//--------------------------------------------------------------------------------
// Copyright 2019-2023 (c) Quanta Sciences, Rama Hoetzlein, ramakarl.com
//
// * Derivative works may append below the copyright notices but should not remove or modify earlier notices.
//
// MIT License:
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
// associated documentation files (the "Software"), to deal in the Software without restriction, including without 
// limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, 
// and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS 
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF 
// OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

// Sample utils
#include "main.h"			// window system 
#include "nv_gui.h"			// gui system
#include "mersenne.h"
#include <GL/glew.h>

class Sample : public Application {
public:
	virtual void startup();
	virtual bool init();
	virtual void display();
	virtual void reshape(int w, int h);
	virtual void motion(AppEnum button, int x, int y, int dx, int dy);
	virtual void keyboard(int keycode, AppEnum action, int mods, int x, int y);
	virtual void mouse(AppEnum button, AppEnum state, int mods, int x, int y);
	virtual void mousewheel(int delta);
	virtual void shutdown();

	bool		intersect_ray_patch (Vector3DF rpos, Vector3DF rdir, Vector3DF q00, Vector3DF q10, Vector3DF q11, Vector3DF q01, float& t, Vector3DF& n ) ;

	Vector3DF	getPatchPnt ( float u, float v );
	void		DrawPatch ();
	void		RaytracePatch ();	
	void		GeneratePatch ();

	void		Resize ();
	void		DrawGrid ();	
	void		UpdateCamera();
	void		MoveCamera ( char t, Vector3DF del );

	Vector3DF	m_p[4];			// input patch

	int			m_numpnts;		
	Vector4DF	m_pnts[16384];	// output intersections
	Vector3DF	m_norms[16384];

	Camera3D*	m_cam;			// camera		
	bool		m_run;
	int			mouse_down;
	Mersenne	m_rand;

};
Sample obj;


bool Sample::init()
{
	addSearchPath(ASSET_PATH);
	init2D( "arial" );
	setText( .02, 1);

	m_run = true;

	m_rand.seed(12);

	// create a camera
	m_cam = new Camera3D;								
	m_cam->setFov ( 60 );
	m_cam->setAspect ( 1 );	
	m_cam->setOrbit ( -20, 20, 0, Vector3DF(0, 2, 0), 30, 1 );	
	Resize ();		
	UpdateCamera();

	// output points
	m_numpnts = 0;

	// create patch
	m_p[0].Set ( -4, 2, +5 );
	m_p[1].Set ( -5, 7, -5 );	
	m_p[2].Set (  4, 4, -3 );	
	m_p[3].Set (  4, 6, +3 );	

	return true;
}

void Sample::Resize ()
{
}

void Sample::DrawGrid ()
{
	for (int i = -10; i <= 10; i++) {
		drawLine3D(float(i),-0.01f, -10.f, float(i), -0.01f, 10.f, .2f, .2f, .2f, 1.f);
		drawLine3D(-10.f,	-0.01f, float(i), 10.f, -0.01f, float(i), .2f, .2f, .2f, 1.f);
	}	
	drawLine3D ( 0,0,0, 1,0,0, 1,0,0,1 );
	drawLine3D ( 0,0,0, 0,1,0, 0,1,0,1 );
	drawLine3D ( 0,0,0, 0,0,1, 0,0,1,1 );

}

void Sample::GeneratePatch ()
{
	// generate a random patch by fixing the x-axis and allowing random y/z.
	// this can possibly generate an overlapping or twisted patch. intersector handles it ok.
	//
	m_p[0] = Vector3DF( -1, m_rand.randF(0,1), m_rand.randF(-1,1) ) * 5.0f;
	m_p[1] = Vector3DF( -1, m_rand.randF(0,1), m_rand.randF(-1,1) ) * 5.0f;
	m_p[2] = Vector3DF(  1, m_rand.randF(0,1), m_rand.randF(-1,1) ) * 5.0f;
	m_p[3] = Vector3DF(  1, m_rand.randF(0,1), m_rand.randF(-1,1) ) * 5.0f;
}

Vector3DF Sample::getPatchPnt ( float u, float v )
{
	Vector3DF a, b, c;
	a = m_p[0] + (m_p[1] - m_p[0]) * u;
	b = m_p[3] + (m_p[2] - m_p[3]) * u;
	c = a + (b - a) * v;
	return c;
}


void Sample::DrawPatch ()
{
	// draw patch u,v lines
	Vector3DF a,b;
	for (float v=0; v <= 1.01; v += 0.1) {
		a = getPatchPnt ( 0, v );
		b = getPatchPnt ( 1, v );
		drawLine3D( a, b, Vector4DF(1,1,0,1));
	}
	for (float u=0; u <= 1.01; u += 0.1) {
		a = getPatchPnt ( u, 0 );
		b = getPatchPnt ( u, 1 );
		drawLine3D( a, b, Vector4DF(1,1,0,1));
	}
	a = getPatchPnt(0,0); drawLine3D ( a, Vector3DF(a.x,0,a.z), Vector4DF(.5,.5,.5,1));
	a = getPatchPnt(1,0); drawLine3D ( a, Vector3DF(a.x,0,a.z), Vector4DF(.5,.5,.5,1));
	a = getPatchPnt(1,1); drawLine3D ( a, Vector3DF(a.x,0,a.z), Vector4DF(.5,.5,.5,1));
	a = getPatchPnt(0,1); drawLine3D ( a, Vector3DF(a.x,0,a.z), Vector4DF(.5,.5,.5,1));
	
	// draw hit points
	Vector3DF n, V = m_cam->getDir ();
	V.Normalize();
		
	Vector4DF clr;
	float rc = (m_run) ? 1 : 0;
	
	for (int p=0; p < m_numpnts; p++) {
		a = Vector3DF ( m_pnts[p] );
		clr = (m_pnts[p].w==1) ? Vector4DF(1, 1, rc, 1) : Vector4DF(1, 0, rc,1);
		drawCircle3D ( a, a+V, 0.02, clr);

		n = m_norms[p];
		drawLine3D ( a, a+n*0.2f, Vector4DF(0,1,1,0.5) );
	}

}

#define fsgn(x)		((x<0) ? -1 : 1)

// ray/bilinear patch intersect
// - the coordinates q should be specified in clockwise order
// 
bool Sample::intersect_ray_patch (Vector3DF rpos, Vector3DF rdir, Vector3DF q00, Vector3DF q10, Vector3DF q11, Vector3DF q01, float& t, Vector3DF& n) 
{
	float t1, v1, t2, v2;
	Vector3DF pa, pb;
	Vector3DF e10 = q10 - q00;
	Vector3DF e11 = q11 - q10;
	Vector3DF e00 = q01 - q00;
	Vector3DF qn = e10.Cross ( (q01-q11) );   // normal of the diagonals
	q00 -= rpos;
	q10 -= rpos;
	float a = q00.Cross ( rdir ).Dot ( e00 );
	float c = qn.Dot ( rdir );
	float b = q10.Cross ( rdir ).Dot ( e11 );
	b -= a + c;
	float det = b*b - 4*a*c;
	if (det < 0 ) return false;
	det = sqrt(det);
	float u1, u2;
	float u=0, v=0;
	t = 1e10;
	if ( c==0 ) {				// trapezid. only one root
		u1 = -a/b; 
		u2 = -1;
	} else {
		u1 = (-b - det*fsgn(b) )/2;
		u2 = a/u1;
		u1 /= c;
	}
	if (0 <= u1 && u1 <=1 ) {
		pa = q00 + (q10-q00)*u1;
		pb = e00 + (e11-e00)*u1;
		n = rdir.Cross ( pb );
		det = n.Dot ( n );
		n = n.Cross ( pa );
		t1 = n.Dot ( pb );
		v1 = n.Dot ( rdir );
		if ( t1 > 0 && 0 <= v1 && v1 <= det ) {
			t = t1/det; u = u1; v = v1/det;
		}
	}
	if ( 0 <= u2 && u2 <= 1 ) {
		pa = q00 + (q10-q00)*u2;
		pb = e00 + (e11-e00)*u2;
		n = rdir.Cross ( pb );
		det = n.Dot ( n );
		n = n.Cross ( pa );
		t2 = n.Dot ( pb ) / det;
		v2 = n.Dot ( rdir );
		if ( 0 <= v2 && v2 <= det && t > t2 && t2 > 0 ) {
			t = t2; u = u2; v = v2/det;
		}
	}
	// surface normal
	Vector3DF du = e10 + ((q11-q01)-e10) * v;
	Vector3DF dv = e00 + (e11-e00) * u;
	n = dv.Cross ( du );	
	n.Normalize ();
	
	// compute shading normal given normals at each vertex (vn)
	// pa = vn[0] + (vn[1]-vn[0])*u;
	// pb = vn[3] + (vn[2]-vn[3])*u;
	// n = pa + (pb-pa)* v;
	// ..OR..
	// n = vn[2]*u*v + vn[1]*(1-u)*v + vn[3]*u*(1-v) + vn[0]*(1-u)*(1-v);  // (slower)

	return true;
}


void Sample::RaytracePatch ()
{
	float t;
	Vector3DF rpos, rdir, n;

	m_numpnts = 0;

	for (int y=0; y < getWidth(); y += 40) {
		for (int x=0; x < getHeight(); x += 40) {

			// Get camera ray
			rpos = m_cam->getPos();
			rdir = m_cam->inverseRay ( x, y, getWidth(), getHeight() );	
			rdir.Normalize();

			Vector3DF V = m_cam->getDir ();
			V.Normalize();

			// Intersect bilinear patch
			if ( intersect_ray_patch ( rpos, rdir, m_p[0], m_p[1], m_p[2], m_p[3], t, n ) ) {

				// compute side relative to ray, store in w of hit pnt
				float side = (n.Dot ( rdir ) > 0 ) ? 0 : 1;

				// record hit
				m_pnts[ m_numpnts ] = Vector4DF( rpos + rdir * t, side );
				m_norms[ m_numpnts ] = n;
				m_numpnts++;
			}
		}
	}
}

void Sample::display()
{

	clearGL();

	if (m_run)
		RaytracePatch ();

	// Draw grid
	start3D(m_cam);	
		DrawGrid();	
	
		DrawPatch ();
	end3D();

	draw3D();

	appPostRedisplay();								// Post redisplay since simulation is continuous
}

void Sample::UpdateCamera()
{
	Vector3DF a, t;
	a = m_cam->getAng();
	t = m_cam->getToPos();
	//dbgprintf ( "angs: %3.4f %3.4f %3.4f, to: %3.4f %3.4f %3.4f\n", a.x,a.y,a.z, t.x,t.y,t.z );

	appPostRedisplay();		// update display
}

void Sample::motion(AppEnum btn, int x, int y, int dx, int dy)
{
	float fine = 0.5;

	switch (mouse_down) {
	case AppEnum::BUTTON_LEFT: {
	} break;

	case AppEnum::BUTTON_MIDDLE: {
		// Adjust target pos		
		m_cam->moveRelative(float(dx) * fine * m_cam->getOrbitDist() / 1000, float(-dy) * fine * m_cam->getOrbitDist() / 1000, 0);
		m_cam->moveRelative(float(dx) * fine * m_cam->getOrbitDist() / 1000, float(-dy) * fine * m_cam->getOrbitDist() / 1000, 0);
		UpdateCamera();
	} break;

	case AppEnum::BUTTON_RIGHT: {

		// Adjust camera orbit 
		Vector3DF angs = m_cam->getAng();
		angs.x += dx * 0.2f * fine;
		angs.y -= dy * 0.2f * fine;
		m_cam->setOrbit(angs, m_cam->getToPos(), m_cam->getOrbitDist(), m_cam->getDolly());		
		m_cam->setOrbit(angs, m_cam->getToPos(), m_cam->getOrbitDist(), m_cam->getDolly());				
		UpdateCamera();
		
		Vector3DF to = m_cam->getToPos();

		//dbgprintf ( "cam: angs %f,%f,%f  to %f,%f,%f  dist %f\n", angs.x, angs.y, angs.z, to.x, to.y, to.z, m_cam->getOrbitDist() );
	} break;
	}
}

void Sample::mouse(AppEnum button, AppEnum state, int mods, int x, int y)
{
	if (guiHandler(button, state, x, y)) return;
	
	mouse_down = (state == AppEnum::BUTTON_PRESS) ? button : -1;		// Track when we are in a mouse drag

	if ( mouse_down==AppEnum::BUTTON_LEFT) {
	}
}

void Sample::mousewheel(int delta)
{
	// Adjust zoom
	float zoomamt = 1.0;
	float dist = m_cam->getOrbitDist();
	float dolly = m_cam->getDolly();
	float zoom = (dist - dolly) * 0.001f;
	dist -= delta * zoom * zoomamt;

	m_cam->setOrbit( m_cam->getAng(), m_cam->getToPos(), dist, dolly);
	m_cam->setOrbit( m_cam->getAng(), m_cam->getToPos(), dist, dolly);

	UpdateCamera();
}


void Sample::MoveCamera ( char t, Vector3DF del )
{
	switch (t) {	
	case 'p': {
		float orbit = m_cam->getOrbitDist() - del.z;
		m_cam->setOrbit( m_cam->getAng(), m_cam->getToPos(), orbit, m_cam->getDolly());
		UpdateCamera();
		} break;
	case 't':
		m_cam->moveRelative(float(del.x) * m_cam->getOrbitDist() / 1000, float(-del.y) * m_cam->getOrbitDist() / 1000, 0);	
		UpdateCamera();
		break;
	}
}

void Sample::keyboard(int keycode, AppEnum action, int mods, int x, int y)
{
	if (action==AppEnum::BUTTON_RELEASE) return;

	float s = (mods==1) ? 10.0 : 1.0;

	switch (keycode) {
	case 'r':
		// generate random patch
		GeneratePatch();
		m_run = true;
		break;
	case 'f': {
		// flip orientation of patch
		Vector3DF t;
		t = m_p[1]; m_p[1] = m_p[2]; m_p[2] = t;
		t = m_p[0]; m_p[0] = m_p[3]; m_p[3] = t;
		m_run = true;
		} break;
	case ' ': 
		m_run = !m_run;
		break;	
	case 'a': case 'A':		MoveCamera('t', Vector3DF(-s, 0, 0));	break;	
	case 'd': case 'D':		MoveCamera('t', Vector3DF(+s, 0, 0));	break;
	case 'q': case 'Q':		MoveCamera('t', Vector3DF(0, -s, 0));	break;
	case 'z': case 'Z':		MoveCamera('t', Vector3DF(0, +s, 0));	break;
	};
}

void Sample::reshape(int w, int h)
{
	glViewport(0, 0, w, h);
	setview2D(w, h);

	m_cam->setSize( w, h );
	m_cam->setAspect(float(w) / float(h));
	m_cam->setOrbit( m_cam->getAng(), m_cam->getToPos(), m_cam->getOrbitDist(), m_cam->getDolly());
	m_cam->updateMatricies();

	appPostRedisplay();
}

void Sample::startup()
{
	int w = 2048, h = 2048;
	appStart("Displacement Map Sphere", "Displacement Map Sphere", w, h, 4, 2, 16, false);
}

void Sample::shutdown()
{
}




