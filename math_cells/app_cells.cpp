//--------------------------------------------------------
// JUST MATH:
// Cell membrane simulation
// 
// Demonstrates a cellular membrane simulation using the physics of colliding circles.
// Cellulose is represented as a collection of plant cells (circles).
// Cells are attracted to a membrane circle, which can be adjusted in size (<,> keys)
// Cells may also respond to temperature ([,] and t keys), with cells migrating across the temperature boundary.
// Rotation may also be induced in the cells (-,+ keys)
// 
//--------------------------------------------------------------------------------
// Copyright 2019-2022 (c) Quanta Sciences, Rama Hoetzlein, ramakarl.com
//
// * Derivative works may append the above copyright notice but should not remove or modify earlier notices.
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

#include <time.h>
#include "main.h"			// window system 
#include "gxlib.h"		// rendering
using namespace glib;

#include "quaternion.h"
#include "mersenne.h"

struct Cell {
	Vec3F	pos;
	Vec3F	vel;
	Vec3F	force;
	float		radius;
	float		temp;
};

class Sample : public Application {
public:
	virtual bool init();
	virtual void startup ();
	virtual void display();
	virtual void reshape(int w, int h);
	virtual void motion (AppEnum button, int x, int y, int dx, int dy);	
	virtual void keyboard(int keycode, AppEnum action, int mods, int x, int y);
	virtual void mouse (AppEnum button, AppEnum state, int mods, int x, int y);	
	virtual void mousewheel(int delta);
	virtual void shutdown();
	
	void		CreateCells();
	void		SimCells();
	void		DivideCells();
	void		InvertTemp();

	void		drawCells();

	void		drawGrid();

	float		m_rotate;
	float		m_stable_radius;
	float		m_temp_radius;
	bool		m_temp;

	std::vector<Cell>	m_cells;

	Vec3F	m_center;

	Mersenne	mt;
	float		m_time;
	bool		m_run;
	Camera3D*	m_cam;
	int			mouse_down;
};
Sample obj;

bool Sample::init ()
{
	int w = getWidth(), h = getHeight();			// window width & height
	m_run = true;

	addSearchPath ( ASSET_PATH );

	init2D ( "arial" );	
	setTextSz ( 16, 1 );		
	
	m_cam = new Camera3D;
	m_cam->SetOrbit ( Vec3F(0, 90,0), Vec3F(0,0,0), 100, 1 );

	mt.seed(164);

	m_rotate = 0;
	m_stable_radius = 5;
	m_temp_radius = 5;
	m_temp = false;
	m_center.Set(0,0,0);

	CreateCells();

	m_time = 0;


	return true;
}

void Sample::CreateCells()
{
	float ring[4];
	float r[4];
	float ang, pct;
	int cnt;
	Cell c;
	ring[0] = 7.0;
	ring[1] = 5.0;
	ring[2] = 3.0;
	
	//r[0] = 1.2; r[1] = 1.2; r[2] = 1.2;			pct = 0.6;
	//r[0] = 1.2; r[1] = 1.1; r[2] = 0.9;			pct = 0.5;
	//r[0] = 0.3; r[1] = 0.4; r[2] = 0.6;				pct = 0.7;
	r[0] = 0.03; r[1] = 0.05; r[2] = 0.3;			pct = 0.95;

	c.pos.Set(0, 0, 0);
	c.vel.Set(0, 0, 0);
	c.force.Set(0, 0, 0);
	c.radius = ring[0];
	c.temp = 0;
	//m_cells.push_back(c);

	
	for (int g=0; g < 3; g++) {
		cnt = int(2.0 * PI * ring[g] / (2 * r[g])) * pct;    //- (g+2)*4;

		for (int i=0; i < cnt; i++) {	
			ang = (i+g*0.3)*360.0/cnt;
			c.pos.Set( cos(ang*DEGtoRAD)*ring[g], 0, sin(ang*DEGtoRAD)*ring[g] );
			c.pos += mt.randV3(-1,1)*0.1f;
			c.pos.y = 0;
		
			c.vel.Set(0,0,0);
			c.force.Set(0, 0, 0);
			c.radius = r[g];
			c.temp = 1; //-g;
			m_cells.push_back ( c );
		}
	}
}

#include "geom_helper.h"

void Sample::SimCells()
{
	float dt = 0.005;
	int coll;

	float dst, r, v;
	float a1, a2, a3, p;
	Vec3F n, f1, f2, ipos;

	// collision detection
	for (int i=0; i < m_cells.size(); i++) {

		// target position
		ipos = m_cells[i].pos + m_cells[i].vel * dt;	 // advance by velocity
		
		// check for collisions
		coll = 0;
		for (int j=0; j < m_cells.size(); j++) {

			if ( i==j ) continue;
			n = ipos - m_cells[j].pos;
			r = m_cells[i].radius + m_cells[j].radius; 
			dst = n.x*n.x + n.z*n.z;

			if ( dst > 0 && dst < r*r ) {
				coll++;
				dst = r-sqrt(dst);						// dst = depth of collision
				n.Normalize();

				// reposition at contact point (and use as new test point)
				ipos = m_cells[j].pos + n * (r+0.001f);	
							
				// note: ideally, find the point along vel where i just touches j, not along the vector n. 
				// using line-circle intersection, where the circle is size R (r1+r2) centered at j.
				
				/* Vec3F v = m_cells[i].vel; v.Normalize();
				double cdn = v.Dot ( m_cells[j].pos - m_cells[i].pos );
				double dd = m_cells[j].pos.Dot ( m_cells[j].pos );
				double cc = m_cells[i].pos.Dot ( m_cells[i].pos );
				double rr = double(r)*r;
				float t1 = cdn - sqrt (cdn*cdn - (dd + cc - rr - 2*m_cells[i].pos.Dot(m_cells[j].pos) ) );
				float t2 = cdn + sqrt (cdn*cdn - (dd + cc - rr - 2*m_cells[i].pos.Dot(m_cells[j].pos) ) );
				ipos = m_cells[i].pos + v * (std::min(t1,t2) - 0.0001f);
				n = ipos - m_cells[j].pos; n.Normalize(); */

				dst = 0.1 + 4.0 * dst * dst;								

				a1 = m_cells[i].vel.Dot ( n );			// reflect around  normal (elastic)
				a2 = m_cells[j].vel.Dot ( n );
				a3 = (m_cells[i].vel.Dot(m_cells[j].vel) - 1) * 0.3;	// directional avoidance
				p = (a3+2*(a1-a2)) / r;					
				p = clamp(p, -20, 20);
				
				f1 = n * m_cells[j].radius * -p * dst;
				f2 = n * m_cells[i].radius * p * dst;
				m_cells[i].force += f1 / dt;
				m_cells[j].force += f2 / dt; 
			}
		}
		// simplified position based dynamics:
		//   0 collisions => advance by velocity
		// 1,2 collisions => reposition at contact point
		//  3+ collisions => best not to move
		if ( coll <= 2 ) {
			m_cells[i].pos = ipos;		
		}
		//dbgprintf("%f %f\n", p, dst);
	}

	// advance
	Vec3F crs;
	Vec3F ctr (0,0,0);

	for (int i = 0; i < m_cells.size(); i++) {
		
		n = m_cells[i].pos - m_center; 
		r = n.Length();
		n *= 1.0/r;

		// vortex		
		crs = crs.Cross ( n, Vec3F(0,1,0) );
		//v = (int(m_time/200.0) % 2)==0 ? 1.0f : -1.0f;		//-- oscillating rotation
		m_cells[i].force += crs * r * m_rotate;
		
		// temperature
		if ( m_cells[i].temp != 0 && m_temp) {
			if (r > m_temp_radius + 0.5) m_cells[i].temp -= 0.021;
			if (r < m_temp_radius - 0.5) m_cells[i].temp += 0.021;
			if ( m_cells[i].temp < -1 ) m_cells[i].temp = -1;
			if (m_cells[i].temp > 1) m_cells[i].temp = 1;
		}

		// convection (inward & outward)
		r -= m_stable_radius;
		if ( m_temp ) m_cells[i].force += n * (r * r) * 1.5f * m_cells[i].temp * (0.5f-m_cells[i].radius);
		if ( m_stable_radius > 0) m_cells[i].force += n * (r * r * r) * -0.5f;

		// update velocity
		m_cells[i].vel += m_cells[i].force * dt;		
		m_cells[i].force = 0;	

		// velocity limiter
		v = m_cells[i].vel.Length();
		if ( v > 3.0 ) m_cells[i].vel *= 0.96f;
			
		ctr += m_cells[i].pos;
	}
	ctr *= 1.0f/m_cells.size();

	m_center = ctr;
	m_center.y = 0;
}

void Sample::InvertTemp()
{
	for (int i = 0; i < m_cells.size(); i++) 
		 m_cells[i].temp = - m_cells[i].temp;
}

void Sample::DivideCells ()
{
	// advance
	float a;
	Vec3F v;
	Cell c;
	int cnt = m_cells.size();
	for (int i = 0; i < cnt; i++) {
		v = m_cells[i].vel;	v.Normalize();
		
		if ( m_cells[i].radius >= 1 ) {
			a = m_cells[i].radius*0.5;
			c.pos =				m_cells[i].pos + v * a;
			m_cells[i].pos =	m_cells[i].pos + v * -a;

			c.radius =			m_cells[i].radius * 0.6;		// half radius
			m_cells[i].radius = m_cells[i].radius * 0.6;
	
			c.vel = v;
			c.force = 0;

			m_cells.push_back(c);
		}
	}
}

void Sample::drawCells()
{
	Vec3F p;
	float r, t;
	Vec4F clr;

	Vec3F N = m_cam->getPos() - m_center;

	drawCircle3D ( m_center, N, m_stable_radius, Vec4F(0.0,0.5,0.0,1));
	if ( m_temp ) {
		drawCircle3D(m_center, N, m_temp_radius+ 0.5, Vec4F(0.0,0.0,0.5, 1));
		drawCircle3D(m_center, N, m_temp_radius- 0.5, Vec4F(0.5, 0.0, 0.0, 1));
	}	

	for (int n=0; n < m_cells.size(); n++) {
		p = m_cells[n].pos;
		r = m_cells[n].radius;
		t = m_cells[n].temp;
		clr = ( t == 0 || m_temp==false ) ? Vec4F(1, 1, 1, 1) : Vec4F( (t+1)*0.5, 0, 1-(t+1)*0.5, 1);
		N = m_cam->getPos() - p;
		drawCircle3D( p, N, r, clr);
	}

}

void Sample::drawGrid()
{
	float o	 = -0.05;		// offset
	for (int n=-100; n <= 100; n+=10 ) {
		drawLine3D ( Vec3F(n, o,-100), Vec3F(n,o,100), Vec4F(.5,.5,.5, .3) );
		drawLine3D ( Vec3F(-100, o, n), Vec3F(100, o, n), Vec4F(.5, .5, .5, .3) );
	}
}


void Sample::display ()
{	
	Vec3F pnt;
	Vec4F clr;

	int w = getWidth();
	int h = getHeight();

	clearGL();
	start2D( w, h );		
		drawText( Vec2F(10, 20), "Input:", Vec4F(1,1,1,1) );
		drawText( Vec2F(10, 65), "  -, +     Induce rotation left/right", Vec4F(1,1,1,1));
		drawText( Vec2F(10, 35), "  <, >     Increase/decrease membrane radius", Vec4F(1,1,1,1));
		drawText( Vec2F(10, 50), "  [, ]     Increase/decrease temperature radius", Vec4F(1,1,1,1));				
		drawText( Vec2F(10, 80), "  t        Turn temperature on/off", Vec4F(1,1,1,1));
	end2D();

	if (m_run) {
		SimCells ();
		m_time += .1;
	}

	start3D( m_cam );
		drawGrid();
	
		drawCells();

	end3D();

	drawAll();

	appPostRedisplay();								// Post redisplay since simulation is continuous
}



void Sample::mouse(AppEnum button, AppEnum state, int mods, int x, int y)
{
	int w = getWidth(), h = getHeight();				// window width & height

	mouse_down = (state == AppEnum::BUTTON_PRESS) ? button : -1;

	if (mouse_down == AppEnum::BUTTON_LEFT) {
		
	}
}


void Sample::motion (AppEnum button, int x, int y, int dx, int dy) 
{
	// Get camera for scene
	bool shift = (getMods() & KMOD_SHIFT);		// Shift-key to modify light
	float fine = 0.5f;
	Vec3F dang; 

	switch ( mouse_down ) {	
	case AppEnum::BUTTON_LEFT:  {	
	
		} break;

	case AppEnum::BUTTON_MIDDLE: {
		// Adjust target pos		
		float zoom = (m_cam->getOrbitDist() - m_cam->getDolly()) * 0.0003f;
		m_cam->moveRelative ( float(dx) * zoom, float(-dy) * zoom, 0 );	
		} break; 

	case AppEnum::BUTTON_RIGHT: {
		// Adjust orbit angles
		Vec3F angs = m_cam->getAng();
		angs.x += dx*0.2f*fine;
		angs.y -= dy*0.2f*fine;				
		m_cam->SetOrbit ( angs, m_cam->getToPos(), m_cam->getOrbitDist(), m_cam->getDolly() );
		} break;	

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

	m_cam->SetOrbit(m_cam->getAng(), m_cam->getToPos(), dist, dolly);		
}




void Sample::keyboard(int keycode, AppEnum action, int mods, int x, int y)
{
	if (action == AppEnum::BUTTON_RELEASE) return;

	switch ( keycode ) {
	case ' ':	m_run = !m_run;	break;
	case ',': case '<':	m_stable_radius -= 0.1; break;
	case '.': case '>':	m_stable_radius += 0.1; break;
	case '[': case '{':	m_temp_radius -= 0.1; break;
	case ']': case '}':	m_temp_radius += 0.1; break;
	case '-': case '_':	m_rotate -= 0.05; break;
	case '=': case '+':	m_rotate += 0.05; break;
	case 't':	m_temp = !m_temp;	break;
	};
}

void Sample::reshape (int w, int h)
{
	glViewport ( 0, 0, w, h );	

	m_cam->setAspect(float(w) / float(h));
	m_cam->SetOrbit(m_cam->getAng(), m_cam->getToPos(), m_cam->getOrbitDist(), m_cam->getDolly());	
		
	appPostRedisplay();	
}

void Sample::startup ()
{
	int w = 1900, h = 1000;
	appStart ( "Cell simulation", "Cell simulation", w, h, 4, 2, 16, false );
}

void Sample::shutdown()
{
}

