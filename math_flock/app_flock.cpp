//--------------------------------------------------------
// JUST MATH:
// Flockv v2
// 
// Rama Hoetzlein, 2023
//
//--------------------------------------------------------------------------------
// Copyright 2019-2023 (c) Quanta Sciences, Rama Hoetzlein, ramakarl.com
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
#include "main.h"					// window system 
#include "nv_gui.h"				// gui system
#include "quaternion.h"
#include "datax.h"
#include "mersenne.h"

struct Bird {

	Vector3DF		pos, vel, accel;
	Vector3DF		lift, thrust, drag, force;	
	Quaternion	orient, ctrl;
	Vector3DF		target;
	float				speed;	
	float				pitch_adv, power, aoa;
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

	void			AddBird ( Vector3DF pos, Vector3DF vel, Vector3DF target, float power );

	void			Reset ();
	void			Advance ();	
	void			CameraToBird ( int b );
	void			CameraToCockpit( int b );
	void			drawGrid( Vector4DF clr );

	int				m_num_birds;
	DataX			m_Birds;

	float			m_DT;
	Vector3DF	m_wind;
	Mersenne  m_rnd;	


	float			m_time, m_max_speed;
	bool			m_run, m_flightcam;
	Camera3D*	m_cam;
	int				mouse_down;
	int			  m_bird_sel;
	bool		  m_cockpit_view;
};

Sample obj;

void Sample::AddBird ( Vector3DF pos, Vector3DF vel, Vector3DF target, float power )
{
	Bird b;
	b.pos = pos;
	b.vel = vel;	
	b.target = target;
	b.orient.fromDirectionAndRoll ( Vector3DF(0,0,1), b.target.x );
	b.power = power;
	b.pitch_adv = 0;
	b.accel.Set(0,0,0);
	
	int ndx = m_Birds.AddElem (0);
	m_Birds.SetElem (0, ndx, &b );
}


bool Sample::init ()
{
	int w = getWidth(), h = getHeight();			// window width &f height
	m_run = true;
	m_flightcam = true;
	m_bird_sel = 0;
	m_cockpit_view = false;
	m_rnd.seed (12);

	addSearchPath ( ASSET_PATH );
	init2D ( "arial" );
	setview2D ( w, h );	
	setText ( 16, 1 );		
	
	m_cam = new Camera3D;
	m_cam->setFov ( 120 );
	m_cam->setNearFar ( 1.0, 100000 );
	m_cam->SetOrbit ( Vector3DF(-30,30,0), Vector3DF(0,5,0), 300, 1 );

	// Initialize birds
	// * birds are placed into a DataX structure to allow
	// for easy sharing between CPU and GPU
	
	m_num_birds = 100;

  m_Birds.AddBuffer ( 0, "bird",		sizeof(Bird),	m_num_birds, DT_CPU | DT_CUMEM | DT_GLVBO );		// cuda-opengl interop

  Reset ();

	m_max_speed = 500.0;		// top speed, 500 m/s = 1800 kph = 1118 mph
	m_DT = 0.001;
	m_wind.Set (0, 0, 0);

	return true;
}

void Sample::Reset ()
{
	Vector3DF pos, vel;

	m_Birds.ClearBuffer (0);

	for (int n=0; n < m_num_birds; n++ ) {
		pos.x = m_rnd.randF( -50, 50 );
		pos.y = m_rnd.randF(  100, 200 );
		pos.z = m_rnd.randF( -50, 50 );
		vel = m_rnd.randV3( -10, 10 );
		AddBird ( pos, vel, Vector3DF(0, 0, 90), 3);
	}
}


void Sample::drawGrid( Vector4DF clr )
{
	Vector3DF a;
	float o = 0.02;

	// center section
	o = -0.02;			// offset
	for (int n=-5000; n <= 5000; n += 50 ) {
		drawLine3D ( Vector3DF(n, o,-5000), Vector3DF(n, o, 5000), Vector4DF(1,1,1,0.3) );
		drawLine3D ( Vector3DF(-5000, o, n), Vector3DF(5000, o, n), Vector4DF(1,1,1,0.3) );
	}
	
	// large sections
	for (int j=-5; j <=5 ; j++) {
		for (int k=-5; k <=5; k++) {
			a = Vector3DF(j, 0, k) * Vector3DF(5000,0,5000);
			if (j==0 && k==0) continue;
			for (int n=0; n <= 5000; n+= 200) {
				drawLine3D ( Vector3DF(a.x,   o, a.z+n), Vector3DF(a.x+5000, o, a.z+n), Vector4DF(1,1,1,0.2) );
				drawLine3D ( Vector3DF(a.x+n,-o, a.z  ), Vector3DF(a.x+n,    o, a.z+5000), Vector4DF(1,1,1,0.2) );
			}
		}
	}

}

float circleDelta (float b, float a)
{
	return fmin ( b-a, 360+a-b );
}

void Sample::Advance ()
{
	Vector3DF fwd, up, right;
	Vector3DF force, torque, vaxis;
	Quaternion ctrl_pitch;
	float airflow, dynamic_pressure;
	float m_LiftFactor = 0.001;
	float m_DragFactor = 0.001;
	float mass = 0.1;							// body mass (kg)
	float CL, L;
	Quaternion ctrlq, tq;
	Vector3DF angs;
	Quaternion angvel;


	Bird* b;
	for (int n=0; n < m_Birds.GetNumElem(0); n++) {

		b = (Bird*) m_Birds.GetElem(0, n);


		// Body orientation
		fwd = Vector3DF(1,0,0) * b->orient;			// X-axis is body forward
		up  = Vector3DF(0,1,0) * b->orient;			// Y-axis is body up
		right = Vector3DF(0,0,1) * b->orient;		// Z-axis is body right

		// Direction of motion
		b->speed = b->vel.Length();
		vaxis = b->vel / b->speed;	
		if ( b->speed < 0 ) b->speed =  0;		// planes dont go in reverse
		if ( b->speed > m_max_speed ) b->speed = m_max_speed;
		if ( b->speed==0) vaxis = fwd;
			
		b->orient.toEuler ( angs );				
		angs.z = fmod (angs.z, 360.0 );

		float reaction_delay = 0.00001;

		// Banking - correlate banking with yaw
		b->target.x = circleDelta(b->target.z, angs.z);

		// Roll - Control input
		// - orient the body by roll
		ctrlq.fromAngleAxis ( (b->target.x - angs.x) * reaction_delay, Vector3DF(1,0,0) * b->orient );
		b->orient *= ctrlq;	b->orient.normalize();
		
		// Pitch & Yaw - Control inputs
		// - apply 'torque' by rotating the velocity vector based on pitch & yaw inputs		
		ctrlq.fromAngleAxis ( circleDelta(b->target.z, angs.z) * reaction_delay, Vector3DF(0,-1,0) * b->orient );
		vaxis *= ctrlq; vaxis.Normalize();	
		ctrlq.fromAngleAxis ( (b->target.y - angs.y) * reaction_delay, Vector3DF(0,0,1) * b->orient );
		vaxis *= ctrlq; vaxis.Normalize();	

		b->vel = vaxis * b->speed;

		b->force = 0;
		torque = 0;

		// Dynamic pressure		
		airflow = b->speed + m_wind.Dot ( fwd*-1.0f );		// airflow = aircraft speed + wind over wing
		float p = 1.225;											// air density, kg/m^3
		float dynamic_pressure = 0.5f * p * airflow * airflow;

		// Lift force
		b->aoa = acos( fwd.Dot( vaxis ) )*RADtoDEG + 1;				// angle-of-attack = angle between velocity and body forward		
 		if (isnan(b->aoa)) b->aoa = 1;
		CL = sin( b->aoa * 0.2);					// CL = coeff of lift, approximate CL curve with sin
		L = CL * dynamic_pressure * m_LiftFactor * 0.5;		// lift equation. L = CL (1/2 p v^2) A
		b->lift = up * L;
		b->force += b->lift;	

		// Drag force	
		b->drag = vaxis * dynamic_pressure * m_DragFactor * -1.0f;			// drag equation. D = Cd (1/2 p v^2) A
		b->force += b->drag; 

		// Thrust force
		b->thrust = fwd * b->power;
		b->force += b->thrust;
	
		// Integrate position		
		b->accel = b->force / mass;				// body forces	
		b->accel += Vector3DF(0,-9.8,0);	// gravity
		b->accel += m_wind * p * 0.1f;		// wind force. Fw = w^2 p * A, where w=wind speed, p=air density, A=frontal area
	
		b->pos += b->vel * m_DT;


		// Level flight
		b->target.y = 0;
		
		// Ground avoidance
		if ( b->pos.y < 10 && fwd.y < 0) {			
			b->target.y = 20;
		}

		// Ground condition
		if (b->pos.y <= 0.00001 ) { 
			// Ground forces
			b->pos.y = 0; b->vel.y = 0; 
			b->accel += Vector3DF(0,9.8,0);	// ground force (upward)
			b->vel *= 0.9999;				// ground friction
			b->orient.fromDirectionAndRoll ( Vector3DF(fwd.x, 0, fwd.z), 0 );	// zero pitch & roll			
		} 
	
		// Integrate velocity
		b->vel += b->accel * m_DT;		

		vaxis = b->vel;	vaxis.Normalize ();

		// Update Orientation
		// Directional stability: airplane will typically reorient toward the velocity vector
		//  see: https://en.wikipedia.org/wiki/Directional_stability
		// this is an assumption yet much simpler/faster than integrating body orientation
		// this way we dont need torque, angular vel, or rotational inertia.
		// stalls are possible but not flat spins or 3D flying		
		angvel.fromRotationFromTo ( fwd, vaxis, .1 );
		if ( !isnan(angvel.X) ) {
			b->orient *= angvel;
			b->orient.normalize();			
		}

	}
}

void Sample::CameraToBird ( int n )
{
	Bird* b = (Bird*) m_Birds.GetElem(0, n);


	m_cam->SetOrbit ( m_cam->getAng(), b->pos, m_cam->getOrbitDist(), m_cam->getDolly() );

}


void Sample::CameraToCockpit(int n )
{
	Bird* b = (Bird*) m_Birds.GetElem(0, n);

	// View direction	
	Vector3DF fwd = b->vel; fwd.Normalize();
	Vector3DF angs;
	b->orient.toEuler ( angs );

	// Set eye level above centerline
	Vector3DF p = b->pos + Vector3DF(0,2,0);	  
	
	m_cam->setDirection ( p, p + fwd, -angs.x );
}

void Sample::display ()
{	
	char msg[2048];
	Vector3DF x,y,z;
	Vector3DF pnt;
	Vector4DF clr;
	int w = getWidth();
	int h = getHeight();
		
	Bird* b;

	if (m_run) { 		
		Advance ();
	}	

	if (m_cockpit_view) {
		CameraToCockpit ( m_bird_sel);
	} else {
		CameraToBird ( m_bird_sel );
  }

	clearGL();
	start2D ();
	setview2D ( w, h );

	b =  (Bird*) m_Birds.GetElem(0, m_bird_sel);
	sprintf ( msg, "%f %f %f, %f\n", b->target.x, b->target.y, b->target.z, b->speed );
	drawText ( 10, 10, msg, 1,1,1,1);



	start3D(m_cam);

		// Draw ground
		drawGrid( (m_flightcam) ? Vector4DF(1,1,1,1) : Vector4DF(0.2,0.2,0.2,1) );


		for (int n=0; n < m_Birds.GetNumElem(0); n++) {

			b = (Bird*) m_Birds.GetElem(0, n);
		
			x = Vector3DF(1,0,0) * b->orient;
			y = Vector3DF(0,1,0) * b->orient;
			z = Vector3DF(0,0,1) * b->orient;

			drawLine3D ( b->pos-z,	b->pos+z,					Vector4DF(1,1,1,0.4) );			// wings
			drawLine3D ( b->pos,	  b->pos+y,					Vector4DF(1,1,0,  1) );			// up
			drawLine3D ( b->pos,	  b->pos+x,					Vector4DF(1,1,0,  1) );			// fwd
			drawLine3D ( b->pos,		b->pos +b->vel*0.1f,		Vector4DF(1,1,1,  1) );			// velocity
			drawLine3D ( b->pos,		b->pos +b->lift + x*0.1f,		Vector4DF(0,1,0,  1) );			// lift (green)
			drawLine3D ( b->pos,		b->pos+b->thrust, Vector4DF(1,0,0,  1) );			// thrust (red)
			drawLine3D ( b->pos,		b->pos+b->drag,		Vector4DF(1,0,1,  1) );			
			drawLine3D ( b->pos,		b->pos+b->force,	Vector4DF(0,1,1,  1) );		

			/*x = Vector3DF(1,0,0) * b->ctrl;			
			z = Vector3DF(0,0,1) * b->ctrl;			
			drawLine3D ( b->pos,	  b->pos+x,					Vector4DF(0,1,1,  1) );			// ctrl fwd
			drawLine3D ( b->pos-z,	b->pos+z,					Vector4DF(0,1,1,  1) );	*/

		}
	end3D();

	draw3D ();
	draw2D (); 	
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
	Vector3DF dang; 

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
		Vector3DF angs = m_cam->getAng();
		angs.x += dx*0.2f;
		angs.y -= dy*0.2f;				
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
	if (action == AppEnum::BUTTON_RELEASE) 
		return;

	switch ( keycode ) {
  case 'c': m_cockpit_view = !m_cockpit_view; break;
	case 'r': Reset(); break;
	case ' ':	m_run = !m_run;	break;	
	case 'z': 
		m_bird_sel--; 
		if (m_bird_sel < 0) m_bird_sel = 0; 
		break;
	case 'x':
		m_bird_sel++; 
		if (m_bird_sel > m_Birds.GetNumElem(0)) m_bird_sel = m_Birds.GetNumElem(0)-1;
		break;
	};
}

void Sample::reshape (int w, int h)
{
	glViewport ( 0, 0, w, h );
	setview2D ( w, h );

	m_cam->setAspect(float(w) / float(h));
	m_cam->SetOrbit(m_cam->getAng(), m_cam->getToPos(), m_cam->getOrbitDist(), m_cam->getDolly());	
		
	appPostRedisplay();	
}

void Sample::startup ()
{
	int w = 1900, h = 1000;
	appStart ( "Flock v2", "Flock v2", w, h, 4, 2, 16, false );
}

void Sample::shutdown()
{
}


