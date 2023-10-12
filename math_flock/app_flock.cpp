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

// Particle data
#define FPOINT		0		
#define FGCELL		1
#define FGNDX			2	

// Acceleration grid data
#define AGRID			0	
#define AGRIDCNT		1
#define	AGRIDOFF		2
#define AAUXARRAY1	3
#define AAUXSCAN1		4
#define AAUXARRAY2	5
#define AAUXSCAN2		6

#define GRID_UNDEF				2147483647			// max int
#define SCAN_BLOCKSIZE		512		

struct Bird {
	Vector3DF		pos, vel, accel;
	Vector3DF		lift, thrust, drag, force;	
	Quaternion	orient, ctrl;
	Vector3DF		target;
	float				speed;	
	float				pitch_adv, power, aoa;
	Vector4DF   clr;

	Vector3DF		ave_pos, ave_vel;
	int				  near_j, near_in, nbr_cnt;
	float			  near_in_ang;
};


struct Accel {
	Vector3DF	bound_min, bound_max;
	float			psmoothradius, sim_scale;
	float			grid_size, grid_density;
	Vector3DF	gridSize, gridDelta, gridMin, gridMax;
	Vector3DI	gridRes, gridScanMax;
	int				gridSrch, gridTotal, gridAdjCnt, gridActive;
	int				gridAdj[64];	
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

	Bird*			AddBird ( Vector3DF pos, Vector3DF vel, Vector3DF target, float power );	
	void			FindNeighbors ();

	void			InitializeGrid ();
	void			InsertIntoGrid ();
	void			PrefixSumGrid ();
	void			DrawAccelGrid ();
	

	void			Reset ();
	void		  Run ();
	void			Advance ();	
	void			CameraToBird ( int b );
	void			CameraToCockpit( int b );
	void			drawGrid( Vector4DF clr );

	int				m_num_birds;
	DataX			m_Birds;
	DataX			m_Grid;
	Accel			m_Accel;

	float			m_DT;	
	Vector3DF	m_wind;
	Mersenne  m_rnd;	


	float			m_time;
	float			m_min_speed, m_max_speed;
	bool			m_run, m_flightcam;
	Camera3D*	m_cam;
	Vector3DF	m_cam_fwd;
	int				mouse_down;
	int			  m_bird_sel;
	bool		  m_cockpit_view;
	bool			m_draw_sphere;
	bool		  m_draw_grid;

	std::vector<Vector4DF>  m_mark;
};

Sample obj;

Bird* Sample::AddBird ( Vector3DF pos, Vector3DF vel, Vector3DF target, float power )
{
	Vector3DF dir, angs;
	Bird b;
	b.pos = pos;
	b.vel = vel;		
	b.target = target;	
	b.power = power;
	b.pitch_adv = 0;
	b.accel.Set(0,0,0);

	dir = b.vel; dir.Normalize();	
	b.orient.fromDirectionAndUp ( dir, Vector3DF(0,1,0) );
	b.orient.normalize();
	b.orient.toEuler ( angs );			
	
	int ndx = m_Birds.AddElem (0);
	m_Birds.SetElem (0, ndx, &b );

	return (Bird*) m_Birds.GetElem(0, ndx);
}

void Sample::Reset ()
{
	Vector3DF pos, vel;
	float h;
	int grp;
	Bird* b;

	int numPoints = m_num_birds;

	m_Birds.DeleteAllBuffers ();
	m_Birds.AddBuffer ( FPOINT, "bird",		sizeof(Bird),	numPoints, DT_CPU | DT_CUMEM );
	m_Birds.AddBuffer ( FGCELL, "gcell",	sizeof(uint),	numPoints, DT_CPU | DT_CUMEM );
	m_Birds.AddBuffer ( FGNDX,  "gndx",		sizeof(uint),	numPoints, DT_CPU | DT_CUMEM );		

	for (int n=0; n < m_num_birds; n++ ) {
		
		/*bool ok = false;
		while (!ok) {
			pos = m_rnd.randV3( -50, 50 );
			if (pos.Length() < 50 ) {				
				grp = (n % 2);
				pos += Vector3DF( 0, 100, grp ? -60 : 60 );
				vel = Vector3DF(  0,   0, grp ?  50 :-50 );
				h = grp ? 90 : -90;
				b = AddBird ( pos, vel, Vector3DF(0, 0, h), 3); 
				b->clr = (grp==0) ? Vector4DF(1,0,0,1) : Vector4DF(0,1,0,1);
				ok = true;
			}
		} */
		 
		pos = m_rnd.randV3( -100, 100 ) + Vector3DF(0, 100, 0);
		vel = m_rnd.randV3( -20, 20 );
		h = m_rnd.randF(-180, 180);
		AddBird ( pos, vel, Vector3DF(0, 0, h), 3); 
	}

	printf ( "Added %d birds.\n", m_num_birds );
}


void Sample::drawGrid( Vector4DF clr )
{
	Vector3DF a;
	float o = 0.02;

	// center section
	o = -0.02;			// offset
	for (int n=-5000; n <= 5000; n += 50 ) {
		drawLine3D ( Vector3DF(n, o,-5000), Vector3DF(n, o, 5000), clr );
		drawLine3D ( Vector3DF(-5000, o, n), Vector3DF(5000, o, n), clr );
	}

}


// Ideal grid cell size (gs) = 2 * smoothing radius = 0.02*2 = 0.04
// Ideal domain size = k * gs / d = k*0.02*2/0.005 = k*8 = {8, 16, 24, 32, 40, 48, ..}
//    (k = number of cells, gs = cell size, d = simulation scale)
//
void Sample::InitializeGrid ()
{
	// Grid size - cell spacing in SPH units
	m_Accel.grid_size = m_Accel.psmoothradius / m_Accel.grid_density;	
																					
	// Grid bounds - one cell beyond fluid domain
	m_Accel.gridMin = m_Accel.bound_min;		m_Accel.gridMin -= float(2.0*(m_Accel.grid_size / m_Accel.sim_scale ));
	m_Accel.gridMax = m_Accel.bound_max;		m_Accel.gridMax += float(2.0*(m_Accel.grid_size / m_Accel.sim_scale ));
	m_Accel.gridSize = m_Accel.gridMax - m_Accel.gridMin;	
	
	float grid_size = m_Accel.grid_size;
	float world_cellsize = grid_size / m_Accel.sim_scale;		// cell spacing in world units
	float sim_scale = m_Accel.sim_scale;

	// Grid res - grid volume uniformly sub-divided by grid size
	m_Accel.gridRes.x = (int) ceil ( m_Accel.gridSize.x / world_cellsize );		// Determine grid resolution
	m_Accel.gridRes.y = (int) ceil ( m_Accel.gridSize.y / world_cellsize );
	m_Accel.gridRes.z = (int) ceil ( m_Accel.gridSize.z / world_cellsize );
	m_Accel.gridSize.x = m_Accel.gridRes.x * world_cellsize;						// Adjust grid size to multiple of cell size
	m_Accel.gridSize.y = m_Accel.gridRes.y * world_cellsize;
	m_Accel.gridSize.z = m_Accel.gridRes.z * world_cellsize;	
	m_Accel.gridDelta = Vector3DF(m_Accel.gridRes) / m_Accel.gridSize;		// delta = translate from world space to cell #	
	
	// Grid total - total number of grid cells
	m_Accel.gridTotal = (int) (m_Accel.gridRes.x * m_Accel.gridRes.y * m_Accel.gridRes.z);

	// Number of cells to search:
	// n = (2r / w) +1,  where n = 1D cell search count, r = search radius, w = world cell width
	//
	m_Accel.gridSrch = (int) (floor(2.0f*(m_Accel.psmoothradius / sim_scale) / world_cellsize) + 1.0f);
	if ( m_Accel.gridSrch < 2 ) m_Accel.gridSrch = 2;
	m_Accel.gridAdjCnt = m_Accel.gridSrch * m_Accel.gridSrch * m_Accel.gridSrch;
	m_Accel.gridScanMax = m_Accel.gridRes - Vector3DI( m_Accel.gridSrch, m_Accel.gridSrch, m_Accel.gridSrch );

	if ( m_Accel.gridSrch > 6 ) {
		dbgprintf ( "ERROR: Neighbor search is n > 6. \n " );
		exit(-1);
	}

	// Auxiliary buffers - prefix sums sizes
	int blockSize = SCAN_BLOCKSIZE << 1;
	int numElem1 = m_Accel.gridTotal;
	int numElem2 = int ( numElem1 / blockSize ) + 1;
	int numElem3 = int ( numElem2 / blockSize ) + 1;

	int numPoints = m_num_birds;

	int mem_usage = DT_CPU | DT_CUMEM;

	// Allocate acceleration
	m_Grid.DeleteAllBuffers ();
	m_Grid.AddBuffer ( AGRID,		  "grid",			sizeof(uint), numPoints,					mem_usage );
	m_Grid.AddBuffer ( AGRIDCNT,	"gridcnt",	sizeof(uint), m_Accel.gridTotal,	mem_usage );
	m_Grid.AddBuffer ( AGRIDOFF,	"gridoff",	sizeof(uint), m_Accel.gridTotal,	mem_usage );
	m_Grid.AddBuffer ( AAUXARRAY1, "aux1",		sizeof(uint), numElem2,						mem_usage );
	m_Grid.AddBuffer ( AAUXSCAN1,  "scan1",		sizeof(uint), numElem2,						mem_usage );
	m_Grid.AddBuffer ( AAUXARRAY2, "aux2",		sizeof(uint), numElem3,						mem_usage );
	m_Grid.AddBuffer ( AAUXSCAN2,  "scan2",		sizeof(uint), numElem3,						mem_usage );

	for (int b=0; b <= AAUXSCAN2; b++)
		m_Grid.SetBufferUsage ( b, DT_UINT );		// for debugging

	// Grid adjacency lookup - stride to access neighboring cells in all 6 directions
	int cell = 0;
	for (int y=0; y < m_Accel.gridSrch; y++ ) 
		for (int z=0; z < m_Accel.gridSrch; z++ ) 
			for (int x=0; x < m_Accel.gridSrch; x++ ) 
				m_Accel.gridAdj [ cell++]  = ( y * m_Accel.gridRes.z+ z ) * m_Accel.gridRes.x +  x ;			

	// Done
	dbgprintf ( "  Accel Grid: %d, Res: %dx%dx%d\n", m_Accel.gridTotal, (int) m_Accel.gridRes.x, (int) m_Accel.gridRes.y, (int) m_Accel.gridRes.z );		
}


void Sample::InsertIntoGrid ()
{
	int numPoints = m_num_birds;

	// Reset all grid cells to empty		
	memset( m_Grid.bufUI(AGRIDCNT),	0,	m_Accel.gridTotal*sizeof(uint));
	memset( m_Grid.bufUI(AGRIDOFF),	0,	m_Accel.gridTotal*sizeof(uint));

	memset( m_Birds.bufUI(FGCELL),	0,	numPoints*sizeof(int));
	memset( m_Birds.bufUI(FGNDX),		0,	numPoints*sizeof(int));

	float poff = m_Accel.psmoothradius / m_Accel.sim_scale;

	// Insert each particle into spatial grid
	Vector3DF gcf;
	Vector3DI gc;
	int gs; 
	Vector3DF ppos;
	uint* pgcell =	  m_Birds.bufUI (FGCELL);
	uint* pgndx =			m_Birds.bufUI (FGNDX);		

	Bird* b;
	
	for ( int n=0; n < m_num_birds; n++ ) {		
		
		b = (Bird*) m_Birds.GetElem(0, n);
		ppos = b->pos;

		gcf = (ppos - m_Accel.gridMin) * m_Accel.gridDelta; 
		gc = Vector3DI( int(gcf.x), int(gcf.y), int(gcf.z) );
		gs = (gc.y * m_Accel.gridRes.z + gc.z)*m_Accel.gridRes.x + gc.x;
	
		if ( gc.x >= 1 && gc.x <= m_Accel.gridScanMax.x && gc.y >= 1 && gc.y <= m_Accel.gridScanMax.y && gc.z >= 1 && gc.z <= m_Accel.gridScanMax.z ) {
			*pgcell = gs;
			*pgndx = *m_Grid.bufUI(AGRIDCNT, gs);
			(*m_Grid.bufUI(AGRIDCNT, gs))++;			
		} else {
			*pgcell = GRID_UNDEF;				
		}					
		pgcell++;
		pgndx++;		
	}


}

void Sample::PrefixSumGrid ()
{
	int numPoints = m_num_birds;
	int numCells = m_Accel.gridTotal;
	uint* mgrid = (uint*) m_Grid.bufI(AGRID);
	uint* mgcnt = (uint*) m_Grid.bufI(AGRIDCNT);
	uint* mgoff = (uint*) m_Grid.bufI(AGRIDOFF);

	// compute prefix sums for offsets
	int sum = 0;	
	for (int n=0; n < numCells; n++) {
		mgoff[n] = sum;
		sum += mgcnt[n];
	}

	// compute master grid list
	uint* pgcell =	  m_Birds.bufUI (FGCELL);
	uint* pgndx =			m_Birds.bufUI (FGNDX);		
	int gs, sort_ndx;
	for (int k=0; k < numPoints; k++) {
    mgrid[k] = GRID_UNDEF;
	}
	for (int j=0; j < numPoints; j++) {

		if ( *pgcell != GRID_UNDEF ) {			
			sort_ndx = mgoff [ *pgcell ] + *pgndx;
			mgrid[ sort_ndx ] = j;			
		} 
		pgcell++;
		pgndx++;
	}
}

void Sample::FindNeighbors ()
{
	// Find neighborhood of each bird to compute:
	// - near_pos - position of nearest bird
	// - ave_pos  - average centroid of neighbor birds
	// - ave_dir  - direction of neighbor birds	
	//

	float d = m_Accel.sim_scale;
	float d2 = d * d;
	float rd2 = (m_Accel.psmoothradius*m_Accel.psmoothradius) / d2;	
	int	nadj = (m_Accel.gridRes.z + 1)*m_Accel.gridRes.x + 1;
	uint j, cell;
	Vector3DF posi, posj, dist;	
	Vector3DF diri, dirj;
	Vector3DF cdir;
	float dsq;
	float nearest, nearest_fwd;
	
	uint*		grid		=	m_Grid.bufUI(AGRID);
	uint*		gridcnt = m_Grid.bufUI(AGRIDCNT);
	uint*   fgc     = m_Grid.bufUI(FGCELL);

	Bird *bi, *bj;
	
	float fov = cos ( 120 * RADtoDEG );
	float fov_fwd = 60;

	float ang, birdang;

	m_mark.clear ();	

	

	// for each bird..
	for (int i=0; i < m_Birds.GetNumElem(0); i++) {

		bi = (Bird*) m_Birds.GetElem(0, i);
		posi = bi->pos;

		// pre-compute for efficiency
		diri = bi->vel;			diri.Normalize();

		// mark for debug draw
		if (i == m_bird_sel ) {
			m_mark.push_back ( Vector4DF( posi, m_Accel.psmoothradius ) );
		}
		
		// clear current bird info
		bi->ave_pos.Set(0,0,0);
		bi->ave_vel.Set(0,0,0);
		bi->near_j = -1;
		bi->near_in = -1;
		bi->nbr_cnt = 0;

		nearest = rd2;
		nearest_fwd = rd2;

		// search neighbors
		int gc = m_Birds.bufUI(FGCELL)[i];
		if ( gc != GRID_UNDEF ) {

			gc -= nadj;

			for (int c=0; c < m_Accel.gridAdjCnt; c++) {
				cell = gc + m_Accel.gridAdj[c];
				int clast = m_Grid.bufUI(AGRIDOFF)[cell] + m_Grid.bufUI(AGRIDCNT)[cell];

				for ( int cndx = m_Grid.bufUI(AGRIDOFF)[cell]; cndx < clast; cndx++ ) {		

					  // get next possible neighbor
					  j = m_Grid.bufUI(AGRID)[cndx];
						if (i==j) continue;
						bj = (Bird*) m_Birds.GetElem(0, j );
						posj = bj->pos;

						dist = posi - posj;
						dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);

						if ( dsq < rd2 ) {
							// neighbor is within radius..
								
							// confirm bird is within forward field-of-view							
							dirj = posj - posi; dirj.Normalize();
							birdang = diri.Dot (dirj);

							if ( birdang > fov ) {

								// find nearest 
								dsq = sqrt(dsq);
								if ( dsq < nearest ) {
									nearest = dsq;
									bi->near_j = j;
								}
								// find nearest incoming (nbr /w opposing velocity)
								Vector3DF v1, v2;
								v1 = bi->vel;	v1.Normalize();
								v2 = bj->vel; v2.Normalize();
								ang = v1.Dot(v2);
								if ( dsq < nearest_fwd && ang < 0 ) {
									nearest_fwd = dsq;
									bi->near_in = j;
									bi->near_in_ang = ang;
								}
								// average neighbors
								bi->ave_pos += posj;
								bi->ave_vel += bj->vel;
								bi->nbr_cnt++;

								// mark for debug draw
								if (i == m_bird_sel ) {
									m_mark.push_back ( Vector4DF( posj, 1 ) );							  
								}
							}
						}
					}
			  }

			}	
		  if (bi->nbr_cnt > 0) {
				bi->ave_pos *= (1.0f / bi->nbr_cnt);
				bi->ave_vel *= (1.0f / bi->nbr_cnt);
			}
	}
}



float circleDelta (float b, float a)
{	
	float d = b-a;
	d = (d > 180) ? d-360 : (d<-180) ? d+360 : d;
	//float q = fmod( fmin ( b-a, 360+b-a ), 180 );
	//printf ("%f %f\n", d, q);
	return d;
}

void Sample::Advance ()
{
	Vector3DF fwd, up, right;
	Vector3DF force, torque, vaxis;
	Vector3DF diri, dirj;
	Quaternion ctrl_pitch;
	float airflow, dynamic_pressure;

	float m_LiftFactor = 0.005;
	float m_DragFactor = 0.003;

	float mass = 0.1;							// body mass (kg)

	float CL, L, dist, cd, vd;
	float pitch, yaw;
	Quaternion ctrlq, tq;
	Vector3DF angs;
	Quaternion angvel;
	Bird *b, *bj, *bfwd;
	bool leader;
	float ave_yaw;

	float safe_radius = 10.0;

	m_time++;

	Vector3DF flock_pnt = m_rnd.randV3(-200, 200);

	//--- Reynold's behaviors	
	//
  for (int n=0; n < m_Birds.GetNumElem(0); n++) {

		b = (Bird*) m_Birds.GetElem(0, n);
		b->clr.Set(0,0,0,0);

		if ( b->nbr_cnt == 0 ) {			
			continue;		
		} 

		diri = b->vel;			diri.Normalize();
		dirj = b->ave_pos - b->pos; dirj.Normalize();		
		leader = (diri.Dot (dirj) < 0);
		dirj *= b->orient.inverse();	
		ave_yaw = atan2 ( dirj.z, dirj.x )*RADtoDEG;

		// Turn entire flock
		/*if ( int(m_time) % 200 == 0 ) {
			if (b->near_in == -1) {				
				b->target.z = fmod( atan2( flock_pnt.z-b->pos.z, flock_pnt.x-b->pos.x ) * RADtoDEG, 180 );
				b->clr.Set(1,1,0,1);			
			}			
		}*/

	
		// Rule 1. Avoidance - avoid nearest bird
		//			
		// 1a. Side neighbor avoidance
		if ( b->near_j != -1) {
			// get nearest bird
			bj = (Bird*) m_Birds.GetElem(0, b->near_j);
			dirj = bj->pos - b->pos;
			dist = dirj.Length();		  

			if ( dist < safe_radius ) {	

				// Angular avoidance			
				dirj = (dirj/dist) * b->orient.inverse();							
				float ang_avoid = 0.20;			
				yaw = atan2( dirj.z, dirj.x )*RADtoDEG;
				pitch = asin( dirj.y )*RADtoDEG;						
				b->target.z -= yaw * ang_avoid / (dist*dist);
				b->target.y -= pitch * ang_avoid / (dist*dist);					

				// Power adjust				
				vd = (b->vel.Length() - bj->vel.Length()) * 8.0;
				float np = 2 - vd;								
				b->power = np;
				/*if ( fabs(b->power-np) > 2 )  {	
					b->target.z -= ave_yaw ;
				}*/
			}			
		}

		// 1b. Incoming neighbor avoidance
		if ( b->near_in != -1 ) {
			// get incoming bird
			bj = (Bird*) m_Birds.GetElem(0, b->near_in);						
			dirj = bj->pos - b->pos;
			dist = dirj.Length();		

			if ( dist < safe_radius ) {
				
				// Direct collision avoidance
				float s = (m_rnd.randF(0,1) < 0.5) ? 1 : -1;			
				b->target.z += (b->near_in_ang < 0 ) ? 2 * s : 0;		

				// Power avoidance					
				vd = (b->vel.Length() + bj->vel.Length()) * 10.0;
				float np = 2 - pow(vd, 3);				
				//if ( fabs(b->power-np) > 2 ) {					
					b->power = np;												
				
			}					
		}

		if (b->power < -40) b->power = -40;
		if (b->power > 40) b->power = 40;	


	  // Rule 2. Alignment - orient toward average direction		
		dirj = b->ave_vel;
		dirj.Normalize();
		dirj *= b->orient.inverse();		// using inverse orient for world-to-local xform		
		yaw = atan2( dirj.z, dirj.x )*RADtoDEG;
		pitch = asin( dirj.y )*RADtoDEG;
		b->target.z += yaw   * 0.026;
		b->target.y += pitch * 0.026;		 

		// Rule 3. Cohesion - steer toward neighbor centroid
		dirj = b->ave_pos - b->pos;
		dirj.Normalize();
		dirj *= b->orient.inverse();		// using inverse orient for world-to-local xform		
		yaw = atan2( dirj.z, dirj.x )*RADtoDEG;
		pitch = asin( dirj.y )*RADtoDEG;
		b->target.z += yaw   * 0.003;
		b->target.y += pitch * 0.003;		
		
	}

	

	//--- Flight model
	//
	for (int n=0; n < m_Birds.GetNumElem(0); n++) {

		b = (Bird*) m_Birds.GetElem(0, n);

		// Body orientation
		fwd = Vector3DF(1,0,0) * b->orient;			// X-axis is body forward
		up  = Vector3DF(0,1,0) * b->orient;			// Y-axis is body up
		right = Vector3DF(0,0,1) * b->orient;		// Z-axis is body right

		// Direction of motion
		b->speed = b->vel.Length();
		vaxis = b->vel / b->speed;	
		if ( b->speed < m_min_speed ) b->speed = m_min_speed;				// birds dont go in reverse
		if ( b->speed > m_max_speed ) b->speed = m_max_speed;
		if ( b->speed==0) vaxis = fwd;
			
		b->orient.toEuler ( angs );			

		// Target corrections
		angs.z = fmod (angs.z, 180.0 );												
		b->target.z = fmod ( b->target.z, 180 );										// yaw -180/180
		b->target.x = circleDelta(b->target.z, angs.z) * 0.5;				// banking
		b->target.y *= 0.999;																				// level out
		if ( b->target.y > 40 ) b->target.y = 40;
		if ( b->target.y < -40 ) b->target.y = -40;
		if ( fabs(b->target.y) < 0.0001) b->target.y = 0;
		
		float reaction_delay = 0.0008;
		//float reaction_delay = 0.0010;

		// Roll - Control input
		// - orient the body by roll
		ctrlq.fromAngleAxis ( (b->target.x - angs.x) * reaction_delay, Vector3DF(1,0,0) * b->orient );
		b->orient *= ctrlq;	b->orient.normalize();
		
		//if (b->pos.x > 200) {
		//	printf ( "%f %f %f\n", b->target.z, angs.z, circleDelta(b->target.z, angs.z)  );
		//}

		// Pitch & Yaw - Control inputs
		// - apply 'torque' by rotating the velocity vector based on pitch & yaw inputs				
		ctrlq.fromAngleAxis ( circleDelta(b->target.z, angs.z) * reaction_delay, Vector3DF(0,-1,0) * b->orient );
		vaxis *= ctrlq; vaxis.Normalize();	
		ctrlq.fromAngleAxis ( (b->target.y - angs.y) * reaction_delay , Vector3DF(0,0,1) * b->orient );
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

		float cd;	

		// Boundaries
		/*if ( b->pos.x < m_Accel.bound_min.x ) b->target.z = 0;
		if ( b->pos.x > m_Accel.bound_max.x ) b->target.z = 179;
		if ( b->pos.z < m_Accel.bound_min.z ) b->target.z = 90;
		if ( b->pos.z > m_Accel.bound_max.z ) b->target.z = -90;*/

		/*Vector3DF pr ( b->pos.x, 0, b->pos.z );
		f ( pr.Length() > m_Accel.bound_max.x ) {
			b->target.z += 0.5;  //= fmod( atan2( -pr.z, -pr.x ) * RADtoDEG, 180 );
		}*/
		
		if ( b->pos.x < m_Accel.bound_min.x ) b->pos.x = m_Accel.bound_max.x;
		if ( b->pos.x > m_Accel.bound_max.x ) b->pos.x = m_Accel.bound_min.x;
		if ( b->pos.z < m_Accel.bound_min.z ) b->pos.z = m_Accel.bound_max.z;
		if ( b->pos.z > m_Accel.bound_max.z ) b->pos.z = m_Accel.bound_min.z;			  

		// Ground avoidance
		cd = b->pos.y - m_Accel.bound_min.y;
		if ( cd < 50 ) {			
			cd = (50-cd)/50.0f;
			b->target.y += cd * 0.5;			
			b->power = 4;
		} 
		
		// Ceiling avoidance
		cd = m_Accel.bound_max.y - b->pos.y;
		if ( cd < 50  ) {	
			cd = (50-cd)/50.0f;									
			b->target.y -= cd * 0.1; 						
		} 

		if (b->nbr_cnt==0) {						
			b->target.z = fmod( atan2( b->pos.z , b->pos.x ) * RADtoDEG + 180, 180 );			
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

void Sample::Run ()
{

	InsertIntoGrid ();

	PrefixSumGrid ();

	FindNeighbors ();

	Advance ();

}


void Sample::DrawAccelGrid ()
{
	Vector3DF r,a,b;
	float v;

	uint* gc = (uint*) m_Grid.bufUI(AGRIDCNT);

	for (r.y=0; r.y < m_Accel.gridRes.y; r.y++) {
		for (r.z=0; r.z < m_Accel.gridRes.z; r.z++) {
			for (r.x=0; r.x < m_Accel.gridRes.x; r.x++) {
				
				a = m_Accel.gridMin + r / m_Accel.gridDelta;
				b = a + (Vector3DF(0.99f,0.99f,0.99f) / m_Accel.gridDelta );								

				v = fmin(1.0, float(*gc)/2.0f);

				drawBox3D ( a, b, v, 1-v, 1-v, 0.02 + v );

				gc++;
			}
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

	m_cam_fwd = m_cam_fwd * 0.99f + fwd * 0.01f;
	m_cam_fwd.Normalize();

	// Set eye level above centerline
	Vector3DF p = b->pos + Vector3DF(0,2,0);	  
	
	m_cam->setDirection ( p, p + m_cam_fwd, -angs.x );
}


bool Sample::init ()
{
	int w = getWidth(), h = getHeight();			// window width &f height
	m_run = true;
	m_flightcam = true;	
	m_cockpit_view = false;
	m_draw_sphere = false;
	m_draw_grid = false;
	m_time = 0;
	m_rnd.seed (12);

	m_bird_sel = 154;

	addSearchPath ( ASSET_PATH );
	init2D ( "arial" );
	setview2D ( w, h );	
	setText ( 16, 1 );		
	
	m_cam = new Camera3D;
	m_cam->setFov ( 120 );
	m_cam->setNearFar ( 1.0, 100000 );
	m_cam->SetOrbit ( Vector3DF(-30,30,0), Vector3DF(0,100,0), 300, 1 );

	// Initialize birds
	// * birds are placed into a DataX structure to allow
	// for easy sharing between CPU and GPU
	
	//m_num_birds = 5000;
	m_num_birds = 1800;
	//m_num_birds = 400;

  Reset ();

	m_Accel.bound_min = Vector3DF(-200,   0, -200);
	m_Accel.bound_max = Vector3DF( 200, 200,  200);
	m_Accel.psmoothradius = 20;
	m_Accel.grid_density = 1.0;
	m_Accel.sim_scale = 1.0;

	InitializeGrid ();

	// speed limits
	m_min_speed = 20.0;
	m_max_speed = 100.0;	

	m_DT = 0.005;
	//m_DT = 0.0025;
	m_wind.Set (0, 0, 0);

	return true;
}


void Sample::display ()
{	
	char msg[2048];
	Vector3DF x,y,z;
	Vector3DF pnt;
	Vector4DF clr;
	int w = getWidth();
	int h = getHeight();
		
	bool draw_vis = true;

	Bird* b;

	if (m_run) { 		
		Run ();
	}	

	/*if (m_cockpit_view) {
		CameraToCockpit ( m_bird_sel);
	} else {
		CameraToBird ( m_bird_sel );
  }*/

	if ( draw_vis ) {
		glClearColor(0,0,0,1);
	} else {
		glClearColor(.8,.8,.9,1);
	}
	clearGL();

	start3D(m_cam);		

		setLight ( S3D, 0, 400, 0, 0.1, 0.1, 0.1 );	
		setMaterial ( S3D, Vector3DF(0,0,0), Vector3DF(.2,.2,.2), 40, 0 );

		// Draw debug marks
		if (m_draw_sphere) {

			Vector3DF cn, p;
			for (int k=0; k < m_mark.size(); k++) {
				p = Vector3DF(m_mark[k]);			
				drawCircle3D ( p, m_cam->getPos(), m_mark[k].w, Vector4DF(1,1,0,1) );
			}
			b = (Bird*) m_Birds.GetElem(0, m_bird_sel );
			drawCircle3D ( b->pos, m_cam->getPos(), 1.1, Vector4DF(0,1,0,1) );
			if (b->near_j != -1 ) {
				b = (Bird*) m_Birds.GetElem(0, b->near_j );
				drawCircle3D ( b->pos, m_cam->getPos(), 1.2, Vector4DF(1,0,0,1) );
			}
			if (b->near_in != -1 ) {
				b = (Bird*) m_Birds.GetElem(0, b->near_in );
				drawCircle3D ( b->pos, m_cam->getPos(), 1.2, Vector4DF(1,0,1,1) );
			}
		}

		// Draw acceleration grid
		if (m_draw_grid) {
			drawBox3D ( m_Accel.bound_min, m_Accel.bound_max, 0,1,1,0.5 );
			DrawAccelGrid ();
		}

		// Draw ground plane
		if ( draw_vis ) {
			drawLine3D ( Vector3DF(0,0,0), Vector3DF(100,0,0), Vector4DF(1,0,0,1));
			drawLine3D ( Vector3DF(0,0,0), Vector3DF(  0,0,100), Vector4DF(0,0,1,1));
			drawGrid( Vector4DF(0.1,0.1,0.1, 1) );
		}

		for (int n=0; n < m_Birds.GetNumElem(0); n++) {

			b = (Bird*) m_Birds.GetElem(0, n);
		
			x = Vector3DF(1,0,0) * b->orient;
			y = Vector3DF(0,1,0) * b->orient;
			z = Vector3DF(0,0,1) * b->orient;


			if ( draw_vis ) {
				
				// visualize velocity
				float v = (b->vel.Length() - m_min_speed) / (m_max_speed-m_min_speed);			
				float v2 = (b->power + 10)/20.f;
				if (b->clr.w==0) {
					drawLine3D ( b->pos,		b->pos + (b->vel*0.05f),	Vector4DF(v2, 1-v2,1-v2,1) );
				} else {
					drawLine3D ( b->pos,		b->pos + (b->vel*0.05f),	b->clr );
				}

			} else {
				// bird dart
				Vector3DF p,q,r,t;
				p = b->pos - z*1;
				q = b->pos + z*1;
				r = b->pos + x*2; 
				t = y;				
				drawTri3D ( p.x,p.y,p.z, q.x,q.y,q.z, r.x,r.y,r.z, t.x,t.y,t.z, 1,1,1,1 );
			}
			

			/*drawLine3D ( b->pos,	  b->pos+x,					Vector4DF(1,1,0,  1) );					// fwd				
			drawLine3D ( b->pos,		b->pos +b->lift + x*0.1f,		Vector4DF(0,1,0,  1) );	// lift (green)
			drawLine3D ( b->pos,		b->pos+b->thrust, Vector4DF(1,0,0,  1) );						// thrust (red)
			drawLine3D ( b->pos,		b->pos+b->drag,		Vector4DF(1,0,1,  1) );			
			drawLine3D ( b->pos,		b->pos+b->force,	Vector4DF(0,1,1,  1) );			*/		

		}
	end3D();

	/*start2D ();
	setview2D ( w, h );

		b =  (Bird*) m_Birds.GetElem(0, m_bird_sel);
		sprintf ( msg, "%f %f %f, %f\n", b->target.x, b->target.y, b->target.z, b->speed );
		drawText ( 10, 10, msg, 1,1,1,1);

	end2D(); */

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
	case 's': m_draw_sphere = !m_draw_sphere; break;
	case 'g': m_draw_grid = !m_draw_grid; break;
  case 'c': 		
		m_cockpit_view = !m_cockpit_view; 
		//m_cam_orient = 
		break;
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
	printf ( "%d \n", m_bird_sel );
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


