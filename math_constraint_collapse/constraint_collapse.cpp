//--------------------------------------------------------------------------------
// JUST MATH:
// Constraint Collapse
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
// Sample utils
#include <algorithm>
#include "mersenne.h"
#include "dataptr.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <vector>
#include <string>
#include <map>

#include "constraint_collapse.h"

void ConstraintCollapse::AllocBuf (int id, int cnt) 
{
  int flags = m_run_cuda ? (DT_CPU | DT_CUMEM) : DT_CPU;
  m_buf[id].Resize( sizeof(float), cnt, 0x0, flags );
  memset( (void *)(m_buf[id].getPtr(0)), 0, sizeof(float)*cnt);
}
void ConstraintCollapse::AllocF (int id, int nbrs, int vals) 
{
  uint64_t cnt = vals * vals * nbrs;    // size B*B*6
  int flags = m_run_cuda ? (DT_CPU | DT_CUMEM) : DT_CPU;
  m_buf[id].Resize( sizeof(float), cnt, 0x0, flags );
  memset( (void *)(m_buf[id].getPtr(0)), 0, sizeof(float)*cnt);
}
void ConstraintCollapse::ZeroBuf (int id)
{
  m_buf[id].FillBuffer ( 0 );
}

int64_t ConstraintCollapse::getVertex (int x, int y, int z) 
{
  return int64_t(z*m_res.y + y)*m_res.x + x;
}

Vector3DI ConstraintCollapse::getVertexPos(int64_t j) {
  Vector3DI p;
  p.z = j / (m_res.x*m_res.y);  j -= p.z * (m_res.x*m_res.y);
  p.y = j / m_res.x;            j -= p.y * m_res.x;
  p.x = j;
  return p;
}

int64_t ConstraintCollapse::getNeighbor( uint64_t j, int nbr ) 
{

  Vector3DI jp = getVertexPos(j);

  // 3D spatial neighbor function
  //
  switch (nbr) {
  case 0:    return (jp.x < m_res.x-1) ?  j+1 : -1;
  case 1:    return (jp.x > 0) ?          j-1 : -1;
  case 2:    return (jp.y < m_res.y-1) ?  j+m_res.x : -1;
  case 3:    return (jp.y > 0) ?          j-m_res.x : -1;
  case 4:    return (jp.z < m_res.z-1) ?  j+(m_res.x*m_res.y) : -1;
  case 5:    return (jp.z > 0) ?          j-(m_res.x*m_res.y) : -1;
  };
  return -1;
}

int64_t ConstraintCollapse::getFace (int64_t a, int nbr )
{
    Vector3DI ai = getVertexPos(a);
    int rx = m_res.x+1;
    int ry = m_res.y+1;    
    int64_t b = int64_t(ai.z*ry + ai.y)*rx + ai.x;

    int64_t f = -1;
    switch (nbr) {
    case 0: f = (b+1)*3 + 0;        break;
    case 1: f = (b+0)*3 + 0;        break;
    case 2: f = (b+rx)*3 + 1;       break;
    case 3: f = (b+0      )*3 + 1;  break;
    case 4: f = (b+(rx*ry))*3 + 2;  break;
    case 5: f = (b+0      )*3 + 2; break;
    default:
        printf ("ERROR on getFace\n" );
        exit(-8);
        break;
    };
    return f;
}


// GetVertexConstraints
// - this function returns a bitmask for the 6x faces of the current vertex
int ConstraintCollapse::GetVertexConstraints ( int64_t p )
{
    return getVal( BUF_C, p );    
}

// this could be much faster with lookup table
int ConstraintCollapse::CountBits ( int mask )
{
    int c = 0;
    c += (mask   ) & 0x1;
    c += (mask>>1) & 0x1;
    c += (mask>>2) & 0x1;
    c += (mask>>3) & 0x1;
    c += (mask>>4) & 0x1;
    c += (mask>>5) & 0x1;
    return c;
}

void ConstraintCollapse::start ()
{
    m_rand.seed( m_seed );

    randomize ();
    
    write_admissible ();  

    m_resolve = 1.0;
    m_cnt = 0;
    m_step = 0;
    m_temp = 0.5;

    for (int n=0; n < 512; n++) 
        m_clr[n] = m_rand.randV3();
}

void ConstraintCollapse::info()
{
    dbgprintf ( "step %d, constraints %d, noise %f, stuck %d\n", m_step, m_cnt, m_temp, m_stuck_cnt);
}

int ConstraintCollapse::single_step ()
{
    m_temp = float(m_cnt/6) / m_num_verts;

    // generate noise
    generate_noise ();
   
    // fix constraints    
    int cnt = check_constraints ();

    fix_constraints ( m_stuck_cnt >= 3 );    
    
    cnt = check_constraints ();
    if (cnt >= m_cnt) {
        m_stuck_cnt++;
        if (m_stuck_cnt > 3 ) {
            m_stuck_cnt = 0;
        }
    } 

    m_cnt = cnt;
    m_step++;

    return cnt;
}

void ConstraintCollapse::generate_noise ()
{
    // constraint-based noise

    float c, r, rd, dist;
    int t, vnbr;
    Vector3DI ki, kv, kd, kmin, kmax;
    Vector3DF dir;
    int64_t k;
    
    float r_noise = 0.05;         // confidence noise (decrease)
    float flip_noise = 0.07;      // flip noise

    int p;
    int pmin = (m_step % 2) ? 0 : getNumVerts()-1;
    int pmax = (m_step % 2) ? getNumVerts() : -1;
    int pdel = (m_step % 2) ? 1 : -1;
    
    for (p = pmin; p != pmax; p += pdel) {

        r = getVal(BUF_R, p );   

        c = GetVertexConstraints( p );        

        // if cell has constraints..
        if (c > 0) {

            kv = getVertexPos( p );

            rd = fmax(1, fmin(5, m_temp * 20));            

            // bias noise toward center of map
            // based on direction & temp (to handle large maps)
            dir = Vector3DF(m_res)*0.5f - kv;
            dir.z = 0;
            //dist = m_rand.randF( 1, dir.Length() );
            dist = m_rand.randF( 1, rd);            
            dir.Normalize();
            dir *= dist;
            dir.z = 0;
            
            ki = kv + dir;

            if (ki.x >= 0 && ki.x < m_res.x && ki.y >= 0 && ki.y < m_res.y ) {

                k = getVertex(ki.x, ki.y, ki.z);
                r = getVal(BUF_R, k);
                r = r - m_rand.randF() * r_noise; // * (dist / (rd*rd));  

                // if cell has uncertainty..
                if (m_rand.randF() < flip_noise * (1-r) ) {
                    t = m_rand.randI( m_num_values );       
                    r += r_noise * 2.0;
                    if (t != 0 && t!=87 ) {
                        SetVal (BUF_T, k, t );
                    }
                }                    
                                       
                r = (r<0) ? 0 : (r>1) ? 1 : r;   
                SetVal(BUF_R, k, r );
            }

            // generate noise in 2x2 region
            // - necessary for unconstrained neighbors to respond to constrainted cells
            // - bias region toward center
            if ( kv.x > m_res.x/2 ) {
                kmin.x = -1;
                kmax.x =  0;
            } else {
                kmin.x =  0;
                kmax.x =  1;
            }
            if ( kv.y > m_res.y/2 ) {
                kmin.y = -1;
                kmax.y =  0;
            } else {
                kmin.y =  0;
                kmax.y =  1;
            }

            // top-left biased noise generation
            // *presses* the errors out toward the edges and top
            for (kd.y = kmin.y; kd.y <= kmax.y; kd.y++) {
              for (kd.x = kmin.x; kd.x <= kmax.x; kd.x++) {

                ki = kv + kd;
                dist = (kd.x*kd.x + kd.y*kd.y);

                if (ki.x >= 0 && ki.x < m_res.x && ki.y >= 0 && ki.y < m_res.y ) {

                    k = getVertex(ki.x, ki.y, ki.z);
                    r = getVal(BUF_R, k);
                    r = r - m_rand.randF() * r_noise; // * (dist / (rd*rd));  

                    // if cell has uncertainty..
                    if (m_rand.randF() < flip_noise * (1-r) ) {
                        t = m_rand.randI( m_num_values );       
                        r += r_noise * 2.0;
                        if (t != 0 && t!=87 ) {
                            SetVal (BUF_T, k, t );
                        }
                    }                    
                                       
                    r = (r<0) ? 0 : (r>1) ? 1 : r;   
                    SetVal(BUF_R, k, r );
                }
              }
            }  

        }
    }
}



void ConstraintCollapse::decimate ()
{
   float p;
   float r, r2;
   int t;

   dbgprintf ( "decimate.\n");

   for (int64_t v=0; v< getNumVerts(); v++) {
       
       r = getVal(BUF_R, v );       
       
       r = (r<0) ? 0 : pow(r, 1.1);       
       //r = r - 0.1;
       //if (r<0) r= 0;
       SetVal(BUF_R, v, r );

       // probability of switching this tile       
       /*if ( m_rand.randF() > pow(r, 1) ) {
           t = m_rand.randF() * m_num_values;
           SetVal(BUF_T, v, t);
       }*/
   }
}

void ConstraintCollapse::check_r ()
{
   float r, r2;
   
   float* mem = ((float*) m_buf[BUF_R].mCpu) + 128;
   
   for (int64_t p=0; p< getNumVerts(); p++) {
       r = getVal(BUF_R, p );
       if (isnan(r)) {
           dbgprintf ("ERROR: %d, step: %d\n", int(p), m_step);
       }
   }
}

void ConstraintCollapse::randomize ()
{
    float t;
    for (int64_t p=0; p< getNumVerts(); p++) {
        t = m_rand.randF() * m_num_values;
        SetVal ( BUF_T, p, t );         // set random tile
    }
}

void ConstraintCollapse::write_admissible ()
{
    std::vector< int32_t > v32;
    int64_t p;

    for (int n=0; n < m_admissible.size(); n++) {
        
        //if (n==0 || n==5) {
            v32 = m_admissible[n];
            p = getVertex( v32[0], v32[1], v32[2] );
            SetVal ( BUF_T, p, v32[3] );    // tile value
            SetVal ( BUF_R, p, 1.0 );       // fully resolved 
        //}
    }
}
 
int ConstraintCollapse::check_constraints ()
{
    int cnt = 0;
    for (int64_t p=0; p< getNumVerts(); p++) {
        cnt += check_constraints ( p );
    }    
    return cnt;
}        

/* -- OLD CODE, vert to face indexing
Vector3DI ai = getVertexPos(p);
int rx = m_res.x+1;
int ry = m_res.y+1;
int64_t b = int64_t(ai.z*ry + ai.y)*rx + ai.x;    */


int ConstraintCollapse::check_constraints ( int64_t p )
{
    int a, b, cnt, msk;
    float rule;
    int64_t pnbr, f;
    Vector3DI pi;

    // tile value at p
    a = getVal( BUF_T, p ); 

    cnt = 0;
    msk = 0;

    // constraint mask
    // bit order:    
    // 3 = nbr(0), x+1 - high bit
    // 2 = nbr(1), x-1 
    // 1 = nbr(2), y+1
    // 0 = nbr(3), y-1 - low bit  

    for (int nbr=0; nbr < getNumNeighbors(); nbr++) {

        // tile value at neighbor of p
        pnbr = getNeighbor(p, nbr);
        if ( pnbr != -1) {
            // get neighbor tile
            b = getVal( BUF_T, pnbr);        
        } else {
            b = 0;
        }
        rule = getValF ( BUF_F, a, b, nbr );   // rule for b->a                    
        
        // constraint in this direction (0 weight = disallowed)
        if (rule==0) {
            cnt++;
            msk |= 1;
        }        
        msk <<= 1;   
    }
    msk >>= 1;
 
    SetVal( BUF_C, p, msk );        

    return cnt;
}

void ConstraintCollapse::GetMaxResolved ( int64_t p, int prev, int next, float r_edge, float& r_fixed, float& r_max)
{
    int fixed;
    float r[6];
    Vector3DI pi = getVertexPos( p );

    // constraint mask
    // bit order:    
    // 3 = nbr(0), x+1 - high bit
    // 2 = nbr(1), x-1 
    // 1 = nbr(2), y+1
    // 0 = nbr(3), y-1 - low bit

    // get resolved values of neighbors
    // (reversed order because of bitmask)
    r[0] = (pi.x < m_res.x-1) ? getVal(BUF_R, p+1) : r_edge;
    r[1] = (pi.x > 0)         ? getVal(BUF_R, p-1) : r_edge;
    r[2] = (pi.y < m_res.y-1) ? getVal(BUF_R, p+m_res.x) : r_edge;
    r[3] = (pi.y > 0 )        ? getVal(BUF_R, p-m_res.x) : r_edge;
    // r[4] = (pi.z < m_res.z-1) ? getVal(BUF_R, p+(m_res.x*m_res.y)) : r_edge;
    // r[5] = (pi.z > 0 )        ? getVal(BUF_R, p-(m_res.x*m_res.y)) : r_edge;
    
    // identify maximal resolved neighbor that 
    // has eliminated a constraint
    r_fixed = 0;
    r_max = r[0];       
    for (int j=0; j < getNumNeighbors(); j++) {
      fixed = ((prev >> j) & 0x1) - ((next >> j) & 0x1);
      if ( fixed==1 ) {
          r_fixed += r[j];
          if ( r[j] > r_max) r_max = r[j];          
      }
    }    
}

void ConstraintCollapse::fix_constraints (bool stuck)
{
    bool noise, empty;
    float c[6];
    float r0, r, new_r, best_r, fix_r, pb;
    int64_t p;
    Vector3DI v;
    int chk, fixed;
    int save_t, save_mask, start_t, tr;
    int new_t, new_mask, new_cnt;
    int best_t, best_mask, best_cnt;
    bool border;

    int pmin = (m_step % 2) ? 0 : getNumVerts()-1;
    int pmax = (m_step % 2) ? getNumVerts() : -1;
    int pdel = (m_step % 2) ? 1 : -1;

    for (p = pmin; p != pmax; p += pdel) {

        v = getVertexPos(p);

        border = (v.x==0 || v.y==0 || v.x==m_res.x-1 || v.y==m_res.y-1 );
        
        // process on checkerboard grid
        // chk = (v.x + v.y + v.z + m_flip) % 2;  
        // if (chk==0) continue;

        // check if fully resolved
        r = getVal( BUF_R, p );   
        r0 = r;
        if ( r >= 1.0 ) continue;

        // count unresolved face constraints at vertex        
        best_mask = GetVertexConstraints( p );    
        best_cnt = CountBits(best_mask); 
        
        if ( best_cnt==0 ) {
            r += 0.00001;
            r = (r<0) ? 0 : (r>1) ? 1 : r;
            SetVal ( BUF_R, p, r ); 
            continue;
        }

        best_t = getVal( BUF_T, p);
        best_r = r;
        save_mask = best_mask;  // save current configuration
        save_t = best_t;

        // MC - basic monte-carlo. test tile values at random.
        start_t = save_t + m_rand.randI( m_num_values );

        // run through tile values and find one that solves the most constraints               
        for (int t=1; t < m_num_values ; t++) { 

            new_t = (start_t + t) % m_num_values;     // <-- bug. this biases for earlier (X+) tiles.
            //new_t = m_rand.randI( m_num_values ) ;

            if (new_t == save_t) continue;

            if (new_t == 0 || new_t==87 ) continue;            

            // tile weight biasing
            //pb = getVal(BUF_G, new_t) * m_num_values;
            //if ( m_rand.randF() > pb ) continue;

            // test a different tile value
            SetVal ( BUF_T, p, new_t );  

            // check face constraints at this vertex
            check_constraints ( p );            
            new_mask = GetVertexConstraints( p );
            new_cnt = CountBits(new_mask);

            // count number of constraints resolved            
            // fix   = number of new face constraints resolved            
            int fix = best_cnt - new_cnt;

            GetMaxResolved ( p, save_mask, new_mask, 0.01, fix_r, new_r );           
            
            // MCMC - Markov Chain Monte Carlo - next step should fix more constraints than current.            
            // Accept tile.. if:
            // - (fix>=1) ==> tile value solves 1 or more faces, accept immediately
            // - (fix==0 && fix_r>best_r) ==> tile value don't solve any self-constraints, but neighbors (fix_r) are resolved          

            if (border) new_r = 1.0;

            bool noise = m_rand.randF(0,1) > r;
            
            if ( fix >= 1 || (fix==0 && new_r >= best_r) ) {            
            //if ( fix >= 1 || (fix==0 && new_r > best_r) || (fix==-1 && stuck && noise) ) {
                best_mask = new_mask;
                best_cnt = new_cnt;
                best_t = new_t;
                best_r = new_r;
            }            
        }    

        if ( best_t != save_t ) {  

            // get best resolved neighbor    
            //fixed = GetMaxResolved ( p, save_mask, best_mask, fix_r, best_r );            
            r = best_r * 0.85;
            r = (r<0) ? 0 : (r>1) ? 1 : r;

            SetVal( BUF_R, p, r );

            SetVal( BUF_T, p, best_t );           
           
        } else {            
            SetVal( BUF_T, p, save_t ); 
        }

        check_constraints ( p );
    }

    m_flip = 1 - m_flip;
}

Vector4DF ConstraintCollapse::getSample ( int buf, int64_t v )
{
    Vector4DF s;
    int t, c, f[6];
    float x, a;

    Vector3DI vi = getVertexPos(v);

    switch (buf) {
    case BUF_T: 
        t = getVal(buf, v);
        a = pow( getVal(BUF_R, v), 0.1 );
        if (t==0) a=0;
        s = Vector4DF( m_clr[t], a);
        break;
    case BUF_C:
        t = getVal(BUF_T,v);
        c = CountBits( GetVertexConstraints( v ) );
        a = getVal(BUF_R, v);
        x = float(c); // 6.0;
        s = Vector4DF(x,x,x,1);
        break;
    case BUF_E:
        x = getVal(buf, v);        
        s = Vector4DF(x,x,x,x);
        break;
    case BUF_R:
        x = getVal(buf, v);
        //x = pow(x, 0.1 );
        s = Vector4DF(x,x,x,x);
        break;
   };
   return s;
}


//----
int ConstraintCollapse::init (int Rx, int Ry, int Rz, std::string &name_fn, std::string &rule_fn) 
{
  int i, j, ret = 0;

  m_dir_inv[0] = 1;
  m_dir_inv[1] = 0;
  m_dir_inv[2] = 3;
  m_dir_inv[3] = 2;
  m_dir_inv[4] = 5;
  m_dir_inv[5] = 4;

  m_num_values = read_F_CSV(rule_fn, name_fn);
  if (m_num_values <= 0 ) {
     printf ( "ERROR: num tiles values = %d\n", m_num_values );
     exit(-7);
  }
  printf(">>>> # of tiles: %i\n", (int) m_num_values );

  assert ( m_num_values == m_tile_name.size() );

  //---

  m_rand.seed ( m_seed++ );

  m_res.Set ( Rx, Ry, Rz );
  m_num_verts = m_res.x * m_res.y * m_res.z;
  m_num_face = (m_res.x+1) * (m_res.y+1) * (m_res.z+1) * 3;

  m_res.Set ( Rx, Ry, Rz );
  m_flip = 0;

  //-- Allocation
  AllocBuf ( BUF_T, m_num_verts );  
  AllocBuf ( BUF_R, m_num_verts );  
  AllocBuf ( BUF_C, m_num_verts ); 
  AllocBuf ( BUF_E, m_num_verts );  

  ZeroBuf ( BUF_T );
  ZeroBuf ( BUF_R );
  ZeroBuf ( BUF_C );
  ZeroBuf ( BUF_E );

  // options
  m_run_cuda  = false;

  return 0;
}

int ConstraintCollapse::read_F_CSV(std::string &rule_fn, std::string &name_fn) {
  int i, ret;
  int b, maxb=-1, B;

  std::vector< std::vector<float> > tile_rule;
  std::vector< float> tile_weight;

  ret = _read_name_csv(name_fn, m_tile_name, tile_weight );
  if (ret < 0) { return ret; }
  ret = _read_rule_csv(rule_fn, tile_rule);
  if (ret < 0) { return ret; }

  // use tile_rule to determine maximum number of tiles
  //
  for (i=0; i<tile_rule.size(); i++) {
    b = (int)tile_rule[i][0];
    if (b>maxb) { maxb = b; }
    b = (int)tile_rule[i][1];
    if (b>maxb) { maxb = b; }
  }
 
  B = maxb+1;

  printf(">>>> # of tile values (%i)\n", (int) B );

  //---

  // F
  //
  // tile_rule[i][0] = from tile value 
  // tile_rule[i][1] = to tile value
  // tile_rule[i][2] = direction
  // tile_rule[i][3] = weight
  
  AllocF ( BUF_F, 6, B );
  memset( m_buf[BUF_F].getData(), 0, 6*B*B*sizeof(float) );  

  m_num_values = B;      // must be set here to SetValF function to work!
  for (i=0; i<tile_rule.size(); i++) {
    SetValF( BUF_F, tile_rule[i][0], tile_rule[i][1], tile_rule[i][2], tile_rule[i][3] );
  }

  AllocBuf ( BUF_G, B );
  for (i=0; i<m_num_values; i++) {
    if (i < tile_weight.size()) {
      SetVal( BUF_G, i, tile_weight[i] );
    }
  }

  return B;
}




//----

int _read_line(FILE *fp, std::string &line) {
  int ch=0, count=0;

  while (!feof(fp)) {
    ch = fgetc(fp);
    if (ch == '\n') { break; }
    if (ch == EOF) { break; }
    line += (char)ch;
    count++;
  }
  return count;
}

int _read_name_csv(std::string &fn, std::vector<std::string> &name, std::vector< float > &weight) {
  int i, idx;
  FILE *fp;
  std::string line, tok, _s;
  std::vector<std::string> toks;
  float w;

  fp = fopen(fn.c_str(), "r");
  if (!fp) { return -1; }

  while (!feof(fp)) {
    line.clear();
    _read_line(fp, line);

    if (line.size()==0) { continue; }
    if (line[0] == '#') { continue; }
    if (line[0] == ' ') { continue; }

    toks.clear();
    tok.clear();
    for (i=0; i<line.size(); i++) {
      if (line[i]==',') {
        toks.push_back(tok);
        tok.clear();
        continue;
      }
      tok += line[i];
    }
    toks.push_back(tok);

    idx = atoi(toks[0].c_str());
    if (idx <= name.size()) {
        for (i=name.size(); i<=idx; i++) {
            _s.clear();
            name.push_back(_s);
            weight.push_back(0.0);
        }
    }    

    name[idx] = toks[1];
    
    w = 1.0;
    if (toks.size() > 2) { w = atof(toks[2].c_str()); }
    weight[idx] = w;
  }

  fclose(fp);
    
  w = 0.0;
  for (i=0; i<weight.size(); i++) { w += weight[i]; }
  for (i=0; i<weight.size(); i++) { weight[i] /= w; }

  return 0;
}


int _read_rule_csv(std::string &fn, std::vector< std::vector<float> > &rule) {
  int i;
  float val, _weight;
  FILE *fp;
  std::string line, tok;
  std::vector<std::string> toks;
  std::vector<float> v;

  float _tile_src, _tile_dst;

  fp = fopen(fn.c_str(), "r");
  if (!fp) { return -1; }

  while (!feof(fp)) {
    line.clear();
    _read_line(fp, line);

    if (line.size()==0) { continue; }
    if (line[0] == '#') { continue; }
    if (line[0] == ' ') { continue; }

    toks.clear();
    tok.clear();
    for (i=0; i<line.size(); i++) {
      if (line[i]==',') {
        toks.push_back(tok);
        tok.clear();
        continue;
      }
      tok += line[i];
    }
    toks.push_back(tok);

    if ((toks.size() < 3) || (toks.size() > 4)) { continue; }

    _tile_src = atof(toks[0].c_str());
    _tile_dst = atof(toks[1].c_str());
    _weight = 1;

    if ((toks.size() >= 4) &&
        (toks[3].size() != 0) &&
        (toks[3][0] != 'u')) {
      _weight = atof(toks[3].c_str());
    }

    // direction wild card
    //
    if ((toks[2].size()==0) ||
        (toks[2][0] == '*')) {
      v.clear();
      v.push_back(0.0);
      v.push_back(0.0);
      v.push_back(0.0);
      v.push_back(0.0);
      for (i=0; i<6; i++) {
        v[0] = _tile_src;
        v[1] = _tile_dst;
        v[2] = (float)i;
        v[3] = _weight ;
        rule.push_back(v);
      }
    }

    // explicit entry
    //
    else {
      v.clear();
      v.push_back(_tile_src);
      v.push_back(_tile_dst);
      v.push_back(atof(toks[2].c_str()));
      v.push_back(_weight);
      rule.push_back(v);
    }

  }

  fclose(fp);

  return 0;
}


// constraint file format:
//
// <x>,<y>,<z>,<tile_id>
//
int ConstraintCollapse::read_constraints (std::string &fn) {
  int i;
  float val, _weight;
  FILE *fp;
  std::string line, tok;
  std::vector<std::string> toks;
  std::vector<float> v;
  std::map<int64_t, int> histogram;
  std::map<int64_t, int>::iterator it;
  std::vector<std::vector<int32_t>> temp_admit;
  int64_t p;

  float _tile_src, _tile_dst;
  int x, y, z, tileid;

  std::vector< int32_t > v32;

  fp = fopen(fn.c_str(), "r");
  if (!fp) { return -1; }
  
  int tmax=0;

  while (!feof(fp)) {
    line.clear();
    _read_line(fp, line);

    toks.clear();
    tok.clear();
    for (i=0; i<line.size(); i++) {
      if (line[i]==',') {
        toks.push_back(tok);
        tok.clear();
        continue;
      }
      tok += line[i];
    }
    toks.push_back(tok);

    if ((toks.size() < 3) || (toks.size() > 4)) { continue; }

    x = atoi(toks[0].c_str());
    y = atoi(toks[1].c_str());
    z = atoi(toks[2].c_str());
    tileid = atoi(toks[3].c_str());

    if ( x >=0 && x < m_res.x && y >=0 && y < m_res.y && z >=0 && z < m_res.z) {
        // histogram admissable tiles
        p = getVertex (x,y,z);
        it = histogram.find ( p );
        if (it != histogram.end() ) {
            it->second++; 
            if (it->second > tmax) tmax = it->second;
        } else {
            histogram.insert ( std::pair<int64_t,int>( p, 1 ) );
        }

        v32.clear();
        v32.push_back(x);
        v32.push_back(y);
        v32.push_back(z);
        v32.push_back(tileid);

        temp_admit.push_back(v32); 
    }
  }
  fclose (fp);

  // filter out those which allow all tiles
  //
  Vector3DI vec;
  m_admissible.clear();
  dbgprintf ( " max admit: %d\n", tmax );
  for (int n=0; n < temp_admit.size(); n++) {
    v32 = temp_admit[n];
    p = getVertex( v32[0], v32[1], v32[2] );
    it = histogram.find( p );
    if ( it->second != tmax ) {
        vec = getVertexPos(p);
        dbgprintf ( " admissible: %d %d %d => %d\n", vec.x, vec.y, vec.z, v32[3] );    
        m_admissible.push_back ( v32 );
    }
  }
  

  return 0;
}

