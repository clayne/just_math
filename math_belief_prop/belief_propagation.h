//--------------------------------------------------------------------------------
// JUST MATH:
// Belief Propagation - on a 3D grid domain
//
// Demonstration of Sum-Product Belief Propagation on a 3D spatial domain.
// Computes:
//     mu_{i,j}[b] = SUM f_{i,j}[a,b] g_i[a] PROD mu_{k,i}[b]
// The message function 'mu' is stored sparsely for neighboring cells in 3D, with size 6*R^3*B,
// where R is the grid resolution, B is the number of discrete values, and 6 is number of neighbors.
//
// To render the result, the belief is estimated at each vertex (voxel), and
// raytraced as a density volume where value probabilities are mapped to color.
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

#ifndef DEF_BELIEF_PROPAGATION
#define DEF_BELIEF_PROPAGATION

#include <algorithm>
#include "mersenne.h"
#include "dataptr.h"

//extern "C" {
//#include "lib/svdlib.h"
//}

#include <Eigen/SVD>

/*
#ifdef USE_OPENGL
  #include <GL/glew.h>
#endif
#ifdef USE_CUDA
  #include "common_cuda.h"
#endif
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <vector>
#include <string>

#define BELIEF_PROPAGATION_VERSION "0.5.0"

#define OPT_PTRS
#define OPT_MUPTR
#define OPT_FH
#define OPT_MUBOUND

#define MU_NOCOPY 0
#define MU_COPY 1


#define VIZ_NONE        0
#define VIZ_MU          1
#define VIZ_DMU         2
#define VIZ_BELIEF      3
#define VIZ_CONSTRAINT  4
#define VIZ_TILECOUNT   5
#define VIZ_ENTROPY     6
#define VIZ_CHANGE      7
#define VIZ_RESPICK     8

#define ALG_CELL_ANY            32
#define ALG_CELL_MIN_ENTROPY    33

#define ALG_TILE_MAX_BELIEF     34

#define ALG_RUN_VANILLA         35
#define ALG_RUN_RESIDUAL        36

// memory locations: cpu + (z*mUseRY + y)*mUseRX + x

// static buffers (input)                                                                                               // Allocation (B=num_vals)
#define BUF_G           1     // tile weights,  weight, all values - beliefprop, G(a) vector                            // <B, 1, 1>
#define BUF_F           2     // tile rules,    rule-to-rule - beliefprop, F(a,b) vector - <src, dest, direction>       // <B, B, 6>

// dynamic buffers
#define BUF_MU          3     // message prob,  6*B,        all verts - beliefprop, mu{i,j}(a,b) vector (B=# values)    // <6, B, num_vert>
#define BUF_MU_NXT      4     // message prob', 6*B,        all verts - beliefprop, mu'{i,j}(a,b) vector (B=# values)   // <6, B, num_vert>
#define BUF_TILE_IDX    5     // tile indexes,  val list,   all verts - beliefprof                                      // <B, num_vert, 1>
#define BUF_TILE_IDX_N  6     // # of tile idx, 1x int,     all verts - beliefprof                                      // <num_vert, 1, 1>

// scratch buffers
#define BUF_H           7     // temporary,     val list,   single vert - beliefprop                                    // <B, 1, 1>
#define BUF_BELIEF      8     // temporary,     val list,   single vert - beliefprop                                    // <B, 1, 1>
#define BUF_VISITED     9     // temporary,     1x int,     all verts - beliefprop                                      // <num_vert, 1, 1>
#define BUF_NOTE        10    // which cells,   2x int,     all verts (of size cell count (i32))                        // <num_vert, 2, 1>
#define BUF_VIZ         11    // vizualization, 1x float,   all verts                                                   // <num_vert, 1, 1>
#define BUF_TILES       12    // max belief ,   1x int,     all verts                                                   // <num_vert, 1, 1>
#define BUF_C           13    // num constraint,1x int,     all verts                                                   // <num_vert, 1, 1>

// svd & residual bp
#define BUF_SVD_U       14    //                                                                                        // <B,  B*, 6>
#define BUF_SVD_Vt      15    //                                                                                        // <B*, B,  6>
#define BUF_SVD_VEC     16    //                                                                                        // <B*, 1,  1>


// auxiliary buffers for residual belief propagaion
//
// BUF_RESIDUE_HEAP         : heap of absolute differences of mu and mu_nxt (float)
// BUF_RESIDUE_HEAP_CELL_BP : back pointer of heap value location in CELL_HEAP (in64_t)
// BUF_RESIDUE_CELL_HEAP    : mapping of cell (and direction, value) to heap position (in64_t).
//                            That is, mapping of mu index to position in heap
//
// All sizes should be (Vol) * 6 * B.
// That is, {volume} x {#neighbors} x {#values} : (dim[0]*dim[1]*dim[2] * 6 * B).
//
// All these structures are for book keeping so that we can get the maximum difference of
// mu and mu_nxt in addition to allowing arbitrary updates on other cells.
//
// The basic residual bp update step will fetch the maximum difference, copy the mu value
// from the mu_nxt buffer to the mu buffer, then update neighboring cells by updating their
// mu_nxt values, updating the residue_heap along the way.
//
#define BUF_RESIDUE_HEAP          19    //                                                                               // <6*B*num_vert, 1, 1>
#define BUF_RESIDUE_HEAP_CELL_BP  20    //                                                                               // <6*B*num_vert, 1, 1>
#define BUF_RESIDUE_CELL_HEAP     21    //                                                                               // <6*B*num_vert, 1, 1>

class BeliefPropagation {
public:
  BeliefPropagation() {
    m_seed = 17;
    m_verbose = 0;
    m_eps_converge = (1.0/(1024.0));
    //m_eps_converge = (1.0/(1024.0*1024.0));
    //m_eps_zero = (1.0/(1024.0*1024.0));
    m_eps_zero = (1.0/(1024.0*1024.0*1024.0*1024.0));
    m_max_iteration = 1024;

    //m_step_cb = 10;
    m_step_cb = 1;
    m_state_info_d = -1;
    m_state_info_iter = 0;

    m_rate = 0.98;

    m_use_svd = 0;
    m_use_checkerboard = 0;

    m_index_heap_size = 0;

    m_stat_enabled = 1;
    m_stat_avg_iter = 0.0;

    // unused...
    m_stat_second_moment_iter = 0.0;

    m_stat_cur_iter = 0;
    m_stat_max_iter = 0;
    m_stat_num_culled = 0;
    m_stat_num_collapsed = 0;

    // unused...
    m_stat_num_chosen = 0;

    m_eps_converge_beg = m_eps_converge;
    m_eps_converge_end = m_eps_converge;

    m_viz_opt = VIZ_NONE;

    m_alg_cell_opt = ALG_CELL_MIN_ENTROPY;
    m_alg_tile_opt = ALG_TILE_MAX_BELIEF;
    m_alg_run_opt = ALG_RUN_VANILLA;

    m_run_iter = 0;
    m_step_iter = 0;

  };

  //------------------------ high level API

  int       start();

  int       RealizePre();
  int       RealizeRun();
  int       RealizeStep();
  int       RealizePost();
  int       Realize();

  int       CheckConstraints ( int64_t p );
  int       CheckConstraints ();

  void      SetVis (int viz_opt);

  //------------------------ belief propagation, mid-level API



  int   init( int, int, int,
              std::vector< std::string  >           tile_name_list,
              std::vector< float >                  tile_weight_list,
              std::vector< std::vector < float > >  rule_list );

  int   init_SVD(void);

  int   realize();

  int   filter_constraint(std::vector< std::vector< int32_t > > &constraint_list);

  void  gp_state_print();

  // legacy
  int   single_realize (int64_t it);
  int   single_realize_cb (int64_t it, void (*cb)(void *));
  //int   single_realize_lest_belief_cb (int64_t it, void (*cb)(void *));

  // core methods
  int   RAMA_single_realize_max_belief_cb(int64_t it, void (*cb)(void *));

  int   single_realize_max_belief_cb(int64_t it, void (*cb)(void *));
  int   single_realize_min_entropy_max_belief_cb(int64_t it, void (*cb)(void *));

  // min belief algorithm
  int   single_realize_min_belief_cb(int64_t it, void (*cb)(void *));

  // experimental
  int   single_realize_min_entropy_min_belief_cb(int64_t it, void (*cb)(void *));
  int   single_realize_residue_cb(int64_t it, void (*cb)(void *));

  int   _pick_tile(int64_t anch_cell, int64_t *max_cell, int32_t *max_tile, int32_t *max_tile_idx, float *max_belief);
  int   _pick_tile_max_belief(int64_t anch_cell, int64_t *max_cell, int32_t *max_tile, int32_t *max_tile_idx, float *max_belief);
  int   _pick_tile_min_belief(int64_t anch_cell, int64_t *min_cell, int32_t *min_tile, int32_t *min_tile_idx, float *min_belief);
  int   _pick_tile_pdf(int64_t anch_cell, int64_t *max_cell, int32_t *max_tile, int32_t *max_tile_idx, float *max_belief);

  void  init_dir_desc();
  float  step(int update_mu);
  float  step_residue(int32_t idir, int64_t cell, int32_t tile);


  //---------------------------------------- residual belief propagation
  //
  int64_t getMuIdx( int32_t idir, int64_t cell, int32_t tile );
  int64_t getMuPos( int64_t idx, int32_t *idir, int64_t *cell, int32_t *tile );

  void    indexHeap_init(void);
  void    indexHeap_swap(int64_t heap_idx_a, int64_t heap_idx_b);
  int32_t indexHeap_push(float val);

  void    indexHeap_update(int64_t heap_idx, float val);
  void    indexHeap_update_mu_idx(int64_t mu_idx, float val);
  void    indexHeap_update_mu_pos(int32_t idir, int64_t cell, int32_t tile, float val);

  int64_t indexHeap_peek(int64_t *mu_idx, float *val);
  int64_t indexHeap_peek_mu_pos(int32_t *idir, int64_t *cell, int32_t *tile_val, float *val);


  int32_t indexHeap_consistency(void);
  int32_t indexHeap_mu_consistency(void);

  void    indexHeap_debug_print(void);


  //----------------------- belief propagation (low-level)

  void  ConstructStaticBufs ();
  void  ConstructDynamicBufs ();
  void  ConstructTempBufs ();
  void  ConstructConstraintBufs();
  void  ConstructSVDBufs ();

  float BeliefProp();
  float BeliefProp_svd ();

  float BeliefProp_cell_residue(int64_t);
  float BeliefProp_cell_residue_svd(int64_t);

  void  UpdateMU ();

  float getVertexBelief ( uint64_t j );
  float _getVertexBelief ( uint64_t j );

  int   getMaxBeliefTile ( uint64_t j );

  void  cellUpdateBelief(int64_t anch_cell);
  int   chooseMaxBelief(int64_t *max_cell, int32_t *max_tile, int32_t *max_tile_idx, float *max_belief);
  int   chooseMinBelief(int64_t *min_cell, int32_t *min_tile, int32_t *min_tile_idx, float *min_belief);

  int   chooseMaxEntropy(int64_t *max_cell, int32_t *max_tile, int32_t *max_tile_idx, float *max_belief);
  int   chooseMinEntropyMaxBelief(int64_t *max_cell, int32_t *max_tile, int32_t *max_tile_idx, float *max_belief);
  int   chooseMinEntropyMinBelief(int64_t *min_cell, int32_t *min_tile, int32_t *min_tile_idx, float *min_belief);

  void  WriteBoundaryMU ();
  void  WriteBoundaryMUbuf(int buf_id);
  void  TransferBoundaryMU (int src_id, int dst_id);
  float MaxDiffMU();
  float MaxDiffMUCellTile(float *max_diff, int64_t *max_cell, int64_t *max_tile_idx, int64_t *max_dir_idx);

  void  RandomizeMU ();

  void  NormalizeMU ();
  void  NormalizeMU (int id);
  void  NormalizeMU_cell_residue (int buf_id, int64_t cell);

  void  filterKeep(uint64_t pos, std::vector<int32_t> &tile_id);
  void  filterDiscard(uint64_t pos, std::vector<int32_t> &tile_id);
  int32_t tileName2ID (std::string &tile_name);
  int32_t tileName2ID (char *);

  // used for visualization
  void  ComputeDiffMUField ();
  void  ComputeBeliefField ();


  // non "strict" bp functions but helpful still
  //
  int   CullBoundary();
  int   cellConstraintPropagate();
  void  cellFillAccessed(uint64_t vtx, int32_t note_idx);
  int   cellFillSingle(uint64_t vtx, int32_t note_idx);

  int   tileIdxCollapse(uint64_t pos, int32_t tile_idx);
  int   tileIdxRemove(uint64_t pos, int32_t tile_idx);

  // note_idx is the 'plane' of BUF_NOTE to unwind
  //
  void  unfillAccessed(int32_t note_idx);
  int   removeTileIdx(int64_t anch_cell, int32_t anch_tile_idx);
  int   sanityAccessed();

  //----------------------- visualization
  Vector4DF getVisSample ( int64_t v );



  //------------------------ memory management    
  
  void          AllocBuf (int id, char dt, uint64_t cntx=1, uint64_t cnty=1, uint64_t cntz=1 );     // new function
  void          ZeroBuf (int id);

  int64_t       getNeighbor(uint64_t j, int nbr);        // 3D spatial neighbor function
  int64_t       getNeighbor(uint64_t j, Vector3DI jp, int nbr);        // 3D spatial neighbor function
  Vector3DI     getVertexPos(int64_t j);
  int64_t       getVertex(int x, int y, int z);
  int           getTilesAtVertex ( int64_t vtx );
  int           getOppositeDir(int nbr)  { return m_dir_inv[nbr]; }
  
  //----------------------- new accessor functions

  inline void*  getPtr(int id, int x=0, int y=0, int z=0)     {return (void*) m_buf[id].getPtr (x, y, z);}     // caller does type casting

  inline int32_t getValI(int id, int x=0, int y=0, int z=0)            {return *(int32_t*) m_buf[id].getPtr (x, y, z);}
  inline int64_t getValL(int id, int x=0, int y=0, int z=0)            {return *(int64_t*) m_buf[id].getPtr (x, y, z);}
  inline float   getValF(int id, int x=0, int y=0, int z=0)            {return *(float*) m_buf[id].getPtr (x, y, z);}

  inline void   SetValI(int id, int32_t val, int x, int y=0, int z=0)     {*(int32_t*) m_buf[id].getPtr(x, y, z) = val;}
  inline void   SetValL(int id, int64_t val, int x, int y=0, int z=0)     {*(int64_t*) m_buf[id].getPtr(x, y, z) = val;}
  inline void   SetValF(int id, float val, int x, int y=0, int z=0)     {*(float*)   m_buf[id].getPtr(x, y, z) = val;}

  inline int    getNumNeighbors(int j)        {return 6;}
  inline int    getNumValues(int j)          {return m_num_values;}
  inline int    getNumVerts()            {return m_num_verts;}


  //----------------------- LEGACY accessor functions

  //  belief prop residue access functions (int64_t)
  //  index heap - ih
  //
  /* inline int64_t* getPtr_ih(int32_t id, int32_t n, int64_t j, int32_t a)                { return  (int64_t*) m_buf[id].getPtr ( uint64_t(n*m_num_verts + j)*m_num_values + a ); }
  inline int64_t  getVal_ih(int32_t id, int32_t n, int64_t j, int32_t a)                { return *(int64_t*) m_buf[id].getPtr ( uint64_t(n*m_num_verts + j)*m_num_values + a ); }
  inline void     SetVal_ih(int32_t id, int32_t n, int64_t j, int32_t a, int64_t val )  { *(int64_t*) m_buf[id].getPtr ( uint64_t(n*m_num_verts + j)*m_num_values + a ) = val; }

  inline int64_t* getPtr_ih(int32_t id, int64_t idx)                { return  (int64_t*) m_buf[id].getPtr ( idx ); }
  inline int64_t  getVal_ih(int32_t id, int64_t idx)                { return *(int64_t*) m_buf[id].getPtr ( idx ); }
  inline void     SetVal_ih(int32_t id, int64_t idx, int64_t val )  { *(int64_t*) m_buf[id].getPtr ( idx ) = val; }

  inline float* getPtr_ihf(int32_t id, int64_t idx)             { return  (float*) m_buf[id].getPtr ( idx ); }
  inline float  getVal_ihf(int32_t id, int64_t idx)             { return *(float*) m_buf[id].getPtr ( idx ); }
  inline void   SetVal_ihf(int32_t id, int64_t idx, float val ) { *(float*) m_buf[id].getPtr ( idx ) = val; } */

/* inline int32_t getValNote(int id, int i, int a)                { return *(int32_t *) m_buf[id].getPtr ( uint64_t(i*m_num_verts+ a) ); }
   inline void    SetValNote(int id, int i, int a, int32_t val)   { *(int32_t*) m_buf[id].getPtr ( (i*m_num_verts + a) ) = val; }  */


  //-------------------------- wave function collapse
  int   wfc();
  int   wfc_start();
  int   wfc_step(int64_t it);


  //-------------------------- debugging functions

  // helper arrays and functions for ease of testing and simple use
  //
  void  debugPrint();
  void  debugPrintC();
  void  debugPrintS();
  void  debugPrintMU();

  // run time statistics and other information
  //
  void    UpdateRunTimeStat(int64_t num_step);
  int32_t m_stat_enabled;
  double  m_stat_avg_iter,
          m_stat_second_moment_iter;
  int64_t m_stat_cur_iter,
          m_stat_max_iter,
          m_stat_num_culled,
          m_stat_num_collapsed,
          m_stat_num_chosen;


  //------------------------- member variables

  // primary data stored in buffers
  DataPtr   m_buf[128];

  // problem size
  int64_t   m_num_verts;    // Xi = 0..X (graph domain)
  int64_t   m_num_values;   //  B = 0..Bm-1 (value domain)
  Vector3DI m_bpres;        // 3D spatial belief prop res

  Vector3DI m_res;          // volume res

  std::vector< std::string > m_tile_name;
  std::vector< std::string > m_dir_desc;
  int m_dir_inv[6];

  // algorithm state
  int32_t     m_run_opt;

  int32_t     m_viz_opt;
  int32_t     m_alg_cell_opt;
  int32_t     m_alg_tile_opt;
  int32_t     m_alg_run_opt;

  int64_t     m_run_iter;


  // SVD number of non singular values in each direction
  //
  int       m_use_svd;
  int64_t   m_svd_nsv[6];
  int       m_use_checkerboard;

  bool      m_run_cuda=0;
  int       m_seed;
  Mersenne  m_rand;


  uint64_t m_note_n[2];

  int64_t m_grid_note_idx;

  float m_rate;



  int m_verbose;

  float m_eps_converge;

  float m_eps_converge_beg,
        m_eps_converge_end;

  float m_eps_zero;


  int64_t   m_step_cb;
  float     m_state_info_d;
  int64_t   m_state_info_iter;

  int64_t   m_step_iter;
  int64_t   m_max_iteration;

  int64_t   m_index_heap_size;



};


#endif
