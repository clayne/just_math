
#ifndef MAIN_BELIEF_PROPAGATION_HPP
#define MAIN_BELIEF_PROPAGATION_HPP

typedef struct _opt_t {
  float alpha;
  int alg_idx;
  std::string fn_name;
  std::string fn_rule;

  std::string tileset_fn,
              tilemap_fn,
              tileobj_fn;
  int32_t tileset_stride_x,
          tileset_stride_y;
  int32_t tileset_margin,
          tileset_spacing;
  int32_t tileset_width,
          tileset_height;

  int tiled_reverse_y;

} opt_t;

extern opt_t g_opt;

typedef struct constraint_op_type {
  char op;
  std::vector< int > dim_range;
  std::vector< int > tile_range;
} constraint_op_t;

int parse_range(std::vector<int> &range, std::string &s, std::vector<int> &dim);
int parse_bracket_range(std::vector<int> &range, std::string &s, std::vector<int> &dim);
int parse_constraint_dsl(std::vector< constraint_op_t > &op_list, std::string &s, std::vector< int > dim, std::vector< std::string > name);




// tests..
//
int test0();
int test1();
int test2();
int test3();
int test4_();
int test4();
int test5();
int test5_1();
int test_cull0();
int test_cull1();
int test_cull2();
int test_cull3();
int test_cull4();
int test6();
int test_realize0();
int test_realize1();
int test_realize2(int x, int y, int z);
int test_wfc0(int x, int y, int z);

int test_step0();
int test_step1();

int run_test(int test_num);

#endif

