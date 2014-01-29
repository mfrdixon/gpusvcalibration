  // file: Error_Func.c
  #ifdef __cplusplus
  extern "C" {
  #endif

  #include <R.h>
  #include <Rinternals.h>
  #include <R_ext/Rdynload.h>

  // prototype for host function

  extern double error_func (double *, int);
  extern double copy_data ( double *, double *, double *, double *, char *, double, double, double, int);
  extern double dealloc_data ();
  extern double set_model ( char *);
  extern double set_block_size(int);

  // the function callable from R
  SEXP Error_Func(SEXP p0) {

    // object for returning the result
    SEXP out;

    double *_p0, _out;

    R_len_t n;
	
    _p0 = REAL(p0);
    n = length(p0);
    PROTECT(out = allocVector(REALSXP,1));
	
    _out = error_func (_p0, n);
    REAL(out)[0] = _out;

    UNPROTECT(1);

    return out;
  } // fast_sum
  // the function callable from R
  SEXP Copy_Data(SEXP K, SEXP T, SEXP W, SEXP OP, SEXP TYPES, SEXP S0, SEXP R0, SEXP Q0) {

    // object for returning the result
    SEXP out;

    double *_K,*_T,*_W,*_OP,_out;
    double _s0, _r0, _q0;
    char* _types;

    R_len_t n;
	
    _K = REAL(K);
    _T = REAL(T);
    _W = REAL(W);
    _OP= REAL(OP);
    _types = CHAR(STRING_ELT(TYPES,0));
    _s0 = REAL(S0)[0];
    _r0 = REAL(R0)[0];
    _q0 = REAL(Q0)[0];
    n = length(K);

    PROTECT(out = allocVector(REALSXP,1));
	
    _out = copy_data (_K,_T,_W,_OP,_types,_s0, _r0,_q0,n);
    REAL(out)[0] = _out;

    UNPROTECT(1);

    return out;
  } 

  SEXP Dealloc_Data() {
    SEXP out;
    double _out;
    PROTECT(out = allocVector(REALSXP,1));
    _out = dealloc_data ();
    REAL(out)[0] = _out;

    UNPROTECT(1);
    return out;
  }

  SEXP Set_Model(SEXP MODEL) {
    SEXP out;
    double _out;
    char * _model;
    _model = CHAR(STRING_ELT(MODEL,0));
    PROTECT(out = allocVector(REALSXP,1));
    _out = set_model (_model);
    REAL(out)[0] = _out;

    UNPROTECT(1);
    return out;
  }

  SEXP Set_Block_Size(SEXP BLKSZ) {
    SEXP out;
    double _out;
    int _blk_size;
    _blk_size = INTEGER(BLKSZ)[0];
    PROTECT(out = allocVector(REALSXP,1));
    _out = set_block_size (_blk_size);
    REAL(out)[0] = _out;

    UNPROTECT(1);
    return out;
  }

  // DLL scaffolding for R
  R_CallMethodDef callMethods[] = {
    {"Error_Func", (DL_FUNC)&Error_Func, 2},
    {"Copy_Data", (DL_FUNC)&Copy_Data, 9},
    {"Dealloc_Data", (DL_FUNC)&Dealloc_Data, 0},
    {"Set_Model", (DL_FUNC)&Set_Model, 1},
    {"Set_Block_Size", (DL_FUNC)&Set_Block_Size, 1},
    {NULL, NULL, 0}
  };

  void R_init_myLib (DllInfo *info) {
    R_registerRoutines (info, NULL, callMethods,
                        NULL, NULL);
  }

  #ifdef __cplusplus
  } // closing brace for extern "C"
  #endif
