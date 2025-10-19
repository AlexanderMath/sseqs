/* ------------------------------------------------------------------ */
/*  handy half2 helpers                                               */
/* ------------------------------------------------------------------ */
__device__ inline __half2 h2(float r, float i)
{
    return __halves2half2(__float2half_rn(r), __float2half_rn(i));
}
__device__ inline float  real(const __half2 v) { 
    return __half2float(__low2half(v));  
    
    }
__device__ inline float  imag(const __half2 v) { 
    return __half2float(__high2half(v)); 
    
    }
__device__ inline __half2 add(const __half2 a, const __half2 b) { return __hadd2(a, b); }
__device__ inline __half2 sub(const __half2 a, const __half2 b) { return __hsub2(a, b); }

__device__ inline __half2 conjh(const __half2 a)
{//return a; // 0.5ms
    return h2(real(a), -imag(a));
}

__device__ inline __half2 cmul_slow(const __half2 a, const __half2 b)
{
    float ar = real(a), ai = imag(a);
    float br = real(b), bi = imag(b);
    return h2(ar * br - ai * bi,
              ar * bi + ai * br);
}

__device__ inline __half2 cmul(const __half2 a, const __half2 b)
{ // return a; here makes it 8.37ms -> 3.44ms wtf ;; hmm, appears this just cuts away with o3 compiler..
// return __hadd2(a,b); // 8.8->7.76
    // Split the components once (no casts)
    __half ar = __low2half(a);     // ar
    __half ai = __high2half(a);    // ai
    __half br = __low2half(b);     // br
    __half bi = __high2half(b);    // bi
    //return a;  // 4.25ms

    // ar*br – ai*bi  (real part)
    __half real = __hsub(__hmul(ar, br), __hmul(ai, bi));
    //return __halves2half2(real,real);  // 7.2ms 

    // ar*bi + ai*br  (imag part)
    __half imag = __hadd(__hmul(ar, bi), __hmul(ai, br));

    return __halves2half2(real, imag);   // pack back into one register
}


__device__ inline __half2 cmul_fast(const __half2 a, const __half2 b) {
   return __hcmadd(a, b, __float2half2_rn(0.f));
}


/*__device__ __forceinline__ __half2 cmul_fast(const __half2 a, const __half2 b){
    const float2 af = __half22float2(a), bf = __half22float2(b);
    return __floats2half2_rn(af.x * bf.x - af.y * bf.y,  af.x * bf.y + af.y * bf.x);
}*/