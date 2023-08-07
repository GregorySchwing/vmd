/***************************************************************************
 *cr                                                                       
 *cr            (C) Copyright 1995-2019 The Board of Trustees of the           
 *cr                        University of Illinois                       
 *cr                         All Rights Reserved                        
 *cr                                                                   
 ***************************************************************************/

/***************************************************************************
 * RCS INFORMATION:
 *
 *	$RCSfile: util_simd_AVX2.C,v $
 *	$Author: johns $	$Locker:  $		$State: Exp $
 *	$Revision: 1.1 $	$Date: 2020/11/01 05:00:11 $
 *
 ***************************************************************************
 * DESCRIPTION:
 *
 * Hand-coded SIMD loops using compiler provided intrinsics, or inline
 * assembly code to generate highly optimized machine code for time-critical
 * loops that crop up commonly used features of VMD. 
 *
 ***************************************************************************/

#include "WKFThreads.h" // CPU capability flags

#if defined(__SSE2__) || (defined(_MSC_VER) && (_MSC_VER >= 1916))
#define VMDUSESSE 1     // enable SSE in combination with target macros
#endif

#if defined(__AVX__) || (defined(_MSC_VER) && (_MSC_VER >= 1916))
#define VMDUSEAVX 1     // enable AVX with runtime dispatch 
#endif

#if defined(__AVX2__) || (defined(_MSC_VER) && (_MSC_VER >= 1916))
#define VMDUSEAVX2 1     // enable AVX with runtime dispatch 
#endif

#if defined(VMDUSESSE) && defined(VMDUSEAVX) 

#if defined(VMDUSESSE)
#include <emmintrin.h>
#endif
#if defined(VMDUSEAVX)
#include <immintrin.h>
#endif

// #include <string.h>
// #include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(_MSC_VER)
#include <windows.h>
#include <conio.h>
#else
#include <unistd.h>
#endif // _MSC_VER


//
// Helper routine for use when coping with unaligned
// buffers returned by malloc() on many GNU systems:
//   http://gcc.gnu.org/bugzilla/show_bug.cgi?id=24261
//   http://www.sourceware.org/bugzilla/show_bug.cgi?id=206
//
// XXX until all compilers support uintptr_t, we have to do 
//     dangerous and ugly things with pointer casting here...
//
#if defined(_WIN64) /* sizeof(size_t) == sizeof(void*) */
#define myintptrtype size_t
#elif 1 /* sizeof(unsigned long) == sizeof(void*) */
#define myintptrtype unsigned long
#else /* C99 */
#define myintptrtype uintptr_t
#endif


//
// Aligment test routines for vector instructions
//
static int test_alignment_Nbyte_powertwo(const void *ptr, unsigned int alignsz) {
  unsigned int alignmask = alignsz - 1;
  return (((myintptrtype) ptr) == (((myintptrtype) ptr) & (~alignmask)));
}


//
// Small inlinable AVX helper routines to make code easier to read
//
static int hadd_m256i(__m256i sum8) {
  int tmp[sizeof(__m256i)/sizeof(int)];
  _mm256_store_si256((__m256i *)tmp, sum8);
  int sum = tmp[0]+tmp[1]+tmp[2]+tmp[3]+tmp[4]+tmp[5]+tmp[6]+tmp[7];
  return sum;
}


#if defined(VMDUSEAVX)
// Find the first selected atom, the last selected atom,
// and the total number of selected atoms.
int analyze_selection_aligned_avx2(int n, const int *on,
                                   int *firstsel, int *lastsel, int *selected) {
  int sel   = 0;   // set count to zero in case of early-exit
  int first = 0;   // if we early-exit, firstsel is 0
  int last  = -1;  // and lastsel is -1
  int i;

  // assign potential early-exit outcomes
  if (selected != NULL)
    *selected = sel;

  if (firstsel != NULL)
    *firstsel = first;

  if (lastsel != NULL)
    *lastsel = last;

  // find the first selected atom, if any
  if ((firstsel != NULL) || (selected != NULL)) {
    // find the first selected atom, if any
    // roll up to the first 32-byte-aligned array index
    for (i=0; ((i<n) && !test_alignment_Nbyte_powertwo(&on[i], 32)); i++) {
      if (on[i]) {
        first = i; // found first selected atom
        goto foundfirstsel;
      }
    }

    // AVX vectorized search loop
    for (; i<(n-7); i+=8) {
      // aligned load of 8 selection flags
      __m256i on8 = _mm256_load_si256((__m256i*) &on[i]);
      if (!_mm256_testz_si256(on8, on8))
        break; // found a block containing the first selected atom
    }

    for (; i<n; i++) {
      if (on[i]) {
        first = i; // found first selected atom
        goto foundfirstsel;
      }
    }

    // prevent x86 AVX-SSE transition performance loss due to CPU state 
    // transition penalties or false dependence on upper register state
    _mm256_zeroupper();
    return -1; // indicate that no selection was found
  }
foundfirstsel:

  // find the last selected atom, if any
  if ((lastsel != NULL) || (selected != NULL)) {
    // AVX vectorized search loop
    // Roll down to next 32-byte boundary
    for (i=n-1; i>=0; i--) {
      if (on[i]) {
        last = i; // found last selected atom
        goto foundlastsel;
      }

      // drop out of the alignment loop once we hit a 32-byte boundary
      if (test_alignment_Nbyte_powertwo(&on[i], 32))
        break;
    }

    for (i-=8; i>=0; i-=8) {
      // aligned load of 8 selection flags
      __m256i on8 = _mm256_load_si256((__m256i*) &on[i]);
      if (!_mm256_testz_si256(on8, on8))
        break; // found a block containing the last selected atom
    }

    int last8=i;
    for (i=last8+7; i>=last8; i--) {
      if (on[i]) {
        last = i; // found last selected atom
        goto foundlastsel;
      }
    }

    // prevent x86 AVX-SSE transition performance loss due to CPU state 
    // transition penalties or false dependence on upper register state
    _mm256_zeroupper();
    return -1; // indicate that no selection was found
  }
foundlastsel:

  // count the number of selected atoms (there are only 0s and 1s)
  // and determine the index of the last selected atom
  if (selected != NULL) {
    // AVX2 vectorized search loop
    // Roll up to next 32-byte boundary
    for (i=first; ((i<=last) && (!test_alignment_Nbyte_powertwo(&on[i], 32))); i++) {
      sel += on[i];
    }

    // Process groups of 8 flags at a time
    __m256i sum8 = _mm256_setzero_si256();
    for (; i<=(last-7); i+=8) {
      // aligned load of four selection flags
      __m256i on8 = _mm256_load_si256((__m256i*) &on[i]);

      // sum selected atoms vertically
      sum8 = _mm256_add_epi32(sum8, on8);
    }
    sel += hadd_m256i(sum8); // sum horizontally to finalize count

    // check the very end of the array (non-divisible by four)
    for (; i<=last; i++) {
      sel += on[i];
    }
  }

  if (selected != NULL)
    *selected = sel;

  if (firstsel != NULL)
    *firstsel = first;

  if (lastsel != NULL)
    *lastsel = last;

  // prevent x86 AVX-SSE transition performance loss due to CPU state 
  // transition penalties or false dependence on upper register state
  _mm256_zeroupper();
  return 0;
}
#endif

#endif // SSE+AVX

