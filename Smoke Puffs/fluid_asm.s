//
//  fluid_asm.s
//  Smoke Puffs
//
//  Created by Matt Stanton on 8/12/12.
//
//

#import "TargetConditionals.h"

#if !TARGET_IPHONE_SIMULATOR

//void JacobiARM(float* pressure,
//float* div,
//float* inv_count,
//int w,
//int h,
//float* new_pressure) {

// pressure / prev_pressure: r0
// div: r1
// inv_count: r2
// w / l:  r3
// h: r4
// new_pressure / new_pressure_start: r5
// pressure_start: r6
// next_pressure: r8

.globl _arm7_jacobi_iteration
_arm7_jacobi_iteration:

push     {r4-r7, lr}           // save LR, R7, R4-R6
add      r7, sp, #12           // adjust R7 to point to saved R7
push     {r8}
//vpush    {d8-d15}              // save VFP/Advanced SIMD registers D8
                               //  (aka S16-S31, Q4-Q7)
//sub      sp, sp, #36           // allocate space for local storage

ldr     r4, [sp, #+(4*6)]
ldr     r5, [sp, #+(4*7)]
sub     r3, r3, #2
muls    r3, r3, r4
asr     r3, r3, #1

lsls    r4, r4, #2

add     r6, r0, r4
vldr.f64    d3, [r6, #-4]  // initialize trailing y cells
add     r1, r1, r4
add     r2, r2, r4
add     r8, r6, r4
add     r5, r5, r4
add     r6, r6, #-4

// calculate first result
vld1.32     {d0}, [r0]!
vld1.32     {d1}, [r8]!
vldm.f32    r6, {d2-d3}
vadd.f32    d0, d0, d1
vld1.32     {d4}, [r1]!
vadd.f32    d2, d2, d3
vld1.32     {d5}, [r2]!
vadd.f32    d0, d0, d4
subs        r3, r3, #1
add         r6, r6, #8
vadd.f32    d6, d2, d0

.align 4

Ljacobi:

//vldr.f64    d0, [r0]
//vldr.f64    d1, [r8]
//vldr.f64    d2, [r6, #-4]
//vldr.f64    d3, [r6, #+4]
//vldr.f64    d4, [r1]       // load div
//vldr.f64    d5, [r2]       // load inv_count
//subs        r3, r3, #1     // subtract from length
//add         r6, r6, #8
//add         r0, r0, #8
//add         r8, r8, #8
//add         r2, r2, #8
//add         r1, r1, #8
//vadd.f32    d6, d0, d1
//vadd.f32    d7, d2, d3
//vadd.f32    d6, d6, d4
//vadd.f32    d7, d7, d6
//vmul.f32    d7, d7, d5
//vstr.f64    d7, [r5]
//add         r5, r5, #8

vld1.32     {d0}, [r0,:64]!
vld1.32     {d1}, [r8,:64]!
vmul.f32    d6, d6, d5
vldm.f32    r6, {d2-d3}
vadd.f32    d0, d0, d1
vld1.32     {d4}, [r1,:64]!
vadd.f32    d2, d2, d3
vld1.32     {d5}, [r2,:64]!
vadd.f32    d0, d0, d4
subs        r3, r3, #1
add         r6, r6, #8
vst1.32     {d6}, [r5,:64]
add         r5, r5, #8
vadd.f32    d6, d2, d0

bne         Ljacobi

// multiply and store last item
vmul.f32    d6, d6, d5
vst1.32     {d6}, [r5,:64]

//add      sp, sp, #36           // deallocate space for local storage
//vpop     {d8-d15}         // restore VFP/Advanced SIMD registers
pop      {r8}
pop      {r4-r7, pc}           // restore R4-R6, saved R7, return to saved LR

#endif