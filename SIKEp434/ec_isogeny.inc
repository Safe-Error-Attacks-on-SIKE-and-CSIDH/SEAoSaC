/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: elliptic curve and isogeny functions
*********************************************************************************************/

#include "P434_internal.h"
#include "hal.h"

void xDBL(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t C24)
{ // Doubling of a Montgomery point in projective coordinates (X:Z).
  // Input: projective Montgomery x-coordinates P = (X1:Z1), where x1=X1/Z1 and Montgomery curve constants A+2C and 4C.
  // Output: projective Montgomery x-coordinates Q = 2*P = (X2:Z2).
    f2elm_t t0, t1;

    fp2sub(P->X, P->Z, t0);                         // t0 = X1-Z1
    fp2add(P->X, P->Z, t1);                         // t1 = X1+Z1
    fp2sqr_mont(t0, t0);                            // t0 = (X1-Z1)^2
    fp2sqr_mont(t1, t1);                            // t1 = (X1+Z1)^2
    fp2mul_mont(C24, t0, Q->Z);                     // Z2 = C24*(X1-Z1)^2
    fp2mul_mont(t1, Q->Z, Q->X);                    // X2 = C24*(X1-Z1)^2*(X1+Z1)^2
    fp2sub(t1, t0, t1);                             // t1 = (X1+Z1)^2-(X1-Z1)^2
    fp2mul_mont(A24plus, t1, t0);                   // t0 = A24plus*[(X1+Z1)^2-(X1-Z1)^2]
    fp2add(Q->Z, t0, Q->Z);                         // Z2 = A24plus*[(X1+Z1)^2-(X1-Z1)^2] + C24*(X1-Z1)^2
    fp2mul_mont(Q->Z, t1, Q->Z);                    // Z2 = [A24plus*[(X1+Z1)^2-(X1-Z1)^2] + C24*(X1-Z1)^2]*[(X1+Z1)^2-(X1-Z1)^2]
}


void xDBLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24plus, const f2elm_t C24, const int e)
{ // Computes [2^e](X:Z) on Montgomery curve with projective constant via e repeated doublings.
  // Input: projective Montgomery x-coordinates P = (XP:ZP), such that xP=XP/ZP and Montgomery curve constants A+2C and 4C.
  // Output: projective Montgomery x-coordinates Q <- (2^e)*P.
    int i;

    copy_words((digit_t*)P, (digit_t*)Q, 2*2*NWORDS_FIELD);

    for (i = 0; i < e; i++) {
        xDBL(Q, Q, A24plus, C24);
    }
}

#if (OALICE_BITS % 2 == 1)

void get_2_isog(const point_proj_t P, f2elm_t A, f2elm_t C)
{ // Computes the corresponding 2-isogeny of a projective Montgomery point (X2:Z2) of order 2.
  // Input:  projective point of order two P = (X2:Z2).
  // Output: the 2-isogenous Montgomery curve with projective coefficients A/C.

    fp2sqr_mont(P->X, A);                           // A = X2^2
    fp2sqr_mont(P->Z, C);                           // C = Z2^2
    fp2sub(C, A, A);                                // A = Z2^2 - X2^2
}


void eval_2_isog(point_proj_t P, point_proj_t Q)
{ // Evaluates the isogeny at the point (X:Z) in the domain of the isogeny, given a 2-isogeny phi.
  // Inputs: the projective point P = (X:Z) and the 2-isogeny kernel projetive point Q = (X2:Z2).
  // Output: the projective point P = phi(P) = (X:Z) in the codomain.
    f2elm_t t0, t1, t2, t3;

    fp2add(Q->X, Q->Z, t0);                         // t0 = X2+Z2
    fp2sub(Q->X, Q->Z, t1);                         // t1 = X2-Z2
    fp2add(P->X, P->Z, t2);                         // t2 = X+Z
    fp2sub(P->X, P->Z, t3);                         // t3 = X-Z
    fp2mul_mont(t0, t3, t0);                        // t0 = (X2+Z2)*(X-Z)
    fp2mul_mont(t1, t2, t1);                        // t1 = (X2-Z2)*(X+Z)
    fp2add(t0, t1, t2);                             // t2 = (X2+Z2)*(X-Z) + (X2-Z2)*(X+Z)
    fp2sub(t0, t1, t3);                             // t3 = (X2+Z2)*(X-Z) - (X2-Z2)*(X+Z)
    fp2mul_mont(P->X, t2, P->X);                    // Xfinal
    fp2mul_mont(P->Z, t3, P->Z);                    // Zfinal
}

#endif

void get_4_isog(const point_proj_t P, f2elm_t A24plus, f2elm_t C24, f2elm_t* coeff)
{ // Computes the corresponding 4-isogeny of a projective Montgomery point (X4:Z4) of order 4.
  // Input:  projective point of order four P = (X4:Z4).
  // Output: the 4-isogenous Montgomery curve with projective coefficients A+2C/4C and the 3 coefficients
  //         that are used to evaluate the isogeny at a point in eval_4_isog().

    fp2sub(P->X, P->Z, coeff[1]);                   // coeff[1] = X4-Z4
    fp2add(P->X, P->Z, coeff[2]);                   // coeff[2] = X4+Z4
    fp2sqr_mont(P->Z, coeff[0]);                    // coeff[0] = Z4^2
    fp2add(coeff[0], coeff[0], coeff[0]);           // coeff[0] = 2*Z4^2
    fp2sqr_mont(coeff[0], C24);                     // C24 = 4*Z4^4
    fp2add(coeff[0], coeff[0], coeff[0]);           // coeff[0] = 4*Z4^2
    fp2sqr_mont(P->X, A24plus);                     // A24plus = X4^2
    fp2add(A24plus, A24plus, A24plus);              // A24plus = 2*X4^2
    fp2sqr_mont(A24plus, A24plus);                  // A24plus = 4*X4^4
}


void eval_4_isog(point_proj_t P, f2elm_t* coeff)
{ // Evaluates the isogeny at the point (X:Z) in the domain of the isogeny, given a 4-isogeny phi defined
  // by the 3 coefficients in coeff (computed in the function get_4_isog()).
  // Inputs: the coefficients defining the isogeny, and the projective point P = (X:Z).
  // Output: the projective point P = phi(P) = (X:Z) in the codomain.
    f2elm_t t0, t1;

    fp2add(P->X, P->Z, t0);                         // t0 = X+Z
    fp2sub(P->X, P->Z, t1);                         // t1 = X-Z
    fp2mul_mont(t0, coeff[1], P->X);                // X = (X+Z)*coeff[1]
    fp2mul_mont(t1, coeff[2], P->Z);                // Z = (X-Z)*coeff[2]
    fp2mul_mont(t0, t1, t0);                        // t0 = (X+Z)*(X-Z)
    fp2mul_mont(t0, coeff[0], t0);                  // t0 = coeff[0]*(X+Z)*(X-Z)
    fp2add(P->X, P->Z, t1);                         // t1 = (X-Z)*coeff[2] + (X+Z)*coeff[1]
    fp2sub(P->X, P->Z, P->Z);                       // Z = (X-Z)*coeff[2] - (X+Z)*coeff[1]
    fp2sqr_mont(t1, t1);                            // t1 = [(X-Z)*coeff[2] + (X+Z)*coeff[1]]^2
    fp2sqr_mont(P->Z, P->Z);                        // Z = [(X-Z)*coeff[2] - (X+Z)*coeff[1]]^2
    fp2add(t1, t0, P->X);                           // X = coeff[0]*(X+Z)*(X-Z) + [(X-Z)*coeff[2] + (X+Z)*coeff[1]]^2
    fp2sub(P->Z, t0, t0);                           // t0 = [(X-Z)*coeff[2] - (X+Z)*coeff[1]]^2 - coeff[0]*(X+Z)*(X-Z)
    fp2mul_mont(P->X, t1, P->X);                    // Xfinal
    fp2mul_mont(P->Z, t0, P->Z);                    // Zfinal
}


void xTPL(const point_proj_t P, point_proj_t Q, const f2elm_t A24minus, const f2elm_t A24plus)
{ // Tripling of a Montgomery point in projective coordinates (X:Z).
  // Input: projective Montgomery x-coordinates P = (X:Z), where x=X/Z and Montgomery curve constants A24plus = A+2C and A24minus = A-2C.
  // Output: projective Montgomery x-coordinates Q = 3*P = (X3:Z3).
    f2elm_t t0, t1, t2, t3, t4, t5, t6;

    fp2sub(P->X, P->Z, t0);                         // t0 = X-Z
    fp2sqr_mont(t0, t2);                            // t2 = (X-Z)^2
    fp2add(P->X, P->Z, t1);                         // t1 = X+Z
    fp2sqr_mont(t1, t3);                            // t3 = (X+Z)^2
    fp2add(t0, t1, t4);                             // t4 = 2*X
    fp2sub(t1, t0, t0);                             // t0 = 2*Z
    fp2sqr_mont(t4, t1);                            // t1 = 4*X^2
    fp2sub(t1, t3, t1);                             // t1 = 4*X^2 - (X+Z)^2
    fp2sub(t1, t2, t1);                             // t1 = 4*X^2 - (X+Z)^2 - (X-Z)^2
    fp2mul_mont(t3, A24plus, t5);                   // t5 = A24plus*(X+Z)^2
    fp2mul_mont(t3, t5, t3);                        // t3 = A24plus*(X+Z)^3
    fp2mul_mont(A24minus, t2, t6);                  // t6 = A24minus*(X-Z)^2
    fp2mul_mont(t2, t6, t2);                        // t2 = A24minus*(X-Z)^3
    fp2sub(t2, t3, t3);                             // t3 = A24minus*(X-Z)^3 - coeff*(X+Z)^3
    fp2sub(t5, t6, t2);                             // t2 = A24plus*(X+Z)^2 - A24minus*(X-Z)^2
    fp2mul_mont(t1, t2, t1);                        // t1 = [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2]
    fp2add(t3, t1, t2);                             // t2 = [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2] + A24minus*(X-Z)^3 - coeff*(X+Z)^3
    fp2sqr_mont(t2, t2);                            // t2 = t2^2
    fp2mul_mont(t4, t2, Q->X);                      // X3 = 2*X*t2
    fp2sub(t3, t1, t1);                             // t1 = A24minus*(X-Z)^3 - A24plus*(X+Z)^3 - [4*X^2 - (X+Z)^2 - (X-Z)^2]*[A24plus*(X+Z)^2 - A24minus*(X-Z)^2]
    fp2sqr_mont(t1, t1);                            // t1 = t1^2
    fp2mul_mont(t0, t1, Q->Z);                      // Z3 = 2*Z*t1
}


void xTPLe(const point_proj_t P, point_proj_t Q, const f2elm_t A24minus, const f2elm_t A24plus, const int e)
{ // Computes [3^e](X:Z) on Montgomery curve with projective constant via e repeated triplings.
  // Input: projective Montgomery x-coordinates P = (XP:ZP), such that xP=XP/ZP and Montgomery curve constants A24plus = A+2C and A24minus = A-2C.
  // Output: projective Montgomery x-coordinates Q <- (3^e)*P.
    int i;

    copy_words((digit_t*)P, (digit_t*)Q, 2*2*NWORDS_FIELD);

    for (i = 0; i < e; i++) {
        xTPL(Q, Q, A24minus, A24plus);
    }
}


void get_3_isog(const point_proj_t P, f2elm_t A24minus, f2elm_t A24plus, f2elm_t* coeff)
{ // Computes the corresponding 3-isogeny of a projective Montgomery point (X3:Z3) of order 3.
  // Input:  projective point of order three P = (X3:Z3).
  // Output: the 3-isogenous Montgomery curve with projective coefficient A/C.
    f2elm_t t0, t1, t2, t3, t4;

    fp2sub(P->X, P->Z, coeff[0]);                   // coeff0 = X-Z
    fp2sqr_mont(coeff[0], t0);                      // t0 = (X-Z)^2
    fp2add(P->X, P->Z, coeff[1]);                   // coeff1 = X+Z
    fp2sqr_mont(coeff[1], t1);                      // t1 = (X+Z)^2
    fp2add(t0, t1, t2);                             // t2 = (X+Z)^2 + (X-Z)^2
    fp2add(coeff[0], coeff[1], t3);                 // t3 = 2*X
    fp2sqr_mont(t3, t3);                            // t3 = 4*X^2
    fp2sub(t3, t2, t3);                             // t3 = 4*X^2 - (X+Z)^2 - (X-Z)^2
    fp2add(t1, t3, t2);                             // t2 = 4*X^2 - (X-Z)^2
    fp2add(t3, t0, t3);                             // t3 = 4*X^2 - (X+Z)^2
    fp2add(t0, t3, t4);                             // t4 = 4*X^2 - (X+Z)^2 + (X-Z)^2
    fp2add(t4, t4, t4);                             // t4 = 2(4*X^2 - (X+Z)^2 + (X-Z)^2)
    fp2add(t1, t4, t4);                             // t4 = 8*X^2 - (X+Z)^2 + 2*(X-Z)^2
    fp2mul_mont(t2, t4, A24minus);                  // A24minus = [4*X^2 - (X-Z)^2]*[8*X^2 - (X+Z)^2 + 2*(X-Z)^2]
    fp2add(t1, t2, t4);                             // t4 = 4*X^2 + (X+Z)^2 - (X-Z)^2
    fp2add(t4, t4, t4);                             // t4 = 2(4*X^2 + (X+Z)^2 - (X-Z)^2)
    fp2add(t0, t4, t4);                             // t4 = 8*X^2 + 2*(X+Z)^2 - (X-Z)^2
    fp2mul_mont(t3, t4, A24plus);                   // A24plus = [4*X^2 - (X+Z)^2]*[8*X^2 + 2*(X+Z)^2 - (X-Z)^2]
}


void eval_3_isog(point_proj_t Q, const f2elm_t* coeff)
{ // Computes the 3-isogeny R=phi(X:Z), given projective point (X3:Z3) of order 3 on a Montgomery curve and
  // a point P with 2 coefficients in coeff (computed in the function get_3_isog()).
  // Inputs: projective points P = (X3:Z3) and Q = (X:Z).
  // Output: the projective point Q <- phi(Q) = (X3:Z3).
    f2elm_t t0, t1, t2;

    fp2add(Q->X, Q->Z, t0);                       // t0 = X+Z
    fp2sub(Q->X, Q->Z, t1);                       // t1 = X-Z
    fp2mul_mont(t0, coeff[0], t0);                // t0 = coeff0*(X+Z)
    fp2mul_mont(t1, coeff[1], t1);                // t1 = coeff1*(X-Z)
    fp2add(t0, t1, t2);                           // t2 = coeff0*(X+Z) + coeff1*(X-Z)
    fp2sub(t1, t0, t0);                           // t0 = coeff1*(X-Z) - coeff0*(X+Z)
    fp2sqr_mont(t2, t2);                          // t2 = [coeff0*(X+Z) + coeff1*(X-Z)]^2
    fp2sqr_mont(t0, t0);                          // t0 = [coeff1*(X-Z) - coeff0*(X+Z)]^2
    fp2mul_mont(Q->X, t2, Q->X);                  // X3final = X*[coeff0*(X+Z) + coeff1*(X-Z)]^2
    fp2mul_mont(Q->Z, t0, Q->Z);                  // Z3final = Z*[coeff1*(X-Z) - coeff0*(X+Z)]^2
}


void inv_3_way(f2elm_t z1, f2elm_t z2, f2elm_t z3)
{ // 3-way simultaneous inversion
  // Input:  z1,z2,z3
  // Output: 1/z1,1/z2,1/z3 (override inputs).
    f2elm_t t0, t1, t2, t3;

    fp2mul_mont(z1, z2, t0);                      // t0 = z1*z2
    fp2mul_mont(z3, t0, t1);                      // t1 = z1*z2*z3
    fp2inv_mont(t1);                              // t1 = 1/(z1*z2*z3)
    fp2mul_mont(z3, t1, t2);                      // t2 = 1/(z1*z2)
    fp2mul_mont(t2, z2, t3);                      // t3 = 1/z1
    fp2mul_mont(t2, z1, z2);                      // z2 = 1/z2
    fp2mul_mont(t0, t1, z3);                      // z3 = 1/z3
    fp2copy(t3, z1);                              // z1 = 1/z1
}


void get_A(const f2elm_t xP, const f2elm_t xQ, const f2elm_t xR, f2elm_t A)
{ // Given the x-coordinates of P, Q, and R, returns the value A corresponding to the Montgomery curve E_A: y^2=x^3+A*x^2+x such that R=Q-P on E_A.
  // Input:  the x-coordinates xP, xQ, and xR of the points P, Q and R.
  // Output: the coefficient A corresponding to the curve E_A: y^2=x^3+A*x^2+x.
    f2elm_t t0, t1, one = {0};

    fpcopy((digit_t*)&Montgomery_one, one[0]);
    fp2add(xP, xQ, t1);                           // t1 = xP+xQ
    fp2mul_mont(xP, xQ, t0);                      // t0 = xP*xQ
    fp2mul_mont(xR, t1, A);                       // A = xR*t1
    fp2add(t0, A, A);                             // A = A+t0
    fp2mul_mont(t0, xR, t0);                      // t0 = t0*xR
    fp2sub(A, one, A);                            // A = A-1
    fp2add(t0, t0, t0);                           // t0 = t0+t0
    fp2add(t1, xR, t1);                           // t1 = t1+xR
    fp2add(t0, t0, t0);                           // t0 = t0+t0
    fp2sqr_mont(A, A);                            // A = A^2
    fp2inv_mont(t0);                              // t0 = 1/t0
    fp2mul_mont(A, t0, A);                        // A = A*t0
    fp2sub(A, t1, A);                             // Afinal = A-t1
}


void j_inv(const f2elm_t A, const f2elm_t C, f2elm_t jinv)
{ // Computes the j-invariant of a Montgomery curve with projective constant.
  // Input: A,C in GF(p^2).
  // Output: j=256*(A^2-3*C^2)^3/(C^4*(A^2-4*C^2)), which is the j-invariant of the Montgomery curve B*y^2=x^3+(A/C)*x^2+x or (equivalently) j-invariant of B'*y^2=C*x^3+A*x^2+C*x.
    f2elm_t t0, t1;

    fp2sqr_mont(A, jinv);                           // jinv = A^2
    fp2sqr_mont(C, t1);                             // t1 = C^2
    fp2add(t1, t1, t0);                             // t0 = t1+t1
    fp2sub(jinv, t0, t0);                           // t0 = jinv-t0
    fp2sub(t0, t1, t0);                             // t0 = t0-t1
    fp2sub(t0, t1, jinv);                           // jinv = t0-t1
    fp2sqr_mont(t1, t1);                            // t1 = t1^2
    fp2mul_mont(jinv, t1, jinv);                    // jinv = jinv*t1
    fp2add(t0, t0, t0);                             // t0 = t0+t0
    fp2add(t0, t0, t0);                             // t0 = t0+t0
    fp2sqr_mont(t0, t1);                            // t1 = t0^2
    fp2mul_mont(t0, t1, t0);                        // t0 = t0*t1
    fp2add(t0, t0, t0);                             // t0 = t0+t0
    fp2add(t0, t0, t0);                             // t0 = t0+t0
    fp2inv_mont(jinv);                              // jinv = 1/jinv
    fp2mul_mont(jinv, t0, jinv);                    // jinv = t0*jinv
}


void xDBLADD(point_proj_t P, point_proj_t Q, const f2elm_t xPQ, const f2elm_t A24)
{ // Simultaneous doubling and differential addition.
  // Input: projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, affine difference xPQ=x(P-Q) and Montgomery curve constant A24=(A+2)/4.
  // Output: projective Montgomery points P <- 2*P = (X2P:Z2P) such that x(2P)=X2P/Z2P, and Q <- P+Q = (XQP:ZQP) such that = x(Q+P)=XQP/ZQP.
    f2elm_t t0, t1, t2;

    fp2add(P->X, P->Z, t0);                         // t0 = XP+ZP
    fp2sub(P->X, P->Z, t1);                         // t1 = XP-ZP
    fp2sqr_mont(t0, P->X);                          // XP = (XP+ZP)^2
    fp2sub(Q->X, Q->Z, t2);                         // t2 = XQ-ZQ
    fp2correction(t2);
    fp2add(Q->X, Q->Z, Q->X);                       // XQ = XQ+ZQ
    fp2mul_mont(t0, t2, t0);                        // t0 = (XP+ZP)*(XQ-ZQ)
    fp2sqr_mont(t1, P->Z);                          // ZP = (XP-ZP)^2
    fp2mul_mont(t1, Q->X, t1);                      // t1 = (XP-ZP)*(XQ+ZQ)
    fp2sub(P->X, P->Z, t2);                         // t2 = (XP+ZP)^2-(XP-ZP)^2
    fp2mul_mont(P->X, P->Z, P->X);                  // XP = (XP+ZP)^2*(XP-ZP)^2
    fp2mul_mont(t2, A24, Q->X);                     // XQ = A24*[(XP+ZP)^2-(XP-ZP)^2]
    fp2sub(t0, t1, Q->Z);                           // ZQ = (XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)
    fp2add(Q->X, P->Z, P->Z);                       // ZP = A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2
    fp2add(t0, t1, Q->X);                           // XQ = (XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)
    fp2mul_mont(P->Z, t2, P->Z);                    // ZP = [A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2]*[(XP+ZP)^2-(XP-ZP)^2]
    fp2sqr_mont(Q->Z, Q->Z);                        // ZQ = [(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
    fp2sqr_mont(Q->X, Q->X);                        // XQ = [(XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)]^2
    fp2mul_mont(Q->Z, xPQ, Q->Z);                   // ZQ = xPQ*[(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
}

void xDBLADD_SIM(point_proj_t P, point_proj_t Q, const f2elm_t xPQ, const f2elm_t A24, int swap)
{ // Simultaneous doubling and differential addition.
  // Input: projective Montgomery points P=(XP:ZP) and Q=(XQ:ZQ) such that xP=XP/ZP and xQ=XQ/ZQ, affine difference xPQ=x(P-Q) and Montgomery curve constant A24=(A+2)/4.
  // Output: projective Montgomery points P <- 2*P = (X2P:Z2P) such that x(2P)=X2P/Z2P, and Q <- P+Q = (XQP:ZQP) such that = x(Q+P)=XQP/ZQP.
    f2elm_t t0, t1, t2;

    fp2add(P->X, P->Z, t0);                         // t0 = XP+ZP
    fp2sub(P->X, P->Z, t1);                         // t1 = XP-ZP
    fp2sqr_mont(t0, P->X);                          // XP = (XP+ZP)^2
    fp2sub(Q->X, Q->Z, t2);                         // t2 = XQ-ZQ
    fp2correction(t2);
    fp2add(Q->X, Q->Z, Q->X);                       // XQ = XQ+ZQ
    fp2mul_mont(t0, t2, t0);                        // t0 = (XP+ZP)*(XQ-ZQ)
    fp2sqr_mont(t1, P->Z);                          // ZP = (XP-ZP)^2
    fp2mul_mont(t1, Q->X, t1);                      // t1 = (XP-ZP)*(XQ+ZQ)
    fp2sub(P->X, P->Z, t2);                         // t2 = (XP+ZP)^2-(XP-ZP)^2
    fp2mul_mont(P->X, P->Z, P->X);                  // XP = (XP+ZP)^2*(XP-ZP)^2
    fp2mul_mont(t2, A24, Q->X);                     // XQ = A24*[(XP+ZP)^2-(XP-ZP)^2]
    fp2sub(t0, t1, Q->Z);                           // ZQ = (XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)
    fp2add(Q->X, P->Z, P->Z);                       // ZP = A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2
    fp2add(t0, t1, Q->X);                           // XQ = (XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)
    fp2mul_mont(P->Z, t2, P->Z);                    // ZP = [A24*[(XP+ZP)^2-(XP-ZP)^2]+(XP-ZP)^2]*[(XP+ZP)^2-(XP-ZP)^2]
    fp2sqr_mont(Q->Z, Q->Z);                        // ZQ = [(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2
    fp2sqr_mont(Q->X, Q->X);                        // XQ = [(XP+ZP)*(XQ-ZQ)+(XP-ZP)*(XQ+ZQ)]^2
    fp2mul_mont(Q->Z, xPQ, Q->Z);                   // ZQ = xPQ*[(XP+ZP)*(XQ-ZQ)-(XP-ZP)*(XQ+ZQ)]^2



}


static void swap_points(point_proj_t P, point_proj_t Q, const digit_t option)
{ // Swap points.
  // If option = 0 then P <- P and Q <- Q, else if option = 0xFF...FF then P <- Q and Q <- P
    digit_t temp;
    unsigned int i;

    for (i = 0; i < NWORDS_FIELD; i++) {
        temp = option & (P->X[0][i] ^ Q->X[0][i]);
        P->X[0][i] = temp ^ P->X[0][i];
        Q->X[0][i] = temp ^ Q->X[0][i];
        temp = option & (P->Z[0][i] ^ Q->Z[0][i]);
        P->Z[0][i] = temp ^ P->Z[0][i];
        Q->Z[0][i] = temp ^ Q->Z[0][i];
        temp = option & (P->X[1][i] ^ Q->X[1][i]);
        P->X[1][i] = temp ^ P->X[1][i];
        Q->X[1][i] = temp ^ Q->X[1][i];
        temp = option & (P->Z[1][i] ^ Q->Z[1][i]);
        P->Z[1][i] = temp ^ P->Z[1][i];
        Q->Z[1][i] = temp ^ Q->Z[1][i];
    }
}

static void LADDER3PT(const f2elm_t xP, const f2elm_t xQ, const f2elm_t xPQ, const digit_t* m, const unsigned int AliceOrBob, point_proj_t R, const f2elm_t A)
{
    point_proj_t R0 = {0}, R2 = {0};
    f2elm_t A24 = {0};
    digit_t mask;
    int i, nbits, bit, swap, prevbit = 0;

    if (AliceOrBob == ALICE) {
        nbits = OALICE_BITS;
    } else {
        nbits = OBOB_BITS - 1;
    }

    // Initializing constant
    fpcopy((digit_t*)&Montgomery_one, A24[0]);
    fp2add(A24, A24, A24);
    fp2add(A, A24, A24);
    fp2div2(A24, A24);
    fp2div2(A24, A24);  // A24 = (A+2)/4

    // Initializing points
    fp2copy(xQ, R0->X);
    fpcopy((digit_t*)&Montgomery_one, (digit_t*)R0->Z);
    fp2copy(xPQ, R2->X);
    fpcopy((digit_t*)&Montgomery_one, (digit_t*)R2->Z);
    fp2copy(xP, R->X);
    fpcopy((digit_t*)&Montgomery_one, (digit_t*)R->Z);
    fpzero((digit_t*)(R->Z)[1]);

    // Main loop
    for (i = 0; i < nbits; i++) {

        bit = (m[i >> LOG2RADIX] >> (i & (RADIX-1))) & 1;

        swap = bit ^ prevbit;
        prevbit = bit;
        mask = 0 - (digit_t)swap;

        swap_points(R, R2, mask);
        xDBLADD(R0, R2, R->X, A24);        

        fp2mul_mont(R2->X, R->Z, R2->X);
    }
    swap = 0 ^ prevbit;
    mask = 0 - (digit_t)swap;
    swap_points(R, R2, mask);

    // trigger_high();
    // trigger_low();
}

/* loop unrolled */
static void LADDER3PT_ATTACK(const f2elm_t xP, const f2elm_t xQ, const f2elm_t xPQ, const digit_t* m, const unsigned int AliceOrBob, point_proj_t R, const f2elm_t A)
{
    point_proj_t R0 = {0}, R2 = {0};
    f2elm_t A24 = {0};
    digit_t mask;
    int i, nbits, bit, swap, prevbit = 0;
    (void) i;
    (void) nbits;

    if (AliceOrBob == ALICE) {
        nbits = OALICE_BITS;
    } else {
        nbits = OBOB_BITS - 1;
    }

    // Initializing constant
    fpcopy((digit_t*)&Montgomery_one, A24[0]);
    fp2add(A24, A24, A24);
    fp2add(A, A24, A24);
    fp2div2(A24, A24);
    fp2div2(A24, A24);  // A24 = (A+2)/4

    // Initializing points
    fp2copy(xQ, R0->X);
    fpcopy((digit_t*)&Montgomery_one, (digit_t*)R0->Z);
    fp2copy(xPQ, R2->X);
    fpcopy((digit_t*)&Montgomery_one, (digit_t*)R2->Z);
    fp2copy(xP, R->X);
    fpcopy((digit_t*)&Montgomery_one, (digit_t*)R->Z);
    fpzero((digit_t*)(R->Z)[1]);

    /*################################################
    // in order to avoid compilation effects, 
    //  we unrolled the whole loop.
    //################################################*/
    
//------------------------------------
bit = (m[0 >> LOG2RADIX] >> (0 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_0
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[1 >> LOG2RADIX] >> (1 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_1
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[2 >> LOG2RADIX] >> (2 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_2
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[3 >> LOG2RADIX] >> (3 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_3
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[4 >> LOG2RADIX] >> (4 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_4
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[5 >> LOG2RADIX] >> (5 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_5
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[6 >> LOG2RADIX] >> (6 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_6
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[7 >> LOG2RADIX] >> (7 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_7
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[8 >> LOG2RADIX] >> (8 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_8
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[9 >> LOG2RADIX] >> (9 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_9
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[10 >> LOG2RADIX] >> (10 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_10
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[11 >> LOG2RADIX] >> (11 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_11
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[12 >> LOG2RADIX] >> (12 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_12
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[13 >> LOG2RADIX] >> (13 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_13
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[14 >> LOG2RADIX] >> (14 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_14
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[15 >> LOG2RADIX] >> (15 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_15
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[16 >> LOG2RADIX] >> (16 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_16
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[17 >> LOG2RADIX] >> (17 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_17
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[18 >> LOG2RADIX] >> (18 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_18
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[19 >> LOG2RADIX] >> (19 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_19
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[20 >> LOG2RADIX] >> (20 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_20
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[21 >> LOG2RADIX] >> (21 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_21
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[22 >> LOG2RADIX] >> (22 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_22
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[23 >> LOG2RADIX] >> (23 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_23
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[24 >> LOG2RADIX] >> (24 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_24
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[25 >> LOG2RADIX] >> (25 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_25
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[26 >> LOG2RADIX] >> (26 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_26
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[27 >> LOG2RADIX] >> (27 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_27
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[28 >> LOG2RADIX] >> (28 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_28
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[29 >> LOG2RADIX] >> (29 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_29
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[30 >> LOG2RADIX] >> (30 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_30
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[31 >> LOG2RADIX] >> (31 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_31
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[32 >> LOG2RADIX] >> (32 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_32
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[33 >> LOG2RADIX] >> (33 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_33
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[34 >> LOG2RADIX] >> (34 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_34
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[35 >> LOG2RADIX] >> (35 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_35
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[36 >> LOG2RADIX] >> (36 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_36
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[37 >> LOG2RADIX] >> (37 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_37
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[38 >> LOG2RADIX] >> (38 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_38
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[39 >> LOG2RADIX] >> (39 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_39
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[40 >> LOG2RADIX] >> (40 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_40
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[41 >> LOG2RADIX] >> (41 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_41
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[42 >> LOG2RADIX] >> (42 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_42
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[43 >> LOG2RADIX] >> (43 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_43
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[44 >> LOG2RADIX] >> (44 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_44
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[45 >> LOG2RADIX] >> (45 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_45
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[46 >> LOG2RADIX] >> (46 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_46
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[47 >> LOG2RADIX] >> (47 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_47
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[48 >> LOG2RADIX] >> (48 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_48
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[49 >> LOG2RADIX] >> (49 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_49
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[50 >> LOG2RADIX] >> (50 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_50
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[51 >> LOG2RADIX] >> (51 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_51
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[52 >> LOG2RADIX] >> (52 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_52
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[53 >> LOG2RADIX] >> (53 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_53
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[54 >> LOG2RADIX] >> (54 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_54
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[55 >> LOG2RADIX] >> (55 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_55
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[56 >> LOG2RADIX] >> (56 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_56
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[57 >> LOG2RADIX] >> (57 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_57
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[58 >> LOG2RADIX] >> (58 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_58
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[59 >> LOG2RADIX] >> (59 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_59
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[60 >> LOG2RADIX] >> (60 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_60
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[61 >> LOG2RADIX] >> (61 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_61
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[62 >> LOG2RADIX] >> (62 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_62
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[63 >> LOG2RADIX] >> (63 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_63
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[64 >> LOG2RADIX] >> (64 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_64
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[65 >> LOG2RADIX] >> (65 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_65
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[66 >> LOG2RADIX] >> (66 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_66
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[67 >> LOG2RADIX] >> (67 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_67
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[68 >> LOG2RADIX] >> (68 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_68
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[69 >> LOG2RADIX] >> (69 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_69
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[70 >> LOG2RADIX] >> (70 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_70
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[71 >> LOG2RADIX] >> (71 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_71
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[72 >> LOG2RADIX] >> (72 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_72
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[73 >> LOG2RADIX] >> (73 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_73
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[74 >> LOG2RADIX] >> (74 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_74
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[75 >> LOG2RADIX] >> (75 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_75
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[76 >> LOG2RADIX] >> (76 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_76
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[77 >> LOG2RADIX] >> (77 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_77
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[78 >> LOG2RADIX] >> (78 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_78
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[79 >> LOG2RADIX] >> (79 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_79
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[80 >> LOG2RADIX] >> (80 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_80
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[81 >> LOG2RADIX] >> (81 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_81
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[82 >> LOG2RADIX] >> (82 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_82
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[83 >> LOG2RADIX] >> (83 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_83
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[84 >> LOG2RADIX] >> (84 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_84
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[85 >> LOG2RADIX] >> (85 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_85
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[86 >> LOG2RADIX] >> (86 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_86
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[87 >> LOG2RADIX] >> (87 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_87
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[88 >> LOG2RADIX] >> (88 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_88
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[89 >> LOG2RADIX] >> (89 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_89
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[90 >> LOG2RADIX] >> (90 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_90
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[91 >> LOG2RADIX] >> (91 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_91
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[92 >> LOG2RADIX] >> (92 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_92
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[93 >> LOG2RADIX] >> (93 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_93
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[94 >> LOG2RADIX] >> (94 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_94
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[95 >> LOG2RADIX] >> (95 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_95
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[96 >> LOG2RADIX] >> (96 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_96
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[97 >> LOG2RADIX] >> (97 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_97
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[98 >> LOG2RADIX] >> (98 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_98
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[99 >> LOG2RADIX] >> (99 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_99
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[100 >> LOG2RADIX] >> (100 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_100
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[101 >> LOG2RADIX] >> (101 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_101
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[102 >> LOG2RADIX] >> (102 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_102
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[103 >> LOG2RADIX] >> (103 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_103
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[104 >> LOG2RADIX] >> (104 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_104
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[105 >> LOG2RADIX] >> (105 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_105
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[106 >> LOG2RADIX] >> (106 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_106
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[107 >> LOG2RADIX] >> (107 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_107
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[108 >> LOG2RADIX] >> (108 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_108
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[109 >> LOG2RADIX] >> (109 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_109
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[110 >> LOG2RADIX] >> (110 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_110
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[111 >> LOG2RADIX] >> (111 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_111
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[112 >> LOG2RADIX] >> (112 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_112
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[113 >> LOG2RADIX] >> (113 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_113
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[114 >> LOG2RADIX] >> (114 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_114
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[115 >> LOG2RADIX] >> (115 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_115
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[116 >> LOG2RADIX] >> (116 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_116
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[117 >> LOG2RADIX] >> (117 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_117
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[118 >> LOG2RADIX] >> (118 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_118
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[119 >> LOG2RADIX] >> (119 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_119
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[120 >> LOG2RADIX] >> (120 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_120
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[121 >> LOG2RADIX] >> (121 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_121
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[122 >> LOG2RADIX] >> (122 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_122
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[123 >> LOG2RADIX] >> (123 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_123
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[124 >> LOG2RADIX] >> (124 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_124
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[125 >> LOG2RADIX] >> (125 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_125
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[126 >> LOG2RADIX] >> (126 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_126
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[127 >> LOG2RADIX] >> (127 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_127
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[128 >> LOG2RADIX] >> (128 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_128
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[129 >> LOG2RADIX] >> (129 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_129
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[130 >> LOG2RADIX] >> (130 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_130
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[131 >> LOG2RADIX] >> (131 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_131
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[132 >> LOG2RADIX] >> (132 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_132
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[133 >> LOG2RADIX] >> (133 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_133
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[134 >> LOG2RADIX] >> (134 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_134
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[135 >> LOG2RADIX] >> (135 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_135
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[136 >> LOG2RADIX] >> (136 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_136
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[137 >> LOG2RADIX] >> (137 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_137
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[138 >> LOG2RADIX] >> (138 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_138
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[139 >> LOG2RADIX] >> (139 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_139
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[140 >> LOG2RADIX] >> (140 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_140
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[141 >> LOG2RADIX] >> (141 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_141
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[142 >> LOG2RADIX] >> (142 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_142
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[143 >> LOG2RADIX] >> (143 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_143
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[144 >> LOG2RADIX] >> (144 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_144
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[145 >> LOG2RADIX] >> (145 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_145
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[146 >> LOG2RADIX] >> (146 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_146
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[147 >> LOG2RADIX] >> (147 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_147
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[148 >> LOG2RADIX] >> (148 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_148
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[149 >> LOG2RADIX] >> (149 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_149
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[150 >> LOG2RADIX] >> (150 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_150
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[151 >> LOG2RADIX] >> (151 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_151
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[152 >> LOG2RADIX] >> (152 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_152
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[153 >> LOG2RADIX] >> (153 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_153
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[154 >> LOG2RADIX] >> (154 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_154
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[155 >> LOG2RADIX] >> (155 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_155
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[156 >> LOG2RADIX] >> (156 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_156
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[157 >> LOG2RADIX] >> (157 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_157
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[158 >> LOG2RADIX] >> (158 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_158
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[159 >> LOG2RADIX] >> (159 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_159
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[160 >> LOG2RADIX] >> (160 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_160
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[161 >> LOG2RADIX] >> (161 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_161
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[162 >> LOG2RADIX] >> (162 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_162
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[163 >> LOG2RADIX] >> (163 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_163
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[164 >> LOG2RADIX] >> (164 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_164
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[165 >> LOG2RADIX] >> (165 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_165
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[166 >> LOG2RADIX] >> (166 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_166
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[167 >> LOG2RADIX] >> (167 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_167
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[168 >> LOG2RADIX] >> (168 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_168
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[169 >> LOG2RADIX] >> (169 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_169
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[170 >> LOG2RADIX] >> (170 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_170
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[171 >> LOG2RADIX] >> (171 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_171
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[172 >> LOG2RADIX] >> (172 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_172
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[173 >> LOG2RADIX] >> (173 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_173
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[174 >> LOG2RADIX] >> (174 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_174
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[175 >> LOG2RADIX] >> (175 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_175
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[176 >> LOG2RADIX] >> (176 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_176
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[177 >> LOG2RADIX] >> (177 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_177
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[178 >> LOG2RADIX] >> (178 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_178
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[179 >> LOG2RADIX] >> (179 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_179
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[180 >> LOG2RADIX] >> (180 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_180
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[181 >> LOG2RADIX] >> (181 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_181
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[182 >> LOG2RADIX] >> (182 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_182
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[183 >> LOG2RADIX] >> (183 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_183
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[184 >> LOG2RADIX] >> (184 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_184
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[185 >> LOG2RADIX] >> (185 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_185
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[186 >> LOG2RADIX] >> (186 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_186
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[187 >> LOG2RADIX] >> (187 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_187
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[188 >> LOG2RADIX] >> (188 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_188
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[189 >> LOG2RADIX] >> (189 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_189
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[190 >> LOG2RADIX] >> (190 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_190
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[191 >> LOG2RADIX] >> (191 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_191
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[192 >> LOG2RADIX] >> (192 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_192
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[193 >> LOG2RADIX] >> (193 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_193
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[194 >> LOG2RADIX] >> (194 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_194
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[195 >> LOG2RADIX] >> (195 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_195
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[196 >> LOG2RADIX] >> (196 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_196
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[197 >> LOG2RADIX] >> (197 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_197
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[198 >> LOG2RADIX] >> (198 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_198
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[199 >> LOG2RADIX] >> (199 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_199
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[200 >> LOG2RADIX] >> (200 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_200
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[201 >> LOG2RADIX] >> (201 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_201
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[202 >> LOG2RADIX] >> (202 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_202
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[203 >> LOG2RADIX] >> (203 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_203
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[204 >> LOG2RADIX] >> (204 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_204
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[205 >> LOG2RADIX] >> (205 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_205
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[206 >> LOG2RADIX] >> (206 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_206
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[207 >> LOG2RADIX] >> (207 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_207
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[208 >> LOG2RADIX] >> (208 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_208
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[209 >> LOG2RADIX] >> (209 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_209
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[210 >> LOG2RADIX] >> (210 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_210
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[211 >> LOG2RADIX] >> (211 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_211
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[212 >> LOG2RADIX] >> (212 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_212
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[213 >> LOG2RADIX] >> (213 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_213
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[214 >> LOG2RADIX] >> (214 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_214
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[215 >> LOG2RADIX] >> (215 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_215
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[216 >> LOG2RADIX] >> (216 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_216
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

//------------------------------------
bit = (m[217 >> LOG2RADIX] >> (217 & (RADIX-1))) & 1;

#ifdef ATTACK_BIT_217
    trigger_high();
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
    trigger_low();
#else
    swap = bit ^ prevbit;
    prevbit = bit;
    mask = 0 - (digit_t)swap;
#endif

swap_points(R, R2, mask);
xDBLADD(R0, R2, R->X, A24);
fp2mul_mont(R2->X, R->Z, R2->X);

    swap = 0 ^ prevbit;
    mask = 0 - (digit_t)swap;
    swap_points(R, R2, mask);

}
