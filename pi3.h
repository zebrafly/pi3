#include "lhep.h"
/*! \file pi3.h
    \brief File describing the pi3 scheme for verifying the computation

	The input x  of pi3 is not encrypted,and the function F is encrypted as F` by LHE.ENC.So the client pass (x,F`)
  to the server and the server computeF`*x by LHE.EVAL and return an out y to the client.Before passing the
  (x,F`) ,the client will generate a verification key to verify if LHE.DEC(y)=Fx after
  the computation.

*/


/*! \struct lhe_f
 * the function F is encrypted as two parts p0 and p1
 */
typedef struct
{
  fq_poly_t p0;
  fq_poly_t p1;
} lhe_f;


/*! \struct hash_k
 * a and b is used to generate private hash key and public hash key in different
 * ways.
 */
typedef struct
{
  bn_t *a;
  bn_t *b;
  g2_t **K;
} hash_k;


/*! \struct vc_k
 * m and d are the row and column of the F passing in and the vck->F going to generate.
 * F is the encrypted function
 * G and T is used to verify
 * s is the secret key
 * hk is the hash key
 */
typedef struct
{
  int m;
  int d;
  lhe_f** F;
  g1_t *G;
  g2_t *T;
  fq_poly_t s;
  hash_k *hk;
} vc_k;


/*! \struct vc_p
 * c is the encrypted message. in pi3 c is x(no encryption)
 */
typedef struct
{
  bn_t *c;
} vc_p;


/**
 * Initializing the pi3 scheme
 */
int pi3_init();


/**
 * stop the program
 */
int pi3_close();


/**
 * multify two bn_t numbers and mod the result by an bn_t number
 * @param[in] in1        a bn_t number
 * @param[in] in2        a bn_t number
 * @param[in] modulus    the bn_t number used to mod
 * @param[out]out        the out number
 */
int bn_mul_modz(bn_t out, bn_t in1, bn_t in2, bn_t modulus);


/**
 * generate a random fq_t polynomial
 * @param[in] deg        the highest degree
 * @param[in] ctx        a parameter used in flint lib
 * @param[in] bound      the bound
 * @param[out]out        the out polynomial
 */
int fq_poly_rand(fq_poly_t pol, int deg, fq_ctx_t ctx, bn_t bound);


/**
 * generate a LHE key
 * @param[in] par        some parameter in flint lib
 * @param[out]s          the secret key get from this function
 */
int lhe_keygen(fq_poly_t s, lhe_par *par);


/**
 * encrype the information
 * @param[in] s          the secret key
 * @param[in] m          the information wanted to be encryped
 * @param[in] par        some parameter in flint lib
 * @param[out]c          the out polynomial
 */
int lhe_enc(lhe_f *c, fq_poly_t s, fq_poly_t m,lhe_par *par);


/**
 * decode the information
 * @param[in] s          the secret key
 * @param[in] c          the encryped information
 * @param[in] par        some parameter in flint lib
 * @param[out]mi         the decoded result
 */
int lhe_dec(fq_poly_t m1, fq_poly_t s, lhe_f *c, lhe_par *par);


/**
 * compute
 * @param[in] L          one line of the function matrix with size 1*d
 * @param[in] C          the input message
 * @param[in] d          the column of the L and the row of the C
 * @param[in] par        some parameter in flint lib
 * @param[out]c          the computation result
 */
int lhe_eval(lhe_f *c, lhe_f*L, bn_t *C,int d, lhe_par *par);


/**
 * generate the hash key
 * @param[in] par        some parameter in flint lib
 * @param[out]hk         the hash key generated by this function
 */
int hash_keygen(hash_k *hk, lhe_par *par);


/**
 * hash the input
 * @param[in] c          the input wanted to be hashed
 * @param[in] hk         the hashkey
 * @param[in] mode       the hash mode(0:using privatekey  1:using publickey)
 * @param[in] par        some parameter in flint lib
 * @param[out]hv         the hash result
 */
int hash_H(g2_t hv, lhe_f *c, hash_k *hk, int mode, lhe_par *par);


/**
 * generate the verification key and encrype the function
 * @param[in] F 	      the input function F
 * @param[in] par       some useful parameters
 * @param[out]vck       output the encryped function vck->F and the verification key: vck->G ,vck ->T lhe key K
 */
int vc_keygen(vc_k *vck, fq_mat_t F, lhe_par *par);


/**
 * generate a random matrix F stored fq_t
 * @param[in] m        the row of the random matrix
 * @param[in] d        the column of the random matrix
 * @param[in] par      some parameter in flint lib
 * @param[out]F        the random matrix
 */
int fq_mat_randz(fq_mat_t F, int m, int d, lhe_par *par);


/**
 * compute the function matrix multiply the message marix. nv=vck->F*vcp->x
 * @param[in] vck     the vck->F is needed here to compute
 * @param[in] vcp     the vcp->x is needed here to compute
 * @param[in] par     some useful parameters
 * @param[out]nv      the outcome of vck->F*vcp->x with size m*1
 */
int vc_comp(lhe_f *nv, vc_k *vck, vc_p *vcp, lhe_par *par);


/**
 * verify the outcome. output 0 if the outcome is true, otherwise 2
 * @param[in] vck   the vck->G and the vck->T and the key K is needed here to verify
 * @param[in] vcp   the input vcp->x is needed here to verify
 * @param[in] nv    the output nv is needed here to verify
 * @param[in] par   some useful parameters
 */
int vc_vrfy(vc_k *vck, vc_p *vcp, lhe_f *nv, lhe_par *par);
