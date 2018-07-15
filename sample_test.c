#include "pi3.h"

//the ciphertext structure of the lhe scheme
typedef struct
{
  fq_poly_t p0;
  fq_poly_t p1;
} lhe_f;



typedef struct
{
  bn_t *a;
  bn_t *b;
  g2_t **K;
} hash_k;



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



typedef struct
{
  bn_t *c;
} vc_p;

//



int bn_mul_modz(bn_t out, bn_t in1, bn_t in2, bn_t modulus)
{
   bn_mul_basic(out,in1,in2);
   bn_mod_basic(out,out,modulus);
   return 0;
}

//this function generate a random polynomial that has small coefficients
int fq_poly_rand(fq_poly_t pol, int deg, fq_ctx_t ctx, bn_t bound)
{

  fq_poly_init(pol,ctx);

  bn_t z1;
  bn_null(z1);

  fmpz_t z2;
  fmpz_init(z2);

  fq_t z3;
  fq_init(z3,ctx);

  for (int n=0;n<deg;n++)
  {
      bn_rand_mod(z1,bound);
      bn2fmpz(z2,z1);
      fq_set_fmpz(z3,z2,ctx);
      fq_poly_set_coeff(pol,n,z3,ctx);
  }


  // release the memory
  bn_free(z1);
  fmpz_clear(z2);
  fq_clear(z3,ctx);
  return 0;

}


//this function reduce the coefficients of a polynomial
// int msg_modp(fq_poly_t m, lhe_par *par)
// {
//
//    long len=fq_poly_length(m,par->ctx);
//
//
//
//    fq_t zq;
//    fq_init(zq,par->ctx);
//
//    fmpz_t zf;
//    fmpz_init(zf);
//
//    for (int i=0;i<len;i++)
//    {
//
//      fq_poly_get_coeff(zq,m,i,par->ctx);
//      fmpz_set_str(zf,fq_get_str_pretty(zq,par->ctx), 10);
//      fmpz_mod(zf,zf,par->pf);
//      fq_set_fmpz(zq,zf,par->ctx);
//      fq_poly_set_coeff(m,i,zq,par->ctx);
//    }
//
//     fq_clear(zq,par->ctx);
//     fmpz_clear(zf);
//
//    return 0;
//
// }


//the key generation algorithm of the linearly homomorphic encryption: LHE.KeyGen
int lhe_keygen(fq_poly_t s, lhe_par *par)
{
  fq_poly_rand(s,par->ni,par->ctx,par->r);
  return 0;
}

//the encryption algorithm of the linearly homomorphic encryption: LHE.Enc
int lhe_enc(lhe_f *c, fq_poly_t s, fq_poly_t m,lhe_par *par)
{
   fq_poly_t a, e;
   fq_poly_rand(a,par->ni,par->ctx,par->q);
   fq_poly_rand(e,par->ni,par->ctx,par->r);

   fq_poly_init(c->p0,par->ctx);
   fq_poly_init(c->p1,par->ctx);

  //set the value of p1
  fq_t u1,u2;
  fq_init(u1,par->ctx);
  fq_zero(u1,par->ctx);
  fq_init(u2,par->ctx);
  fq_one(u2,par->ctx);
  fq_sub(u1,u1,u2,par->ctx);

  //fq_print(u1,par->ctx);
  fq_poly_scalar_mul_fq(c->p1,a,u1,par->ctx);

  //fq_poly_print(c->p1,par->ctx);

  //set the value of c0
  fq_poly_t v;
  fq_poly_init(v,par->ctx);
  //p0=p0+as
  fq_poly_mul(v,a,s,par->ctx);
  fq_poly_add(c->p0,c->p0,v,par->ctx);
  //p0=p0+pe
  fq_set_fmpz(u1,par->pf,par->ctx);
  //fmpz_print(par->pf);
  fq_poly_scalar_mul_fq(v,e,u1,par->ctx);
  fq_poly_add(c->p0,c->p0,v,par->ctx);
  //p0=p0+m
  fq_poly_add(c->p0,c->p0,m,par->ctx);

  //reduction modulo x^n+1
  fq_poly_t Q;
  fq_poly_init(Q,par->ctx);
  fq_poly_divrem(Q,c->p0,c->p0,par->modf,par->ctx);

  //release the memory
  fq_poly_clear(a,par->ctx);
  fq_poly_clear(e,par->ctx);
  fq_clear(u1,par->ctx);
  fq_clear(u2,par->ctx);
  fq_poly_clear(v,par->ctx);
  fq_poly_clear(Q,par->ctx);

  return 0;

}

//the decryption function of the linearly homomorphic encryption scheme: LEH.Dec
int lhe_dec(fq_poly_t m1, fq_poly_t s, lhe_f *c, lhe_par *par)
{
    fq_poly_init(m1,par->ctx);
    fq_poly_mul(m1,s,c->p1,par->ctx);
    fq_poly_add(m1,m1,c->p0,par->ctx);

    //reduction modulo x^n+1
    fq_poly_t Q;
    fq_poly_init(Q,par->ctx);
    fq_poly_divrem(Q,m1,m1,par->modf,par->ctx);

    //release the memory
    fq_poly_clear(Q,par->ctx);

    return 0;
}


// the eval function of the linearly homomorphic encryption scheme: LHE.Eval
int lhe_eval(lhe_f *c, lhe_f*L, bn_t *C,int d, lhe_par *par)
{
   fq_poly_init(c->p0,par->ctx);
   fq_poly_init(c->p1,par->ctx);

   fq_poly_t z;
   fq_poly_init(z,par->ctx);

   fq_poly_t Q;
   fq_poly_init(Q,par->ctx);

   for (int i=0; i<d; i++)
   {
       fq_t tempfq;
       fq_init(tempfq,par->ctx);
       bn2fq(tempfq,C[i],par->ctx);
       fq_poly_scalar_mul_fq(z,L[i].p0,tempfq,par->ctx);
       fq_poly_add(c->p0,c->p0,z,par->ctx);
       fq_poly_divrem(Q,c->p0,c->p0,par->modf,par->ctx);
       fq_poly_scalar_mul_fq(z,L[i].p1,tempfq,par->ctx);
       fq_poly_add(c->p1,c->p1,z,par->ctx);
       fq_poly_divrem(Q,c->p1,c->p1,par->modf,par->ctx);
   }

   //release the memory
   fq_poly_clear(z,par->ctx);
   fq_poly_clear(Q,par->ctx);
   return 0;

}


// the hash key generation algorithm
int hash_keygen(hash_k *hk, lhe_par *par)
{
   bn_t zt;
   bn_null(zt);


   hk->a=malloc(sizeof(bn_t)*2);
   hk->b=malloc(sizeof(bn_t)*par->ni);
   hk->K=malloc(sizeof(g2_t*));

   for (int i=0;i<2;i++)
   {
      bn_rand_mod(zt,par->q);
      bn_copy(hk->a[i],zt);
   }


   for (int j=0;j<par->ni;j++)
   {
      bn_rand_mod(zt,par->q);
      bn_copy(hk->b[j],zt);
   }


   g2_t zg2;
   g2_new(zg2);
   g2_set_infty(zg2);

   for (int i=0; i<2; i++)
   {
       hk->K[i]=malloc(sizeof(g2_t)*par->ni);
       for (int j=0;j<par->ni;j++)
       {
          bn_mul_modz(zt,hk->a[i],hk->b[j],par->q);
          g2_mul(zg2,par->h,zt);
          g2_copy((hk->K[i])[j],zg2);
       }
   }

   bn_free(zt);
   g2_free(zg2);
   return 0;

}

//compute the hash digest
int hash_H(g2_t hv, lhe_f *c, hash_k *hk, int mode, lhe_par *par)
{
    g2_t zg2;
    g2_new(zg2);
    g2_set_infty(zg2);

    fq_t zq;
    fq_init(zq,par->ctx);

    bn_t s, E;
    bn_new(s);

    bn_new(E);
    bn_zero(s);
    bn_zero(E);

    if (mode==0)
    {
      for (int j=0;j<par->ni;j++)
      {
         fq_poly_get_coeff(zq,c->p0,j,par->ctx);
         bn_read_str(E,fq_get_str_pretty(zq,par->ctx),strlen(fq_get_str_pretty(zq,par->ctx)),  10);
         bn_mul_modz(E,E,hk->a[0],par->q);
         bn_mul_modz(E,E,hk->b[j],par->q);
         bn_add_mod(s,s,E,par->q);
         //printf("\n\n this is E0: \n\n");
         //bn_print(s);
      }

      for (int j=0;j<par->ni;j++)
      {
        fq_init(zq,par->ctx);
        bn_new(E);
          bn_zero(E);
         fq_poly_get_coeff(zq,c->p1,j,par->ctx);
         bn_read_str(E,fq_get_str_pretty(zq,par->ctx),strlen(fq_get_str_pretty(zq,par->ctx)),  10);
         bn_mul_modz(E,E,hk->a[1],par->q);
         bn_mul_modz(E,E,hk->b[j],par->q);
         //printf("\n\n this is E1: \n\n");
         //bn_print(E);
         bn_add_mod(s,s,E,par->q);
      }
      //printf("\n\n this is s: \n\n");
      //bn_print(s);
      g2_mul(zg2,par->h,s);
      g2_copy(hv,zg2);
    }

    else
    {

      for (int j=0;j<par->ni;j++)
      {
        fq_init(zq,par->ctx);
        bn_new(E);
          bn_zero(E);
         fq_poly_get_coeff(zq,c->p0,j,par->ctx);
         bn_read_str(E,fq_get_str_pretty(zq,par->ctx), strlen(fq_get_str_pretty(zq,par->ctx)),  10);
         g2_mul(zg2,hk->K[0][j],E);
         g2_add(hv,hv,zg2);
         g2_norm(hv,hv);
      }
      printf("%s\n","fuck" );
      for (int j=0;j<par->ni;j++)
      {
         fq_poly_get_coeff(zq,c->p1,j,par->ctx);
         bn_read_str(E,fq_get_str_pretty(zq,par->ctx), strlen(fq_get_str_pretty(zq,par->ctx)), 10);
         g2_mul(zg2,hk->K[1][j],E);
         g2_add(hv,hv,zg2);
         g2_norm(hv,hv);
      }
    }

    return 0;
}



//the homomorphic evaluation algorithm of the hash function
// int hash_eval(g2_t hv, g2_t *HV, fq_t *L, int d, lhe_par *par)
// {
//    g2_t zg2;
//    g2_new(zg2);
//    g2_set_infty(zg2);
//
//    bn_t zb;
//    bn_new(zb);
//    bn_read_str(zb,fq_get_str_pretty(L[0],par->ctx), strlen(fq_get_str_pretty(L[0],par->ctx)),  10);
//
//
//    g2_mul(hv,HV[0],zb);
//    for (int i=1;i<d;i++)
//    {
//       bn_read_str(zb,fq_get_str_pretty(L[i],par->ctx), strlen(fq_get_str_pretty(L[i],par->ctx)),  10);
//       g2_mul(zg2,HV[i],zb);
//       g2_add(hv,hv,zg2);
//    }
//
//    bn_free(zb);
//    g2_free(zg2);
//
//    return 0;
// }





//the vc_keygen algorithm: generating the evaluation key and public verification key of the scheme
int vc_keygen(vc_k *vck, fq_mat_t F, lhe_par *par)
{

  //choose secret key for the lhe scheme
  lhe_keygen(vck->s,par);

  //encrypting the function F
  fq_t zq;
  fq_init(zq,par->ctx);
  fq_poly_t zqp;
  fq_poly_init(zqp,par->ctx);

  vck->F=(lhe_f**)malloc(sizeof(lhe_f*)*vck->m);
  for(int i=0;i<vck->m;i++)
  {
    vck->F[i]=(lhe_f*)malloc(sizeof(lhe_f)*vck->d);
  }


  for(int i=0;i<vck->m;i++)
  {
    printf("encrypting the function,line:%d\n",i+1);
    for(int j=0;j<vck->d;j++)
    {
      fq_set(zq,fq_mat_entry(F,i,j),par->ctx);
      fq_poly_set_fq(zqp,zq,par->ctx);
      lhe_enc(&vck->F[i][j],vck->s,zqp,par);
    }
  }
  printf("%s\n","the function has been encryped" );

  // generate the hash key
  hash_k *zhk;
  vck->hk=malloc(sizeof(*zhk));
  hash_keygen(vck->hk,par);
  g2_t ** hvx;
  hvx=(g2_t**)malloc(sizeof(g2_t*)*vck->m);
  for(int i=0;i<vck->m;i++)
  {
    hvx[i]=(g2_t*)malloc(sizeof(g2_t)*vck->d);
  }
  for(int i=0;i<vck->m;i++)
  {
    for(int j=0;j<vck->d;j++)
    {
      hash_H(hvx[i][j],&vck->F[i][j],vck->hk,0,par);
    }
  }


   vck->G=malloc(sizeof(g1_t)*vck->m);
   vck->T=malloc(sizeof(g2_t)*vck->d);
   fq_mat_t R;
   fq_mat_init(R,1,vck->m,par->ctx);
   bn_t zb;
   bn_new(zb);
   fmpz_t zf;
   fmpz_init(zf);
   fq_init(zq,par->ctx);


   for (int i=0;i<vck->m;i++)
   {
      bn_rand_mod(zb,par->q);
      bn2fmpz(zf,zb);

      fq_set_fmpz(zq,zf,par->ctx);
      fq_mat_entry_set(R,0,i,zq,par->ctx);
      g1_new(vck->G[i]);
      g1_mul(vck->G[i],par->g,zb);
   }


   printf("\n\n this is the random vector R:\n\n");
   fq_mat_print_pretty(R,par->ctx);

   for (int j=0;j<vck->d;j++)
   {
      g2_new(vck->T[j]);
      for(int i=0;i<vck->m;i++)
      {
        bn_t tempbn;
        bn_new(tempbn);
        g2_t tempg2;
        g2_new(tempg2);
        g2_set_infty(tempg2);
        fq2bn(tempbn,fq_mat_entry(R,0,i),par->ctx);
        g2_mul(tempg2,hvx[i][j],tempbn);
        g2_add(vck->T[j],vck->T[j],tempg2);
      }

   }





   //release the memory
   fq_mat_clear(R,par->ctx);
   bn_free(zb);
   bn_free(tempbn);
   fmpz_clear(zf);
   fq_clear(zq,par->ctx);
   return 0;
}


//generate a random matrix over fq
int fq_mat_randz(fq_mat_t F, int m, int d, lhe_par *par)
{
   bn_t zb;
   bn_new(zb);
   fmpz_t zf;
   fmpz_init(zf);
   fq_t zq;
   fq_init(zq,par->ctx);

   for (int i=0;i<m;i++)
   {
     for (int j=0;j<d;j++)
     {
         bn_rand_mod(zb,par->p);
         bn2fmpz(zf,zb);
         fq_set_fmpz(zq,zf,par->ctx);
         fq_mat_entry_set(F,i,j,zq,par->ctx);
     }
   }
   return 0;
}


//the vc_comp algorithm
int vc_comp(lhe_f *nv, vc_k *vck, vc_p *vcp, lhe_par *par)
{
  lhe_f *L=malloc(sizeof(lhe_f)*vck->d);
  for (int i=0;i<vck->m;i++)
  {
     for (int j=0;j<vck->d;j++)
     {
       fq_poly_init(L[j].p0,par->ctx);
       fq_poly_init(L[j].p1,par->ctx);
       fq_poly_set(L[j].p0,vck->F[i][j].p0,par->ctx);
       fq_poly_set(L[j].p1,vck->F[i][j].p1,par->ctx);
     }
     lhe_eval(&nv[i],L,vcp->c,vck->d,par);
  }
  return 0;
}

//the vc_vrfy algorithm
int vc_vrfy(vc_k *vck, vc_p *vcp, lhe_f *nv, lhe_par *par)
{
  int flag=-1;
  g2_t hv;
  g2_new(hv);
  g2_set_infty(hv);

  gt_t zgt;
  gt_new(zgt);
  gt_set_unity(zgt);


  gt_t left;
  gt_new(left);
  gt_set_unity(left);

  gt_t right;
  gt_new(right);
  gt_set_unity(right);


  for (int i=0;i<vck->m;i++)
  {

     gt_new(zgt);
     gt_set_unity(zgt);
     g2_new(hv);
     g2_set_infty(hv); // this looks completely strange, why?

     hash_H(hv,&nv[i],vck->hk,1,par);

     g2_norm(hv,hv);
     printf("\n\n this is hv of nv[%d]\n\n",i);
     g2_print(hv);
     printf("\n\n");
     pc_map(zgt,vck->G[i],hv);

     gt_mul(left,left,zgt);

  }

  printf("\n\n the verification tag left is: \n\n");
  gt_print(left);


  for (int i=0;i<vck->d;i++)
  {
     gt_new(zgt);
     gt_set_unity(zgt);
     pc_map(zgt,par->g,vck->T[i]);
     gt_exp(zgt,zgt,vcp->c[i]);
     gt_mul(right,right,zgt);
  }

  printf("\n\n the verification tag right is: \n\n");
  gt_print(right);

  flag=gt_cmp(left,right);



  if (flag ==2)
  {
     return 2;
  }

  return 0;
}



int main(int argc, char **argv)
{

//initialize the project /////////////////////////////////////////////////
        if (pi3_init())
        {
                printf("Testing FAILED\n");
                printf("Problem initializing the library\n");
               return 1;
        }


	lhe_par *par=malloc(sizeof(*par));
	lhep_new(par);
printf("\n--------------------------------------begin test-------------------------\n\n\n");

        //test the rand matrix algorithm
        int mat_m=100;
        int mat_d=100;
        fq_mat_t F;
        fq_mat_init(F,mat_m,mat_d,par->ctx);
        fq_mat_randz(F,mat_m,mat_d,par);

        //print
        printf("\n\n the matrix F is:\n\n");
        fq_mat_print_pretty(F,par->ctx);
        printf("\n" );

        //test the vc_keygen algorithm
        vc_k *vck=malloc(sizeof(*vck));
        vck->m=mat_m;
        vck->d=mat_d;

        vc_keygen(vck,F,par);
        //print
        printf("\n\n the encryped function F->p0 is:\n\n" );
        fq_poly_print_pretty(vck->F[0][0].p0,"x",par->ctx);
        printf("\n\n the integer vck->m=%d\n\n",vck->m);
        printf("\n\n the integer vck->d=%d\n\n",vck->d);
        printf("\n\n the vck->G vector is:\n\n");
        for (int i=0;i<vck->m;i++)
        {
         g1_print(vck->G[i]);
         printf("\n\n");
        }
        printf("\n\n the vck->T vector is:\n\n");
        for (int j=0;j<vck->d;j++)
        {
         g2_print(vck->T[j]);
         printf("\n\n");
	      }



        //test the vc_pgen algorithm
        vc_p *vcp=malloc(sizeof(*vcp));
        //generate a random message vector of d elements
        bn_t x[vck->d];
        for (int j=0;j<vck->d;j++)
          bn_rand_mod(x[j],par->p);
        // compute the problem instance
        vcp->c=x;
        //print
        printf("\n\n the encryption secret key vck->s is\n\n");
        fq_poly_print_pretty(vck->s,"x",par->ctx);


        printf("\n\n the input x is:\n\n");
        for (int j=0;j<vck->d;j++)
        {
           printf("\n\n this is x[%d]:\n\n",j);
           bn_print(x[j]);
           printf("\n\n");
        }


        printf("\n\n the encoded input is:\n\n");
        for (int j=0;j<vck->d;j++)
        {

           printf("\n\n this is vcp->c[%d]:\n\n",j);
           bn_print(vcp->c[j]);
           printf("\n\n\n");
        }


        //test the vc_comp algorithm
        lhe_f *nv=malloc(sizeof(lhe_f)*vck->m);
        fq_poly_t m1;
        fq_poly_init(m1,par->ctx);

        vc_comp(nv,vck,vcp,par);
        //print
        printf("\n\n the computation result is:\n\n");
        for (int i=0;i<vck->m;i++)
        {
           printf("\n\n this is nv[%d]->p0:\n\n",i);
           fq_poly_print_pretty(nv[i].p0,"x",par->ctx);
           printf("\n\n");

           printf("\n\n this is nv[%d]->p1:\n\n",i);
           fq_poly_print_pretty(nv[i].p1,"x",par->ctx);
           printf("\n\n");

           printf("\n\n the decryption of nv[%d] is:\n\n",i);
           lhe_dec(m1,vck->s,&nv[i],par);
           fq_poly_print_pretty(m1,"x",par->ctx);
           printf("\n\n\n");
        }


       // test the vc_vrfy algorithm


       int flag=vc_vrfy(vck,vcp,nv,par);
       printf("\n\n the verification result is flag=%d\n\n",flag);


printf("\n\n\n--------------------------------------end  test------------------------\n\n\n");

/*
        //release the memory
	bn_free(order);
	fq_clear(u,par->ctx);
	g2_free(zg);
	bn_free(zt);
	fq_poly_clear(s,par->ctx);
	fq_poly_clear(m,par->ctx);
	fq_poly_clear(m1,par->ctx);
	fq_poly_clear(m2,par->ctx);
*/


        return 0;
}
