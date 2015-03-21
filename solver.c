#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "image.h"
#include "solver.h"



/*THIS IS A SLOW VERSION BUT READABLE VERSION OF THE MAIN FUNCTION */
void sor_coupled_slow_but_readable(image_t *du, image_t *dv, const image_t *a11, const image_t *a12, const image_t *a22, const image_t *b1, const image_t *b2, const image_t *dpsis_horiz, const image_t *dpsis_vert, int iterations, float omega)
{
  int i,j,iter;
  float sigma_u,sigma_v,sum_dpsis,A11,A22,A12,B1,B2,det;
  int s = du->stride;
  for( iter = 0 ; iter < iterations ; iter++ )
    {
      for(j=0 ; j<du->height ; j++ )
	{
	  for( i=0 ; i<du->width ; i++ )
	    {
	      int o = j*s+i;
	      if( i==0 && j==0 )
		{
		  sum_dpsis = dpsis_horiz->data[o]               + dpsis_vert->data[o]              ;
		  sigma_u   = dpsis_horiz->data[o]*du->data[o+1] + dpsis_vert->data[o]*du->data[o+s];
		  sigma_v   = dpsis_horiz->data[o]*dv->data[o+1] + dpsis_vert->data[o]*dv->data[o+s];
		}
	      else if( i==du->width-1 && j==0 )
		{
		  sum_dpsis = dpsis_horiz->data[o-1]               + dpsis_vert->data[o]              ;
		  sigma_u   = dpsis_horiz->data[o-1]*du->data[o-1] + dpsis_vert->data[o]*du->data[o+s];
		  sigma_v   = dpsis_horiz->data[o-1]*dv->data[o-1] + dpsis_vert->data[o]*dv->data[o+s];
		}
	      else if( j==0 ) 
		{
		  sum_dpsis = dpsis_horiz->data[o-1]               + dpsis_horiz->data[o]               + dpsis_vert->data[o]              ;
		  sigma_u   = dpsis_horiz->data[o-1]*du->data[o-1] + dpsis_horiz->data[o]*du->data[o+1] + dpsis_vert->data[o]*du->data[o+s];
		  sigma_v   = dpsis_horiz->data[o-1]*dv->data[o-1] + dpsis_horiz->data[o]*dv->data[o+1] + dpsis_vert->data[o]*dv->data[o+s];
		}
	      else if ( i== 0 && j==du->height-1 )
		{
		  sum_dpsis = dpsis_horiz->data[o]               + dpsis_vert->data[o-s]              ;
		  sigma_u   = dpsis_horiz->data[o]*du->data[o+1] + dpsis_vert->data[o-s]*du->data[o-s];
		  sigma_v   = dpsis_horiz->data[o]*dv->data[o+1] + dpsis_vert->data[o-s]*dv->data[o-s];
		}
	      else if ( i==du->width-1 && j==du->height-1)
		{
		  sum_dpsis = dpsis_horiz->data[o-1]               + dpsis_vert->data[o-s]              ;
		  sigma_u   = dpsis_horiz->data[o-1]*du->data[o-1] + dpsis_vert->data[o-s]*du->data[o-s];
		  sigma_v   = dpsis_horiz->data[o-1]*dv->data[o-1] + dpsis_vert->data[o-s]*dv->data[o-s];		}
	      else if( j==du->height-1)
		{
		  sum_dpsis = dpsis_horiz->data[o-1]               + dpsis_horiz->data[o]               + dpsis_vert->data[o-s]              ;
		  sigma_u   = dpsis_horiz->data[o-1]*du->data[o-1] + dpsis_horiz->data[o]*du->data[o+1] + dpsis_vert->data[o-s]*du->data[o-s];
		  sigma_v   = dpsis_horiz->data[o-1]*dv->data[o-1] + dpsis_horiz->data[o]*dv->data[o+1] + dpsis_vert->data[o-s]*dv->data[o-s];
		}
	      else if ( i==0 )
		{
		  sum_dpsis = dpsis_horiz->data[o]               + dpsis_vert->data[o-s]               + dpsis_vert->data[o]              ;
		  sigma_u   = dpsis_horiz->data[o]*du->data[o+1] + dpsis_vert->data[o-s]*du->data[o-s] + dpsis_vert->data[o]*du->data[o+s];
		  sigma_v   = dpsis_horiz->data[o]*dv->data[o+1] + dpsis_vert->data[o-s]*dv->data[o-s] + dpsis_vert->data[o]*dv->data[o+s];
		}
	      else if( i==du->width-1)
		{
		  sum_dpsis = dpsis_horiz->data[o-1]               + dpsis_vert->data[o-s]               + dpsis_vert->data[o]              ;
		  sigma_u   = dpsis_horiz->data[o-1]*du->data[o-1] + dpsis_vert->data[o-s]*du->data[o-s] + dpsis_vert->data[o]*du->data[o+s];
		  sigma_v   = dpsis_horiz->data[o-1]*dv->data[o-1] + dpsis_vert->data[o-s]*dv->data[o-s] + dpsis_vert->data[o]*dv->data[o+s];
		}
	      else
		{
		  sum_dpsis = dpsis_horiz->data[o-1]               + dpsis_horiz->data[o]               + dpsis_vert->data[o-s]               + dpsis_vert->data[o]              ;
		  sigma_u   = dpsis_horiz->data[o-1]*du->data[o-1] + dpsis_horiz->data[o]*du->data[o+1] + dpsis_vert->data[o-s]*du->data[o-s] + dpsis_vert->data[o]*du->data[o+s];
		  sigma_v   = dpsis_horiz->data[o-1]*dv->data[o-1] + dpsis_horiz->data[o]*dv->data[o+1] + dpsis_vert->data[o-s]*dv->data[o-s] + dpsis_vert->data[o]*dv->data[o+s];
		}
	      A22 = a11->data[o]+sum_dpsis;
	      A11 = a22->data[o]+sum_dpsis;
	      A12 = -a12->data[o];
	      det = A11*A22-A12*A12;
	      A22 /= det;
	      A11 /= det;
	      A12 /= det;
	      B1 = b1->data[o]+sigma_u;
	      B2 = b2->data[o]+sigma_v;
	      du->data[o] += omega*( A11*B1 + A12*B2-du->data[o]);
	      dv->data[o] += omega*( A12*B1 + A22*B2-dv->data[o]);
	    }
	}
    }
}

/* THIS IS A FASTER VERSION BUT UNREADABLE
the first iteration is done separately from the other to compute the inverse of the 2x2 block diagonal
the loop over pixels is split in different sections: the first line, the middle lines and the last line are split
for each line, the main for loop over columns is done 4 by 4, with the first and last one done independently
only work if width>=2 & height>=2 & iterations>=1*/
void sor_coupled(image_t *du, image_t *dv, const image_t *a11, const image_t *a12, const image_t *a22, const image_t *b1, const image_t *b2, const image_t *dpsis_horiz, const image_t *dpsis_vert, int iterations, float omega)
{
  if(du->width<2 || du->height<2 || iterations < 1)
    sor_coupled_slow_but_readable(du,dv,a11,a12,a22,b1,b2,dpsis_horiz,dpsis_vert,iterations, omega);

  int i,                                         // index for rows
    j,                                           // index columns
    iter,                                        // index of iteration
    incr_line = du->stride - du->width + 1,      // increment to pass from the last column to the first column at the next row
    ibefore,                                     // index of columns for the first part of the iteration
    nbefore = (du->width-2)%4;                   // number of columns to add to have a multiple of 4 (minus 2 for the first and last colum)
  int ifst = du->width-2-nbefore,                // first value of i when decreasing, do not count first and last column, and column to have a multiple of 4
    jfst = du->height-2;                         // first value of j when decreasing without counting first and last line
  
  // to avoid compute them many times
  int stride = du->stride;
  int stride1 = stride+1;
  int stride2 = stride1+1;
  int stride3 = stride2+1;
  int stride_ = -stride;
  int stride_1 = stride_+1;
  int stride_2 = stride_1+1;
  int stride_3 = stride_2+1;

  // [A11 A12 ; A12 A22] = inv([a11 a12 ; a12 a22]) including the dpsis component
  image_t *A11 = image_new(du->width,du->height), 
    *A12 = image_new(du->width,du->height), 
    *A22 = image_new(du->width,du->height);

  float sigma_u,sigma_v,       // contains the sum of the dpsis multiply by u or v coefficient in the line except the diagonal one
    sum_dpsis,                 // sum of the dpsis coefficient in a line
    det,                       // local variable to compute determinant
    B1,B2,                     // local variable

    // Next variables as used to move along the images
    *du_ptr = du->data, *dv_ptr = dv->data, 
    *a11_ptr = a11->data, *a12_ptr = a12->data, *a22_ptr = a22->data,
    *b1_ptr = b1->data, *b2_ptr = b2->data,
    *dpsis_horiz_ptr = dpsis_horiz->data, *dpsis_vert_ptr = dpsis_vert->data,
    *A11_ptr = A11->data, *A12_ptr = A12->data, *A22_ptr = A22->data;

  // ---------------- FIRST ITERATION ----------------- //
  // reminder: inv([a b ; b c]) = [c -b ; -b a]/(ac-bb)
  //
  // for each pixel, compute sum_dpsis, sigma_u and sigma_v
  // add sum_dpsis to a11 and a22 and compute the determinant
  // deduce A11 A12 A22
  // compute B1 and B2
  // update du and dv
  // update pointer

  // ------------ first line, first column
  sum_dpsis = dpsis_horiz_ptr[0]           + dpsis_vert_ptr[0]                    ;
  sigma_u   = dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[0]*du_ptr[stride] ;
  sigma_v   = dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[0]*dv_ptr[stride] ;
  A22_ptr[0] = a11_ptr[0]+sum_dpsis;
  A11_ptr[0] = a22_ptr[0]+sum_dpsis;
  A12_ptr[0] = -a12_ptr[0];
  det = A11_ptr[0]*A22_ptr[0] - A12_ptr[0]*A12_ptr[0];
  A11_ptr[0] /= det;
  A22_ptr[0] /= det;
  A12_ptr[0] /= det;
  B1 = b1_ptr[0]+sigma_u;
  B2 = b2_ptr[0]+sigma_v;
  du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
  dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
  du_ptr++; dv_ptr++; 
  a11_ptr++; a12_ptr++; a22_ptr++;
  A11_ptr++; A12_ptr++; A22_ptr++;
  b1_ptr++; b2_ptr++;
  dpsis_horiz_ptr++; dpsis_vert_ptr++;

  // ------------ first line, column just after the first one to have a multiple of 4
  for(ibefore = nbefore ; ibefore-- ; ) // faster than for(ibefore = 0 ; ibefore < nbefore ; ibefore--)
    {
      sum_dpsis = dpsis_horiz_ptr[-1]            + dpsis_horiz_ptr[0]           + dpsis_vert_ptr[0]                    ;
      sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[0]*du_ptr[stride] ;
      sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[0]*dv_ptr[stride] ;
      A22_ptr[0] = a11_ptr[0]+sum_dpsis;
      A11_ptr[0] = a22_ptr[0]+sum_dpsis;
      A12_ptr[0] = -a12_ptr[0];
      det = A11_ptr[0]*A22_ptr[0] - A12_ptr[0]*A12_ptr[0];
      A11_ptr[0] /= det;
      A22_ptr[0] /= det;
      A12_ptr[0] /= det;
      B1 = b1_ptr[0]+sigma_u;
      B2 = b2_ptr[0]+sigma_v;
      du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
      dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
      du_ptr++; dv_ptr++; 
      a11_ptr++; a12_ptr++; a22_ptr++;
      A11_ptr++; A12_ptr++; A22_ptr++;
      b1_ptr++; b2_ptr++;
      dpsis_horiz_ptr++; dpsis_vert_ptr++;
    }

  // ------------ first line, other columns by 4
  for(i = ifst ; i ; i-=4)
    {
      // 1
      sum_dpsis = dpsis_horiz_ptr[-1]            + dpsis_horiz_ptr[0]           + dpsis_vert_ptr[0]                    ;
      sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[0]*du_ptr[stride] ;
      sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[0]*dv_ptr[stride] ;
      A22_ptr[0] = a11_ptr[0]+sum_dpsis;
      A11_ptr[0] = a22_ptr[0]+sum_dpsis;
      A12_ptr[0] = -a12_ptr[0];
      det = A11_ptr[0]*A22_ptr[0] - A12_ptr[0]*A12_ptr[0];
      A11_ptr[0] /= det;
      A22_ptr[0] /= det;
      A12_ptr[0] /= det;
      B1 = b1_ptr[0]+sigma_u;
      B2 = b2_ptr[0]+sigma_v;
      du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
      dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
      // 2
      sum_dpsis = dpsis_horiz_ptr[0]           + dpsis_horiz_ptr[1]           + dpsis_vert_ptr[1]                      ;
      sigma_u   = dpsis_horiz_ptr[0]*du_ptr[0] + dpsis_horiz_ptr[1]*du_ptr[2] + dpsis_vert_ptr[1]*du_ptr[stride1] ;
      sigma_v   = dpsis_horiz_ptr[0]*dv_ptr[0] + dpsis_horiz_ptr[1]*dv_ptr[2] + dpsis_vert_ptr[1]*dv_ptr[stride1] ;
      A22_ptr[1] = a11_ptr[1]+sum_dpsis;
      A11_ptr[1] = a22_ptr[1]+sum_dpsis;
      A12_ptr[1] = -a12_ptr[1];
      det = A11_ptr[1]*A22_ptr[1] - A12_ptr[1]*A12_ptr[1];
      A11_ptr[1] /= det;
      A22_ptr[1] /= det;
      A12_ptr[1] /= det;
      B1 = b1_ptr[1]+sigma_u;
      B2 = b2_ptr[1]+sigma_v;
      du_ptr[1] += omega*( A11_ptr[1]*B1 + A12_ptr[1]*B2 - du_ptr[1] );
      dv_ptr[1] += omega*( A12_ptr[1]*B1 + A22_ptr[1]*B2 - dv_ptr[1] );
      // 3
      sum_dpsis = dpsis_horiz_ptr[1]           + dpsis_horiz_ptr[2]           + dpsis_vert_ptr[2]                      ;
      sigma_u   = dpsis_horiz_ptr[1]*du_ptr[1] + dpsis_horiz_ptr[2]*du_ptr[3] + dpsis_vert_ptr[2]*du_ptr[stride2] ;
      sigma_v   = dpsis_horiz_ptr[1]*dv_ptr[1] + dpsis_horiz_ptr[2]*dv_ptr[3] + dpsis_vert_ptr[2]*dv_ptr[stride2] ;
      A22_ptr[2] = a11_ptr[2]+sum_dpsis;
      A11_ptr[2] = a22_ptr[2]+sum_dpsis;
      A12_ptr[2] = -a12_ptr[2];
      det = A11_ptr[2]*A22_ptr[2] - A12_ptr[2]*A12_ptr[2];
      A11_ptr[2] /= det;
      A22_ptr[2] /= det;
      A12_ptr[2] /= det;
      B1 = b1_ptr[2]+sigma_u;
      B2 = b2_ptr[2]+sigma_v;
      du_ptr[2] += omega*( A11_ptr[2]*B1 + A12_ptr[2]*B2 - du_ptr[2] );
      dv_ptr[2] += omega*( A12_ptr[2]*B1 + A22_ptr[2]*B2 - dv_ptr[2] );
      // 4
      sum_dpsis = dpsis_horiz_ptr[2]           + dpsis_horiz_ptr[3]           + dpsis_vert_ptr[3]                      ;
      sigma_u   = dpsis_horiz_ptr[2]*du_ptr[2] + dpsis_horiz_ptr[3]*du_ptr[4] + dpsis_vert_ptr[3]*du_ptr[stride3] ;
      sigma_v   = dpsis_horiz_ptr[2]*dv_ptr[2] + dpsis_horiz_ptr[3]*dv_ptr[4] + dpsis_vert_ptr[3]*dv_ptr[stride3] ;
      A22_ptr[3] = a11_ptr[3]+sum_dpsis;
      A11_ptr[3] = a22_ptr[3]+sum_dpsis;
      A12_ptr[3] = -a12_ptr[3];
      det = A11_ptr[3]*A22_ptr[3] - A12_ptr[3]*A12_ptr[3];
      A11_ptr[3] /= det;
      A22_ptr[3] /= det;
      A12_ptr[3] /= det;
      B1 = b1_ptr[3]+sigma_u;
      B2 = b2_ptr[3]+sigma_v;
      du_ptr[3] += omega*( A11_ptr[3]*B1 + A12_ptr[3]*B2 - du_ptr[3] );
      dv_ptr[3] += omega*( A12_ptr[3]*B1 + A22_ptr[3]*B2 - dv_ptr[3] );
      // increment pointer
      du_ptr += 4; dv_ptr += 4; 
      a11_ptr += 4; a12_ptr += 4; a22_ptr += 4;
      A11_ptr += 4; A12_ptr += 4; A22_ptr += 4;
      b1_ptr += 4; b2_ptr += 4;
      dpsis_horiz_ptr += 4; dpsis_vert_ptr += 4;
    }

  // ------------ first line, last column
  sum_dpsis = dpsis_horiz_ptr[-1]            + dpsis_vert_ptr[0]                    ;
  sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_vert_ptr[0]*du_ptr[stride] ;
  sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_vert_ptr[0]*dv_ptr[stride] ;
  A22_ptr[0] = a11_ptr[0]+sum_dpsis;
  A11_ptr[0] = a22_ptr[0]+sum_dpsis;
  A12_ptr[0] = -a12_ptr[0];
  det = A11_ptr[0]*A22_ptr[0] - A12_ptr[0]*A12_ptr[0];
  A11_ptr[0] /= det;
  A22_ptr[0] /= det;
  A12_ptr[0] /= det;
  B1 = b1_ptr[0]+sigma_u;
  B2 = b2_ptr[0]+sigma_v;
  du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
  dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
  // increment pointer to the next line
  du_ptr += incr_line; dv_ptr += incr_line; 
  a11_ptr += incr_line; a12_ptr += incr_line; a22_ptr += incr_line;
  A11_ptr += incr_line; A12_ptr += incr_line; A22_ptr += incr_line;
  b1_ptr += incr_line; b2_ptr += incr_line;
  dpsis_horiz_ptr += incr_line; dpsis_vert_ptr += incr_line;

  // ------------ line in the middle
  for(j = jfst ; j-- ; )    // fast than for(j=1 ; j<du->height-1 ; j--)
    {

      // ------------ line in the middle, first column
      sum_dpsis = dpsis_horiz_ptr[0]           + dpsis_vert_ptr[stride_]                     + dpsis_vert_ptr[0]                    ;
      sigma_u   = dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[stride_]*du_ptr[stride_] + dpsis_vert_ptr[0]*du_ptr[stride] ;
      sigma_v   = dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[stride_]*dv_ptr[stride_] + dpsis_vert_ptr[0]*dv_ptr[stride] ;
      A22_ptr[0] = a11_ptr[0]+sum_dpsis;
      A11_ptr[0] = a22_ptr[0]+sum_dpsis;
      A12_ptr[0] = -a12_ptr[0];
      det = A11_ptr[0]*A22_ptr[0] - A12_ptr[0]*A12_ptr[0];
      A11_ptr[0] /= det;
      A22_ptr[0] /= det;
      A12_ptr[0] /= det;
      B1 = b1_ptr[0]+sigma_u;
      B2 = b2_ptr[0]+sigma_v;
      du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
      dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
      du_ptr++; dv_ptr++; 
      a11_ptr++; a12_ptr++; a22_ptr++;
      A11_ptr++; A12_ptr++; A22_ptr++;
      b1_ptr++; b2_ptr++;
      dpsis_horiz_ptr++; dpsis_vert_ptr++;

      // ------------ line in the middle, column just after the first one to have a multiple of 4
      for(ibefore = nbefore ; ibefore-- ; ) // faster than for(ibefore = 0 ; ibefore < nbefore ; ibefore--)
	{
	  sum_dpsis = dpsis_horiz_ptr[-1]            + dpsis_horiz_ptr[0]           + dpsis_vert_ptr[stride_]                     + dpsis_vert_ptr[0]                    ;
	  sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[stride_]*du_ptr[stride_] + dpsis_vert_ptr[0]*du_ptr[stride] ;
	  sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[stride_]*dv_ptr[stride_] + dpsis_vert_ptr[0]*dv_ptr[stride] ;
	  A22_ptr[0] = a11_ptr[0]+sum_dpsis;
	  A11_ptr[0] = a22_ptr[0]+sum_dpsis;
	  A12_ptr[0] = -a12_ptr[0];
	  det = A11_ptr[0]*A22_ptr[0] - A12_ptr[0]*A12_ptr[0];
	  A11_ptr[0] /= det;
	  A22_ptr[0] /= det;
	  A12_ptr[0] /= det;
	  B1 = b1_ptr[0]+sigma_u;
	  B2 = b2_ptr[0]+sigma_v;
	  du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
	  dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
	  du_ptr++; dv_ptr++; 
	  a11_ptr++; a12_ptr++; a22_ptr++;
	  A11_ptr++; A12_ptr++; A22_ptr++;
	  b1_ptr++; b2_ptr++;
	  dpsis_horiz_ptr++; dpsis_vert_ptr++;
	}

      // ------------ line in the middle, other columns by 4
      for(i = ifst ; i ; i-=4)
	{
	  // 1
	  sum_dpsis = dpsis_horiz_ptr[-1]            + dpsis_horiz_ptr[0]           + dpsis_vert_ptr[stride_]                     + dpsis_vert_ptr[0]                    ;
	  sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[stride_]*du_ptr[stride_] + dpsis_vert_ptr[0]*du_ptr[stride] ;
	  sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[stride_]*dv_ptr[stride_] + dpsis_vert_ptr[0]*dv_ptr[stride] ;
	  A22_ptr[0] = a11_ptr[0]+sum_dpsis;
	  A11_ptr[0] = a22_ptr[0]+sum_dpsis;
	  A12_ptr[0] = -a12_ptr[0];
	  det = A11_ptr[0]*A22_ptr[0] - A12_ptr[0]*A12_ptr[0];
	  A11_ptr[0] /= det;
	  A22_ptr[0] /= det;
	  A12_ptr[0] /= det;
	  B1 = b1_ptr[0]+sigma_u;
	  B2 = b2_ptr[0]+sigma_v;
	  du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
	  dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
	  // 2
	  sum_dpsis = dpsis_horiz_ptr[0]           + dpsis_horiz_ptr[1]           + dpsis_vert_ptr[stride_1]                      + dpsis_vert_ptr[1]                      ;
	  sigma_u   = dpsis_horiz_ptr[0]*du_ptr[0] + dpsis_horiz_ptr[1]*du_ptr[2] + dpsis_vert_ptr[stride_1]*du_ptr[stride_1] + dpsis_vert_ptr[1]*du_ptr[stride1] ;
	  sigma_v   = dpsis_horiz_ptr[0]*dv_ptr[0] + dpsis_horiz_ptr[1]*dv_ptr[2] + dpsis_vert_ptr[stride_1]*dv_ptr[stride_1] + dpsis_vert_ptr[1]*dv_ptr[stride1] ;
	  A22_ptr[1] = a11_ptr[1]+sum_dpsis;
	  A11_ptr[1] = a22_ptr[1]+sum_dpsis;
	  A12_ptr[1] = -a12_ptr[1];
	  det = A11_ptr[1]*A22_ptr[1] - A12_ptr[1]*A12_ptr[1];
	  A11_ptr[1] /= det;
	  A22_ptr[1] /= det;
	  A12_ptr[1] /= det;
	  B1 = b1_ptr[1]+sigma_u;
	  B2 = b2_ptr[1]+sigma_v;
	  du_ptr[1] += omega*( A11_ptr[1]*B1 + A12_ptr[1]*B2 - du_ptr[1] );
	  dv_ptr[1] += omega*( A12_ptr[1]*B1 + A22_ptr[1]*B2 - dv_ptr[1] );
	  // 3
	  sum_dpsis = dpsis_horiz_ptr[1]           + dpsis_horiz_ptr[2]           + dpsis_vert_ptr[stride_2]                      + dpsis_vert_ptr[2]                      ;
	  sigma_u   = dpsis_horiz_ptr[1]*du_ptr[1] + dpsis_horiz_ptr[2]*du_ptr[3] + dpsis_vert_ptr[stride_2]*du_ptr[stride_2] + dpsis_vert_ptr[2]*du_ptr[stride2] ;
	  sigma_v   = dpsis_horiz_ptr[1]*dv_ptr[1] + dpsis_horiz_ptr[2]*dv_ptr[3] + dpsis_vert_ptr[stride_2]*dv_ptr[stride_2] + dpsis_vert_ptr[2]*dv_ptr[stride2] ;
	  A22_ptr[2] = a11_ptr[2]+sum_dpsis;
	  A11_ptr[2] = a22_ptr[2]+sum_dpsis;
	  A12_ptr[2] = -a12_ptr[2];
	  det = A11_ptr[2]*A22_ptr[2] - A12_ptr[2]*A12_ptr[2];
	  A11_ptr[2] /= det;
	  A22_ptr[2] /= det;
	  A12_ptr[2] /= det;
	  B1 = b1_ptr[2]+sigma_u;
	  B2 = b2_ptr[2]+sigma_v;
	  du_ptr[2] += omega*( A11_ptr[2]*B1 + A12_ptr[2]*B2 - du_ptr[2] );
	  dv_ptr[2] += omega*( A12_ptr[2]*B1 + A22_ptr[2]*B2 - dv_ptr[2] );
	  // 4
	  sum_dpsis = dpsis_horiz_ptr[2]           + dpsis_horiz_ptr[3]           + dpsis_vert_ptr[stride_3]                      + dpsis_vert_ptr[3]                      ;
	  sigma_u   = dpsis_horiz_ptr[2]*du_ptr[2] + dpsis_horiz_ptr[3]*du_ptr[4] + dpsis_vert_ptr[stride_3]*du_ptr[stride_3] + dpsis_vert_ptr[3]*du_ptr[stride3] ;
	  sigma_v   = dpsis_horiz_ptr[2]*dv_ptr[2] + dpsis_horiz_ptr[3]*dv_ptr[4] + dpsis_vert_ptr[stride_3]*dv_ptr[stride_3] + dpsis_vert_ptr[3]*dv_ptr[stride3] ;
	  A22_ptr[3] = a11_ptr[3]+sum_dpsis;
	  A11_ptr[3] = a22_ptr[3]+sum_dpsis;
	  A12_ptr[3] = -a12_ptr[3];
	  det = A11_ptr[3]*A22_ptr[3] - A12_ptr[3]*A12_ptr[3];
	  A11_ptr[3] /= det;
	  A22_ptr[3] /= det;
	  A12_ptr[3] /= det;
	  B1 = b1_ptr[3]+sigma_u;
	  B2 = b2_ptr[3]+sigma_v;
	  du_ptr[3] += omega*( A11_ptr[3]*B1 + A12_ptr[3]*B2 - du_ptr[3] );
	  dv_ptr[3] += omega*( A12_ptr[3]*B1 + A22_ptr[3]*B2 - dv_ptr[3] );
	  // increment pointer
	  du_ptr += 4; dv_ptr += 4; 
	  a11_ptr += 4; a12_ptr += 4; a22_ptr += 4;
	  A11_ptr += 4; A12_ptr += 4; A22_ptr += 4;
	  b1_ptr += 4; b2_ptr += 4;
	  dpsis_horiz_ptr += 4; dpsis_vert_ptr += 4;
	}
      
      // ------------ line in the middle, last column
      sum_dpsis = dpsis_horiz_ptr[-1]            + dpsis_vert_ptr[stride_]                     + dpsis_vert_ptr[0]                    ;
      sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_vert_ptr[stride_]*du_ptr[stride_] + dpsis_vert_ptr[0]*du_ptr[stride] ;
      sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_vert_ptr[stride_]*dv_ptr[stride_] + dpsis_vert_ptr[0]*dv_ptr[stride] ;
      A22_ptr[0] = a11_ptr[0]+sum_dpsis;
      A11_ptr[0] = a22_ptr[0]+sum_dpsis;
      A12_ptr[0] = -a12_ptr[0];
      det = A11_ptr[0]*A22_ptr[0] - A12_ptr[0]*A12_ptr[0];
      A11_ptr[0] /= det;
      A22_ptr[0] /= det;
      A12_ptr[0] /= det;
      B1 = b1_ptr[0]+sigma_u;
      B2 = b2_ptr[0]+sigma_v;
      du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
      dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
      // increment pointer to the next line
      du_ptr += incr_line; dv_ptr += incr_line; 
      a11_ptr += incr_line; a12_ptr += incr_line; a22_ptr += incr_line;
      A11_ptr += incr_line; A12_ptr += incr_line; A22_ptr += incr_line;
      b1_ptr += incr_line; b2_ptr += incr_line;
      dpsis_horiz_ptr += incr_line; dpsis_vert_ptr += incr_line;  

    }

  // ------------ last line, first column
  sum_dpsis = dpsis_horiz_ptr[0]           + dpsis_vert_ptr[stride_]                     ;
  sigma_u   = dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[stride_]*du_ptr[stride_] ;
  sigma_v   = dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[stride_]*dv_ptr[stride_] ;
  A22_ptr[0] = a11_ptr[0]+sum_dpsis;
  A11_ptr[0] = a22_ptr[0]+sum_dpsis;
  A12_ptr[0] = -a12_ptr[0];
  det = A11_ptr[0]*A22_ptr[0] - A12_ptr[0]*A12_ptr[0];
  A11_ptr[0] /= det;
  A22_ptr[0] /= det;
  A12_ptr[0] /= det;
  B1 = b1_ptr[0]+sigma_u;
  B2 = b2_ptr[0]+sigma_v;
  du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
  dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
  du_ptr++; dv_ptr++; 
  a11_ptr++; a12_ptr++; a22_ptr++;
  A11_ptr++; A12_ptr++; A22_ptr++;
  b1_ptr++; b2_ptr++;
  dpsis_horiz_ptr++; dpsis_vert_ptr++;

  // ------------ last line, column just after the first one to have a multiple of 4
  for(ibefore = nbefore ; ibefore-- ; ) // faster than for(ibefore = 0 ; ibefore < nbefore ; ibefore--)
    {
      sum_dpsis = dpsis_horiz_ptr[-1]            + dpsis_horiz_ptr[0]           + dpsis_vert_ptr[stride_]                     ;
      sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[stride_]*du_ptr[stride_] ;
      sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[stride_]*dv_ptr[stride_] ;
      A22_ptr[0] = a11_ptr[0]+sum_dpsis;
      A11_ptr[0] = a22_ptr[0]+sum_dpsis;
      A12_ptr[0] = -a12_ptr[0];
      det = A11_ptr[0]*A22_ptr[0] - A12_ptr[0]*A12_ptr[0];
      A11_ptr[0] /= det;
      A22_ptr[0] /= det;
      A12_ptr[0] /= det;
      B1 = b1_ptr[0]+sigma_u;
      B2 = b2_ptr[0]+sigma_v;
      du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
      dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
      du_ptr++; dv_ptr++; 
      a11_ptr++; a12_ptr++; a22_ptr++;
      A11_ptr++; A12_ptr++; A22_ptr++;
      b1_ptr++; b2_ptr++;
      dpsis_horiz_ptr++; dpsis_vert_ptr++;
    }

  // ------------ last line, other columns by 4
  for(i = ifst ; i ; i-=4)
    {
      // 1
      sum_dpsis = dpsis_horiz_ptr[-1]            + dpsis_horiz_ptr[0]           + dpsis_vert_ptr[stride_]                     ;
      sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[stride_]*du_ptr[stride_] ;
      sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[stride_]*dv_ptr[stride_] ;
      A22_ptr[0] = a11_ptr[0]+sum_dpsis;
      A11_ptr[0] = a22_ptr[0]+sum_dpsis;
      A12_ptr[0] = -a12_ptr[0];
      det = A11_ptr[0]*A22_ptr[0] - A12_ptr[0]*A12_ptr[0];
      A11_ptr[0] /= det;
      A22_ptr[0] /= det;
      A12_ptr[0] /= det;
      B1 = b1_ptr[0]+sigma_u;
      B2 = b2_ptr[0]+sigma_v;
      du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
      dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
      // 2
      sum_dpsis = dpsis_horiz_ptr[0]           + dpsis_horiz_ptr[1]           + dpsis_vert_ptr[stride_1]                      ;
      sigma_u   = dpsis_horiz_ptr[0]*du_ptr[0] + dpsis_horiz_ptr[1]*du_ptr[2] + dpsis_vert_ptr[stride_1]*du_ptr[stride_1] ;
      sigma_v   = dpsis_horiz_ptr[0]*dv_ptr[0] + dpsis_horiz_ptr[1]*dv_ptr[2] + dpsis_vert_ptr[stride_1]*dv_ptr[stride_1] ;
      A22_ptr[1] = a11_ptr[1]+sum_dpsis;
      A11_ptr[1] = a22_ptr[1]+sum_dpsis;
      A12_ptr[1] = -a12_ptr[1];
      det = A11_ptr[1]*A22_ptr[1] - A12_ptr[1]*A12_ptr[1];
      A11_ptr[1] /= det;
      A22_ptr[1] /= det;
      A12_ptr[1] /= det;
      B1 = b1_ptr[1]+sigma_u;
      B2 = b2_ptr[1]+sigma_v;
      du_ptr[1] += omega*( A11_ptr[1]*B1 + A12_ptr[1]*B2 - du_ptr[1] );
      dv_ptr[1] += omega*( A12_ptr[1]*B1 + A22_ptr[1]*B2 - dv_ptr[1] );
      // 3
      sum_dpsis = dpsis_horiz_ptr[1]           + dpsis_horiz_ptr[2]           + dpsis_vert_ptr[stride_2]                      ;
      sigma_u   = dpsis_horiz_ptr[1]*du_ptr[1] + dpsis_horiz_ptr[2]*du_ptr[3] + dpsis_vert_ptr[stride_2]*du_ptr[stride_2] ;
      sigma_v   = dpsis_horiz_ptr[1]*dv_ptr[1] + dpsis_horiz_ptr[2]*dv_ptr[3] + dpsis_vert_ptr[stride_2]*dv_ptr[stride_2] ;
      A22_ptr[2] = a11_ptr[2]+sum_dpsis;
      A11_ptr[2] = a22_ptr[2]+sum_dpsis;
      A12_ptr[2] = -a12_ptr[2];
      det = A11_ptr[2]*A22_ptr[2] - A12_ptr[2]*A12_ptr[2];
      A11_ptr[2] /= det;
      A22_ptr[2] /= det;
      A12_ptr[2] /= det;
      B1 = b1_ptr[2]+sigma_u;
      B2 = b2_ptr[2]+sigma_v;
      du_ptr[2] += omega*( A11_ptr[2]*B1 + A12_ptr[2]*B2 - du_ptr[2] );
      dv_ptr[2] += omega*( A12_ptr[2]*B1 + A22_ptr[2]*B2 - dv_ptr[2] );
      // 4
      sum_dpsis = dpsis_horiz_ptr[2]           + dpsis_horiz_ptr[3]           + dpsis_vert_ptr[stride_3]                      ;
      sigma_u   = dpsis_horiz_ptr[2]*du_ptr[2] + dpsis_horiz_ptr[3]*du_ptr[4] + dpsis_vert_ptr[stride_3]*du_ptr[stride_3] ;
      sigma_v   = dpsis_horiz_ptr[2]*dv_ptr[2] + dpsis_horiz_ptr[3]*dv_ptr[4] + dpsis_vert_ptr[stride_3]*dv_ptr[stride_3] ;
      A22_ptr[3] = a11_ptr[3]+sum_dpsis;
      A11_ptr[3] = a22_ptr[3]+sum_dpsis;
      A12_ptr[3] = -a12_ptr[3];
      det = A11_ptr[3]*A22_ptr[3] - A12_ptr[3]*A12_ptr[3];
      A11_ptr[3] /= det;
      A22_ptr[3] /= det;
      A12_ptr[3] /= det;
      B1 = b1_ptr[3]+sigma_u;
      B2 = b2_ptr[3]+sigma_v;
      du_ptr[3] += omega*( A11_ptr[3]*B1 + A12_ptr[3]*B2 - du_ptr[3] );
      dv_ptr[3] += omega*( A12_ptr[3]*B1 + A22_ptr[3]*B2 - dv_ptr[3] );
      // increment pointer
      du_ptr += 4; dv_ptr += 4; 
      a11_ptr += 4; a12_ptr += 4; a22_ptr += 4;
      A11_ptr += 4; A12_ptr += 4; A22_ptr += 4;
      b1_ptr += 4; b2_ptr += 4;
      dpsis_horiz_ptr += 4; dpsis_vert_ptr += 4;
    }

  // ------------ last line, last column
  sum_dpsis = dpsis_horiz_ptr[-1]            + dpsis_vert_ptr[stride_]                     ;
  sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_vert_ptr[stride_]*du_ptr[stride_] ;
  sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_vert_ptr[stride_]*dv_ptr[stride_] ;
  A22_ptr[0] = a11_ptr[0]+sum_dpsis;
  A11_ptr[0] = a22_ptr[0]+sum_dpsis;
  A12_ptr[0] = -a12_ptr[0];
  det = A11_ptr[0]*A22_ptr[0] - A12_ptr[0]*A12_ptr[0];
  A11_ptr[0] /= det;
  A22_ptr[0] /= det;
  A12_ptr[0] /= det;
  B1 = b1_ptr[0]+sigma_u;
  B2 = b2_ptr[0]+sigma_v;
  du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
  dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
  // useless to increment here

  // ---------------- OTHER ITERATION ----------------- //
  // for each pixel, compute sigma_u and sigma_v
  // compute B1 and B2
  // update du and dv
  // update pointer

  for(iter = iterations ; --iter ; ) // faster than for(iter = 1 ; iter<iterations ; iter++)
    {

      // set pointer to the beginning
      du_ptr = du->data; dv_ptr = dv->data;
      A11_ptr = A11->data; A12_ptr = A12->data; A22_ptr = A22->data;
      b1_ptr = b1->data; b2_ptr = b2->data;
      dpsis_horiz_ptr = dpsis_horiz->data; dpsis_vert_ptr = dpsis_vert->data;

      // process all elements as before
      
      // ------------ first line, first column
      sigma_u   = dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[0]*du_ptr[stride] ;
      sigma_v   = dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[0]*dv_ptr[stride] ;
      B1 = b1_ptr[0]+sigma_u;
      B2 = b2_ptr[0]+sigma_v;
      du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
      dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
      du_ptr++; dv_ptr++; 
      A11_ptr++; A12_ptr++; A22_ptr++;
      b1_ptr++; b2_ptr++;
      dpsis_horiz_ptr++; dpsis_vert_ptr++;

      // ------------ first line, column just after the first one to have a multiple of 4
      for(ibefore = nbefore ; ibefore-- ; ) // faster than for(ibefore = 0 ; ibefore < nbefore ; ibefore--)
	{
	  sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[0]*du_ptr[stride] ;
	  sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[0]*dv_ptr[stride] ;
	  B1 = b1_ptr[0]+sigma_u;
	  B2 = b2_ptr[0]+sigma_v;
	  du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
	  dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
	  du_ptr++; dv_ptr++; 
	  A11_ptr++; A12_ptr++; A22_ptr++;
	  b1_ptr++; b2_ptr++;
	  dpsis_horiz_ptr++; dpsis_vert_ptr++;
	}

      // ------------ first line, other columns by 4
      for(i = ifst ; i ; i-=4)
	{
	  // 1
	  sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[0]*du_ptr[stride] ;
	  sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[0]*dv_ptr[stride] ;
	  B1 = b1_ptr[0]+sigma_u;
	  B2 = b2_ptr[0]+sigma_v;
	  du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
	  dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
	  // 2
	  sigma_u   = dpsis_horiz_ptr[0]*du_ptr[0] + dpsis_horiz_ptr[1]*du_ptr[2] + dpsis_vert_ptr[1]*du_ptr[stride1] ;
	  sigma_v   = dpsis_horiz_ptr[0]*dv_ptr[0] + dpsis_horiz_ptr[1]*dv_ptr[2] + dpsis_vert_ptr[1]*dv_ptr[stride1] ;
	  B1 = b1_ptr[1]+sigma_u;
	  B2 = b2_ptr[1]+sigma_v;
	  du_ptr[1] += omega*( A11_ptr[1]*B1 + A12_ptr[1]*B2 - du_ptr[1] );
	  dv_ptr[1] += omega*( A12_ptr[1]*B1 + A22_ptr[1]*B2 - dv_ptr[1] );
	  // 3
	  sigma_u   = dpsis_horiz_ptr[1]*du_ptr[1] + dpsis_horiz_ptr[2]*du_ptr[3] + dpsis_vert_ptr[2]*du_ptr[stride2] ;
	  sigma_v   = dpsis_horiz_ptr[1]*dv_ptr[1] + dpsis_horiz_ptr[2]*dv_ptr[3] + dpsis_vert_ptr[2]*dv_ptr[stride2] ;
	  B1 = b1_ptr[2]+sigma_u;
	  B2 = b2_ptr[2]+sigma_v;
	  du_ptr[2] += omega*( A11_ptr[2]*B1 + A12_ptr[2]*B2 - du_ptr[2] );
	  dv_ptr[2] += omega*( A12_ptr[2]*B1 + A22_ptr[2]*B2 - dv_ptr[2] );
	  // 4
	  sigma_u   = dpsis_horiz_ptr[2]*du_ptr[2] + dpsis_horiz_ptr[3]*du_ptr[4] + dpsis_vert_ptr[3]*du_ptr[stride3] ;
	  sigma_v   = dpsis_horiz_ptr[2]*dv_ptr[2] + dpsis_horiz_ptr[3]*dv_ptr[4] + dpsis_vert_ptr[3]*dv_ptr[stride3] ;
	  B1 = b1_ptr[3]+sigma_u;
	  B2 = b2_ptr[3]+sigma_v;
	  du_ptr[3] += omega*( A11_ptr[3]*B1 + A12_ptr[3]*B2 - du_ptr[3] );
	  dv_ptr[3] += omega*( A12_ptr[3]*B1 + A22_ptr[3]*B2 - dv_ptr[3] );
	  // increment pointer
	  du_ptr += 4; dv_ptr += 4; 
	  A11_ptr += 4; A12_ptr += 4; A22_ptr += 4;
	  b1_ptr += 4; b2_ptr += 4;
	  dpsis_horiz_ptr += 4; dpsis_vert_ptr += 4;
	}
      
      // ------------ first line, last column
      sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_vert_ptr[0]*du_ptr[stride] ;
      sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_vert_ptr[0]*dv_ptr[stride] ;
      B1 = b1_ptr[0]+sigma_u;
      B2 = b2_ptr[0]+sigma_v;
      du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
      dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
      // increment pointer to the next line
      du_ptr += incr_line; dv_ptr += incr_line; 
      A11_ptr += incr_line; A12_ptr += incr_line; A22_ptr += incr_line;
      b1_ptr += incr_line; b2_ptr += incr_line;
      dpsis_horiz_ptr += incr_line; dpsis_vert_ptr += incr_line;
      
      // ------------ line in the middle
      for(j = jfst ; j-- ; )    // fast than for(j=1 ; j<du->height-1 ; j--)
	{
	  
	  // ------------ line in the middle, first column
	  sigma_u   = dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[stride_]*du_ptr[stride_] + dpsis_vert_ptr[0]*du_ptr[stride] ;
	  sigma_v   = dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[stride_]*dv_ptr[stride_] + dpsis_vert_ptr[0]*dv_ptr[stride] ;
	  B1 = b1_ptr[0]+sigma_u;
	  B2 = b2_ptr[0]+sigma_v;
	  du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
	  dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
	  du_ptr++; dv_ptr++; 
	  A11_ptr++; A12_ptr++; A22_ptr++;
	  b1_ptr++; b2_ptr++;
	  dpsis_horiz_ptr++; dpsis_vert_ptr++;
	  
	  // ------------ line in the middle, column just after the first one to have a multiple of 4
	  for(ibefore = nbefore ; ibefore-- ; ) // faster than for(ibefore = 0 ; ibefore < nbefore ; ibefore--)
	    {
	      sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[stride_]*du_ptr[stride_] + dpsis_vert_ptr[0]*du_ptr[stride] ;
	      sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[stride_]*dv_ptr[stride_] + dpsis_vert_ptr[0]*dv_ptr[stride] ;
	      B1 = b1_ptr[0]+sigma_u;
	      B2 = b2_ptr[0]+sigma_v;
	      du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
	      dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
	      du_ptr++; dv_ptr++; 
	      A11_ptr++; A12_ptr++; A22_ptr++;
	      b1_ptr++; b2_ptr++;
	      dpsis_horiz_ptr++; dpsis_vert_ptr++;
	    }
	  
	  // ------------ line in the middle, other columns by 4
	  for(i = ifst ; i ; i-=4)
	    {
	      // 1
	      sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[stride_]*du_ptr[stride_] + dpsis_vert_ptr[0]*du_ptr[stride] ;
	      sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[stride_]*dv_ptr[stride_] + dpsis_vert_ptr[0]*dv_ptr[stride] ;
	      B1 = b1_ptr[0]+sigma_u;
	      B2 = b2_ptr[0]+sigma_v;
	      du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
	      dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
	      // 2
	      sigma_u   = dpsis_horiz_ptr[0]*du_ptr[0] + dpsis_horiz_ptr[1]*du_ptr[2] + dpsis_vert_ptr[stride_1]*du_ptr[stride_1] + dpsis_vert_ptr[1]*du_ptr[stride1] ;
	      sigma_v   = dpsis_horiz_ptr[0]*dv_ptr[0] + dpsis_horiz_ptr[1]*dv_ptr[2] + dpsis_vert_ptr[stride_1]*dv_ptr[stride_1] + dpsis_vert_ptr[1]*dv_ptr[stride1] ;
	      B1 = b1_ptr[1]+sigma_u;
	      B2 = b2_ptr[1]+sigma_v;
	      du_ptr[1] += omega*( A11_ptr[1]*B1 + A12_ptr[1]*B2 - du_ptr[1] );
	      dv_ptr[1] += omega*( A12_ptr[1]*B1 + A22_ptr[1]*B2 - dv_ptr[1] );
	      // 3
	      sigma_u   = dpsis_horiz_ptr[1]*du_ptr[1] + dpsis_horiz_ptr[2]*du_ptr[3] + dpsis_vert_ptr[stride_2]*du_ptr[stride_2] + dpsis_vert_ptr[2]*du_ptr[stride2] ;
	      sigma_v   = dpsis_horiz_ptr[1]*dv_ptr[1] + dpsis_horiz_ptr[2]*dv_ptr[3] + dpsis_vert_ptr[stride_2]*dv_ptr[stride_2] + dpsis_vert_ptr[2]*dv_ptr[stride2] ;
	      B1 = b1_ptr[2]+sigma_u;
	      B2 = b2_ptr[2]+sigma_v;
	      du_ptr[2] += omega*( A11_ptr[2]*B1 + A12_ptr[2]*B2 - du_ptr[2] );
	      dv_ptr[2] += omega*( A12_ptr[2]*B1 + A22_ptr[2]*B2 - dv_ptr[2] );
	      // 4
	      sigma_u   = dpsis_horiz_ptr[2]*du_ptr[2] + dpsis_horiz_ptr[3]*du_ptr[4] + dpsis_vert_ptr[stride_3]*du_ptr[stride_3] + dpsis_vert_ptr[3]*du_ptr[stride3] ;
	      sigma_v   = dpsis_horiz_ptr[2]*dv_ptr[2] + dpsis_horiz_ptr[3]*dv_ptr[4] + dpsis_vert_ptr[stride_3]*dv_ptr[stride_3] + dpsis_vert_ptr[3]*dv_ptr[stride3] ;
	      B1 = b1_ptr[3]+sigma_u;
	      B2 = b2_ptr[3]+sigma_v;
	      du_ptr[3] += omega*( A11_ptr[3]*B1 + A12_ptr[3]*B2 - du_ptr[3] );
	      dv_ptr[3] += omega*( A12_ptr[3]*B1 + A22_ptr[3]*B2 - dv_ptr[3] );
	      // increment pointer
	      du_ptr += 4; dv_ptr += 4; 
	      A11_ptr += 4; A12_ptr += 4; A22_ptr += 4;
	      b1_ptr += 4; b2_ptr += 4;
	      dpsis_horiz_ptr += 4; dpsis_vert_ptr += 4;
	    }
	  
	  // ------------ line in the middle, last column
	  sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_vert_ptr[stride_]*du_ptr[stride_] + dpsis_vert_ptr[0]*du_ptr[stride] ;
	  sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_vert_ptr[stride_]*dv_ptr[stride_] + dpsis_vert_ptr[0]*dv_ptr[stride] ;
	  B1 = b1_ptr[0]+sigma_u;
	  B2 = b2_ptr[0]+sigma_v;
	  du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
	  dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
	  // increment pointer to the next line
	  du_ptr += incr_line; dv_ptr += incr_line; 
	  A11_ptr += incr_line; A12_ptr += incr_line; A22_ptr += incr_line;
	  b1_ptr += incr_line; b2_ptr += incr_line;
	  dpsis_horiz_ptr += incr_line; dpsis_vert_ptr += incr_line;  
	  
	}
      
      // ------------ last line, first column
      sigma_u   = dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[stride_]*du_ptr[stride_] ;
      sigma_v   = dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[stride_]*dv_ptr[stride_] ;
      B1 = b1_ptr[0]+sigma_u;
      B2 = b2_ptr[0]+sigma_v;
      du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
      dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
      du_ptr++; dv_ptr++; 
      A11_ptr++; A12_ptr++; A22_ptr++;
      b1_ptr++; b2_ptr++;
      dpsis_horiz_ptr++; dpsis_vert_ptr++;
      
      // ------------ last line, column just after the first one to have a multiple of 4
      for(ibefore = nbefore ; ibefore-- ; ) // faster than for(ibefore = 0 ; ibefore < nbefore ; ibefore--)
	{
	  sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[stride_]*du_ptr[stride_] ;
	  sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[stride_]*dv_ptr[stride_] ;
	  B1 = b1_ptr[0]+sigma_u;
	  B2 = b2_ptr[0]+sigma_v;
	  du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
	  dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
	  du_ptr++; dv_ptr++; 
	  A11_ptr++; A12_ptr++; A22_ptr++;
	  b1_ptr++; b2_ptr++;
	  dpsis_horiz_ptr++; dpsis_vert_ptr++;
	}
      
      // ------------ last line, other columns by 4
      for(i = ifst ; i ; i-=4)
	{
	  // 1
	  sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_horiz_ptr[0]*du_ptr[1] + dpsis_vert_ptr[stride_]*du_ptr[stride_] ;
	  sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_horiz_ptr[0]*dv_ptr[1] + dpsis_vert_ptr[stride_]*dv_ptr[stride_] ;
	  B1 = b1_ptr[0]+sigma_u;
	  B2 = b2_ptr[0]+sigma_v;
	  du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
	  dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
	  // 2
	  sigma_u   = dpsis_horiz_ptr[0]*du_ptr[0] + dpsis_horiz_ptr[1]*du_ptr[2] + dpsis_vert_ptr[stride_1]*du_ptr[stride_1] ;
	  sigma_v   = dpsis_horiz_ptr[0]*dv_ptr[0] + dpsis_horiz_ptr[1]*dv_ptr[2] + dpsis_vert_ptr[stride_1]*dv_ptr[stride_1] ;
	  B1 = b1_ptr[1]+sigma_u;
	  B2 = b2_ptr[1]+sigma_v;
	  du_ptr[1] += omega*( A11_ptr[1]*B1 + A12_ptr[1]*B2 - du_ptr[1] );
	  dv_ptr[1] += omega*( A12_ptr[1]*B1 + A22_ptr[1]*B2 - dv_ptr[1] );
	  // 3
	  sigma_u   = dpsis_horiz_ptr[1]*du_ptr[1] + dpsis_horiz_ptr[2]*du_ptr[3] + dpsis_vert_ptr[stride_2]*du_ptr[stride_2] ;
	  sigma_v   = dpsis_horiz_ptr[1]*dv_ptr[1] + dpsis_horiz_ptr[2]*dv_ptr[3] + dpsis_vert_ptr[stride_2]*dv_ptr[stride_2] ;
	  B1 = b1_ptr[2]+sigma_u;
	  B2 = b2_ptr[2]+sigma_v;
	  du_ptr[2] += omega*( A11_ptr[2]*B1 + A12_ptr[2]*B2 - du_ptr[2] );
	  dv_ptr[2] += omega*( A12_ptr[2]*B1 + A22_ptr[2]*B2 - dv_ptr[2] );
	  // 4
	  sigma_u   = dpsis_horiz_ptr[2]*du_ptr[2] + dpsis_horiz_ptr[3]*du_ptr[4] + dpsis_vert_ptr[stride_3]*du_ptr[stride_3] ;
	  sigma_v   = dpsis_horiz_ptr[2]*dv_ptr[2] + dpsis_horiz_ptr[3]*dv_ptr[4] + dpsis_vert_ptr[stride_3]*dv_ptr[stride_3] ;
	  B1 = b1_ptr[3]+sigma_u;
	  B2 = b2_ptr[3]+sigma_v;
	  du_ptr[3] += omega*( A11_ptr[3]*B1 + A12_ptr[3]*B2 - du_ptr[3] );
	  dv_ptr[3] += omega*( A12_ptr[3]*B1 + A22_ptr[3]*B2 - dv_ptr[3] );
	  // increment pointer
	  du_ptr += 4; dv_ptr += 4; 
	  A11_ptr += 4; A12_ptr += 4; A22_ptr += 4;
	  b1_ptr += 4; b2_ptr += 4;
	  dpsis_horiz_ptr += 4; dpsis_vert_ptr += 4;
	}
      
      // ------------ last line, last column
      sigma_u   = dpsis_horiz_ptr[-1]*du_ptr[-1] + dpsis_vert_ptr[stride_]*du_ptr[stride_] ;
      sigma_v   = dpsis_horiz_ptr[-1]*dv_ptr[-1] + dpsis_vert_ptr[stride_]*dv_ptr[stride_] ;
      B1 = b1_ptr[0]+sigma_u;
      B2 = b2_ptr[0]+sigma_v;
      du_ptr[0] += omega*( A11_ptr[0]*B1 + A12_ptr[0]*B2 - du_ptr[0] );
      dv_ptr[0] += omega*( A12_ptr[0]*B1 + A22_ptr[0]*B2 - dv_ptr[0] );
      // useless to increment here
    }
  // delete allocated images
  image_delete(A11); image_delete(A12); image_delete(A22);
}

