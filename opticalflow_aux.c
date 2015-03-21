#include <stdlib.h>
#include <math.h>

#include "opticalflow_aux.h"

#define datanorm 0.1f*0.1f//0.01f // square of the normalization factor
#define epsilon_color (0.001f*0.001f)//0.000001f
#define epsilon_grad (0.001f*0.001f)//0.000001f
#define epsilon_desc (0.001f*0.001f)//0.000001f
#define epsilon_smooth (0.001f*0.001f)//0.000001f

/* warp a color image according to a flow. src is the input image, wx and wy, the input flow. dst is the warped image and mask contains 0 or 1 if the pixels goes outside/inside image boundaries */
void image_warp(color_image_t *dst, image_t *mask, const color_image_t *src, const image_t *wx, const image_t *wy)
{
  int i, j, offset, offset_color, incr_line = mask->stride-mask->width, x, y, x1, x2, y1, y2;
  float xx, yy, dx, dy;
  for(j=0,offset_color=0,offset=0 ; j<src->height ; j++)
    {
      for(i=0 ; i<src->width ; i++,offset++,offset_color++)
	{
	  xx = i+wx->data[offset];
	  yy = j+wy->data[offset];
	  x = floor(xx);
	  y = floor(yy);
	  dx = xx-x;
	  dy = yy-y;
	  mask->data[offset] = (xx>=0 && xx<=src->width-1 && yy>=0 && yy<=src->height-1);
	  x1 = MINMAX(x,src->width);
	  x2 = MINMAX(x+1,src->width);
	  y1 = MINMAX(y,src->height);
	  y2 = MINMAX(y+1,src->height);
	  dst->c1[offset_color] = 
	    src->c1[y1*src->width+x1]*(1.0f-dx)*(1.0f-dy) +
	    src->c1[y1*src->width+x2]*dx*(1.0f-dy) +
	    src->c1[y2*src->width+x1]*(1.0f-dx)*dy +
	    src->c1[y2*src->width+x2]*dx*dy;
	  dst->c2[offset_color] = 
	    src->c2[y1*src->width+x1]*(1.0f-dx)*(1.0f-dy) +
	    src->c2[y1*src->width+x2]*dx*(1.0f-dy) +
	    src->c2[y2*src->width+x1]*(1.0f-dx)*dy +
	    src->c2[y2*src->width+x2]*dx*dy;
	  dst->c3[offset_color] = 
	    src->c3[y1*src->width+x1]*(1.0f-dx)*(1.0f-dy) +
	    src->c3[y1*src->width+x2]*dx*(1.0f-dy) +
	    src->c3[y2*src->width+x1]*(1.0f-dx)*dy +
	    src->c3[y2*src->width+x2]*dx*dy;
	}
      offset += incr_line;
    }
}

/* compute image first and second order spatio-temporal derivatives of a color image */
void get_derivatives(const color_image_t *im1, const color_image_t *im2, const convolution_t *deriv,
		     color_image_t *dx, color_image_t *dy, color_image_t *dt, 
		     color_image_t *dxx, color_image_t *dxy, color_image_t *dyy, color_image_t *dxt, color_image_t *dyt)
{
  // derivatives are computed on the mean of the first image and the warped second image
  color_image_t *tmp_im2 = color_image_new(im2->width,im2->height); 
  int i=0;
  for(i=0 ; i<im1->height*im1->width ; i++)
    {
      tmp_im2->c1[i] = 0.5f * (im2->c1[i]+im1->c1[i]);
      tmp_im2->c2[i] = 0.5f * (im2->c2[i]+im1->c2[i]);
      tmp_im2->c3[i] = 0.5f * (im2->c3[i]+im1->c3[i]);	  
      dt->c1[i] = im2->c1[i]-im1->c1[i];
      dt->c2[i] = im2->c2[i]-im1->c2[i];
      dt->c3[i] = im2->c3[i]-im1->c3[i];
    }   
  // compute all other derivatives
  color_image_convolve_hv(dx, tmp_im2, deriv, NULL);
  color_image_convolve_hv(dy, tmp_im2, NULL, deriv);
  color_image_convolve_hv(dxx, dx, deriv, NULL);
  color_image_convolve_hv(dxy, dx, NULL, deriv);
  color_image_convolve_hv(dyy, dy, NULL, deriv);
  color_image_convolve_hv(dxt, dt, deriv, NULL);
  color_image_convolve_hv(dyt, dt, NULL, deriv);
  // free memory
  color_image_delete(tmp_im2);
}

/* compute the smoothness term */
/* It is represented as two images, the first one for horizontal smoothness, the second for vertical
   in dst_horiz, the pixel i,j represents the smoothness weight between pixel i,j and i,j+1
   in dst_vert, the pixel i,j represents the smoothness weight between pixel i,j and i+1,j */
void compute_smoothness(image_t *dst_horiz, image_t *dst_vert, const image_t *uu, const image_t *vv, const convolution_t *deriv_flow, float quarter_alpha)
{
  int width = uu->width, height = vv->height, stride = uu->stride, j;
  image_t *ux = image_new(width,height), *vx = image_new(width,height), *uy = image_new(width,height), *vy = image_new(width,height), *smoothness = image_new(width,height);
  // compute derivatives [-0.5 0 0.5]
  convolve_horiz(ux, uu, deriv_flow);
  convolve_horiz(vx, vv, deriv_flow);
  convolve_vert(uy, uu, deriv_flow);
  convolve_vert(vy, vv, deriv_flow);
  // compute smoothness
  for(j=0 ; j<height ; j++)
    {
      int i, offset = j*stride;
      for(i=0 ; i<width ; i++,offset++)
	{
	  float tmp;
	  tmp  = ux->data[offset]*ux->data[offset];
	  tmp += uy->data[offset]*uy->data[offset];
	  tmp += vx->data[offset]*vx->data[offset];
	  tmp += vy->data[offset]*vy->data[offset];
	  smoothness->data[offset] = quarter_alpha / sqrt(tmp+epsilon_smooth);
	}
    }
  image_delete(ux); image_delete(uy); image_delete(vx); image_delete(vy); 
  // compute dst_horiz
  for(j=0 ; j<height ; j++)
    {
      int i, offset = j*stride;
      for(i=0 ; i<width-1 ; i++,offset++)
	{ 
	  dst_horiz->data[offset] = smoothness->data[offset]+smoothness->data[offset+1];
	}
    }
  // compute dst_vert
  for(j=0 ; j<height-1 ; j++)
    {
      int i, offset = j*stride;
      for(i=0 ; i<width ; i++,offset++)
	{ 
	  dst_vert->data[offset] = smoothness->data[offset]+smoothness->data[offset+stride];
	}
    }
  image_delete(smoothness);
}

/* sub the laplacian (smoothness term) to the right-hand term */
void sub_laplacian(image_t *dst, const image_t *src, const image_t *weight_horiz, const image_t *weight_vert)
{
  float tmp;
  int i,j,offsetline = src->stride-src->width;
  float *src_ptr = src->data, *dst_ptr = dst->data, *weight_horiz_ptr = weight_horiz->data, *weight_vert_ptr = weight_vert->data;
  // horizontal filtering
  j = src->height;
  while(j--) // faster than for(j=0;j<src->height;j++)
    {
      i = src->width-1;
      while(i--)
	{
	  tmp = (*weight_horiz_ptr)*((*(src_ptr+1))-(*src_ptr));
	  *dst_ptr += tmp;
	  *(dst_ptr+1) -= tmp;
	  dst_ptr++;
	  src_ptr++;
	  weight_horiz_ptr++;
	}
      dst_ptr += offsetline+1;
      src_ptr += offsetline+1;
      weight_horiz_ptr += offsetline+1;
    }
  
  src_ptr = src->data;
  dst_ptr = dst->data;
  // vertical filtering
  j = src->height-1;
  while(j--)
    {
      i = src->width;
      while(i--)
	{
	  tmp = (*weight_vert_ptr)*((*(src_ptr+src->stride))-(*src_ptr));
	  *dst_ptr += tmp;
	  *(dst_ptr+src->stride) -= tmp;
	  dst_ptr++;
	  src_ptr++;
	  weight_vert_ptr++;
	}
      dst_ptr += offsetline;
      src_ptr += offsetline;
      weight_vert_ptr += offsetline;
    }
}

/* compute the dataterm and the matching term
   a11 a12 a22 represents the 2x2 diagonal matrix, b1 and b2 the right hand side
   other (color) images are input */
void compute_data_and_match(image_t *a11, image_t *a12, image_t *a22, image_t *b1, image_t *b2, image_t *mask, image_t *wx, image_t *wy, image_t *du, image_t *dv, image_t *uu, image_t *vv, color_image_t *Ix, color_image_t *Iy, color_image_t *Iz, color_image_t *Ixx, color_image_t *Ixy, color_image_t *Iyy, color_image_t *Ixz, color_image_t *Iyz, image_t *desc_weight, image_t *desc_flow_x, image_t *desc_flow_y, float half_delta_over3, float half_beta, float half_gamma_over3)
{
  float *du_ptr, *dv_ptr, 
    *uu_ptr, *vv_ptr,
    *wx_ptr, *wy_ptr,
    *mask_ptr, 
    *a11_ptr, *a12_ptr, *a22_ptr, 
    *b1_ptr, *b2_ptr, 
    *ix_c1_ptr, *iy_c1_ptr, *iz_c1_ptr, *ixx_c1_ptr, *ixy_c1_ptr, *iyy_c1_ptr, *ixz_c1_ptr, *iyz_c1_ptr,
    *ix_c2_ptr, *iy_c2_ptr, *iz_c2_ptr, *ixx_c2_ptr, *ixy_c2_ptr, *iyy_c2_ptr, *ixz_c2_ptr, *iyz_c2_ptr,
    *ix_c3_ptr, *iy_c3_ptr, *iz_c3_ptr, *ixx_c3_ptr, *ixy_c3_ptr, *iyy_c3_ptr, *ixz_c3_ptr, *iyz_c3_ptr,
    *desc_weight_ptr, *desc_flow_x_ptr, *desc_flow_y_ptr;
  int offsetline = uu->stride-uu->width, i, j;
  float n1, n2, n3, n4, n5, n6;

  double tmp, tmp2, tmp3, tmp4, tmp5, tmp6;
  desc_weight_ptr = desc_weight->data; desc_flow_x_ptr = desc_flow_x->data; desc_flow_y_ptr = desc_flow_y->data;
  du_ptr = du->data; dv_ptr = dv->data;
  uu_ptr = uu->data; vv_ptr = vv->data;
  wx_ptr = wx->data; wy_ptr = wy->data;
  mask_ptr = mask->data; 
  a11_ptr = a11->data; a12_ptr = a12->data; a22_ptr = a22->data;
  b1_ptr = b1->data; b2_ptr = b2->data;
  ix_c1_ptr = Ix->c1; iy_c1_ptr = Iy->c1; iz_c1_ptr = Iz->c1; ixx_c1_ptr = Ixx->c1; ixy_c1_ptr = Ixy->c1; iyy_c1_ptr = Iyy->c1; ixz_c1_ptr = Ixz->c1; iyz_c1_ptr = Iyz->c1;
  ix_c2_ptr = Ix->c2; iy_c2_ptr = Iy->c2; iz_c2_ptr = Iz->c2; ixx_c2_ptr = Ixx->c2; ixy_c2_ptr = Ixy->c2; iyy_c2_ptr = Iyy->c2; ixz_c2_ptr = Ixz->c2; iyz_c2_ptr = Iyz->c2;
  ix_c3_ptr = Ix->c3; iy_c3_ptr = Iy->c3; iz_c3_ptr = Iz->c3; ixx_c3_ptr = Ixx->c3; ixy_c3_ptr = Ixy->c3; iyy_c3_ptr = Iyy->c3; ixz_c3_ptr = Ixz->c3; iyz_c3_ptr = Iyz->c3;

  j = uu->height;
  while(j--) // faster than for (j=0 ; j<uu->height ; j++)
    {
      i = uu->width;
      while(i--)
	{
	  if(*mask_ptr)
	    {
	      // color
	      tmp  = *iz_c1_ptr + (*ix_c1_ptr)*(*du_ptr) + (*iy_c1_ptr)*(*dv_ptr);
	      n1 = (*ix_c1_ptr) * (*ix_c1_ptr) + (*iy_c1_ptr) * (*iy_c1_ptr) + datanorm;
	      tmp2 = *iz_c2_ptr + (*ix_c2_ptr)*(*du_ptr) + (*iy_c2_ptr)*(*dv_ptr);
	      n2 = (*ix_c2_ptr) * (*ix_c2_ptr) + (*iy_c2_ptr) * (*iy_c2_ptr) + datanorm;
	      tmp3 = *iz_c3_ptr + (*ix_c3_ptr)*(*du_ptr) + (*iy_c3_ptr)*(*dv_ptr);
	      n3 = (*ix_c3_ptr) * (*ix_c3_ptr) + (*iy_c3_ptr) * (*iy_c3_ptr) + datanorm;
	      tmp = half_delta_over3 / sqrt(tmp*tmp/n1 + tmp2*tmp2/n2 + tmp3*tmp3/n3 + epsilon_color);
	      *a11_ptr = tmp * (*ix_c1_ptr) * (*ix_c1_ptr) / n1;
	      *a12_ptr = tmp * (*ix_c1_ptr) * (*iy_c1_ptr) / n1;
	      *a22_ptr = tmp * (*iy_c1_ptr) * (*iy_c1_ptr) / n1;
	      *b1_ptr = -tmp * (*iz_c1_ptr) * (*ix_c1_ptr) / n1;
	      *b2_ptr = -tmp * (*iz_c1_ptr) * (*iy_c1_ptr) / n1;
	      *a11_ptr += tmp * (*ix_c2_ptr) * (*ix_c2_ptr) / n2;
	      *a12_ptr += tmp * (*ix_c2_ptr) * (*iy_c2_ptr) / n2;
	      *a22_ptr += tmp * (*iy_c2_ptr) * (*iy_c2_ptr) / n2;
	      *b1_ptr -= tmp * (*iz_c2_ptr) * (*ix_c2_ptr) / n2;
	      *b2_ptr -= tmp * (*iz_c2_ptr) * (*iy_c2_ptr) / n2 ;
	      *a11_ptr += tmp * (*ix_c3_ptr) * (*ix_c3_ptr) / n3;
	      *a12_ptr += tmp * (*ix_c3_ptr) * (*iy_c3_ptr) / n3;
	      *a22_ptr += tmp * (*iy_c3_ptr) * (*iy_c3_ptr) / n3;
	      *b1_ptr -= tmp * (*iz_c3_ptr) * (*ix_c3_ptr) / n3;
	      *b2_ptr -= tmp * (*iz_c3_ptr) * (*iy_c3_ptr) / n3;
	      // gradient 
	      n1 = (*ixx_c1_ptr) * (*ixx_c1_ptr) + (*ixy_c1_ptr) * (*ixy_c1_ptr) + datanorm;
	      n2 = (*iyy_c1_ptr) * (*iyy_c1_ptr) + (*ixy_c1_ptr) * (*ixy_c1_ptr) + datanorm;
	      tmp  = *ixz_c1_ptr + (*ixx_c1_ptr) * (*du_ptr) + (*ixy_c1_ptr) * (*dv_ptr);
	      tmp2 = *iyz_c1_ptr + (*ixy_c1_ptr) * (*du_ptr) + (*iyy_c1_ptr) * (*dv_ptr);
	      n3 = (*ixx_c2_ptr) * (*ixx_c2_ptr) + (*ixy_c2_ptr) * (*ixy_c2_ptr) + datanorm;
	      n4 = (*iyy_c2_ptr) * (*iyy_c2_ptr) + (*ixy_c2_ptr) * (*ixy_c2_ptr) + datanorm;
	      tmp3 = *ixz_c2_ptr + (*ixx_c2_ptr) * (*du_ptr) + (*ixy_c2_ptr) * (*dv_ptr);
	      tmp4 = *iyz_c2_ptr + (*ixy_c2_ptr) * (*du_ptr) + (*iyy_c2_ptr) * (*dv_ptr);
	      n5 = (*ixx_c3_ptr) * (*ixx_c3_ptr) + (*ixy_c3_ptr) * (*ixy_c3_ptr) + datanorm;
	      n6 = (*iyy_c3_ptr) * (*iyy_c3_ptr) + (*ixy_c3_ptr) * (*ixy_c3_ptr) + datanorm;
	      tmp5 = *ixz_c3_ptr + (*ixx_c3_ptr) * (*du_ptr) + (*ixy_c3_ptr) * (*dv_ptr);
	      tmp6 = *iyz_c3_ptr + (*ixy_c3_ptr) * (*du_ptr) + (*iyy_c3_ptr) * (*dv_ptr);
	      tmp = half_gamma_over3 / sqrt(tmp*tmp/n1 + tmp2*tmp2/n2 + tmp3*tmp3/n3 + tmp4*tmp4/n4 + tmp5*tmp5/n5 + tmp6*tmp6/n6 + epsilon_grad);
	      *a11_ptr += tmp * ( (*ixx_c1_ptr)*(*ixx_c1_ptr)/n1 + (*ixy_c1_ptr)*(*ixy_c1_ptr)/n2 );
	      *a12_ptr += tmp * ( (*ixx_c1_ptr)*(*ixy_c1_ptr)/n1 + (*ixy_c1_ptr)*(*iyy_c1_ptr)/n2 );
	      *a22_ptr += tmp * ( (*iyy_c1_ptr)*(*iyy_c1_ptr)/n2 + (*ixy_c1_ptr)*(*ixy_c1_ptr)/n1 );
	      *b1_ptr -= tmp * ( (*ixx_c1_ptr)*(*ixz_c1_ptr)/n1 + (*ixy_c1_ptr)*(*iyz_c1_ptr)/n2 );
	      *b2_ptr -= tmp * ( (*iyy_c1_ptr)*(*iyz_c1_ptr)/n2 + (*ixy_c1_ptr)*(*ixz_c1_ptr/n1) );
	      *a11_ptr += tmp * ( (*ixx_c2_ptr)*(*ixx_c2_ptr)/n3 + (*ixy_c2_ptr)*(*ixy_c2_ptr)/n4 );
	      *a12_ptr += tmp * ( (*ixx_c2_ptr)*(*ixy_c2_ptr)/n3 + (*ixy_c2_ptr)*(*iyy_c2_ptr)/n4 );
	      *a22_ptr += tmp * ( (*iyy_c2_ptr)*(*iyy_c2_ptr)/n4 + (*ixy_c2_ptr)*(*ixy_c2_ptr)/n3 );
	      *b1_ptr -= tmp * ( (*ixx_c2_ptr)*(*ixz_c2_ptr)/n3 + (*ixy_c2_ptr)*(*iyz_c2_ptr)/n4 );
	      *b2_ptr -= tmp * ( (*iyy_c2_ptr)*(*iyz_c2_ptr)/n4 + (*ixy_c2_ptr)*(*ixz_c2_ptr)/n3 );
	      *a11_ptr += tmp * ( (*ixx_c3_ptr)*(*ixx_c3_ptr)/n5 + (*ixy_c3_ptr)*(*ixy_c3_ptr)/n6 );
	      *a12_ptr += tmp * ( (*ixx_c3_ptr)*(*ixy_c3_ptr)/n5 + (*ixy_c3_ptr)*(*iyy_c3_ptr)/n6 );
	      *a22_ptr += tmp * ( (*iyy_c3_ptr)*(*iyy_c3_ptr)/n6 + (*ixy_c3_ptr)*(*ixy_c3_ptr)/n5 );
	      *b1_ptr -= tmp * ( (*ixx_c3_ptr)*(*ixz_c3_ptr)/n5 + (*ixy_c3_ptr)*(*iyz_c3_ptr)/n6 );
	      *b2_ptr -= tmp * ( (*iyy_c3_ptr)*(*iyz_c3_ptr)/n6 + (*ixy_c3_ptr)*(*ixz_c3_ptr)/n5 );     
	    }
	  else
	    {
	      *a11_ptr = 0.0f;
	      *a12_ptr = 0.0f;
	      *a22_ptr = 0.0f;
	      *b1_ptr = 0.0f;
	      *b2_ptr = 0.0f;
	    }
	  if(half_beta)
	    {
	      // descriptor
	      tmp  = *uu_ptr - (*desc_flow_x_ptr);
	      tmp2 = *vv_ptr - (*desc_flow_y_ptr);
	      tmp = half_beta*(*desc_weight_ptr)/sqrt(tmp*tmp+tmp2*tmp2+epsilon_desc);
	      *a11_ptr += tmp;
	      *a22_ptr += tmp;
	      *b1_ptr -= tmp*((*wx_ptr)-(*desc_flow_x_ptr));
	      *b2_ptr -= tmp*((*wy_ptr)-(*desc_flow_y_ptr));
	      uu_ptr++; vv_ptr++;
	      wx_ptr++; wy_ptr++;
	      desc_flow_x_ptr++; desc_flow_y_ptr++; desc_weight_ptr++;
	    }
	  du_ptr++; dv_ptr++;
	  mask_ptr++; 
	  a11_ptr++; a12_ptr++; a22_ptr++; 
	  b1_ptr++; b2_ptr++; 
	  ix_c1_ptr++; iy_c1_ptr++; iz_c1_ptr++; ixx_c1_ptr++; ixy_c1_ptr++; iyy_c1_ptr++; ixz_c1_ptr++; iyz_c1_ptr++;
	  ix_c2_ptr++; iy_c2_ptr++; iz_c2_ptr++; ixx_c2_ptr++; ixy_c2_ptr++; iyy_c2_ptr++; ixz_c2_ptr++; iyz_c2_ptr++;
	  ix_c3_ptr++; iy_c3_ptr++; iz_c3_ptr++; ixx_c3_ptr++; ixy_c3_ptr++; iyy_c3_ptr++; ixz_c3_ptr++; iyz_c3_ptr++;
	}
      du_ptr += offsetline; dv_ptr += offsetline; 
      uu_ptr += offsetline; vv_ptr += offsetline; 
      wx_ptr += offsetline; wy_ptr += offsetline; 
      mask_ptr += offsetline; 
      a11_ptr += offsetline; a12_ptr += offsetline; a22_ptr += offsetline;
      b1_ptr += offsetline; b2_ptr += offsetline;
      desc_weight_ptr += offsetline; desc_flow_x_ptr += offsetline; desc_flow_y_ptr += offsetline;
    }
}


/* resize the descriptors to the new size using a weighted mean */
void descflow_resize(image_t *dst_flow_x, image_t *dst_flow_y, image_t *dst_weight, const image_t *src_flow_x, const image_t *src_flow_y, const image_t *src_weight)
{
  int src_width = src_flow_x->width, src_height = src_flow_x->height, src_stride = src_flow_x->stride,
    dst_width = dst_flow_x->width, dst_height = dst_flow_x->height, dst_stride = dst_flow_x->stride;
  float scale_x = ((float)dst_width-1)/((float)src_width-1), scale_y = ((float)dst_height-1)/((float)src_height-1);
  image_erase(dst_flow_x); image_erase(dst_flow_y); image_erase(dst_weight);
  int j;
  for( j=0 ; j<src_height ; j++)
    {
      float yy = ((float)j)*scale_y;
      float yyf = floor(yy);
      float dy = yy-yyf;
      int y1 = MINMAX( (int) yyf   , dst_height);
      int y2 = MINMAX( (int) yyf+1 , dst_height);
      int i;
      for( i=0 ; i<src_width ; i++ )
	{
	  float weight = src_weight->data[j*src_stride+i];
	  if( weight<0.0000000001f )
	    continue;
	  float xx = ((float)i)*scale_x;
	  float xxf = floor(xx);
	  float dx = xx-xxf;
	  int x1 = MINMAX( (int) xxf   , dst_width);
	  int x2 = MINMAX( (int) xxf+1 , dst_width);
	  float weightxy, newweight;
	  if( dx )
	    {
	      if( dy )
		{
		  weightxy = weight*dx*dy;
		  newweight = dst_weight->data[y2*dst_stride+x2] + weightxy;
		  dst_flow_x->data[y2*dst_stride+x2] = (dst_flow_x->data[y2*dst_stride+x2]*dst_weight->data[y2*dst_stride+x2] + src_flow_x->data[j*src_stride+i]*weightxy*scale_x)/newweight;
		  dst_flow_y->data[y2*dst_stride+x2] = (dst_flow_y->data[y2*dst_stride+x2]*dst_weight->data[y2*dst_stride+x2] + src_flow_y->data[j*src_stride+i]*weightxy*scale_y)/newweight;
		  dst_weight->data[y2*dst_stride+x2] = newweight;
		}
	      weightxy = weight*dx*(1.0f-dy);
	      newweight = dst_weight->data[y1*dst_stride+x2] + weightxy;
	      dst_flow_x->data[y1*dst_stride+x2] = (dst_flow_x->data[y1*dst_stride+x2]*dst_weight->data[y1*dst_stride+x2] + src_flow_x->data[j*src_stride+i]*weightxy*scale_x)/newweight;
	      dst_flow_y->data[y1*dst_stride+x2] = (dst_flow_y->data[y1*dst_stride+x2]*dst_weight->data[y1*dst_stride+x2] + src_flow_y->data[j*src_stride+i]*weightxy*scale_y)/newweight;
	      dst_weight->data[y1*dst_stride+x2] = newweight;
	    }
	  if( dy )
	    {
	      weightxy = weight*(1.0f-dx)*dy;
	      newweight = dst_weight->data[y2*dst_stride+x1] + weightxy;
	      dst_flow_x->data[y2*dst_stride+x1] = (dst_flow_x->data[y2*dst_stride+x1]*dst_weight->data[y2*dst_stride+x1] + src_flow_x->data[j*src_stride+i]*weightxy*scale_x)/newweight;
	      dst_flow_y->data[y2*dst_stride+x1] = (dst_flow_y->data[y2*dst_stride+x1]*dst_weight->data[y2*dst_stride+x1] + src_flow_y->data[j*src_stride+i]*weightxy*scale_y)/newweight;
	      dst_weight->data[y2*dst_stride+x1] = newweight;
	    }
	  weightxy = weight*(1.0f-dx)*(1.0f-dy);
	  newweight = dst_weight->data[y1*dst_stride+x1] + weightxy;
	  dst_flow_x->data[y1*dst_stride+x1] = (dst_flow_x->data[y1*dst_stride+x1]*dst_weight->data[y1*dst_stride+x1] + src_flow_x->data[j*src_stride+i]*weightxy*scale_x)/newweight;
	  dst_flow_y->data[y1*dst_stride+x1] = (dst_flow_y->data[y1*dst_stride+x1]*dst_weight->data[y1*dst_stride+x1] + src_flow_y->data[j*src_stride+i]*weightxy*scale_y)/newweight;
	  dst_weight->data[y1*dst_stride+x1] = newweight;
	}
    }
}

/* resize the descriptors to the new size using a nearest neighbor method while keeping the descriptor with the higher weight at the end */
void descflow_resize_nn(image_t *dst_flow_x, image_t *dst_flow_y, image_t *dst_weight, const image_t *src_flow_x, const image_t *src_flow_y, const image_t *src_weight)
{
  int src_width = src_flow_x->width, src_height = src_flow_x->height, src_stride = src_flow_x->stride,
    dst_width = dst_flow_x->width, dst_height = dst_flow_x->height, dst_stride = dst_flow_x->stride;
  float scale_x = ((float)dst_width-1)/((float)src_width-1), scale_y = ((float)dst_height-1)/((float)src_height-1);
  image_erase(dst_flow_x); image_erase(dst_flow_y); image_erase(dst_weight);
  int j;
  for( j=0 ; j<src_height ; j++)
    {
      float yy = ((float)j)*scale_y;
      int y = (int) 0.5f+yy; // equivalent to round(yy)
      int i;
      for( i=0 ; i<src_width ; i++ )
	{
	  float weight = src_weight->data[j*src_stride+i];
	  if( !weight )
	    continue;
	  float xx = ((float)i)*scale_x;
	  int x = (int) 0.5f+xx; // equivalent to round(xx)
	  if( dst_weight->data[y*dst_stride+x] < weight )
	    {
	      dst_weight->data[y*dst_stride+x] = weight;
	      dst_flow_x->data[y*dst_stride+x] = src_flow_x->data[j*src_stride+i]*scale_x;
	      dst_flow_y->data[y*dst_stride+x] = src_flow_y->data[j*src_stride+i]*scale_y;
	    }
	}
    }
}
