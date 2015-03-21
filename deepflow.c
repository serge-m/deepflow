#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "image.h"
#include "opticalflow.h"
#include "io.h"

void usage()
{
  printf("usage:\n");
  printf("./deepflow image1 image2 outputfile [options] \n");
  printf("Compute the flow between two images and store it into a file\n");
  printf("Images must be in PPM or JPG format");
  printf("\n");
  printf("options:\n"); 
  printf("    -h, --help                                                print this message\n");
  printf("    -a, -alpha              <float>(12.0)                    weight of smoothness terms\n");
  printf("    -b, -beta               <float>(300.0)                   weight of descriptor matching\n");
  printf("    -g, -gamma              <float>(3.0)                     weight of gradient constancy assumption\n");
  printf("    -d, -delta              <float>(2.0)                     weight of color constancy assumption\n");
  printf("    -s, -sigma              <float>(0.8)                     standard deviation of Gaussian presmoothing kernel\n");
  printf("    -e, -eta                <float>(0.95)                    ratio factor for coarse-to-fine scheme\n");
  printf("    -minsize                <interger>(25)                   size of the coarsest level\n");
  printf("    -inner                  <integer>(5)                     number of inner fixed point iterations\n");
  printf("    -iter                   <integer>(25)                    number of iterations for the solver\n");
  printf("    -soromega               <float>(1.6)                     omega parameter of the sor method\n");
  printf("    -bk                     <float>(0.0)                     use decreasing beta i.e. beta(k) = beta*( k / kmax )^betak, if 0, a last iteration is done with beta=0\n");
  printf("\n");
  printf("    -match                                                   '-match filename' reads matches from a file and '-match' from stdin. Matches must be given with the format <int> <int> <int> <int> <float> x1 y1 x2 y2 score at a resolution of 512x256\n");
  printf("    -matchf                                                  same as -match with images at original resolution\n");
  printf("\n");
  printf("    -sintel                                                  set the parameters to the one used in the ICCV paper for MPI-Sintel dataset\n");
  printf("    -middlebury                                              set the parameters to the one used in the ICCV paper for middlebury dataset\n");
  printf("    -kitti                                                   set the parameters to the one used in the ICCV paper for KITTI dataset\n");
  printf("\n");
}

void require_argument(const char *arg)
{
  fprintf(stderr, "Require an argument after option %s\n", arg);
  exit(1);
}

int main(int argc, char ** argv)
{
  image_t *match_x = NULL, *match_y = NULL, *match_z = NULL;
  
  // load images
  if(argc < 4)
    {
      fprintf(stderr,"Wrong command, require at least 3 arguments.\n\n");
      usage();
      exit(1);
    }
  color_image_t *im1 = color_image_load(argv[1]), *im2 = color_image_load(argv[2]);
  if(im1->width != im2->width || im1->height != im2->height)
    {
      fprintf(stderr,"Image dimensions does not match\n");
      exit(1);
    }
  
  // set params to default
  optical_flow_params_t* params = (optical_flow_params_t*) malloc(sizeof(optical_flow_params_t));
  if(!params)
    {
      fprintf(stderr,"error deepflow(): not enough memory\n");
      exit(1);
    }
  optical_flow_params_default(params);
  
  // parse options
  int current_arg = 4;
  while(1)
    {	
      if( current_arg >= argc)
	{
	  break;
	}
      if(!strcmp(argv[current_arg],"-h") || !strcmp(argv[current_arg],"--help") )
	{
	  usage();
	  exit(1);
	}
      else if(!strcmp(argv[current_arg],"-a") || !strcmp(argv[current_arg],"-alpha") )
	{
	  current_arg++;
	  if(current_arg >= argc) require_argument("alpha");
	  float alpha = atof(argv[current_arg++]);
	  if(alpha<0)
	    {
	      fprintf(stderr,"Alpha argument cannot be negative\n");
	      exit(1);
	    }
	  params->alpha = alpha;
	}
      else if(!strcmp(argv[current_arg],"-b") || !strcmp(argv[current_arg],"-beta") )
	{
	  current_arg++;
	  if(current_arg >= argc) require_argument("beta");
	  float beta = atof(argv[current_arg++]);
	  if(beta<0)
	    {
	      fprintf(stderr,"Beta argument cannot be negative\n");
	      exit(1);
	    }
	  params->beta = beta;
	}
      else if(!strcmp(argv[current_arg],"-g") || !strcmp(argv[current_arg],"-gamma") )
	{
	  current_arg++;
	  if(current_arg >= argc) require_argument("gamma");
	  float gamma = atof(argv[current_arg++]);
	  if(gamma<0)
	    {
	      fprintf(stderr,"Gamma argument cannot be negative\n");
	      exit(1);
	    }
	  params->gamma = gamma;
	}
      else if(!strcmp(argv[current_arg],"-d") || !strcmp(argv[current_arg],"-delta") )
	{
	  current_arg++;
	  if(current_arg >= argc) require_argument("delta");
	  float delta = atof(argv[current_arg++]);
	  if(delta<0)
	    {
	      fprintf(stderr,"Delta argument cannot be negative\n");
	      exit(1);
	    }
	  params->delta = delta;
	}
      else if(!strcmp(argv[current_arg],"-s") || !strcmp(argv[current_arg],"-sigma") )
	{
	  current_arg++;
	  if(current_arg >= argc) require_argument("sigma");
	  float sigma = atof(argv[current_arg++]);
	  if(sigma<0)
	    {
	      fprintf(stderr,"Sigma argument is negative\n");
	      exit(1);
	    }
	  params->sigma = sigma;
	} 
      else if(!strcmp(argv[current_arg],"-bk"))
	{
	  current_arg++;
	  if(current_arg >= argc) require_argument("bk");
	  float betak = atof(argv[current_arg++]);
	  if(betak<0.0f)
	    {
	      fprintf(stderr,"Bk argument must be positive\n");
	      exit(1);
	    }
	  params->bk = betak;
	}
      else if(!strcmp(argv[current_arg],"-e") || !strcmp(argv[current_arg],"-eta") )
	{
	  current_arg++;
	  if(current_arg >= argc) require_argument("eta");
	  float eta = atof(argv[current_arg++]);
	  if(eta<0.25 || eta>0.98)
	    {
	      fprintf(stderr,"Eta argument has to be between 0.25 and 0.98\n");
	      exit(1);
	    }
	  params->eta = eta;
	} 
      else if( !strcmp(argv[current_arg],"-minsize") )
	{
	  current_arg++;
	  if(current_arg >= argc) require_argument("minsize");
	  int minsize = atoi(argv[current_arg++]);
	  if(minsize < 10)
	    {
	      fprintf(stderr,"Minsize argument has to be higher than 10\n");
	      exit(1);
	    }
	  params->min_size = minsize;
	} 
      else if(!strcmp(argv[current_arg],"-inner") )
	{
	  current_arg++;
	  if(current_arg >= argc) require_argument("inner");
	  int inner = atoi(argv[current_arg++]);
	  if(inner<=0)
	    {
	      fprintf(stderr,"Inner argument must be strictly positive\n");
	      exit(1);
	    }
	  params->n_inner_iteration = inner;
	} 
      else if(!strcmp(argv[current_arg],"-iter") )
	{
	  current_arg++;
	  if(current_arg >= argc) require_argument("iter");
	  int iter = atoi(argv[current_arg++]);
	  if(iter<=0)
	    {
	      fprintf(stderr,"Iter argument must be strictly positive\n");
	      exit(1);
	    }
	  params->n_solver_iteration = iter;
	}
      else if( !strcmp(argv[current_arg],"-match") || !strcmp(argv[current_arg],"-matchf"))
	{
	  int wm = im1->width, hm = im1->height;
	  if( !strcmp(argv[current_arg++],"-match") )
	    {
	      wm = 512;
	      hm = 256;
	    }
	  image_delete(match_x); image_delete(match_y); image_delete(match_z);
	  match_x = image_new(wm, hm); match_y = image_new(wm, hm); match_z = image_new(wm, hm); 
	  image_erase(match_x); image_erase(match_y); image_erase(match_z);
	  FILE *fid = stdin;

	  if( current_arg<argc && argv[current_arg][0] != '-')
	    {
	      fid = fopen(argv[current_arg++], "r");
	      if(fid==NULL)
		{
		  fprintf(stderr, "Cannot read matches from file %s", argv[current_arg-1]);
		  exit(1);
		}
	    }
	  int x1, x2, y1, y2;
	  float score;
	  while(!feof(fid) && fscanf(fid, "%d %d %d %d %f\n", &x1, &y1, &x2, &y2, &score)==5)
	    {
	      if( x1<0 || y1<0 || x2<0 || y2<0 || x1>=wm || y1>=hm || x2>=wm || y2>=hm)
		{
		  fprintf(stderr, "Error while reading matches %d %d -> %d %d, out of bounds\n", x1, y1, x2, y2);
		  exit(1);
		}
	      match_x->data[ y1*match_x->stride+x1 ] = (float) (x2-x1);
	      match_y->data[ y1*match_x->stride+x1 ] = (float) (y2-y1);
	      match_z->data[ y1*match_x->stride+x1 ] = score;
	    }
	}
      else if ( !strcmp(argv[current_arg],"-sintel") )
	{
		current_arg++;
	  optical_flow_params_sintel(params);
	}
      else if ( !strcmp(argv[current_arg],"-middlebury") )
	{
		current_arg++;
	  optical_flow_params_middlebury(params);
	}
      else if ( !strcmp(argv[current_arg],"-kitti") )
	{
		current_arg++;
	  optical_flow_params_kitti(params);
	}	
      else
	{
	  if(argv[current_arg][0] == '-')
	    {
	      fprintf(stderr,"Unknow options %s\n",argv[current_arg]);
	    }
	  else
	    {
	      fprintf(stderr,"Error while reading options, %s\n",argv[current_arg]);
	    }
	  exit(1);
	}
    }
	
  image_t *wx = image_new(im1->width,im1->height), *wy = image_new(im1->width,im1->height);
  optical_flow(wx, wy, im1, im2, params, match_x, match_y, match_z);
  writeFlowFile(argv[3], wx, wy);
  image_delete(wx);
  image_delete(wy);
  image_delete(match_x); image_delete(match_y); image_delete(match_z);
  color_image_delete(im1); color_image_delete(im2);
  free(params);
  return 0;
  }
