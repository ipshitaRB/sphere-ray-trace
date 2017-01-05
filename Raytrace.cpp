

/*******************************************************************************

  Purpose:	This program renders three spheres and four planes using 
			non-recursive ray tracing.

			The purpose of this program is to demonstrate
			1. shadow, reflection and refraction using ray tracing 
			2. three types of textures- woodgrain , checkerboard and bump map .

			Only the following have been implemented :

			1. Bump Mapping
			2. Wood Grain 
			3. Recursive Reflection '
			4. Shadow 
			5. Checker Board


 Author : Ipshita Roy Burman 


*******************************************************************************/
#include "stdafx.h"
#include <stdlib.h>
#include <cstdlib>
#include <stdio.h>
#include <math.h>
#include <gl/glut.h>

#define DEBUG   0

/*  Define some constants.  */

#define	PI	3.1415926

/* Max image size allowed. */

#define MAX_SIZE 512

/*  Define some structures.  */

#define MAX_DEPTH 4



struct	points	{
		   float   x, y, z;
};

typedef struct	rgb_struct	{
		   float   r, g, b;
} rgb;

//for bump mapping
#define TABLE_SIZE_1 65

//Noise Table
float noise_tabl[TABLE_SIZE_1][TABLE_SIZE_1];

/*  Viewing parameters.  */

struct	points	from, at, up;
float	VXR, VXL, VYB, VYT;
float	ax, ay, az, bx, by, bz, cx, cy, cz;
float	viewangle, angle, tanv2;
float	xinterval, yinterval;

/*  Illumination parameters.  */

float	lsx, lsy, lsz;
rgb	il, ia;
rgb	ka1, kd1, ks1;
rgb	ka2, kd2, ks2;
rgb	ka3, kd3, ks3;
rgb	ka4, kd4, ks4;
rgb	ka5, kd5, ks5;
rgb	ka6, kd6, ks6;
rgb	tka1, tkd1, tka2, tkd2, tka3, tkd3;
rgb	tka4, tkd4, tka5, tkd5, tka6, tkd6;

int	phong1, phong2, phong3, phong4, phong5, phong6;

/*  Image parameters.  */

int		xmax_pixel, ymax_pixel;

/* Image buffer.  A more efficient approach is to use one single array (texture_RGB)
   rather than using the three rays.  */

float *texture_R;
float *texture_G;
float *texture_B;

/*  Object parameters.  */

float	xc1, yc1, zc1, r1, xc2, yc2, zc2, r2, xc3, yc3, zc3, r3, xc4, yc4, zc4, r4;
float	a1, b1, c1, d1, a2, b2, c2, d2, a3, b3, c3, d3;
float	px1_min, px1_max, py1_min, py1_max, pz1_min, pz1_max;
float	px2_min, px2_max, py2_min, py2_max, pz2_min, pz2_max;
float	px3_min, px3_max, py3_min, py3_max, pz3_min, pz3_max;

rgb Reflection(int depth, float fromx, float fromy, float fromz, float tox, float toy, float toz, float *rr, float *rg, float *rb, int object_id);

/* Image output file. */

FILE	*outpfile;

/*******************************************************************************

  Title:	Read_Information

  Purpose:	This function reads in the information about the objects (three spheres
            and three planes) to be rendered.  The information is assumed to be 
			stored in an ASCII text file called "ray.dat".	

*******************************************************************************/

int Read_Information()
{
	char	str[132];
	FILE	*inpfile;

	if ( (inpfile = fopen("ray_IRB.dat","r"))==NULL) {
		printf("ERROR: Could not open ray.dat for read!\n");
		return(0);
	}

/*  Read in viewing information.  */

	fscanf(inpfile,"%s %f %f %f\n",str, &from.x, &from.y, &from.z);
	fscanf(inpfile,"%s %f %f %f\n",str, &at.x, &at.y, &at.z);
	fscanf(inpfile,"%s %f %f %f\n",str, &up.x, &up.y, &up.z);
	fscanf(inpfile,"%s %f\n",str, &viewangle);
	angle = viewangle * PI/180.0;
	tanv2 = tan(angle/2.0);
	printf("tanv2 = %f\n", tanv2);

/* Read in the viewport information. */

	fscanf(inpfile,"%s %f %f\n", str, &VXL, &VXR);
	fscanf(inpfile,"%s %f %f\n", str, &VYB, &VYT);

/* Read in the light vector (L).  Note: To specify a light source position, the light 
   vector will need to be computed in Compute_Color().   */

	fscanf(inpfile,"%s %f %f %f\n",str, &lsx, &lsy, &lsz);

	printf("From: %f %f %f\n", from.x, from.y, from.z);
	printf("At: %f %f %f\n", at.x, at.y, at.z);
	printf("Up: %f %f %f  View Angle: %f \n", up.x, up.y, up.z, viewangle);

/*  Read in spheres' information.  */

	fscanf(inpfile,"%s %f %f %f %f\n",str, &xc1, &yc1, &zc1, &r1);
	fscanf(inpfile,"%s %f %f %f %f\n",str, &xc2, &yc2, &zc2, &r2);
	fscanf(inpfile,"%s %f %f %f %f\n",str, &xc3, &yc3, &zc3, &r3);
	fscanf(inpfile,"%s %f %f %f %f\n",str, &xc4, &yc4, &zc4, &r4);

/*  Read in checker boards' (planes') information.  */

	fscanf(inpfile,"%s %f %f %f %f\n",str, &a1, &b1, &c1, &d1);
	fscanf(inpfile,"%s %f %f %f %f\n",str, &a2, &b2, &c2, &d2);
	fscanf(inpfile,"%s %f %f %f %f\n",str, &a3, &b3, &c3, &d3);

/*  Read in the boundaries of the planes.  */

	fscanf(inpfile,"%s %f %f %f %f %f %f\n", str, &px1_min, &py1_min,
		&pz1_min, &px1_max, &py1_max, &pz1_max);
	fscanf(inpfile,"%s %f %f %f %f %f %f\n", str, &px2_min, &py2_min,
		&pz2_min, &px2_max, &py2_max, &pz2_max);
	fscanf(inpfile,"%s %f %f %f %f %f %f\n", str, &px3_min, &py3_min,
		&pz3_min, &px3_max, &py3_max, &pz3_max);

/*  Read in the image size.  */

	fscanf(inpfile,"%s %d %d\n",str, &xmax_pixel, &ymax_pixel);

	/* Make sure the image size does not exceed MAX_SIZE x MAX_SIZE.  */

	if (xmax_pixel > MAX_SIZE || ymax_pixel > MAX_SIZE) {
		printf("Error: Exceeded max image size %d x %d\n", xmax_pixel, ymax_pixel);
		printf("Reset to max image size: %d x %d\n", MAX_SIZE, MAX_SIZE);
		xmax_pixel = MAX_SIZE-1;
		ymax_pixel = MAX_SIZE - 1;
	}

	fclose(inpfile);

/*  Open an output file to store the intensity values of the output image.  */

	if ((outpfile = fopen("image_IRB.out","wb")) == NULL) {
		printf("ERROR:  cannot open image.out for write.\n");
		return(0);
	}

/*  Allocate memory for the image buffer.  */

	texture_R = new float [xmax_pixel * ymax_pixel ];
	texture_G = new float [xmax_pixel * ymax_pixel ];
	texture_B = new float [xmax_pixel * ymax_pixel ];
	printf("image_buf allocated.  Image size %d x %d\n", xmax_pixel, ymax_pixel);

	return(1);
}


/*******************************************************************************

  Title:	Normalize

  Purpose:	This function normalizes the given vector.

*******************************************************************************/

void Normalize(float *x,float *y,float *z)
{
	float	norm;

	norm = sqrt( *x * *x + *y * *y + *z * *z );
	if (norm != 0.0) {
		*x = *x / norm;
		*y = *y / norm;
		*z = *z / norm;
	}
}

/*******************************************************************************

  Title:	Power

  Purpose:	This function computes the power of the given base and
		    exponent.

*******************************************************************************/

float 	Power(float base,int exp)
{
	int	i;
	float	value;

	value = 1.0;
	for (i=1; i<=exp; i++)
		value *= base;

	return( value );
}


/*******************************************************************************

  Title:	Compute_M

  Purpose:	This function computes the transformation matrix to be used
		    in the perspective viewing model.

*******************************************************************************/

void Compute_M()
{

/*  Compute the line-of-sight vector, c.  */

	cx = at.x - from.x;
	cy = at.y - from.y;
	cz = at.z - from.z;
	Normalize(&cx, &cy, &cz);

/*  Compute the cross product of vector c and the up vector.  */

	ax = cy*up.z - up.y*cz;
	ay = up.x*cz - cx*up.z;
	az = cx*up.y - up.x*cy;
	Normalize(&ax, &ay, &az);

/*  Compute the cross product of vector a and c.  */

	bx = ay*cz - cy*az;
	by = cx*az - ax*cz;
	bz = ax*cy - cx*ay;
}

/*******************************************************************************

  Title:	Setup_Parameters

  Purpose:	This function sets up the necessary parameters for 
		    performing the ray trace.  It first computes the 
		    transformation matrix for the perspective viewing model, then
		    sets up the default illumination parameters.

*******************************************************************************/

void Setup_Parameters()
{

/*  Compute the transformation matrix for converting world coordinates to eye
    coordinates.  */

	Compute_M(); 

/*  Normalized the given directional light vector. 

	Note:  DO NOT normalize this vector, if the position of the light source is given.
	       The light vector (L) would need to computed for each intersction point,
		   then normalized in Compute_Color().   */

	Normalize(&lsx, &lsy, &lsz);
	printf("light %f %f %f\n", lsx, lsy, lsz);

/*  Set up the conversion factors for converting from pixel coordinates to 
    view port coordinates.  */

	xinterval = (VXR - VXL) / xmax_pixel;
	yinterval = (VYT - VYB) / ymax_pixel;

/*  Set up default illumination (Phong lighting) parameters.  */

	il.r = 1.0;	il.g = 1.0;	il.b = 1.0;
	ia.r = 1.0;	ia.g = 1.0;	ia.b = 1.0;

/*  Phone lighting parameters for the three spheres.  */

	ka1.r = 0.3;	ka1.g = 0.0;	ka1.b = 0.0;
	kd1.r = 0.7;	kd1.g = 0.0;	kd1.b = 0.0;
	ks1.r = 1.0;	ks1.g = 1.0;	ks1.b = 1.0;
	tka1.r = 0.2;	tka1.g = 0.0;	tka1.b = 0.0;
	tkd1.r = 0.2;	tkd1.g = 0.0;	tkd1.b = 0.0;
	phong1 = 10;

	ka2.r = 0.0;	ka2.g = 0.3;	ka2.b = 0.0;
	kd2.r = 0.0;	kd2.g = 0.7;	kd2.b = 0.0;
	ks2.r = 1.0;	ks2.g = 1.0;	ks2.b = 1.0;
	tka2.r = 0.0;	tka2.g = 0.6;	tka2.b = 0.0;
	tkd2.r = 0.0;	tkd2.g = 0.6;	tkd2.b = 0.0;
	phong2 = 90;

	ka3.r = 0.0;	ka3.g = 0.0;	ka3.b = 0.3;
	kd3.r = 0.0;	kd3.g = 0.0;	kd3.b = 0.7;
	ks3.r = 1.0;	ks3.g = 1.0;	ks3.b = 1.0;
	tka3.r = 0.0;	tka3.g = 0.0;	tka3.b = 0.2;
	tkd3.r = 0.0;	tkd3.g = 0.0;	tkd3.b = 0.2;
	phong3 = 20;

/*  Phone lighting parameters for the three planes (not shown).  */

	ka4.r = 0.1;	ka4.g = 0.1;	ka4.b = 0.0;
	kd4.r = 0.7;	kd4.g = 0.7;	kd4.b = 0.0;
	ks4.r = 1.0;	ks4.g = 1.0;	ks4.b = 1.0;
	tka4.r = 0.1;	tka4.g = 0.0;	tka4.b = 0.0;
	tkd4.r = 0.7;	tkd4.g = 0.0;	tkd4.b = 0.0;
	phong4 = 120;

	ka5.r = 0.1;	ka5.g = 0.0;	ka5.b = 0.1;
	kd5.r = 0.7;	kd5.g = 0.0;	kd5.b = 0.7;
	ks5.r = 1.0;	ks5.g = 1.0;	ks5.b = 1.0;
	tka5.r = 0.0;	tka5.g = 0.0;	tka5.b = 0.1;
	tkd5.r = 0.0;	tkd5.g = 0.0;	tkd5.b = 0.7;
	phong5 = 120;

	ka6.r = 0.1;	ka6.g = 0.1;	ka6.b = 0.1;
	kd6.r = 0.8;	kd6.g = 0.6;	kd6.b = 0.2;
	ks6.r = 1.0;	ks6.g = 1.0;	ks6.b = 1.0;
	tka6.r = 0.0;	tka6.g = 0.1;	tka6.b = 0.0;
	tkd6.r = 0.0;	tkd6.g = 0.7;	tkd6.b = 0.0;
	phong6 = 120;
}




/*******************************************************************************

  Title:	Check_Sphere

  Purpose:	This function determines if the give ray intercepts the given
		    sphere.

*******************************************************************************/

void Check_Sphere(float px,float py,float pz,float dx,float dy,float dz,float xc,float yc,float zc,float r,float *t1, float *t2)
{
	float	a, b, c, xdiff, ydiff, zdiff, discr;

	xdiff = px-xc;
	ydiff = py-yc;
	zdiff = pz-zc;
	a = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff - r*r;
	b = 2.0*( dx*xdiff + dy*ydiff + dz*zdiff );
	c = dx*dx + dy*dy + dz*dz;

/*  Check if there are any intersections.  */

	discr = b*b - 4.0*a*c;
	if (discr < 0.0) {
		*t1 = -1.0;
		*t2 = -1.0;
	}
	else if (discr == 0.0) {
		*t1 = -b / (2.0*c);
		*t2 = -1.0;
	}
	else {
		discr = sqrt(discr);
		*t1 = (-b + discr) / (2.0*c);
		*t2 = (-b - discr) / (2.0*c);
	}
}


/*******************************************************************************

  Title:	Check_Plane

  Purpose:	This function checks if the given ray intercepts the given
		    plane.

*******************************************************************************/

void Check_Plane(float px,float py,float pz,float dx,float dy,float dz,float a,float b,float c,float d,float *t1)
{
	*t1 = (-a*px - b*py - c*pz - d) / (a*dx + b*dy + c*dz);
}


/*******************************************************************************

  Title:	Compute_Intersection

  Purpose:	This function computes the intersection of ray with an
		    object.  The intersection point is given by a parametric value
		    t, where ray = p + d*t, d = the direction of the ray, and p is
		    the starting point of the ray.

*******************************************************************************/

void Compute_Intersection(float px,float py,float pz,float dx,float dy, float dz,float t,float *newx,float *newy,float *newz)
{
	*newx = px + t*dx;
	*newy = py + t*dy;
	*newz = pz + t*dz;
}

/*******************************************************************************

Title:	Check Shadow

Purpose:	This function checks for intersection to see if a shadow is generated
and returns 1 or 0 (shadow flag value)

*******************************************************************************/


int Check_Shadow(float x, float y, float z, int object_id)
{
	float t1, t2, t_min, ipx, ipy, ipz;
	int i;

	t_min = 999.0;



	/*  Check if the current ray intercepts sphere #1.  */

	if (object_id != 1){
		Check_Sphere(x, y, z, lsx, lsy, lsz, xc1,
			yc1, zc1, r1, &t1, &t2);

		if (t1 >= 0.0 && t1 < t_min) {
			return 1;
		}

		if (t2 >= 0.0 && t2 < t_min) {
			return 1;
		}
	}
	/*  Check if the current ray intercepts sphere #2.  */
	if (object_id != 2){
		Check_Sphere(x, y, z, lsx, lsy, lsz, xc2,
			yc2, zc2, r2, &t1, &t2);

		if (t1 >= 0.0 && t1 < t_min) {
			return 1;
		}

		if (t2 >= 0.0 && t2 < t_min) {
			return 1;
		}
	}
	/*  Check if the current ray intercepts sphere #3.  */
	if (object_id != 3){
		Check_Sphere(x, y, z, lsx, lsy, lsz, xc3,
			yc3, zc3, r3, &t1, &t2);

		if (t1 >= 0.0 && t1 < t_min) {
			return 1;
		}

		if (t2 >= 0.0 && t2 < t_min) {
			return 1;
		}
	}
	/*  Check if the current ray intercepts plane #1.  */
	if (object_id != 4){
		Check_Plane(x, y, z, lsx, lsy, lsz, a1, b1, c1, d1, &t1);

		if (t1 >= 0.0 && t1 < t_min) {
			/*  Check if the intersection point is inside the min/max values. */

			Compute_Intersection(x, y, z, lsx, lsy, lsz, t1, &ipx, &ipy, &ipz);

			if (ipx >= px1_min && ipx <= px1_max &&
				ipy >= py1_min  && ipy <= py1_max &&
				ipz >= pz1_min && ipz <= pz1_max) {
				return 1;
			}
		}
	}
	/*  Check if the current ray intercepts plane #2.  */
	if (object_id != 5){
		Check_Plane(x, y, z, lsx, lsy, lsz, a2, b2, c2, d2, &t1);

		if (t1 >= 0.0 && t1 < t_min) {
			/*  Check if the intersection point is inside the min/max values. */

			Compute_Intersection(x, y, z, lsx, lsy, lsz, t1, &ipx, &ipy, &ipz);

			if (ipx >= px2_min && ipx <= px2_max &&
				ipy >= py2_min  && ipy <= py2_max &&
				ipz >= pz2_min && ipz <= pz2_max) {
				return 1;
			}
		}
	}
	/*  Check if the current ray intercepts plane #3.  */
	if (object_id != 6){
		Check_Plane(x, y, z, lsx, lsy, lsz, a3, b3, c3, d3, &t1);

		if (t1 >= 0.0 && t1 < t_min) {
			/*  Check if the intersection point is inside the min/max values. */

			Compute_Intersection(x, y, z, lsx, lsy, lsz, t1, &ipx, &ipy, &ipz);

			if (ipx >= px3_min && ipx <= px3_max &&
				ipy >= py3_min  && ipy <= py3_max &&
				ipz >= pz3_min && ipz <= pz3_max) {
				return 1;
			}
		}
	}
	return 0;
}

/*******************************************************************************

  Title:	Compute_Color

  Purpose:	This function computes the intensity of the color for the
		    given location based on the Phong lighting model.

*******************************************************************************/

rgb Compute_Color(int shadow_flag, float ipx,float ipy,float  ipz,float  nx,float  ny,float  nz,
				   rgb ia,rgb ka,rgb kd, rgb ks,int n,float *r,float *g, float *b , int obj_num, int depth, float fromx, float fromy, float fromz, float tox, float toy, float toz)
{
	float	vx, vy, vz, rx, ry, rz;
	float	ndotl, vdotr, cosalphapower;
	rgb color;
/*  Compute the view vector.  */

	vx = from.x - ipx;
	vy = from.y - ipy;
	vz = from.z - ipz;
	Normalize(&vx, &vy, &vz);

/*  Compute the R (reflection) vector.  */

	ndotl = nx*lsx + ny*lsy + nz*lsz;
	rx = 2.0*ndotl*nx - lsx;
	ry = 2.0*ndotl*ny - lsy;
	rz = 2.0*ndotl*nz - lsz;

/* Compute the V (view) vector. */

	vdotr = vx*rx + vy*ry + vz*rz;

/* Compute Ia * Ka.  */

	*r = ia.r * ka.r;
	*g = ia.g * ka.g;
	*b = ia.b * ka.b;


/* Compute diffuse reflection. */

	if (ndotl >= 0.0 && shadow_flag==0) { 

		/*  diffuse reflection = kd * N dot L * Il  */

		*r = *r + kd.r*ndotl*il.r;
		*g = *g + kd.g*ndotl*il.g;
		*b = *b + kd.b*ndotl*il.b;

		if (vdotr >= 0.0) {

			/*  specular reflection = ks * cos(alpha)**K^n * Il */

			cosalphapower = Power(vdotr, n);

			*r = *r + ks.r*cosalphapower*il.r;
			*g = *g + ks.g*cosalphapower*il.g;
			*b = *b + ks.b*cosalphapower*il.b;
		}
	}
	else
		if(shadow_flag == 1){

			//only provide the ambient light
			*r = ia.r * ka.r;
			*g = ia.g * ka.g;
			*b = ia.b * ka.b;
		}

/*  Make sure that the color is within range.  */

	if (*r > 1.0) *r = 1.0;
	if (*g > 1.0) *g = 1.0;
	if (*b > 1.0) *b = 1.0;

	if(depth > 0 && depth < MAX_DEPTH && obj_num == 1){
		color =  Reflection(depth++,fromx, fromy, fromz,tox, toy, toz, r, g, b, obj_num); 
		*r += ks.r*color.r;
		*g += ks.g*color.g;
		*b += ks.b*color.b;
	}

	color.r = *r;
	color.g = *g;
	color.b = *b;

	return color;
}

/*******************************************************************************

Title:	Calculate Noise

Purpose:	This function is used to create noise for bump map.

*******************************************************************************/

float calc_noise(float iu, float iv, int direction)
{
	int i, j, x, y, left, right;
	float noise, u, v, w, u1, v1, w1;
	float val;

	//Filling noise table with random noise
	for (i = 0; i < 65; i++)
		for (j = 0; j < 65; j++)
		{

			val = rand() % 255;
			noise_tabl[i][j] = val / 256;
		}

	i = (int)iu;
	j = (int)iv;
	x = i%TABLE_SIZE_1;
	y = j%TABLE_SIZE_1;

	if (direction == 1)
	{
		if (x < 0)
			left = 0;
		else
			left = x - 1;

		if (x >= TABLE_SIZE_1)
			right = TABLE_SIZE_1;
		else
			right = x + 1;
		//Getting the value of noise
		noise = (noise_tabl[right][y] - noise_tabl[left][y]) / 2.0;
	}

	else {
		if (y <= 0)
			left = 0;
		else
			left = y - 1;

		if (y >= TABLE_SIZE_1)
			right = TABLE_SIZE_1;
		else
			right = y + 1;
		//Getting the value of random noise
		noise = (noise_tabl[x][right] - noise_tabl[x][left]) / 2.0;
	}
	return(noise);

}

/*******************************************************************************

Title:	Bump_Map

Purpose:	Perturb the given normal vector based on bump mapping

*******************************************************************************/
void Bump_Map(float x, float y, float z, float xc, float yc, float zc, float r, float *nx, float *ny, float *nz)
{
	//Variables
	float xp, yp, zp, iu, iv, xu, yu, zu, xv, yv, zv;
	float fu, fv, a, dx, dy, dz, u, v, nnx, nny, nnz;

	//translate to origin
	xp = (x - xc) / r;
	yp = (y - yc) / r;
	zp = (z - zc) / r;

	//convert to uv coordinates
	u = asin(zp);
	v = atan2(yp, xp);

	//convert to integer (iu, iv ) coordinates
	iu = (u + PI) / (2 * PI) * (TABLE_SIZE_1);
	iv = (v + (PI / 2)) / PI * (TABLE_SIZE_1);

	//get the partials
	xu = -r*cos(v)*sin(u);
	xv = -r*sin(v)*cos(u);
	yu = r*cos(v)*cos(u);
	yv = -r*sin(v)*sin(u);
	zu = 0.0;
	zv = r*cos(v);

	//calculate the perturbations
	fu = calc_noise(iu, iv, 1);
	fv = calc_noise(iu, iv, 2);

	int noise_method = 0;
	//compute D
	if(noise_method == 1){
	dx = fu*xu + fv*xv;
	dy = fu*yu + fv*yv;
	dz = fu*zu + fv*zv;
	}
	else
	{
	nnx = *nx;
	nny = *ny;
	nnz = *nz;
	Normalize(&nnx, &nny, &nnz);
	


	//computing cross product of Pu x N
	dx = fv* (yu*nnz - nny*zu);
	dy = fv* (nnx*zu - xu*nnz);
	dz = fv* (xu*nny - nnx*yu);

	//computing cross product of Pv x N
	dx = fu* (yv*nnz - nny*zv);
	dy = fu* (nnx*zv - xv*nnz);
	dz = fu* (xv*nny - nnx*yv);
	}

	//normalize and scale D
	Normalize(&dx, &dy, &dz);

	a = sqrt(fu*fu + fv*fv);
	dx *= a;
	dy *= a;
	dz *= a;

	//Finding new normals after displacement
	*nx += dx;


	

	*ny += dy;
	*nz += dz;
}

/*******************************************************************************

Title:	WoodGrain

Purpose:	This function is used to create a woodgrain texture.

*******************************************************************************/

rgb woodgrain(float u, float v, float w, float *r, float *g, float *b)
{

	float radius, angle;
	int grain;
	rgb grain_color ;

	//Creating texture
	radius = sqrt((u*u) + (w*w));
	if (w == 0)
		angle = PI / 2;
	else
		angle = atan2(u, w);

	radius = radius + 2 * sin(20 * angle + v / 5);
	grain = (int)radius % 6;

	if (grain < 4)
	{

		//Setting rgb values
		*r = 0.6;
		*g = 0.4;
		*b = 0.12;
	}

	else
	{
		//Setting rgb values
		*r = 0.4;
		*g = 0.2;
		*b = 0.0;
	}

	grain_color.r = *r;
	grain_color.g = *g;
	grain_color.b = *b;
	
	return grain_color;
					
}

rgb Reflection(int depth, float fromx, float fromy, float fromz, float tox, float toy, float toz, float *rr, float *rg, float *rb, int object_id){
	int	    xp, yp, obj_num, shadow_flag;
	int	    texture, buf_ptr;
	int     num_image = 1;
	float	xv, yv, dx, dy, dz, nx, ny, nz;
	float	t_min, t1, t2, ipx, ipy, ipz;
	float	r, g, b;
	float   u, v;
	rgb color;
	float   interx, intery, interz;
	float   reflectedx, reflectedy, reflectedz;
	float   dot;
	t_min = 999.0;
			obj_num = 0;
			texture = 0;
	Check_Sphere(fromx, fromy, fromz, tox, toy, toz, xc1, 
				yc1, zc1, r1, &t1, &t2);

			if (t1>=0.0) {
				t_min = t1;
				obj_num = 1;
				Compute_Intersection(fromx, fromy, fromz, 
					tox, toy, toz, t1, &ipx, &ipy, &ipz);
			}

			if (t2>=0.0 && t2<t_min) {
				t_min = t2;
				obj_num = 1;
				Compute_Intersection(fromx, fromy, fromz, 
					tox, toy, toz, t2, &ipx, &ipy, &ipz);
			}

/*  Check if the current ray intercepts sphere #2.  */

			Check_Sphere(fromx, fromy, fromz, tox, toy, toz, xc2, 
				yc2, zc2, r2, &t1, &t2);

			if (t1>=0.0 && t1<t_min) {
				t_min = t1;
				obj_num = 2;
				Compute_Intersection(fromx, fromy, fromz, 
					tox, toy, toz, t1, &ipx, &ipy, &ipz);
			}

			if (t2>=0.0 && t2<t_min) {
				t_min = t2;
				obj_num = 2;
				Compute_Intersection(fromx, fromy, fromz, 
					tox, toy, toz, t2, &ipx, &ipy, &ipz);
			}

/*  Check if the current ray intercepts sphere #3.  */

			Check_Sphere(fromx, fromy, fromz, tox, toy, toz, xc3, 
				yc3, zc3, r3, &t1, &t2);

			if (t1>=0.0 && t1<t_min) {
				t_min = t1;
				obj_num = 3;
				Compute_Intersection(fromx, fromy, fromz, 
					tox, toy, toz, t1, &ipx, &ipy, &ipz);
			}

			if (t2>=0.0 && t2<t_min) {
				t_min = t2;
				obj_num = 3;
				Compute_Intersection(fromx, fromy, fromz, 
					tox, toy, toz, t2, &ipx, &ipy, &ipz);
			}

			/*  Check if the current ray intercepts plane #1.  */

			Check_Plane(fromx, fromy, fromz, tox, toy, toz, a1, b1, c1, d1, &t1);

			if (t1 >= 0.0 && t1<t_min) {
				/*  Check if the intersection point is inside the min/max values. */

				Compute_Intersection(fromx, fromy, fromz, 
					tox, toy, toz, t1, &ipx, &ipy, &ipz);

				if (ipx >= px1_min && ipx <= px1_max && 
					ipy >= py1_min  && ipy <= py1_max && 
					ipz >=  pz1_min && ipz <= pz1_max ) { 

						t_min = t1;
						obj_num = 4;
				}
			}

/*  Check if the current ray intercepts plane #2.  */

			Check_Plane(fromx, fromy, fromz, tox, toy, toz, a2, b2, c2, d2, &t1);

			if (t1 >= 0.0 && t1<t_min) {
				/*  Check if the intersection point is inside the min/max values. */

				Compute_Intersection(fromx, fromy, fromz, 
					tox, toy, toz, t1, &ipx, &ipy, &ipz);

				if (ipx >= px2_min && ipx <= px2_max && 
					ipy >= py2_min  && ipy <= py2_max && 
					ipz >=  pz2_min && ipz <= pz2_max ) { 

						t_min = t1;
						obj_num = 5;
				}
			}

			/*  Check if the current ray intercepts plane #3.  */

			Check_Plane(fromx, fromy, fromz, tox, toy, toz, a3, b3, c3, d3, &t1);

			if (t1 >= 0.0 && t1<t_min) {
				/*  Check if the intersection point is inside the min/max values. */

				Compute_Intersection(fromx, fromy, fromz, 
					tox, toy, toz, t1, &ipx, &ipy, &ipz);

				if (ipx >= px3_min && ipx <= px3_max && 
					ipy >= py3_min  && ipy <= py3_max && 
					ipz >=  pz3_min && ipz <= pz3_max ) { 

						t_min = t1;
						obj_num = 6;
				}
			}
			rgb graincolor;

/*  Compute the intensity to use at the current pixel.  */

			switch (obj_num) {

/*  The current ray does not intersect any of the objects.  */

				case 0 : r = 0.0;
					 g = 0.7;
					 b = 0.8;
					 color.r = r;
					 color.g = g;
					 color.b 
						 = b;
					 break;

/*  The current ray intercept sphere #1.  */

				case 1 : 
					 nx = ipx - xc1;
					 ny = ipy - yc1;
					 nz = ipz - zc1;
					 Normalize(&nx, &ny, &nz);
					 
					 shadow_flag = 0;
					 shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
					 texture = 0;

					 //Finding intersection ray
				interx = tox - ipx;
				intery = toy - ipy;
				interz = toz - ipz;
				//Finding dot product and reflected ray
				dot = nx*interx + ny*intery + nz*interz;
				reflectedx = interx - 2 * (nx*dot);
				reflectedy = intery - 2 * (ny*dot);
				reflectedz = interz - 2 * (nz*dot);

					 color = Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka1, kd1, ks1, phong1, &r, &g, &b, obj_num, 0,ipx, ipy, ipz, reflectedx, reflectedy, reflectedz);
					 break;

/*  The current ray intercepts sphere #2.  */

				case 2 : 
					 nx = ipx - xc2;
					 ny = ipy - yc2;
					 nz = ipz - zc2;
					 Normalize(&nx, &ny, &nz);
					 shadow_flag = 0;
					 shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
					
					 texture = 0;
					 Bump_Map(ipx, ipy, ipz, tox, toy, toz, 1, &nx, &ny, &nz);
				 	 color = Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka2, kd2, ks2, phong2, &r, &g, &b, obj_num, 0 , fromx, fromy, fromz, tox, toy, toz);
					 break;

/*  The current ray intercepts sphere #3.  */

				case 3 : 
					 nx = ipx - xc3;
					 ny = ipy - yc3;
					 nz = ipz - zc3;
					 Normalize(&nx, &ny, &nz);
					 shadow_flag = 0;
					 shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);	 
					// texture = 0;

					
					 graincolor = woodgrain(ipx, ipy, ipz, &r, &g, &b );
					
					

					 //if (texture==1)
						// Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka3, tkd3, ks3, phong3, &r, &g, &b);
					 //else
					 color = Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka3, kd3, ks3, phong3, &graincolor.r, &graincolor.g, &graincolor.b, obj_num, 0 ,fromx, fromy, fromz, tox, toy, toz);
					 break;

/*  The current ray intercepts checker board #1.  */

				case 4 : 
					 nx = a1;
					 ny = b1;
					 nz = c1;
					 shadow_flag = 0;
					shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
					woodgrain(ipx, ipy, ipz, &r, &g, &b);

					 if (ipx < 2.0 || (ipx>=4.0 && ipx<6.0)) {
						if ((ipy>=2.0 && ipy<4.0) || (ipy>=6.0))
								texture = 1;
						else
							texture = 0;
					}
					else {
						if ((ipy<2.0) || (ipy>=4.0 && ipy<6.0)) 
							texture = 1;
						else
							texture = 0;
					}
					if (texture==1)
					color = 	Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka4, tkd4, ks4, phong4, &r, &g, &b, obj_num, 0,fromx, fromy, fromz, tox, toy, toz);
					else
					color = 	Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka4, kd4, ks4, phong4, &r, &g, &b, obj_num, 0,fromx, fromy, fromz, tox, toy, toz);
					 break;

/*  The current ray intercepts checker board #2.  */

				case 5 : 
					 nx = a2;
					 ny = b2;
					 nz = c2;
					 shadow_flag = 0;
					 shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
					
					 if (ipz < 2.0 || (ipz>=4.0 && ipz<6.0)) {
						if ((ipy>=2.0 && ipy<4.0) || (ipy>=6.0))
								texture = 1;
						else
								texture = 0;
					}
					else {
						if ((ipy<2.0) || (ipy>=4.0 && ipy<6.0)) 
							texture = 1;
						else
							texture = 0;
					}
					if (texture == 1)
					color = 	 Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka5, tkd5, ks5, phong5, &r, &g, &b, obj_num, 0,fromx, fromy, fromz, tox, toy, toz);
					else
					color = 	 Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka5, kd5, ks5, phong5, &r, &g, &b, obj_num, 0,fromx, fromy, fromz, tox, toy, toz);
					 break;

/*  The current ray intercepts checker board #3.  */

				case 6 : 
					 nx = a3;
					 ny = b3;
					 nz = c3;

					 shadow_flag = 0;
					 shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
	
					 if (ipx < 2.0 || (ipx>=4.0 && ipx<6.0)) {
						if ((ipz>=2.0 && ipz<4.0) || (ipz>=6.0))
							texture = 1;
						else
							texture = 0;
					 }
					 else {
						if ((ipz<2.0) || (ipz>=4.0 && ipz<6.0)) 
							texture = 1;
						else
							texture = 0;
					 }
					
					 if (texture == 1)
					color = 	Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka6, tkd6, ks6, phong6, &r, &g, &b, obj_num, 0 ,fromx, fromy, fromz, tox, toy, toz);
					 else
					color = 	 Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka6, kd6, ks6, phong6, &r, &g, &b, obj_num, 0,fromx, fromy, fromz, tox, toy, toz);
					 break;
			}

			return color;
}


/*******************************************************************************

  Title:	Ray_Trace	

  Purpose:	This function performs simple ray tracing.  This is a non-recursive
            ray tracing without any reflection and refraction rays. 

*******************************************************************************/

void Ray_Trace()
{
	int	    xp, yp, obj_num, shadow_flag;
	int	    texture, buf_ptr;
	int     num_image = 1;
	float	xv, yv, dx, dy, dz, nx, ny, nz;
	float	t_min, t1, t2, ipx, ipy, ipz;
	float	r, g, b;
	float   u, v;
		t_min = 999.0;
			obj_num = 0;
			texture = 0;
			rgb color;
			float   interx, intery, interz;
	float   reflectedx, reflectedy, reflectedz;
	float   dot;

/*  Generate a ray for each pixel in the desired image.  */

	printf("Ray tracing...\n");
	buf_ptr = 0;
	for (xp=0; xp<xmax_pixel; xp++) {
		u = (float)xp/xmax_pixel;

		for (yp=0; yp<ymax_pixel; yp++) {
			v = (float)yp/ymax_pixel;

/*  Compute the corresponding view port coordinates.  */

			xv = VXL + xp * xinterval;
			yv = VYB + yp * yinterval;

/*  Compute the direction of the current ray from the "From" point to the 
    current position on the image.  */

			dx = ax*xv*tanv2 + bx*yv*tanv2 + cx;
			dy = ay*xv*tanv2 + by*yv*tanv2 + cy;
			dz = az*xv*tanv2 + bz*yv*tanv2 + cz;

/*  Check if the current ray intercepts sphere #1.  */

			t_min = 999.0;
			obj_num = 0;
			texture = 0;

			Check_Sphere(from.x, from.y, from.z, dx, dy, dz, xc1, 
				yc1, zc1, r1, &t1, &t2);

			if (t1>=0.0) {
				t_min = t1;
				obj_num = 1;
				Compute_Intersection(from.x, from.y, from.z, 
					dx, dy, dz, t1, &ipx, &ipy, &ipz);
			}

			if (t2>=0.0 && t2<t_min) {
				t_min = t2;
				obj_num = 1;
				Compute_Intersection(from.x, from.y, from.z, 
					dx, dy, dz, t2, &ipx, &ipy, &ipz);
			}

/*  Check if the current ray intercepts sphere #2.  */

			Check_Sphere(from.x, from.y, from.z, dx, dy, dz, xc2, 
				yc2, zc2, r2, &t1, &t2);

			if (t1>=0.0 && t1<t_min) {
				t_min = t1;
				obj_num = 2;
				Compute_Intersection(from.x, from.y, from.z, 
					dx, dy, dz, t1, &ipx, &ipy, &ipz);
			}

			if (t2>=0.0 && t2<t_min) {
				t_min = t2;
				obj_num = 2;
				Compute_Intersection(from.x, from.y, from.z, 
					dx, dy, dz, t2, &ipx, &ipy, &ipz);
			}

/*  Check if the current ray intercepts sphere #3.  */

			Check_Sphere(from.x, from.y, from.z, dx, dy, dz, xc3, 
				yc3, zc3, r3, &t1, &t2);

			if (t1>=0.0 && t1<t_min) {
				t_min = t1;
				obj_num = 3;
				Compute_Intersection(from.x, from.y, from.z, 
					dx, dy, dz, t1, &ipx, &ipy, &ipz);
			}

			if (t2>=0.0 && t2<t_min) {
				t_min = t2;
				obj_num = 3;
				Compute_Intersection(from.x, from.y, from.z, 
					dx, dy, dz, t2, &ipx, &ipy, &ipz);
			}

			/*  Check if the current ray intercepts plane #1.  */

			Check_Plane(from.x, from.y, from.z, dx, dy, dz, a1, b1, c1, d1, &t1);

			if (t1 >= 0.0 && t1<t_min) {
				/*  Check if the intersection point is inside the min/max values. */

				Compute_Intersection(from.x, from.y, from.z, 
					dx, dy, dz, t1, &ipx, &ipy, &ipz);

				if (ipx >= px1_min && ipx <= px1_max && 
					ipy >= py1_min  && ipy <= py1_max && 
					ipz >=  pz1_min && ipz <= pz1_max ) { 

						t_min = t1;
						obj_num = 4;
				}
			}

/*  Check if the current ray intercepts plane #2.  */

			Check_Plane(from.x, from.y, from.z, dx, dy, dz, a2, b2, c2, d2, &t1);

			if (t1 >= 0.0 && t1<t_min) {
				/*  Check if the intersection point is inside the min/max values. */

				Compute_Intersection(from.x, from.y, from.z, 
					dx, dy, dz, t1, &ipx, &ipy, &ipz);

				if (ipx >= px2_min && ipx <= px2_max && 
					ipy >= py2_min  && ipy <= py2_max && 
					ipz >=  pz2_min && ipz <= pz2_max ) { 

						t_min = t1;
						obj_num = 5;
				}
			}

			/*  Check if the current ray intercepts plane #3.  */

			Check_Plane(from.x, from.y, from.z, dx, dy, dz, a3, b3, c3, d3, &t1);

			if (t1 >= 0.0 && t1<t_min) {
				/*  Check if the intersection point is inside the min/max values. */

				Compute_Intersection(from.x, from.y, from.z, 
					dx, dy, dz, t1, &ipx, &ipy, &ipz);

				if (ipx >= px3_min && ipx <= px3_max && 
					ipy >= py3_min  && ipy <= py3_max && 
					ipz >=  pz3_min && ipz <= pz3_max ) { 

						t_min = t1;
						obj_num = 6;
				}
			}
			rgb graincolor;

/*  Compute the intensity to use at the current pixel.  */

			switch (obj_num) {

/*  The current ray does not intersect any of the objects.  */

				case 0 : r = 0.0;
					 g = 0.7;
					 b = 0.8;
					 break;

/*  The current ray intercept sphere #1.  */

				case 1 : 
					 nx = ipx - xc1;
					 ny = ipy - yc1;
					 nz = ipz - zc1;
					 Normalize(&nx, &ny, &nz);
					 
					 shadow_flag = 0;
					 shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
					 texture = 0;
					 //Finding intersection ray
				interx = dx - ipx;
				intery = dy - ipy;
				interz = dz - ipz;
				//Finding dot product and reflected ray
				dot = nx*interx + ny*intery + nz*interz;
				reflectedx = interx - 2 * (nx*dot);
				reflectedy = intery - 2 * (ny*dot);
				reflectedz = interz - 2 * (nz*dot);

				 	color = Reflection(1,ipx, ipy, ipz, reflectedx, reflectedy, reflectedz, &r, &g, &b, obj_num);

					r = color.r;
					g = color.g;
					b = color.b;
					 break;

/*  The current ray intercepts sphere #2.  */

				case 2 : 
					 nx = ipx - xc2;
					 ny = ipy - yc2;
					 nz = ipz - zc2;
					 Normalize(&nx, &ny, &nz);
					 shadow_flag = 0;
					 shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
					
					 texture = 0;
					 Bump_Map(ipx, ipy, ipz, dx, dy, dz, 1, &nx, &ny, &nz);
				 	 Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka2, kd2, ks2, phong2, &r, &g, &b, obj_num, 0,from.x, from.y, from.z, 
					dx, dy, dz);
					 break;

/*  The current ray intercepts sphere #3.  */

				case 3 : 
					 nx = ipx - xc3;
					 ny = ipy - yc3;
					 nz = ipz - zc3;
					 Normalize(&nx, &ny, &nz);
					 shadow_flag = 0;
					 shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);	 
					// texture = 0;

					
					 graincolor = woodgrain(ipx, ipy, ipz, &r, &g, &b );
					
					

					 //if (texture==1)
						// Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka3, tkd3, ks3, phong3, &r, &g, &b);
					 //else
					 Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka3, kd3, ks3, phong3, &graincolor.r, &graincolor.g, &graincolor.b, obj_num, 0,from.x, from.y, from.z, 
					dx, dy, dz);
					 break;

/*  The current ray intercepts checker board #1.  */

				case 4 : 
					 nx = a1;
					 ny = b1;
					 nz = c1;
					 shadow_flag = 0;
					shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
					woodgrain(ipx, ipy, ipz, &r, &g, &b);

					 if (ipx < 2.0 || (ipx>=4.0 && ipx<6.0)) {
						if ((ipy>=2.0 && ipy<4.0) || (ipy>=6.0))
								texture = 1;
						else
							texture = 0;
					}
					else {
						if ((ipy<2.0) || (ipy>=4.0 && ipy<6.0)) 
							texture = 1;
						else
							texture = 0;
					}
					if (texture==1)
						Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka4, tkd4, ks4, phong4, &r, &g, &b, obj_num, 0,from.x, from.y, from.z, 
					dx, dy, dz);
					else
						Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka4, kd4, ks4, phong4, &r, &g, &b, obj_num, 0,from.x, from.y, from.z, 
					dx, dy, dz);
					 break;

/*  The current ray intercepts checker board #2.  */

				case 5 : 
					 nx = a2;
					 ny = b2;
					 nz = c2;
					 shadow_flag = 0;
					 shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
					
					 if (ipz < 2.0 || (ipz>=4.0 && ipz<6.0)) {
						if ((ipy>=2.0 && ipy<4.0) || (ipy>=6.0))
								texture = 1;
						else
								texture = 0;
					}
					else {
						if ((ipy<2.0) || (ipy>=4.0 && ipy<6.0)) 
							texture = 1;
						else
							texture = 0;
					}
					if (texture == 1)
						 Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka5, tkd5, ks5, phong5, &r, &g, &b, obj_num, 0,from.x, from.y, from.z, 
					dx, dy, dz);
					else
						 Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka5, kd5, ks5, phong5, &r, &g, &b, obj_num, 0,from.x, from.y, from.z, 
					dx, dy, dz);
					 break;

/*  The current ray intercepts checker board #3.  */

				case 6 : 
					 nx = a3;
					 ny = b3;
					 nz = c3;

					 shadow_flag = 0;
					 shadow_flag = Check_Shadow(ipx, ipy, ipz, obj_num);
	
					 if (ipx < 2.0 || (ipx>=4.0 && ipx<6.0)) {
						if ((ipz>=2.0 && ipz<4.0) || (ipz>=6.0))
							texture = 1;
						else
							texture = 0;
					 }
					 else {
						if ((ipz<2.0) || (ipz>=4.0 && ipz<6.0)) 
							texture = 1;
						else
							texture = 0;
					 }
					
					 if (texture == 1)
						Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, tka6, tkd6, ks6, phong6, &r, &g, &b, obj_num, 0 ,from.x, from.y, from.z, 
					dx, dy, dz);
					 else
						 Compute_Color(shadow_flag, ipx, ipy, ipz, nx, ny, nz, ia, ka6, kd6, ks6, phong6, &r, &g, &b, obj_num, 0,from.x, from.y, from.z, 
					dx, dy, dz);
					 break;
			}

			/* Save the computed color intensity to the image buffer. */

			texture_R[xp + xmax_pixel * yp] = r;
			texture_G[xp + xmax_pixel * yp] = g;
			texture_B[xp + xmax_pixel * yp] = b;
		}
	}

/*  Write the image to the output file.  */

	printf("Writing to image...\n");
	fwrite(&xmax_pixel, sizeof(int), 1, outpfile);
	fwrite(&ymax_pixel, sizeof(int), 1, outpfile);
	
	fwrite(texture_R, sizeof(float), xmax_pixel*ymax_pixel, outpfile);
	fwrite(texture_G, sizeof(float), xmax_pixel*ymax_pixel, outpfile);
	fwrite(texture_B, sizeof(float), xmax_pixel*ymax_pixel, outpfile);

	fclose(outpfile);
}

/* Initialize the projection matrix.  */

void myinit(void)
{
/* attributes */

      glClearColor(1.0, 1.0, 1.0, 1.0); /* white background */

/* set up viewing */
/* 512 x 512 window with origin lower left */

      glMatrixMode(GL_PROJECTION);
      glLoadIdentity();
      gluOrtho2D(0.0, 512.0, 0.0, 512.0);
      glMatrixMode(GL_MODELVIEW);
}

/* Display the ray traced image.   A more efficient method is 
   to use glDrawPixels(). */

void display( void )
{
	int s, t;
	float  r, g, b;

	glClear(GL_COLOR_BUFFER_BIT);  /*clear the window */

	for(t = 0; t < ymax_pixel; t++) {
		for(s = 0; s < xmax_pixel; s++) {

			r = texture_R[s + xmax_pixel * t];
			g = texture_G[s + xmax_pixel * t];
			b = texture_B[s + xmax_pixel * t];

			glColor3f(r, g, b);
			glBegin(GL_POINTS);
               glVertex2f(s,t); 
			glEnd();
		}
	 }

     glFlush(); /* clear buffers */
}

/*  Main routine.  */

int main(int argc, char**argv)
{

	Read_Information();
	Setup_Parameters();
	Ray_Trace();

/* Standard GLUT initialization */

    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB); /* default, not needed */
    glutInitWindowSize(500,500); /* 500 x 500 pixel window */
    glutInitWindowPosition(0,0); /* place window top left on display */
    glutCreateWindow("Ray Trace"); /* window title */
    glutDisplayFunc(display); /* display callback invoked when window opened */

    myinit(); /* set attributes */
    glutMainLoop(); /* enter event loop */

	return(0);
}


