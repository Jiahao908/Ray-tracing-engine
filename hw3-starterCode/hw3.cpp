/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Jiahao Fu
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <glm/glm.hpp>
#include <iostream>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0
#define FAR -100000.0

using namespace std;

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

double RATIO = WIDTH / HEIGHT;
double xamnt, yamnt; //x and y coordinates in the pixel plane

glm::vec3 sphere_inter_color;

//This function will return a pair, the first vec3 of which is the coordinates of intersection points and
//the second of which is the index of sphere that the ray intersects with.
//This function will return -1 if don't have any intersections.
pair<glm::vec3, int> has_sphere_intersection(glm::vec3 start, glm::vec3 end)
{
  vector< glm::vec3 > sphere_intersection;
  int sphere_i = -1;
  glm::vec3 intersect = glm::vec3(0, 0, FAR);
  int t_min = 100000000;
  double x0 = start[0];
  double y0 = start[1];
  double z0 = start[2];
  double xd = end[0];
  double yd = end[1];
  double zd = end[2];
  for(int i = 0; i < num_spheres; i++) // iterate through every light source and check intersections
  {
    double xc = spheres[i].position[0];
    double yc = spheres[i].position[1];
    double zc = spheres[i].position[2];
    double r = spheres[i].radius;
    double b = 2 * (xd * (x0 - xc) + yd * (y0 - yc) + zd * (z0 - zc));
    double c = pow(x0 - xc, 2) + pow(y0 - yc, 2) + pow(z0 - zc, 2) - pow(r, 2);
    double t0, t1;
    double check = b * b - 4 * c;
    if(check < 0)
    {
      continue;
    }
    else
    {
      t0 = (-b - sqrt(check))/2.0;
      t1 = (-b + sqrt(check))/2.0;
    }
    if(t0 < t_min && t0 > 0.001)
    {
      intersect = glm::vec3(x0 + xd * t0, y0 + yd * t0, z0 + zd * t0);
      sphere_i = i;
      t_min = t0; //update t_min, sphere_i and intersect point
    }
  }
  /*cout << "DONE!" << endl;
  cout << sphere_intersection.size() << endl;*/
  return make_pair(intersect, sphere_i);
}

//calculate alpha, beta and gamma for triangles phong shading
glm::vec3 calculate_abg(glm::vec3 intersect, glm::vec3 Va, glm::vec3 Vb, glm::vec3 Vc)
{
  double alpha, beta, gamma;
  glm::vec3 cp_bca, B_A, C_A;
  B_A = Vb - Va;
  C_A = Vc - Va;
  cp_bca = glm::cross(B_A, C_A);
  double area_abc = 0.5 * sqrt(pow(cp_bca[0],2) + pow(cp_bca[1],2) + pow(cp_bca[2], 2));

  glm::vec3 cp_bci, B_I, C_I;
  B_I = Vb - intersect;
  C_I = Vc - intersect;
  cp_bci = glm::cross(B_I, C_I);
  double area_bci = 0.5 * sqrt(pow(cp_bci[0], 2) + pow(cp_bci[1], 2) + pow(cp_bci[2], 2));
  alpha = area_bci / area_abc;

  glm::vec3 cp_cai, A_I;
  A_I = Va - intersect;
  C_I = Vc - intersect;
  cp_cai = glm::cross(C_I, A_I);
  double area_cai = 0.5 * sqrt(pow(cp_cai[0], 2) + pow(cp_cai[1], 2) + pow(cp_cai[2], 2));
  beta = area_cai / area_abc;

  glm::vec3 cp_abi;
  B_I = Vb - intersect;
  A_I = Va - intersect;
  cp_abi = glm::cross(A_I, B_I);
  double area_abi = 0.5 * sqrt(pow(cp_abi[0], 2) + pow(cp_abi[1], 2) + pow(cp_abi[2], 2));
  gamma = area_abi / area_abc;

  return glm::vec3(alpha, beta, gamma);

}

//3D check whether the intersection point is inside the triangle
bool inside_triangle_3D(glm::vec3 intersect, glm::vec3 Va, glm::vec3 Vb, glm::vec3 Vc)
{
  double alpha, beta, gamma;
  glm::vec3 cp_bca, B_A, C_A;
  B_A = Vb - Va;
  C_A = Vc - Va;
  cp_bca = glm::cross(B_A, C_A);
  double area_abc = 0.5 * sqrt(pow(cp_bca[0],2) + pow(cp_bca[1],2) + pow(cp_bca[2], 2));

  glm::vec3 cp_bci, B_I, C_I;
  B_I = Vb - intersect;
  C_I = Vc - intersect;
  cp_bci = glm::cross(B_I, C_I);
  double area_bci = 0.5 * sqrt(pow(cp_bci[0], 2) + pow(cp_bci[1], 2) + pow(cp_bci[2], 2));
  alpha = area_bci/area_abc;
  if(glm::dot(cp_bca, cp_bci) >= 0)
    alpha = alpha;
  else
  {
    alpha = 0-alpha;
  }
  

  glm::vec3 cp_cai, A_I;
  A_I = Va - intersect;
  C_I = Vc - intersect;
  cp_cai = glm::cross(C_I, A_I);
  double area_cai = 0.5 * sqrt(pow(cp_cai[0], 2) + pow(cp_cai[1], 2) + pow(cp_cai[2], 2));
  beta = area_cai / area_abc;
  if(glm::dot(cp_bca, cp_cai) >= 0)
    beta = beta;
  else
  {
    beta = 0-beta;
  }


  glm::vec3 cp_abi;
  B_I = Vb - intersect;
  A_I = Va - intersect;
  cp_abi = glm::cross(A_I, B_I);
  double area_abi = 0.5 * sqrt(pow(cp_abi[0], 2) + pow(cp_abi[1], 2) + pow(cp_abi[2], 2));
  gamma = area_abi / area_abc;
  if(glm::dot(cp_bca, cp_abi) >= 0)
    gamma = gamma;
  else
  {
    gamma = 0-gamma;
  }
  //cout << "Alpha: " << alpha << " Beta: " << beta << " Gamma: " <<  gamma << endl;
  if(alpha >= 0 && alpha <= 1 && beta >= 0 && beta <= 1 && gamma >= 0 && gamma <=1 && alpha + beta + gamma >=0.99 && alpha + beta + gamma <= 1.01)
  {
    //cout << "Alpha: " << alpha << " Beta: " << beta << " Gamma: " <<  gamma << endl;
    return true;
  }
  else
  {
    return false;
  }
}

//2D check whether th intersection point is inside the triangle
bool inside_triangle_2D(glm::vec3 Vi, glm::vec3 Va, glm::vec3 Vb, glm::vec3 Vc)
{
  double alpha, beta, gamma;
  double area_abc = 0.5 * ((Vb[0] - Va[0]) * (Vc[1] - Va[1]) - (Vc[0] - Va[0]) * (Vb[1] - Va[1]));
  double area_abi = 0.5 * ((Vb[0] - Va[0]) * (Vi[1] - Va[1]) - (Vi[0] - Va[0]) * (Vb[1] - Va[1]));
  double area_ibc = 0.5 * ((Vb[0] - Vi[0]) * (Vc[1] - Vi[1]) - (Vc[0] - Vi[0]) * (Vb[1] - Vi[1]));
  double area_aic = 0.5 * ((Vi[0] - Va[0]) * (Vc[1] - Va[1]) - (Vc[0] - Va[0]) * (Vi[1] - Va[1]));


  alpha = area_ibc / area_abc;
  beta = area_aic / area_abc;
  gamma = area_abi / area_abc;

  if(alpha >= 0 && alpha <= 1 && beta >= 0 && beta <= 1 && gamma >= 0 && gamma <=1 &&  alpha + beta + gamma >= 0.999 && alpha + beta + gamma <= 1.001)
  {
    //cout << "Alpha: " << alpha << " Beta: " << beta << " Gamma:" <<  gamma << endl;
    return true;
  }
  else
  {
    return false;
  }
}

//This function will return a pair, the first vec3 of which is the coordinates of intersection points and
//the second of which is the index of triangle that the ray intersects with. 
//This function will return -1 if don't have any intersections.
pair<glm::vec3, int> has_triangle_intersection(glm::vec3 start, glm::vec3 end)
{
  glm::vec3 inter = glm::vec3(0, 0, FAR);
  int triangle_i = -1;
  int t_min = 10000000;
  double x0 = start[0];
  double y0 = start[1];
  double z0 = start[2];
  double xd = end[0];
  double yd = end[1];
  double zd = end[2];
  for(int i = 0; i < num_triangles; i++)
  {
    double t;
    double xa = triangles[i].v[0].position[0];
    double ya = triangles[i].v[0].position[1];
    double za = triangles[i].v[0].position[2];
    glm::vec3 Va = glm::vec3(xa, ya, za);

    double xb = triangles[i].v[1].position[0];
    double yb = triangles[i].v[1].position[1];
    double zb = triangles[i].v[1].position[2];
    glm::vec3 Vb = glm::vec3(xb, yb, zb);

    double xc = triangles[i].v[2].position[0];
    double yc = triangles[i].v[2].position[1];
    double zc = triangles[i].v[2].position[2];
    glm::vec3 Vc = glm::vec3(xc, yc, zc);

    glm::vec3 Vn = glm::normalize(glm::cross(Vb - Va, Vc - Va));
    double na = Vn[0];
    double nb = Vn[1];
    double nc = Vn[2];
    double nd = 0 - na * xa - nb * ya - nc * za;
    double check = glm::dot(Vn, end);
    if(check == 0)
      continue;
    else
    {
      t = -(na * x0 + nb * y0 + nc * z0 + nd)/check;
      glm::vec3 temp = glm::vec3(xd * t, yd * t, zd * t);
      if(inside_triangle_3D(temp, Va, Vb, Vc) && t < t_min && t > 0.001)
      {
        inter = temp;
        triangle_i = i;
        t_min = t;
      } 
    }
  }
  return make_pair(inter, triangle_i);
}

//Phong shading for sphere intersection
glm::vec3 phong_sphere(glm::vec3 inter_pos, int i)
{
  glm::vec3 pnt_color = glm::vec3(0, 0, 0);
  glm::vec3 sphere_pos = glm::vec3(spheres[i].position[0], spheres[i].position[1], spheres[i].position[2]);
  for(int lgt_i = 0; lgt_i < num_lights; lgt_i++)
  {
    glm::vec3 light_pos = glm::vec3(lights[lgt_i].position[0], lights[lgt_i].position[1], lights[lgt_i].position[2]);
    glm::vec3 light_dir = glm::normalize(light_pos - inter_pos);
    int tri = has_triangle_intersection(inter_pos, light_dir).second;
    int sph = has_sphere_intersection(inter_pos, light_dir).second;
    //cout << "tri: " <<  tri  << " sph: " << sph << endl;
    if(tri != -1 || sph != -1)
    {
      continue;
    }
    else
    {
      double dot_ln, dot_rv;
      glm::vec3 pnt_normal = glm::normalize(inter_pos - sphere_pos);
      glm::vec3 pnt_lgt = glm::normalize(light_pos - inter_pos);
      glm::vec3 pnt_vsn = glm::normalize(glm::vec3(0,0,0)-inter_pos);
      glm::vec3 pnt_rft;
      if(glm::dot(pnt_lgt, pnt_normal) <= 0)
      {
        dot_ln = 0.0;
        pnt_rft = glm::normalize(glm::vec3(0,0,0)-pnt_lgt);
      }
      else
      {
        dot_ln = glm::dot(pnt_lgt, pnt_normal);
        pnt_rft = glm::normalize(2 * glm::dot(pnt_lgt, pnt_normal) * pnt_normal - pnt_lgt);
      }
      if(glm::dot(pnt_vsn, pnt_rft) <= 0)
      {
        dot_rv = 0;
      }
      else
      {
        dot_rv = glm::dot(pnt_vsn, pnt_rft);
      }
      
      double kd_red = spheres[i].color_diffuse[0];
      double kd_green = spheres[i].color_diffuse[1];
      double kd_blue = spheres[i].color_diffuse[2];
      double ks_red = spheres[i].color_specular[0];
      double ks_green = spheres[i].color_specular[1];
      double ks_blue = spheres[i].color_specular[2];
      double shine = spheres[i].shininess;
      
      pnt_color[0] += lights[lgt_i].color[0] * (dot_ln * kd_red + ks_red * pow(dot_rv, shine));
      pnt_color[1] += lights[lgt_i].color[1] * (dot_ln * kd_green + ks_green * pow(dot_rv, shine));
      pnt_color[2] += lights[lgt_i].color[2] * (dot_ln * kd_blue + ks_blue * pow(dot_rv, shine));
      //compute the light separately
    }
  }
  pnt_color[0] += ambient_light[0];
  pnt_color[1] += ambient_light[1];
  pnt_color[2] += ambient_light[2];
  if(pnt_color[0] > 1)
    pnt_color[0] = 1;
  if(pnt_color[1] > 1)
    pnt_color[1] = 1;
  if(pnt_color[2] > 1)
    pnt_color[2] = 1;
  return pnt_color;
}

//Phong shading for triangle intersection
glm::vec3 phong_triangle(glm::vec3 inter_pos, int i)
{
  glm::vec3 pnt_color = glm::vec3(0, 0, 0);
  for(int lgt_i = 0; lgt_i < num_lights; lgt_i++) //iterate through every light source
  {
    glm::vec3 light_pos = glm::vec3(lights[lgt_i].position[0], lights[lgt_i].position[1], lights[lgt_i].position[2]);
    glm::vec3 light_dir = glm::normalize(light_pos - inter_pos);
    int tri = has_triangle_intersection(inter_pos, light_dir).second;
    int sph = has_sphere_intersection(inter_pos, light_dir).second;
    if(tri != -1 || sph != -1)
    {
      continue;
    }
    else
    {
      double xa = triangles[i].v[0].position[0];
      double ya = triangles[i].v[0].position[1];
      double za = triangles[i].v[0].position[2];
      glm::vec3 Va = glm::vec3(xa, ya, za);
      double xb = triangles[i].v[1].position[0];
      double yb = triangles[i].v[1].position[1];
      double zb = triangles[i].v[1].position[2];
      glm::vec3 Vb = glm::vec3(xb, yb, zb);
      double xc = triangles[i].v[2].position[0];
      double yc = triangles[i].v[2].position[1];
      double zc = triangles[i].v[2].position[2];
      glm::vec3 Vc = glm::vec3(xc, yc, zc);

      double alpha, beta, gamma;
      /*double area_abc = 0.5 * ((Vb[0] - Va[0]) * (Vc[1] - Va[1]) - (Vc[0] - Va[0]) * (Vb[1] - Va[1]));
      double area_abi = 0.5 * ((Vb[0] - Va[0]) * (inter_pos[1] - Va[1]) - (inter_pos[0] - Va[0]) * (Vb[1] - Va[1]));
      double area_ibc = 0.5 * ((Vb[0] - inter_pos[0]) * (Vc[1] - inter_pos[1]) - (Vc[0] - inter_pos[0]) * (Vb[1] - inter_pos[1]));
      double area_aic = 0.5 * ((inter_pos[0] - Va[0]) * (Vc[1] - Va[1]) - (Vc[0] - Va[0]) * (inter_pos[1] - Va[1]));


      alpha = area_ibc / area_abc;
      beta = area_aic / area_abc;
      gamma = area_abi / area_abc;*/
    
      glm::vec3 abg = calculate_abg(inter_pos, Va, Vb, Vc);
      alpha = abg[0];
      beta = abg[1];
      gamma = abg[2];

      glm::vec3 normal_A = glm::vec3(triangles[i].v[0].normal[0], triangles[i].v[0].normal[1], triangles[i].v[0].normal[2]);
      glm::vec3 normal_B = glm::vec3(triangles[i].v[1].normal[0], triangles[i].v[1].normal[1], triangles[i].v[1].normal[2]);
      glm::vec3 normal_C = glm::vec3(triangles[i].v[2].normal[0], triangles[i].v[2].normal[1], triangles[i].v[2].normal[2]);
      glm::vec3 normal = glm::normalize((float)alpha * normal_A + (float)beta * normal_B + (float)gamma * normal_C); 

      
      glm::vec3 diffuse_A = glm::vec3(triangles[i].v[0].color_diffuse[0], triangles[i].v[0].color_diffuse[1], triangles[i].v[0].color_diffuse[2]);
      glm::vec3 diffuse_B = glm::vec3(triangles[i].v[1].color_diffuse[0], triangles[i].v[1].color_diffuse[1], triangles[i].v[1].color_diffuse[2]);
      glm::vec3 diffuse_C = glm::vec3(triangles[i].v[2].color_diffuse[0], triangles[i].v[2].color_diffuse[1], triangles[i].v[2].color_diffuse[2]);
      glm::vec3 diffuse = glm::vec3((float)alpha * diffuse_A + (float)beta * diffuse_B + (float)gamma * diffuse_C); 

      glm::vec3 specular_A = glm::vec3(triangles[i].v[0].color_specular[0], triangles[i].v[0].color_specular[1], triangles[i].v[0].color_specular[2]);
      glm::vec3 specular_B = glm::vec3(triangles[i].v[1].color_specular[0], triangles[i].v[1].color_specular[1], triangles[i].v[1].color_specular[2]);
      glm::vec3 specular_C = glm::vec3(triangles[i].v[2].color_specular[0], triangles[i].v[2].color_specular[1], triangles[i].v[2].color_specular[2]);
      glm::vec3 specular = glm::vec3((float)alpha * specular_A + (float)beta * specular_B + (float)gamma * specular_C); 

      
      double shininess_A = triangles[i].v[0].shininess;
      double shininess_B = triangles[i].v[1].shininess;
      double shininess_C = triangles[i].v[2].shininess;
      double shininess = alpha * shininess_A + beta * shininess_B + gamma * shininess_C;

      double dot_ln, dot_rv;
      glm::vec3 lgt = glm::normalize(light_pos - inter_pos);
      glm::vec3 vsn = glm::normalize(glm::vec3(0,0,0)-inter_pos);
      glm::vec3 rft;
      if(glm::dot(lgt, normal) <= 0)
      {
        dot_ln = 0.0;
        rft = glm::normalize(glm::vec3(0,0,0) - lgt);
      }
      else
      {
        dot_ln = glm::dot(lgt, normal);
        rft = glm::normalize((float)(2 * dot_ln) * normal - lgt);
      }
      if(glm::dot(vsn, rft) <= 0)
      {
        dot_rv = 0;
      }
      else
      {
        dot_rv = glm::dot(vsn, rft);
      }

      double kd_red = diffuse[0];
      double kd_green = diffuse[1];
      double kd_blue = diffuse[2];
      double ks_red = specular[0];
      double ks_green  = specular[1];
      double ks_blue = specular[2];
      double shine = shininess;

      pnt_color[0] += lights[lgt_i].color[0] * (dot_ln * kd_red + ks_red * pow(dot_rv, shine));
      pnt_color[1] += lights[lgt_i].color[1] * (dot_ln * kd_green + ks_green * pow(dot_rv, shine)); 
      pnt_color[2] += lights[lgt_i].color[2] * (dot_ln * kd_blue + ks_blue * pow(dot_rv, shine));
    }
  }
  pnt_color[0] += ambient_light[0];
  pnt_color[1] += ambient_light[1];
  pnt_color[2] += ambient_light[2];
  if(pnt_color[0] > 1)
    pnt_color[0] = 1;
  if(pnt_color[1] > 1)
    pnt_color[1] = 1;
  if(pnt_color[2] > 1)
    pnt_color[2] = 1;
  return pnt_color;
}

//MODIFY THIS FUNCTION
void draw_scene()
{
  //a simple test output
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      //xamnt = ((x + 0.5)/WIDTH) * RATIO - (((WIDTH - HEIGHT)/(double)HEIGHT)/2);
      //yamnt = ((HEIGHT - y) + 0.5)/HEIGHT;
      xamnt = (((double)(2 * x) - WIDTH)/ HEIGHT) * 0.57735;
      yamnt = (((double)(2 * y) - HEIGHT)/ HEIGHT) * 0.57735;
      glm::vec3 start_point = glm::vec3(0,0,0);
      glm::vec3 pixel_point = glm::normalize(glm::vec3(xamnt, yamnt, -1));
      glm::vec3 tri_inter = has_triangle_intersection(start_point, pixel_point).first;
      glm::vec3 sph_inter = has_sphere_intersection(start_point, pixel_point).first;
      int tri_i = has_triangle_intersection(start_point, pixel_point).second;
      int sph_i = has_sphere_intersection(start_point, pixel_point).second;
      if(tri_inter[2] > sph_inter[2])//select a intersection with bigger z value
      {
        glm::vec3 pnt_color = phong_triangle(tri_inter, tri_i);
        plot_pixel(x, y, pnt_color[0] * 255, pnt_color[1] * 255, pnt_color[2] * 255);
        //triangle phong shading;
      }
      else if(tri_inter[2] < sph_inter[2])
      {
        //sphere phong shading;
        glm::vec3 pnt_color = phong_sphere(sph_inter, sph_i);
        plot_pixel(x, y, pnt_color[0] * 255, pnt_color[1] * 255, pnt_color[2] * 255);
      }
      else
      {
        //background plot
        plot_pixel(x, y, 255, 255, 255);
      }
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); 
  fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

