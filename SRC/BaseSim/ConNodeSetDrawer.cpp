#include "ConNodeSetDrawer.h"
#include <math.h>
#include <GL/glew.h>
#include <GL/gl.h>
using namespace LSW_SIM;

void ConNodeSetDrawer::draw(const set<pConNodesOfFrame> &node_set){
  
  glDisable(GL_LIGHTING);
  glColor4d(0,0,1,0);
  glPointSize( 8.0 );

  set<pConNodesOfFrame>::iterator it=node_set.begin();
  for (; it != node_set.end(); ++it){
	const VectorXd &points = (*it)->getBarycenter();
	if (points.size() >= 3){
	 
	  glEnableClientState( GL_VERTEX_ARRAY );
	  glVertexPointer(3, GL_DOUBLE, 0, &points[0]);
	  glDrawArrays(GL_POINTS, 0, points.size()/3);
	  glDisableClientState( GL_VERTEX_ARRAY );
	}
  }

  glEnable(GL_LIGHTING);
}

void ConNodeSetDrawer::draw(const set<pConNodesOfFrame> &node_set, const float radius){

  glDisable(GL_LIGHTING);
  glColor4d(0,1.0f,0,1.0f);
  set<pConNodesOfFrame>::iterator it=node_set.begin();
  for (; it != node_set.end(); ++it){
	const VectorXd &points = (*it)->getBarycenter();
	for(int i = 0 ;i*3+2 < points.size(); ++i){
	  glPushMatrix();
	  glTranslatef(points[i*3+0],points[i*3+1],points[i*3+2]);
	  gluSphere(gluNewQuadric(),radius,10,10);
	  glPopMatrix();
	}
  }
  glEnable(GL_LIGHTING);
}
