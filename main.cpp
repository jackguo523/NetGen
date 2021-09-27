//****************************README****************************
// This is a light software for generating network phantoms and saving as OBJ files.
// Author: JIAMING GUO, STIM LAB, UH, 5/30/2017

// STD include
#include <vector>
#include <thread>

// OPENGL include
#include <GL/glut.h>
#include <GL/freeglut.h>

// STIM include
#include <stim/visualization/gl_network.h>
#include <stim/visualization/gl_aaboundingbox.h>
#include <stim/parser/arguments.h>
#include <stim/visualization/camera.h>


//****************************parameters****************************
// hard-coded
float zoom_factor = 10.0f;			// zooming factor
float radius_factor;				// radius increasing/decreasing factor
float z_factor;						// z-axis moving factor
float camera_factor = 2.0f;			// start point of the camera as a function of X and Y size
float orbit_factor = 0.01f;			// degrees per pixel used to orbit the camera
float eps = 10.0f;					// threshold epsilon
GLint subdivision = 20;				// slices and stacks size
float default_radius = 5.0f;		// default radius

// structure
struct vertex {
	stim::vec3<float> c;		// position
	float r = default_radius;	// radius
};
struct edge {
	unsigned v[2];								// start and end vertex indices
};

// global parameter
int vX, vY;				// viewport size
int X, Y;				// workspace size
stim::camera cam;		// camera object
float move_pace;		// camera move pace
std::vector<vertex> V;	// list of vertices
std::vector<edge> E;	// list of edges
unsigned num = 0;		// number of vertices in a new line of edges
int iter = 0;			// iterator indicates index of current vertex
int name = 0;			// output network's main name in sequences
int sub_name = 0;		// output network's sub name in sequences
vertex new_vertex;		// new vertex
vertex tmp_vertex;		// temp vertex
edge new_edge;			// new edge
edge tmp_edge;			// temp edge
std::string stackdir;

// glut events
int mouse_x, mouse_y;	// mouse click position
int mods;				// special keyboard input
bool LTbutton = false;	// true->left button down
stim::vec3<float> A = stim::vec3<float>(FLT_MAX, FLT_MAX, FLT_MAX);		// minimum point in the bounding box
stim::vec3<float> B = stim::vec3<float>(-FLT_MAX, -FLT_MAX, -FLT_MAX);	// maximum point in the bounding box
std::vector<unsigned> color;	// color scheme for each edge

// indexing
unsigned radius_index = UINT_MAX;
unsigned z_index = UINT_MAX;
unsigned color_index = 0;

// flags
bool first_click = true;			// flag indicates first clicking for a new line of edges
bool ortho = true;					// flag indicates in ortho mode or in 3D projection mode
bool render_z = false;				// flag indicates rendering z position value aside the vertex in green
bool render_radius = false;			// flag indicates rendering radius value aside the vertex in green

// colors
#define JACK_CTRL_PTS 11
static float JACKCP[JACK_CTRL_PTS * 3] = { 0.671f, 0.851f, 0.914f,
										0.502f, 0.804f, 0.757f,
										0.651f, 0.851f, 0.416f,
										0.945f, 0.714f, 0.855f,
										0.600f, 0.439f, 0.671f,
										0.914f, 0.761f, 0.490f,
										0.729f, 0.729f, 0.729f,
										0.957f, 0.647f, 0.510f,
										0.996f, 0.878f, 0.565f,
										0.992f, 0.722f, 0.388f,
										0.957f, 0.427f, 0.263f };


//****************************auxiliary functions*********************************
// find the nearest vertex of current click position
inline bool epsilon_vertex(int x, int y, float eps, unsigned& v) {

	float d = FLT_MAX;									// minimum distance between 2 vertices
	float tmp_d = 0.0f;									// temporary stores distance for loop
	unsigned tmp_i = 0;									// temporary stores connection index for loop
	stim::vec3<float> tmp_v;							// temporary stores current loop point
	stim::vec3<float> copy_v;							// copy vertex of the existing vertex (change z position to 0.0f)
	d = FLT_MAX;										// set to max of float number

	for (unsigned i = 0; i < V.size(); i++) {
		tmp_v = stim::vec3<float>((float)x, (float)(vY - y), 0.0f);
		tmp_v[0] = tmp_v[0] * (float)X / vX;
		tmp_v[1] = tmp_v[1] * (float)Y / vY;
		copy_v = stim::vec3<float>(V[i].c[0], V[i].c[1], 0.0f);

		tmp_v = tmp_v - copy_v;							// calculate a vector between two vertices
		tmp_d = tmp_v.len();							// calculate length of that vector
		if (tmp_d < d) {
			d = tmp_d;									// if found a nearer vertex 
			tmp_i = i;									// get the index of that vertex
		}
	}

	if (d < eps) {										// if current click is close to vertex we set before
		// must have at least three point to make a plane or loop 
		if (tmp_i < num && (tmp_i == V.size() - 1 || tmp_i == V.size() - 2) && !first_click && !mods) {
			// do nothing
		}
		else {
			v = tmp_i;									// copy the extant vertex's index to v
		}
		return true;
	}

	return false;
}

// check whether there is a edge between two vertices
inline bool is_edge(unsigned idx) {

	for (unsigned i = 0; i < E.size(); i++) {	// brute force method
		if (E[i].v[0] == new_edge.v[0] && E[i].v[1] == idx)
			return true;
		else if (E[i].v[1] == new_edge.v[0] && E[i].v[0] == idx)
			return true;
	}

	return false;
}

// get the bounding box of current network
inline void get_bounding_box() {
	
	unsigned n = V.size();
	float tmp;
	for (unsigned i = 0; i < n; i++) {
		for (unsigned j = 0; j < 3; j++) {
			tmp = V[i].c[j];
			if (tmp < A[j])
				A[j] = tmp;
			if (tmp > B[j])
				B[j] = tmp;
		}
	}
}


//****************************glut functions****************************
// set up squash transform
void glut_render_ortho_projection() {
	
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	vX = glutGet(GLUT_WINDOW_WIDTH);
	vY = glutGet(GLUT_WINDOW_HEIGHT);
	glViewport(0, 0, vX, vY);
	glOrtho(0.0f, X, 0.0f, Y, -200.0f, 200.0f);		// hard coded 400 z-stack
}

void glut_render_projection() {

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	vX = glutGet(GLUT_WINDOW_WIDTH);
	vY = glutGet(GLUT_WINDOW_HEIGHT);
	glViewport(0, 0, vX, vY);
	float aspect = (float)vX / (float)vY;			// aspect ratio
	gluPerspective(60, aspect, 0.1, 1000000);
}

// translate camera to origin
void glut_render_ortho_modelview() {
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void glut_render_modelview() {
	
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	stim::vec3<float> eye = cam.getPosition();		// camera position
	stim::vec3<float> focus = cam.getLookAt();		// camera focus position
	stim::vec3<float> up = cam.getUp();				// camera UP direction

	gluLookAt(eye[0], eye[1], eye[2], focus[0], focus[1], focus[2], up[0], up[1], up[2]);	// set up camera properties
}

// render vertex as sphere
void glut_render_sphere(GLint subdivision) {
	
	if (ortho)
		glColor3f(0.992f, 0.859f, 0.780f);
	else
		glColor3f(JACKCP[0], JACKCP[1], JACKCP[2]);
	unsigned n = V.size();					// get the number of vertices
	for (unsigned i = 0; i < n; i++) {
		glPushMatrix();
		glTranslatef(V[i].c[0], V[i].c[1], V[i].c[2]);
		glutSolidSphere(V[i].r, subdivision, subdivision);
		glPopMatrix();
	}
}

// render edge as cylinder
void glut_render_cylinder(GLint subdivision) {

	stim::vec3<float> tmp_d;
	stim::vec3<float> center1;
	stim::vec3<float> center2;
	stim::circle<float> tmp_c;
	float r1;
	float r2;
	std::vector<typename stim::vec3<float> > cp1(subdivision + 1);
	std::vector<typename stim::vec3<float> > cp2(subdivision + 1);

	unsigned n = E.size();
	for (unsigned i = 0; i < n; i++) {
		if (ortho) {
			glEnable(GL_BLEND);									// enable color blend
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);	// set blend function
			glColor4f(JACKCP[color[i] * 3 + 0], JACKCP[color[i] * 3 + 1], JACKCP[color[i] * 3 + 2], 0.7f);
		}
		else
			glColor3f(JACKCP[0], JACKCP[1], JACKCP[2]);

		center1 = V[E[i].v[0]].c;
		center2 = V[E[i].v[1]].c;
		r1 = V[E[i].v[0]].r;
		r2 = V[E[i].v[1]].r;

		tmp_d = center2 - center1;
		tmp_d = tmp_d.norm();
		tmp_c.rotate(tmp_d);
		stim::circle<float> c1(center1, r1, tmp_d, tmp_c.U);
		cp1 = c1.glpoints(subdivision);
		stim::circle<float> c2(center2, r2, tmp_d, tmp_c.U);
		cp2 = c2.glpoints(subdivision);

		// render the edge as cylinder
		glBegin(GL_QUAD_STRIP);
		for (unsigned j = 0; j < subdivision + 1; j++) {
			stim::vec3<float> n1 = cp1[j] - center1;
			stim::vec3<float> n2 = cp2[j] - center2;
			n1 = n1.norm();
			n2 = n2.norm();
			glNormal3f(n1[0], n1[1], n1[2]);
			glVertex3f(cp1[j][0], cp1[j][1], cp1[j][2]);
			glNormal3f(n2[0], n2[1], n2[2]);
			glVertex3f(cp2[j][0], cp2[j][1], cp2[j][2]);
		}
		glEnd();

		if (ortho)		// disable BENLD after enabled
			glDisable(GL_BLEND);
	}
	glFlush();
}

// glut light source
void glut_light() {
	
	stim::vec3<float> p1 = (cam.getUp() + cam.getUp().cross(cam.getPosition())) * 100000;
	stim::vec3<float> p2 = (cam.getUp() + cam.getPosition().cross(cam.getUp())) * 100000;

	// light source
	GLfloat global_ambient[] = { 0.4, 0.4, 0.4, 1.0 };
	GLfloat ambient[] = { 0.2, 0.2, 0.2, 1.0 };
	GLfloat diffuse1[] = { 0.4, 0.4, 0.4, 1.0 };
	GLfloat diffuse2[] = { 0.2, 0.2, 0.2, 1.0 };
	GLfloat specular[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat position1[] = { p1[0], p1[1], p1[2], 1.0 };		// upper-right light source
	GLfloat position2[] = { p2[0], p2[1], p2[2], 1.0 };		// lower-left light source
	
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glShadeModel(GL_SMOOTH);

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);

	glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);				// set ambient for light 0
	glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse1);				// set diffuse for light 0
	glLightfv(GL_LIGHT0, GL_SPECULAR, specular);			// set specular for light 0
	glLightfv(GL_LIGHT0, GL_POSITION, position1);			// set position for light 0

	glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);				// set ambient for light 1
	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse2);				// set diffuse for light 1
	glLightfv(GL_LIGHT1, GL_SPECULAR, specular);			// set specular for light 1
	glLightfv(GL_LIGHT1, GL_POSITION, position2);			// set position for light 1
}

// main render function
void glut_render() {
	
	if (!ortho) {
		glut_light();
		glEnable(GL_COLOR_MATERIAL);
	}

	// initialize glut window
	glEnable(GL_SMOOTH);
	if (ortho)
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	else
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	if (ortho) {
		glut_render_ortho_projection();
		glut_render_ortho_modelview();
	}
	else {
		glEnable(GL_DEPTH_TEST);
		glut_render_projection();
		glut_render_modelview();
	}

	glut_render_cylinder(subdivision);
	glut_render_sphere(subdivision);

	if (!ortho)
		glDisable(GL_COLOR_MATERIAL);

	if (!first_click && ortho) {							// render a transparent line to indicate your next click position
		glEnable(GL_BLEND);									// enable color blend
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);	// set blend function
		glColor4f(JACKCP[color_index * 3 + 0], JACKCP[color_index * 3 + 1], JACKCP[color_index * 3 + 2], 0.2f);
		stim::vec3<float> tmp_d;
		stim::circle<float> tmp_c;
		std::vector<typename stim::vec3<float> > cp1(subdivision + 1);
		std::vector<typename stim::vec3<float> > cp2(subdivision + 1);
		tmp_d = tmp_vertex.c - V[tmp_edge.v[0]].c;
		tmp_d = tmp_d.norm();
		tmp_c.rotate(tmp_d);
		stim::circle<float> c1(V[tmp_edge.v[0]].c, V[tmp_edge.v[0]].r, tmp_d, tmp_c.U);
		stim::circle<float> c2(tmp_vertex.c, tmp_vertex.r, tmp_d, tmp_c.U);
		cp1 = c1.glpoints(subdivision);
		cp2 = c2.glpoints(subdivision);
		glBegin(GL_QUAD_STRIP);
		for (unsigned j = 0; j < subdivision + 1; j++) {
			glVertex3f(cp1[j][0], cp1[j][1], cp1[j][2]);
			glVertex3f(cp2[j][0], cp2[j][1], cp2[j][2]);
		}
		glEnd();
		glFlush();
		glDisable(GL_BLEND);
	}

	// render the vertex index aside the vertex
	if (ortho) {
		for (unsigned i = 0; i < V.size(); i++) {
			if (i != radius_index && i != z_index) {				// only render the one that doen't need to display radius
				glColor3f(0.0f, 0.0f, 0.0f);
				glRasterPos3f(V[i].c[0], V[i].c[1], V[i].c[2]);		// mark index right above the vertex
				std::stringstream ss;
				ss << i + 1;
				glutBitmapString(GLUT_BITMAP_HELVETICA_18, (const unsigned char*)(ss.str().c_str()));
			}
		}
	}

	// render z position value
	if (render_z) {
		glColor3f(0.835f, 0.243f, 0.310f);
		glRasterPos3f(V[z_index].c[0], V[z_index].c[1], V[z_index].c[2]);
		std::stringstream ss;
		ss << "z=" << V[z_index].c[2];
		glutBitmapString(GLUT_BITMAP_HELVETICA_18, (const unsigned char*)(ss.str().c_str()));
	}

	// render radius value
	if (render_radius) {
		glColor3f(0.835f, 0.243f, 0.310f);
		glRasterPos3f(V[radius_index].c[0], V[radius_index].c[1], V[radius_index].c[2]);
		std::stringstream ss;
		ss << "r=" << V[radius_index].r;
		glutBitmapString(GLUT_BITMAP_HELVETICA_18, (const unsigned char*)(ss.str().c_str()));
	}

	glutSwapBuffers();
}

// register mouse click events
void glut_mouse(int button, int state, int x, int y) {

	mods = glutGetModifiers();					// get special keyboard input
	
	mouse_x = x;
	mouse_y = y;

	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
		LTbutton = true;
	else if (button == GLUT_LEFT_BUTTON && state == GLUT_UP)
		LTbutton = false;

	// generate a new line of edges by clicking on window console in ortho-projection
	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN && ortho && !mods) {
		unsigned idx = UINT_MAX;				// stores the vertex's index

		if (first_click) {						// first time clicking for a new line of edges
			bool flag = epsilon_vertex(x, y, eps, idx);
			if (flag) {
				new_edge.v[0] = idx;			// store the geometry starting vertex index
				tmp_edge.v[0] = idx;
				num++;							// number of edge increments
			}
			else {
				new_vertex.c = stim::vec3<float>(x, (vY - y), 0);	// make a new vertex
				new_vertex.c[0] = new_vertex.c[0] * (float)X / vX;
				new_vertex.c[1] = new_vertex.c[1] * (float)Y / vY;
				new_edge.v[0] = iter;								// make a new edge and set the starting vertex
				tmp_edge.v[0] = iter;
				V.push_back(new_vertex);							// push a new vertex
				iter++;												// iterator + 1
				num++;												// added a vertex
			}
			first_click = false;
		}
		else {													// following click of one line of edge
			bool flag = epsilon_vertex(x, y, eps, idx);
			if (flag) {
				if (!is_edge(idx)) {							// no edge between two vertices
					if (idx != UINT_MAX) {						// acceptable click
						new_edge.v[1] = idx;
						if (new_edge.v[0] != new_edge.v[1]) {	// simple graph, no self-loop
							E.push_back(new_edge);
							color.push_back(color_index);		// record the color scheme
							first_click = true;
							num = 0;							// start a new line of edges
							color_index = (color_index == JACK_CTRL_PTS - 1) ? 0 : color_index + 1;	// update color scheme for new line of edges
						}
						else {
							std::cout << "Wrong click";
							std::cout.flush();
						}
					}
				}
				else {
					std::cout << "There exists an edge";
					std::cout.flush();
				}
			}
			else {
				new_vertex.c = stim::vec3<float>(x, (vY - y), 0);	// make a new vertex
				new_vertex.c[0] = new_vertex.c[0] * (float)X / vX;
				new_vertex.c[1] = new_vertex.c[1] * (float)Y / vY;
				new_edge.v[1] = iter;								// make a new edge and set the starting vertex to current
				V.push_back(new_vertex);							// push a new vertex
				E.push_back(new_edge);								// push a new edge
				color.push_back(color_index);						// record the color scheme
				new_edge.v[0] = iter;								// make a new edge and set the starting vertex to current
				tmp_edge.v[0] = iter;
				iter++;												// iterator + 1
				num++;												// added a vertex
			}
		}
	}
	// keyboard input z position value for one vertex
	else if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN && ortho && mods == GLUT_ACTIVE_SHIFT) {
		
		bool flag = epsilon_vertex(x, y, eps, z_index);
		if (flag) {
			render_z = true;
			float value = 0.0f;
			std::cout << "Please type in one z-position value for vertex " << z_index + 1 << std::endl;
			std::cin >> value;
			V[z_index].c[2] = value;
		}
		else {
			z_index = UINT_MAX;
			render_z = false;
		}	
	}
}

// register mouse move events
void glut_motion(int x, int y) {
		
	if (!ortho && LTbutton){		// move camera in 3D projection mode
		float theta = orbit_factor * (mouse_x - x);		// determine the number of degrees along the x-axis to rotate
		float phi = orbit_factor * (y - mouse_y);		// number of degrees along the y-axis to rotate

		cam.OrbitFocus(theta, phi);						// rotate the camera around the focal point
	}

	mouse_x = x;
	mouse_y = y;

	glutPostRedisplay();
}

// register mouse passive move events
void glut_passive_motion(int x, int y) {
	
	mods = glutGetModifiers();

	if (ortho) {
		tmp_vertex.c = stim::vec3<float>((float)x, (float)(vY - y), 0.0f);
		tmp_vertex.c[0] = tmp_vertex.c[0] * (float)X / vX;
		tmp_vertex.c[1] = tmp_vertex.c[1] * (float)Y / vY;
	}
	
	if (ortho) {		// display z position
		bool flag = epsilon_vertex(x, y, eps, z_index);
		if (!flag) {
			render_radius = false;
			radius_index = UINT_MAX;
			render_z = false;
			z_index = UINT_MAX;
		}
	}

	mouse_x = x;
	mouse_y = y;

	glutPostRedisplay();
}

// regiser wheel events
void glut_wheel(int wheel, int direction, int x, int y) {

	mods = glutGetModifiers();					// get special keyboard input

	// move camera in 3D projection
	if (!ortho) {
		if (direction > 0)						// if button 3(up), move closer
			move_pace = zoom_factor;
		else
			move_pace = -zoom_factor;

		cam.Push(move_pace);
	}
	else {
		bool flag = epsilon_vertex(x, y, eps, radius_index);
		// increase vertex radius
		if (flag && !mods) {
			render_radius = true;
			render_z = false;
			if (direction > 0)
				V[radius_index].r += radius_factor;	// increase radius
			else
				V[radius_index].r -= radius_factor;
			if (V[radius_index].r <= 0)				// degenerated case
				V[radius_index].r = default_radius;
		}
		else if (!flag && !mods){
			radius_index = UINT_MAX;
			render_radius = false;
		}
		// change vertex z position
		else if (flag && mods == GLUT_ACTIVE_SHIFT) {
			z_index = radius_index;
			render_radius = false;
			render_z = true;

			if (direction > 0)
				V[z_index].c[2] -= z_factor;
			else
				V[z_index].c[2] += z_factor;		// move vertex along positive z-axis
		}
		else if (!flag && mods == GLUT_ACTIVE_SHIFT) {
			z_index = UINT_MAX;
			render_z = false;
		} 
	}

	glutPostRedisplay();
}

// register keyboard inputs
void glut_keyboard(unsigned char key, int x, int y) {
	
	switch (key) {
	// press space to start a new line of edges when pressing SPACE
	case 32:
		first_click = true;
		num = 0;
		color_index = (color_index == JACK_CTRL_PTS - 1) ? 0 : color_index + 1;	// update color scheme for new line of edges
		break;

	// save current network as OBJ file in stackdir directory
	case 's': {
		std::stringstream output_ss;
		output_ss << name << "_" << sub_name << "_net" << ".obj";
		std::string output_filename = output_ss.str();
		std::ofstream output_file;

		output_file.open(output_filename.c_str());
		for (unsigned i = 0; i < V.size(); i++)
			output_file << "v" << " " << V[i].c[0] << " " << V[i].c[1] << " " << V[i].c[2] << std::endl;
		for (unsigned i = 0; i < V.size(); i++)
			output_file << "vt" << " " << V[i].r << std::endl;
		for (unsigned i = 0; i < E.size(); i++)
			output_file << "l" << " " << E[i].v[0] + 1 << "/" << E[i].v[0] + 1 << " " << E[i].v[1] + 1 << "/" << E[i].v[1] + 1 << std::endl;
		output_file.close();
		sub_name++;			// sub name change

		break;
	}

	// undo
	case 'u': {
		
		// first vertex on a new line of edges
		if (num == 1) {
			bool flag = false;						// check whether current vertex belongs to another edge
			for (unsigned i = 0; i < E.size(); i++) {
				if (new_edge.v[0] == E[i].v[0] || new_edge.v[0] == E[i].v[1]) {
					flag = true;
					break;
				}
			}
			if (new_edge.v[0] == V.size() - 1 && !flag) {	// new vertex
				V.pop_back();								// pop back new vertex
				iter--;
			}	
			first_click = true;
			num = 0;
		}
		// not first vertex
		else if(num > 1) {
			new_edge.v[0] = E[E.size() - 1].v[0];
			tmp_edge.v[0] = new_edge.v[0];
			E.pop_back();							// pop back new "things"
			color.pop_back();
			V.pop_back();
			iter--;
			num--;
		}
		break;
	}
	
	// close window and exit application when pressing ESC
	case 27:
		std::exit(1);
	}

	glutPostRedisplay();
}

// dynamically set menu
void glut_set_menu(int i) {
	
	if (i == 2) {
		glutChangeToMenuEntry(1, "3D projection", 3);
		glutAddMenuEntry("new network", 1);
	}
		
	else if (i == 3) {
		glutChangeToMenuEntry(1, "2D projection", 2);
		glutRemoveMenuItem(2);
	}
		
}

// register glut menu options
void glut_menu(int value) {

	if (value == 1) {		// create new network
		ortho = true;
		first_click = true;

		// clear up previous work
		glClear(GL_COLOR_BUFFER_BIT);
		V.clear();
		E.clear();
		iter = 0;
		num = 0;
	}

	if (value == 2) {		// ortho mode
		ortho = true;
		glDisable(GL_DEPTH_TEST);
		glDisable(GL_LIGHTING);
		glDisable(GL_LIGHT0);
		glDisable(GL_LIGHT1);
	}

	if (value == 3) {		// 3D projection mode
		ortho = false;
		render_z = false;

		// initialize camera object
		get_bounding_box();							// get the bounding box of current network
		stim::vec3<float> c = (A + B) * 0.5f;		// get the center of the bounding box
		stim::vec3<float> size = (B - A);			// get the size of the bounding box
		cam.setPosition(c + stim::vec3<float>(0.0f, 0.0f, camera_factor * std::max(size[0], size[1])));
		cam.LookAt(c[0], c[1], c[2]);

		// enable light
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_LIGHT1);
	}

	// set up new menu
	glut_set_menu(value);

	glutPostRedisplay();
}

// window reshape function
void glut_reshape(int x, int y) {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glut_render_ortho_projection();
	glut_render_ortho_modelview();
}

// glut initialization
void glut_initialize() {
	int myargc = 1;
	char* myargv[1];
	myargv[0] = strdup("GENERATE_NETWORK");

	glutInit(&myargc, myargv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
	glutInitWindowPosition(800, 0);
	glutInitWindowSize(1000, 1000);
	glutCreateWindow("Generate Network Phantom");

	glutDisplayFunc(glut_render);
	glutMouseFunc(glut_mouse);
	glutMotionFunc(glut_motion);
	glutPassiveMotionFunc(glut_passive_motion);
	glutMouseWheelFunc(glut_wheel);
	glutKeyboardFunc(glut_keyboard);
	glutReshapeFunc(glut_reshape);

	// initilize menu
	glutCreateMenu(glut_menu);					// create a menu object 
	glutAddMenuEntry("3D projection", 3);
	glutAddMenuEntry("new network", 1);
	glutAttachMenu(GLUT_RIGHT_BUTTON);			// register right mouse to open menu option
}

// output an advertisement for the lab, authors and usage information
void advertise() {
	std::cout << std::endl << std::endl;
	std::cout << " =======================================================================================" << std::endl;
	std::cout << "|				     Thank you for using the network phantom builder                      |" << std::endl;
	std::cout << "|        Scalable Tissue Imaging and Modeling (STIM) Lab, University of Houston         |" << std::endl;
	std::cout << "|                              Developers: Jiaming Guo                                  |" << std::endl;
	std::cout << " =======================================================================================" << std::endl << std::endl;

	std::cout << "Keyboard inputs: " << std::endl;
	std::cout << "SPACE : start/end a line of edges." << std::endl;
	std::cout << "u : undo clicking." << std::endl;
	std::cout << "s : save current network in OBJ file." << std::endl;
}

// main function
int main(int argc, char* argv[]) {

	stim::arglist args;				// create an instance of arglist

	// add arguments
	args.add("help", "print help");
	args.add("radius", "set default radii", "5", "real value > 0");
	args.add("rstep", "radius changing rate", "0.1", "real value > 0");
	args.add("zstep", "step size along z-axis", "0.1", "real value > 0");
	args.add("workspace", "sets the size of the workspace", "400", "real value > 0");
	args.add("stackdir", "set the directory of the output image stack", "", "any existing directory (ex. /home/name/network)");

	args.parse(argc, argv);			// parse the command line

	// print advertise and help if need
	if (args["help"].is_set()) {
		//
		std::cout << args.str();
		std::exit(1);
	}

	// get the default radius
	default_radius = args["radius"].as_float();

	// get the radius changing factor
	radius_factor = args["rstep"].as_float();

	// get the z-step size
	z_factor = args["zstep"].as_float();

	// get the workspace size
	X = Y = args["workspace"].as_float();

	// get the save directory
	if (args["stackdir"].is_set())
		stackdir = args["stackdir"].as_string();

	// glut main loop
	glut_initialize();
	glutMainLoop();

	return 0;
}