#include <igl/copyleft/blender/abf.h>
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>

#include "tutorial_shared_path.h"

Eigen::MatrixXd V_uv;
Eigen::MatrixXd V;
Eigen::MatrixXi F;

bool show_uv = false;

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
	using namespace Eigen;

	if (key == '1')
		show_uv = false;
	else if (key == '2')
		show_uv = true;
	
	if (show_uv)
	{
		viewer.data.set_mesh(V_uv, F);

		MatrixXd V_uv_3(V_uv.rows(), 3);
		V_uv_3 << V_uv, MatrixXd::Zero(V_uv.rows(), 1);

		viewer.core.align_camera_center(V_uv_3, F);
	}
	else
	{
		viewer.data.set_mesh(V, F);
		viewer.core.align_camera_center(V, F);
	}

	viewer.data.compute_normals();

	return false;
}

int main(int argc, char *argv[])
{
	using namespace std;

	printf("Load a mesh in OFF format\n");
	// Load a mesh in OFF format
	//igl::readOFF(TUTORIAL_SHARED_PATH "/camelhead.off", V, F);
	igl::readOFF(TUTORIAL_SHARED_PATH "/own/1triangle.off", V, F);
	//igl::readOFF(TUTORIAL_SHARED_PATH "/own/tetrahedron_bottomless.off", V, F);
	//igl::readOFF(TUTORIAL_SHARED_PATH "/own/bunny_bottomless.off", V, F);

	igl::copyleft::abf_solve(V, F, V_uv, true);

	// Scale UV to make the texture more clear
	V_uv *= 20;

	// Plot the mesh
	igl::viewer::Viewer viewer;
	viewer.data.set_mesh(V, F);
	viewer.data.set_uv(V_uv);
	viewer.callback_key_down = &key_down;

	// set color and width of wireframe
	const Eigen::Vector4f color(1, 0, 0, 1);
	viewer.core.line_color = color;
	viewer.core.line_width = 1.3;

	// Draw checkerboard texture
	viewer.core.show_texture = true;

	// Launch the viewer
	viewer.launch();
}
