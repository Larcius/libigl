#include <igl/copyleft/blender/abf.h>
#include <igl/arap.h>
#include <igl/lscm.h>
#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/map_vertices_to_circle.h>
#include <igl/viewer/Viewer.h>
#include <nanogui/formhelper.h>
#include <igl/file_dialog_open.h>

#include <igl/distortion.h>
#include <iomanip>
#include <string>
#include <locale>
#include <codecvt>

#include "tutorial_shared_path.h"


igl::viewer::Viewer viewer;

Eigen::MatrixXd V;
Eigen::MatrixXi F;

// all kinds of parametrizations
Eigen::MatrixXd& V_uv_grid = Eigen::MatrixXd();
Eigen::MatrixXd& V_uv_lscm = Eigen::MatrixXd();
Eigen::MatrixXd& V_uv_arap = Eigen::MatrixXd();
Eigen::MatrixXd& V_uv_abf = Eigen::MatrixXd();
Eigen::MatrixXd& initial_guess = Eigen::MatrixXd();

// be sure that these vectors aggree in size and order
std::vector<std::reference_wrapper<Eigen::MatrixXd>> methods = {
	V_uv_grid,
	V_uv_lscm,
	initial_guess,
	V_uv_arap,
	V_uv_abf
};
std::vector<std::string> methodNames = {
	"Simple grid param",
	"Least squares conformal map",
	"Harmonic map",
	"As-rigid-as-possible",
	"Direct angle based flattening ++"
};
std::vector<std::string> methodAbbreviations = {
	"grid param",
	"LSCM",
	"Harmonic map",
	"ARAP",
	"dABF++"
};

// be sure that this vector is synchronized with the table
std::vector<std::string> norms = {
	"stretch",
	"shear",
	"angle",
	"area"
};

#ifdef IGL_VIEWER_WITH_NANOGUI
nanogui::ComboBox *parameterizationComboBox;
#endif

bool showUvMapping = false;
int selectedMethod = 0;

void showMesh() {
	using namespace Eigen;

#ifdef IGL_VIEWER_WITH_NANOGUI
	// temporary disable callback
	std::function<void(bool)>& tempCallback = parameterizationComboBox->changeCallback();
	// to update selected index
	parameterizationComboBox->setSelectedIndex(selectedMethod);
	// and finally enable callback
	parameterizationComboBox->setChangeCallback(tempCallback);
#endif

	viewer.data.set_uv(methods[selectedMethod], F);

	if (showUvMapping)
	{
		viewer.data.set_mesh(viewer.data.V_uv, F);

		MatrixXd V_uv_3(V.rows(), 3);
		V_uv_3 << viewer.data.V_uv, MatrixXd::Zero(V.rows(), 1);

		//viewer.core.model << Eigen::Matrix4f::Identity(4, 4);
		viewer.core.align_camera_center(V_uv_3, F);
	}
	else
	{
		viewer.data.set_mesh(V, F);
		viewer.core.align_camera_center(V, F);
	}

	viewer.data.compute_normals();
}

void showSelectedMethod() {
	assert(selectedMethod >= 0 && selectedMethod < methods.size() && "Invalid method selection");

	viewer.data.set_uv(methods[selectedMethod]);

	showMesh();
}

void load(std::string fname) {
	using namespace std;
	using namespace Eigen;

	if (fname.length() == 0)
		return;

	cout << "loading " << fname << endl;

	// clear everything
	viewer.data.clear();
	V = viewer.data.V;
	F = viewer.data.F;
	V_uv_grid = viewer.data.V_uv;
	initial_guess = viewer.data.V_uv;
	V_uv_arap = viewer.data.V_uv;
	V_uv_lscm = viewer.data.V_uv;

	//try {
	//	...
	//}
	//catch (const exception& e) {
	//	cerr << e.what() << endl;
	//}

	if (!viewer.load_mesh_from_file(fname.c_str())) {
		viewer.data.clear();
		throw exception("could not load object");
	}

	cout << "loaded" << endl;

	V = viewer.data.V;
	F = viewer.data.F;
	V_uv_grid = viewer.data.V_uv;

	showUvMapping = false;

	viewer.data.set_mesh(V, F);
	viewer.data.set_uv(V_uv_grid);
	viewer.data.compute_normals();
	viewer.core.align_camera_center(V, F);

	cout << "computing boundary..." << endl;

	VectorXi bnd;

	igl::boundary_loop(F, bnd);

	if (bnd.size() == 0) {
		throw exception("object has no boundary (e.g. a sphere) and thus cannot be flattened");
	}

	cout << "boundary computed" << endl;

	cout << "computing parametrizations..." << endl;

	// Fix two points on the boundary
	MatrixXd bc(2, 2);
	bc << 0, 0, 1, 0;
	VectorXi b(2, 1);
	b(0) = bnd(0);
	b(1) = bnd(round(bnd.size() / 2));

	// LSCM parametrization
	igl::lscm(V, F, b, bc, V_uv_lscm);
	MatrixXd mirrorX(2, 2);
	mirrorX << -1, 0, 0, 1;
	V_uv_lscm *= mirrorX;


	// Compute the initial solution for ARAP (harmonic parametrization)
	MatrixXd bnd_uv;
	igl::map_vertices_to_circle(V, bnd, bnd_uv);

	igl::harmonic(V, F, bnd, bnd_uv, 1, initial_guess);


	// Add dynamic regularization to avoid to specify boundary conditions
	igl::ARAPData arap_data;
	arap_data.with_dynamics = true;
	b = VectorXi::Zero(0);
	bc = MatrixXd::Zero(0, 0);

	// Initialize ARAP
	arap_data.max_iter = 100;
	// 2 means that we're going to *solve* in 2d
	arap_precomputation(V, F, 2, b, arap_data);

	// Solve arap using the harmonic map as initial guess
	V_uv_arap = initial_guess;

	arap_solve(bc, arap_data, V_uv_arap);

	igl::copyleft::abf_solve(V, F, V_uv_abf, true);

	// Scale UV for better comparison
	double scaleFactor = 16;
	for (unsigned int i = 0; i < methods.size(); i++) {
		MatrixXd &V_uv = methods[i];		
		V_uv *= scaleFactor / sqrt(igl::compute_area(V_uv, F));
	}

	cout << "parametrizations computed" << endl;

	cout << endl << "distortions:" << endl;

	string methodHeadline = "method";
	int methodNameMaxLength = methodHeadline.length();
	for (unsigned int i = 0; i < methods.size(); i++) {
		string methodName = methodAbbreviations[i];

		if (methodName.length() > methodNameMaxLength) {
			methodNameMaxLength = methodName.length();
		}
	}

	int precision = 3; // some value greater than 0
	int normsLength = precision + 6;

	for (unsigned int i = 0; i < norms.size(); i++) {
		if (norms[i].length() > normsLength) {
			normsLength = norms[i].length();
		}
	}

	cout.setf(ios::scientific);
	cout.precision(precision);

	cout << (char)0xC9 << string(methodNameMaxLength + 2, (char)0xCD) << (char)0xCB << string(normsLength + 2, (char)0xCD) << (char)0xCB << string(normsLength + 2, (char)0xCD) << (char)0xCB << string(normsLength + 2, (char)0xCD) << (char)0xCB << string(normsLength + 2, (char)0xCD) << (char)0xBB << endl;

	cout << (char)0xBA << " " << methodHeadline << string(methodNameMaxLength - methodHeadline.length(), ' ') << " " << (char)0xBA;
	for (unsigned int i = 0; i < norms.size(); i++) {
		if (i > 0) {
			cout << (char)0xB3;
		}

		string norm = norms[i];

		cout << " " << norm << string(normsLength - norm.length(), ' ') << " ";
	}
	cout << (char)0xBA << endl;

	cout << (char)0xCC << string(methodNameMaxLength + 2, (char)0xCD) << (char)0xCE << string(normsLength + 2, (char)0xCD) << (char)0xCE << string(normsLength + 2, (char)0xCD) << (char)0xCE << string(normsLength + 2, (char)0xCD) << (char)0xCE << string(normsLength + 2, (char)0xCD) << (char)0xB9 << endl;

	for (unsigned int i = 0; i < methods.size(); i++) {
		if (i > 0) {
			cout << (char)0xCC << string(methodNameMaxLength + 2, (char)0xC4) << (char)0xCE << string(normsLength + 2, (char)0xC4) << (char)0xC5 << string(normsLength + 2, (char)0xC4) << (char)0xC5 << string(normsLength + 2, (char)0xC4) << (char)0xC5 << string(normsLength + 2, (char)0xC4) << (char)0xB9 << endl;
		}

		string methodName = methodAbbreviations[i];
		MatrixXd &V_uv = methods[i];

		cout << (char)0xBA << " " << methodName << string(methodNameMaxLength - methodName.length(), ' ') << " " << (char)0xBA;
		cout << " " << igl::compute_stretch_distortion(V_uv, V, F) << " " << (char)0xB3;
		cout << " " << igl::compute_shear_distortion(V_uv, V, F) << " " << (char)0xB3;
		cout << " " << igl::compute_angle_distortion(V_uv, V, F) << " " << (char)0xB3;
		cout << " " << igl::compute_area_distortion(V_uv, V, F) << " " << (char)0xBA << endl;
	}

	cout << (char)0xC8 << string(methodNameMaxLength + 2, (char)0xCD) << (char)0xCA << string(normsLength + 2, (char)0xCD) << (char)0xCA << string(normsLength + 2, (char)0xCD) << (char)0xCA << string(normsLength + 2, (char)0xCD) << (char)0xCA << string(normsLength + 2, (char)0xCD) << (char)0xBC << endl;

	showSelectedMethod();
}

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier)
{
	using namespace std;
	using namespace Eigen;

	key = tolower(key);


	if (key == 'q') {
		showUvMapping = false;
		showMesh();
	}
	else if (key == 'w') {
		showUvMapping = true;
		showMesh();
	}
	else if (key >= '1' && key <= '9') {
		int selection = key - '1';

		if (selection < methods.size()) {
			selectedMethod = selection;

			showMesh();
			return true;
		}
	}
	else if (key == 'p') {
		load(igl::file_dialog_open());
		return true;
	}
	else {
		return false;
	}
}

void addCustomGUI(igl::viewer::Viewer& viewer) {
#ifdef IGL_VIEWER_WITH_NANOGUI
	using namespace nanogui;

	// Extend viewer menu
	viewer.callback_init = [&](igl::viewer::Viewer& viewer)
	{
		FormHelper *ngui = viewer.ngui;

		Window *window = ngui->addWindow(Eigen::Vector2i(230, 10), "ParamCompare");

		// load/save buttons
		ngui->addGroup("File");

		ngui->addButton("Load", [&]() {
			load(igl::file_dialog_open());
		});

		ngui->addButton("Save UV", [&]() {
			showUvMapping = true;
			showMesh();

			viewer.open_dialog_save_mesh();
		});

		// Parameterization
		Label *parameterizationGroup = ngui->addGroup("Parameterization");

		parameterizationComboBox = new ComboBox(window, methodNames, methodAbbreviations);
		parameterizationComboBox->setCallback([&](int selected) {
			selectedMethod = selected;
			showMesh();
		});

		ngui->addWidget("Method", parameterizationComboBox);
		ngui->addVariable<bool>("show UV mapping", [&](bool checked)
		{
			showUvMapping = checked;
			showMesh();
		}, [&]()
		{
			return showUvMapping;
		});

		// call to generate menu
		viewer.screen->performLayout();
		return false;
	};

#endif
}

int main(int argc, char *argv[])
{
	printf("tutorial::599_ParamCompare usage:\n");
	printf("  Q,q     Show 3D model\n");
	printf("  W,w     Show UV map\n");
	printf("  P,p     Open file dialog to load a model\n");
	printf(" params:\n");

	for (unsigned int i = 0; i < methods.size(); i++) {
		printf("  %d       %s\n", i + 1, methodAbbreviations[i]);
	}
	printf("\n\n");

	viewer.callback_key_down = &key_down;

	// set color and width of wireframe
	const Eigen::Vector4f color(1, 0, 0, 1);
	viewer.core.line_color = color;
	viewer.core.line_width = 1.3;

	// Draw checkerboard texture
	viewer.core.show_texture = true;

	addCustomGUI(viewer);

	// Launch the viewer
	viewer.launch();
}