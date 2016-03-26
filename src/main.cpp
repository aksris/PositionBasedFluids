#include "main.hpp"

using namespace std;


int main()
{
	Scene scene;
    char* filename = "src/resources/scene.json";

    scene.parseScene(filename, scene);
    int width = 1024;
    int height = 768;
	Viewer *view = new Viewer(width, height, scene);
    view->display();
    return 0;
}

