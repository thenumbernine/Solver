#include <string>

void test_discreteLaplacian();
void test_smallDense();

int main(int argc, char** argv) {
	std::string test = "discreteLaplacian";
	if (argc > 1) test = argv[1];
	if (test == "smallDense") {
		test_smallDense();
	} else {
		test_discreteLaplacian();
	}
}
