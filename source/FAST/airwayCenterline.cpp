#include "FAST/Exporters/VTKMeshFileExporter.hpp"
#include "FAST/Importers/ImageFileImporter.hpp"
#include "FAST/Algorithms/CenterlineExtraction/CenterlineExtraction.hpp"
#include "FAST/Tools/CommandLineParser.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

using namespace fast;

int main(int argc, char** argv) {
  Reporter::setGlobalReportMethod(Reporter::COUT); // TODO remove
	CommandLineParser parser("Centerline extraction");
	parser.addVariable("mask_path", true, "segmentation mask .mhd");
	parser.addVariable("output_path", true, "path to output .vtk");
	parser.parse(argc, argv);

	std::string volPath = parser.get("mask_path");
	std::ifstream maskHeader(volPath);

	Vector3f offset;
	Vector3f spacing;
	Vector3i imgDim;

	while (!maskHeader.eof()) {
		std::string line;
		getline(maskHeader, line);
		std::string attr = line.substr(0, line.find(" "));

		if (attr != "Offset" && attr != "ElementSpacing" && attr != "DimSize") {
			continue;
		}

		std::string value = line.substr(line.find("=") + 2, std::string::npos);
		std::string elmParts[3];

		for (int i = 0; i < 3; ++i) {
			size_t pos = value.find(" ");
			if (pos == std::string::npos) {
				elmParts[i] = value;
			} else {
				elmParts[i] = value.substr(0, pos);
				value.erase(0, pos + 1);
			}
		}

		if (attr == "Offset") {
			offset = Vector3f(std::stof(elmParts[0]), std::stof(elmParts[1]), std::stof(elmParts[2]));
		} else if (attr == "ElementSpacing") {
			spacing = Vector3f(std::stof(elmParts[0]), std::stof(elmParts[1]), std::stof(elmParts[2]));
		} else if (attr == "DimSize") {
			imgDim = Vector3i(std::stoi(elmParts[0]), std::stoi(elmParts[1]), std::stoi(elmParts[2]));
		}
	}
	maskHeader.close();

	// Import CT data
	auto importer = ImageFileImporter::New();
	importer->setFilename(volPath);

	// Extract centerline from segmentation
	auto centerline = CenterlineExtraction::New();
	centerline->setInputConnection(importer->getOutputPort());
	centerline->setOffset(offset);

  // Export centerline
	auto vtkExporter = VTKMeshFileExporter::New();
	vtkExporter->setFilename(parser.get("output_path"));
	vtkExporter->setInputConnection(centerline->getOutputPort());
	vtkExporter->update();
}
