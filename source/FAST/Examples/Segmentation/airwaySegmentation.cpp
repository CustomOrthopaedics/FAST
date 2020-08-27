/**
 * Examples/Segmentation/airwaySegmentation.cpp
 *
 * If you edit this example, please also update the wiki and source code file in the repository.
 */
#include <FAST/Exporters/MetaImageExporter.hpp>
#include <FAST/Exporters/VTKMeshFileExporter.hpp>
#include "FAST/Algorithms/AirwaySegmentation/AirwaySegmentation.hpp"
#include "FAST/Importers/ImageFileImporter.hpp"
#include "FAST/Algorithms/SurfaceExtraction/SurfaceExtraction.hpp"
#include "FAST/Algorithms/CenterlineExtraction/CenterlineExtraction.hpp"
#include "FAST/Tools/CommandLineParser.hpp"
#include "nlohmann/json.hpp"
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

using namespace fast;

Vector3i mmToVx(Vector3f seed, Vector3f origin, Vector3f spacing) {
	int x = std::round(((seed[0] - origin[0]) / spacing[0]) + 0.5);
	int y = std::round(((seed[1] - origin[1]) / spacing[1]) + 0.5);
	int z = std::round(((seed[2] - origin[2]) / spacing[2]) + 0.5);

	return Vector3i(x, y, z);
}

Vector3f vxToMm(Vector3i seed, Vector3f origin, Vector3f spacing) {
	float x = (float(seed[0]) * spacing[0]) + origin[0];
	float y = (float(seed[1]) * spacing[1]) + origin[1];
	float z = (float(seed[2]) * spacing[2]) + origin[2];

	return Vector3f(x, y, z);
}

int main(int argc, char** argv) {
    Reporter::setGlobalReportMethod(Reporter::COUT); // TODO remove
	CommandLineParser parser("Airway segmentation");
	parser.addVariable("scan", false, "scan.mhd file");
    parser.addVariable("smoothing", "0.5", "How much smoothing to apply before segmentation");
    parser.addVariable("seed", false, "Manual seed point coordinate (--seed x,y,z)");
    parser.addVariable("export_segmentation", false, "Filename to export segmentation volume to. Example: segmentation.mhd");
    parser.addVariable("export_centerline", false, "Filename to export centerline to. Example: centerline.vtk");
	parser.addVariable("seg_config", true, "seg config .json");
	parser.parse(argc, argv);

	std::ifstream configFile(parser.get("seg_config"));
	nlohmann::json segData;
	configFile >> segData;
	configFile.close();

	std::string volPath = segData["vol_path"];
	std::ifstream scanHeader(volPath);

	Vector3f offset;
	Vector3f spacing;
	Vector3i imgDim;

	while (!scanHeader.eof()) {
		std::string line;
		getline(scanHeader, line);
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
	scanHeader.close();

	// Import CT data
	auto importer = ImageFileImporter::New();
	importer->setFilename(volPath);

	// Do airway segmentation
	auto segmentation = AirwaySegmentation::New();
	segmentation->setInputConnection(importer->getOutputPort());
	segmentation->setSmoothing(parser.get<float>("smoothing"));
	segmentation->setVoxSpacing(spacing);

	auto initialSeeds = segData["tracheal_points_mm"]["initial"];
	auto additionalSeeds = segData["tracheal_points_mm"]["additional"];
	for (int i = 0; i < initialSeeds.size(); ++i) {
		auto seed = initialSeeds[i];
		Vector3i seedVx = mmToVx(Vector3f(seed[0], seed[1], seed[2]), offset, spacing);
		segmentation->addSeedPoint(seedVx);
	}
	for (int i = 0; i < additionalSeeds.size(); ++i) {
		auto seed = additionalSeeds[i];
		Vector3i seedVx = mmToVx(Vector3f(seed[0], seed[1], seed[2]), offset, spacing);
		segmentation->addSeedPoint(seedVx);
	}

	// Extract centerline from segmentation
	auto centerline = CenterlineExtraction::New();
	centerline->setInputConnection(segmentation->getOutputPort());
	centerline->setOffset(offset);

	auto exporter = MetaImageExporter::New();
	std::string airwayPath = (std::string)segData["results_path"] + "/seg_" + std::to_string((int)segData["seg_id"]) + ".mhd";
	exporter->setFilename(airwayPath);
	exporter->setInputConnection(segmentation->getOutputPort());
	exporter->update();

	auto vtkExporter = VTKMeshFileExporter::New();
	std::string centerlinePath = (std::string)segData["results_path"] + "/seg_" + std::to_string((int)segData["seg_id"]) + ".centerline.vtk";
	vtkExporter->setFilename(centerlinePath);
	vtkExporter->setInputConnection(centerline->getOutputPort());
	vtkExporter->update();

    // Export data if user has specified to do so
    if(parser.gotValue("export_segmentation")) {
        auto exporter = MetaImageExporter::New();
        exporter->setFilename(parser.get("export_segmentation"));
        exporter->setInputConnection(segmentation->getOutputPort());
        exporter->update();
    }
    if(parser.gotValue("export_centerline")) {
        auto exporter = VTKMeshFileExporter::New();
        exporter->setFilename(parser.get("export_centerline"));
        exporter->setInputConnection(centerline->getOutputPort());
        exporter->update();
		
    }
// 	if(parser.gotValue("export_segmentation") || parser.gotValue("export_centerline"))
// 		return 0; // Dont render if exporting..

// 	// Extract surface from segmentation
// 	auto extraction = SurfaceExtraction::New();
// 	extraction->setInputConnection(segmentation->getOutputPort());

// 	// Set up renderers and window
// 	auto renderer = TriangleRenderer::New();
// 	renderer->addInputConnection(extraction->getOutputPort());

// 	auto lineRenderer = LineRenderer::New();
// 	lineRenderer->addInputConnection(centerline->getOutputPort());
// 	lineRenderer->setDefaultDrawOnTop(true);
// 	lineRenderer->setDefaultColor(Color::Blue());

// 	auto window = SimpleWindow::New();
// 	window->addRenderer(renderer);
// 	window->addRenderer(lineRenderer);
// #ifdef FAST_CONTINUOUS_INTEGRATION
// 	// This will automatically close the window after 5 seconds, used for CI testing
//     window->setTimeout(5*1000);
// #endif
// 	window->start();

	if (initialSeeds.size() + additionalSeeds.size() == 0) {
		Vector3i autoSeed = segmentation->autoSeed;
		if (autoSeed == Vector3i::Zero()) {
			// no seed found. catch error and handle here?
		}

		Vector3f seedVx = vxToMm(autoSeed, offset, spacing);

		segData["tracheal_points_mm"]["initial"] = { {seedVx[0], seedVx[1], seedVx[2]} };
		segData["tracheal_points_mm"]["additional"] = std::vector<float>{};
	} else {
		segData["tracheal_points_mm"]["initial"] = std::vector<float>{};
		segData["tracheal_points_mm"]["additional"] = std::vector<float>{};
	}

	std::string outPath = (std::string)segData["results_path"] + "/seg_" + std::to_string((int)segData["seg_id"]) + ".tracheal_points.json";
	std::ofstream algOut(outPath);
	algOut << std::setw(4) << segData << std::endl;
	algOut.close();
}
