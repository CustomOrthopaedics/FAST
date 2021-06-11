#include "FAST/Exporters/MetaImageExporter.hpp"
#include "FAST/Exporters/VTKMeshFileExporter.hpp"
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

Vector3i mmToVx(Vector3f seed, Vector3f origin, Vector3f spacing);
Vector3f vxToMm(Vector3i seed, Vector3f origin, Vector3f spacing);
double toThreeDecimalPlaces(double d);

int main(int argc, char** argv) {
  Reporter::setGlobalReportMethod(Reporter::COUT); // TODO remove
	CommandLineParser parser("Airway segmentation");
	parser.addVariable("seg_config", true, "seg config .json");
	parser.addVariable("debug", "false", "whether or not to output voxData.json file");
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

	bool debug = parser.get("debug") == "true" ? true : false;

	// Do airway segmentation
	auto segmentation = AirwaySegmentation::New();
	segmentation->setInputConnection(importer->getOutputPort());
	segmentation->setVoxSpacing(spacing);
	segmentation->setDebug(debug);

	if (segData.contains("sensitivity")) {
		segmentation->setSensitivity(segData["sensitivity"]);
	} else {
		segmentation->setSensitivity(5);
	}

	auto newSeeds = segData["tracheal_points_vx"]["new"];
	for (auto seed : newSeeds) {
		segmentation->addSeedPoint(Vector3i(seed[0], seed[1], seed[2]));
	}

	auto existingSeeds = segData["tracheal_points_vx"]["existing"];
	for (auto seed : existingSeeds) {
		segmentation->addSeedPoint(Vector3i(seed[0], seed[1], seed[2]));
	}

	if (segData.contains("seg_bb_min_mm") && segData.contains("seg_bb_max_mm") && segData["seg_bb_min_mm"].size() == 3 && segData["seg_bb_max_mm"].size() == 3) {
		Vector3i bbMin(mmToVx(Vector3f(segData["seg_bb_min_mm"][0], segData["seg_bb_min_mm"][1], segData["seg_bb_min_mm"][2]), offset, spacing));
		Vector3i bbMax(mmToVx(Vector3f(segData["seg_bb_max_mm"][0], segData["seg_bb_max_mm"][1], segData["seg_bb_max_mm"][2]), offset, spacing));

		segmentation->setBB(bbMin, bbMax);
	}

	auto exporter = MetaImageExporter::New();
	std::string airwayPath = (std::string)segData["results_path"] + "/seg_" + std::to_string((int)segData["seg_id"]) + ".mhd";
	exporter->setFilename(airwayPath);
	exporter->setInputConnection(segmentation->getOutputPort());
	exporter->update();

	if (!debug) {
		return 0;
	}

	nlohmann::json voxelJSON;

	for (Voxel vox : segmentation->maskVoxels) {
		char key[12];
		sprintf(key, "%03d,%03d,%03d", vox.point.x(), vox.point.y(), vox.point.z());
		Vector3f posMM = vxToMm(vox.point, offset, spacing);
		voxelJSON[key]["posMM"] = std::vector<double>{posMM.x(), posMM.y(), posMM.z()};
		voxelJSON[key]["cent"] = toThreeDecimalPlaces(vox.centricity);
		voxelJSON[key]["minR"] = toThreeDecimalPlaces(vox.minRadius);
		voxelJSON[key]["meanR"] = toThreeDecimalPlaces(vox.meanRadii);
		voxelJSON[key]["sHU"] = vox.rays[0].startHU;
		voxelJSON[key]["mIdx"] = vox.maskIdx;
		voxelJSON[key]["pIdx"] = vox.pathVoxelIdx;

		// for (int i = 0; i < vox.rays.size(); ++i) {
		// 	voxelJSON[key]["rays"][i]["dir"] = std::vector<double>{toThreeDecimalPlaces(vox.rays[i].direction.x()), toThreeDecimalPlaces(vox.rays[i].direction.y()), toThreeDecimalPlaces(vox.rays[i].direction.z())};
		// 	voxelJSON[key]["rays"][i]["len"] = toThreeDecimalPlaces(vox.rays[i].length);
		// 	voxelJSON[key]["rays"][i]["eHU"] = vox.rays[i].endHU;
		// }
	}

	std::ofstream voxFile((std::string)segData["results_path"] + "/voxData.json");

	voxFile << voxelJSON << std::endl;

	voxFile.close();
}

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

// https://github.com/nlohmann/json/issues/1191#issuecomment-412688086
double toThreeDecimalPlaces(double d) {

  int i;

  if (d >= 0)

    i = static_cast<int>(d * 1000 + 0.5);

  else

    i = static_cast<int>(d * 1000 - 0.5);


  return (i / 1000.0);

}
