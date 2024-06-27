#pragma once
#ifndef OBJ_SEQUENCE_RENDER_H
#define OBJ_SEQUENCE_RENDER_H

#include "BaseSimulator.h"
#include "Model\tiny_obj_loader.h"

class ObjSequenceRender: public BaseSimulator
{
public:
public:
	ObjSequenceRender(std::string simName = "render_obj_sequence", RunPlatform runPlatform = RunPlatform::CPU) :
		BaseSimulator(simName, runPlatform)
	{}

	bool readSequence(const std::string config)
	{
		std::ifstream configFin;
		configFin.open(config);
		if (!configFin.is_open()) return false;
		std::string line;
		std::getline(configFin, dir);
		std::getline(configFin, prefix);
		configFin >> frameNum;
		configFin.close();
		for (int i = 0; i < frameNum; i++)
		{
			std::string filename;
			QString i_str;
			i_str.setNum(i);
			if (i < 10)
				filename = dir + prefix + "000" + i_str.toStdString() + ".obj";
			else if(i < 100)
				filename = dir + prefix + "00" + i_str.toStdString() + ".obj";
			else if (i < 1000)
				filename = dir + prefix + "0" + i_str.toStdString() + ".obj";
			else
				filename = dir + prefix + i_str.toStdString() + ".obj";

			if (i == 0)
			{
				BaseModel* model = new BaseModel();
				model->setValid(model->readMeshFromObjFormat(filename, this));
				model->init(false, false);
				models.push_back(model);

			}
			else
			{
				std::ifstream fin(filename.c_str());
				if (!fin.is_open())
					return false;

				tinyobj::attrib_t attrib;
				std::vector<tinyobj::shape_t> shapes;
				std::vector<tinyobj::material_t> groupMaterials;

				std::string err;
				std::string base_dir = getPathDir(filename);
				if (base_dir.empty())
				{
					base_dir = ".";
				}

				bool ret = tinyobj::LoadObj(&attrib, &shapes, &groupMaterials, &err, filename.c_str(), base_dir.c_str(), true);
				if (!err.empty())
				{
					std::cerr << err << std::endl;
				}
				if (!ret) {
					std::cerr << "Error: Failed to load " << filename << " !" << std::endl;
					fin.close();
					return false;
				}
				int pointsNum = (int)(attrib.vertices.size()) / 3;
				VectorX seq(totalPointsNum * 3);
				assert(pointsNum == totalPointsNum);
				for (int c = 0; c < attrib.vertices.size(); c++)
					seq[c] = attrib.vertices[c];
				fin.close();
				frames.push_back(seq);
			}
			std::cout << "read " << filename << std::endl;
		}
		return true;
	}

	bool readSequence()
	{
		for (int i = 0; i < frameNum; i++)
		{
			std::string filename;
			QString i_str;
			i_str.setNum(i);
			if (i < 10)
				filename = dir + prefix + "000" + i_str.toStdString() + ".obj";
			else if (i < 100)
				filename = dir + prefix + "00" + i_str.toStdString() + ".obj";
			else if (i < 1000)
				filename = dir + prefix + "0" + i_str.toStdString() + ".obj";
			else
				filename = dir + prefix + i_str.toStdString() + ".obj";

			if (i == 0)
			{
				BaseModel* model = new BaseModel();
				model->setValid(model->readMeshFromObjFormat(filename, this));
				model->init(false, false);
				models.push_back(model);

			}
			else
			{
				std::ifstream fin(filename.c_str());
				if (!fin.is_open())
					return false;

				tinyobj::attrib_t attrib;
				std::vector<tinyobj::shape_t> shapes;
				std::vector<tinyobj::material_t> groupMaterials;

				std::string err;
				std::string base_dir = getPathDir(filename);
				if (base_dir.empty())
				{
					base_dir = ".";
				}

				bool ret = tinyobj::LoadObj(&attrib, &shapes, &groupMaterials, &err, filename.c_str(), base_dir.c_str(), true);
				if (!err.empty())
				{
					std::cerr << err << std::endl;
				}
				if (!ret) {
					std::cerr << "Error: Failed to load " << filename << " !" << std::endl;
					fin.close();
					return false;
				}
				int pointsNum = (int)(attrib.vertices.size()) / 3;
				VectorX seq(totalPointsNum * 3);
				assert(pointsNum == totalPointsNum);
				for (int c = 0; c < attrib.vertices.size(); c++)
					seq[c] = attrib.vertices[c];
				fin.close();
				frames.push_back(seq);
			}
			std::cout << "read " << filename << std::endl;
		}
		return true;
	}

	virtual bool readSimulatorFromConfigFile(const std::string filename, TiXmlElement * item)
	{
		std::string itemName = item->Value();
		if (itemName != std::string("simulator"))
			return false;

		TiXmlElement* subItem = item->FirstChildElement();
		while (subItem)
		{
			std::string subItemName = subItem->Value();
			if (subItemName == std::string("dir"))
			{
				dir = subItem->GetText();
			}
			else if (subItemName == std::string("prefix"))
			{
				prefix = subItem->GetText();
			}
			else if (subItemName == std::string("frames"))
			{
				std::string text = subItem->GetText();
				std::strstream ss;
				ss << text;
				ss >> frameNum;
			}
			subItem = subItem->NextSiblingElement();
		}
		return readSequence();
	}



	virtual void run(int frame = 0)
	{
		if (frame >= frames.size())
			return;
		std::copy(frames[frame].data(), frames[frame].data() + 3 * totalPointsNum, models[0]->points.buffer);
	}

	std::string dir;
	std::string prefix;
	int frameNum;
	std::vector<VectorX> frames;
};




#endif