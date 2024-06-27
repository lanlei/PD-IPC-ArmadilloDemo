#include "BaseRenderMaterial.h"
#include <iostream>
#include <string>
#include "Commom\DataConversion.h"
BaseRenderMaterial::BaseRenderMaterial()
{
//	ambient = QColor(64, 64, 166);
	ambient = QColor(240, 128, 128);	
	diffuse = QColor(76, 76, 76);
	specular = QColor(76, 76, 76);
	shinness = 6.0;
	
	for (int i = 0; i < 4; i++)
	{
		textureMapFilename[i] = std::string("");
		useTextureMap[i] = 0;
		_textureMap[i] = nullptr;
	}
	_transparent = false;
}

bool BaseRenderMaterial::readTextureMap(const QString filename, int mapIndex)
{
	std::string img_filename = filename.toStdString();
	QImage image = QImage(img_filename.c_str());
	
	if (image.isNull())
	{
		std::cout << "The texture " << img_filename << " don't exist !" << std::endl;
		return false;
	}

	if (_textureMap[mapIndex] != nullptr)
		delete _textureMap[mapIndex];
	_textureMap[mapIndex] = new QOpenGLTexture(QImage(img_filename.c_str()).mirrored());

	if (_textureMap[mapIndex]->isCreated())
	{
		_textureMap[mapIndex]->setMinificationFilter(QOpenGLTexture::Nearest);
		_textureMap[mapIndex]->setMagnificationFilter(QOpenGLTexture::Linear);
		_textureMap[mapIndex]->setWrapMode(QOpenGLTexture::Repeat);
		textureMapFilename[mapIndex] = img_filename;
		useTextureMap[mapIndex] = 1;
	}else useTextureMap[mapIndex] = 0;
	return true;
}

bool BaseRenderMaterial::readCubeTextureMap(std::vector<QString>& imgFilename, int mapIndex)
{
	if (_textureMap[mapIndex] != nullptr)
		delete _textureMap[mapIndex];
	_textureMap[mapIndex] = new QOpenGLTexture(QOpenGLTexture::TargetCubeMap);
	
	const QImage posx = QImage(imgFilename[0]).convertToFormat(QImage::Format_RGBA8888);
	const QImage negx = QImage(imgFilename[1]).convertToFormat(QImage::Format_RGBA8888);

	const QImage posy = QImage(imgFilename[2]).convertToFormat(QImage::Format_RGBA8888);
	const QImage negy = QImage(imgFilename[3]).convertToFormat(QImage::Format_RGBA8888);

	const QImage posz = QImage(imgFilename[4]).convertToFormat(QImage::Format_RGBA8888);
	const QImage negz = QImage(imgFilename[5]).convertToFormat(QImage::Format_RGBA8888);

	if (!_textureMap[mapIndex]->isCreated())
		_textureMap[mapIndex]->create();

	_textureMap[mapIndex]->setSize(posx.width(), posx.height(), posx.depth());
	_textureMap[mapIndex]->setFormat(QOpenGLTexture::RGBA8_UNorm);
	_textureMap[mapIndex]->allocateStorage();

	_textureMap[mapIndex]->setData(0, 0, QOpenGLTexture::CubeMapPositiveX,
		QOpenGLTexture::RGBA, QOpenGLTexture::UInt8,
		posx.constBits(), Q_NULLPTR);

	_textureMap[mapIndex]->setData(0, 0, QOpenGLTexture::CubeMapPositiveY,
		QOpenGLTexture::RGBA, QOpenGLTexture::UInt8,
		posy.constBits(), Q_NULLPTR);

	_textureMap[mapIndex]->setData(0, 0, QOpenGLTexture::CubeMapPositiveZ,
		QOpenGLTexture::RGBA, QOpenGLTexture::UInt8,
		posz.constBits(), Q_NULLPTR);

	_textureMap[mapIndex]->setData(0, 0, QOpenGLTexture::CubeMapNegativeX,
		QOpenGLTexture::RGBA, QOpenGLTexture::UInt8,
		negx.constBits(), Q_NULLPTR);

	_textureMap[mapIndex]->setData(0, 0, QOpenGLTexture::CubeMapNegativeY,
		QOpenGLTexture::RGBA, QOpenGLTexture::UInt8,
		negy.constBits(), Q_NULLPTR);

	_textureMap[mapIndex]->setData(0, 0, QOpenGLTexture::CubeMapNegativeZ,
		QOpenGLTexture::RGBA, QOpenGLTexture::UInt8,
		negz.constBits(), Q_NULLPTR);

	_textureMap[mapIndex]->setWrapMode(QOpenGLTexture::ClampToEdge);
	_textureMap[mapIndex]->setMinificationFilter(QOpenGLTexture::LinearMipMapLinear);
	_textureMap[mapIndex]->setMagnificationFilter(QOpenGLTexture::LinearMipMapLinear);
	
	return false;
}

void BaseRenderMaterial::transferToShader(QOpenGLShaderProgram * program)
{
	program->setUniformValue("material.ambient", getOpenglRGBA(ambient));
	program->setUniformValue("material.diffuse", getOpenglRGBA(diffuse));
	program->setUniformValue("material.specular", getOpenglRGBA(specular));
	program->setUniformValue("material.shininess", getOpenglRGBA(shinness));

	if (useTextureMap[AmbientMapIndex] == 1 && _textureMap[AmbientMapIndex]->isCreated())
	{
		program->setUniformValue("material.useAmbientMap", useTextureMap[AmbientMapIndex]);
		getTextureMap(AmbientMapIndex)->bind(1);
	}else program->setUniformValue("material.useAmbientMap", 0);
		
	if (useTextureMap[DiffuseMapIndex] == 1 && _textureMap[DiffuseMapIndex]->isCreated())
	{
		program->setUniformValue("material.useDiffuseMap", useTextureMap[DiffuseMapIndex]);
		getTextureMap(DiffuseMapIndex)->bind(2);
	}else program->setUniformValue("material.useDiffuseMap", 0);
		
	if (useTextureMap[SpecularMapIndex] == 1 && _textureMap[SpecularMapIndex]->isCreated())
	{
		program->setUniformValue("material.useSpecularMap", useTextureMap[SpecularMapIndex]);
		getTextureMap(SpecularMapIndex)->bind(3);
	}else program->setUniformValue("material.useSpecularMap", 0);
		
	if (useTextureMap[BumpMapIndex] == 1 && _textureMap[BumpMapIndex]->isCreated())
	{
		program->setUniformValue("material.useBumpMap", useTextureMap[BumpMapIndex]);
		getTextureMap(BumpMapIndex)->bind(4);
	}else program->setUniformValue("material.useBumpMap", 0);
		
}
