#ifndef _VTKWRITER_H_
#define _VTKWRITER_H_

//采用xml格式输出vtk格式文件

#include "../domain/scalarfield2D.h"
#include "../domain/vectorfield2D.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

class VtkWriter {
private:
	
	std::ofstream vtkFile;
	
	std::vector<void*> scalars, vectors, 
		scalarsSpecial, vectorsSpecial;//后两个是为了访问特殊的场，例如相场，它的数据类型不是double
	
	std::vector<std::string> scalarFieldName, vectorFieldName, 
		scalarSpecialFieldName, vectorSpecialFieldName;

public:

	VtkWriter(){}

	//1
	template <typename T>
	void add(LBM::ScalarField2D<T>& scalarField, const std::string& fieldName) {
		if (sizeof(T) == 8) {	//做这样的判断是为了区分一般的double类型的场与其他类型的场
			scalars.push_back(&scalarField);
			scalarFieldName.push_back(fieldName);
		}
		else {
			scalarsSpecial.push_back(&scalarField);
			scalarSpecialFieldName.push_back(fieldName);
		}
	}
	//2
	template <typename T>
	void add(LBM::VectorField2D<T>& vectorField, const std::string& fieldName) {
		if (sizeof(T) == 8) {	//做这样的判断是为了区分一般的double类型的场与其他类型的场
			vectors.push_back(&vectorField);
			vectorFieldName.push_back(fieldName);
		}
		else {
			vectorsSpecial.push_back(&vectorField);
			vectorSpecialFieldName.push_back(fieldName);
		}
	}

	//一般情况下 x0=1, x1=X, y0=1, y1=Y
	template <typename T1 = int, typename T2 = double>	
	void ascii(const std::string& fileName, const int& x0, const int& x1, const int& y0, const int& y1) {

		//一般T2不用改，因为是用于密度速度场的，直接double类型
		//如果有其他的场，例如相场，可以改T1，也可以不改，只要类型字节大小与int相同即可
		
		LBM::ScalarField2D<T2>* pScalarField = nullptr;	//用于访问一般标量场
		LBM::VectorField2D<T2>* pVectorField = nullptr;	//用于访问一般向量场
		LBM::ScalarField2D<T1>* pScalarSpecialField = nullptr;	//用于访问特殊标量场
		LBM::VectorField2D<T1>* pVectorSpecialField = nullptr;	//用于访问特殊向量场
		
		vtkFile.open(fileName);
		//第1行
		vtkFile << "<?xml version=\"1.0\"?>" << std::endl;
		//第2行
		vtkFile << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
		//第3行
		vtkFile
			<< "<ImageData WholeExtent= \""
			<< x0
			<< " "
			<< x1
			<< " "
			<< y0
			<< " "
			<< y1
			<< " 0 0\" Origin=\"0 0 0\" Spacing=\"1 1 1\">" << std::endl;
		//第4行
		vtkFile
			<< "<Piece Extent = \""
			<< x0 
			<< " "
			<< x1
			<< " "
			<< y0
			<< " "
			<< y1
			<< " 0 0\">" << std::endl;
		//第5行
		vtkFile
			<< "<PointData";

		//保证vectorField[0]可以在paraview中以向量形式呈现,注意 只有一个向量场可以以向量形式呈现,vectorField[1]就不行了
		if (!vectorFieldName.empty()) {
			vtkFile << " Vectors=\"" << vectorFieldName[0] << "\"";
		}
		vtkFile << ">" << std::endl;
		
		/**************************************** DataArray ****************************************/
		for (int n = 0; n < scalars.size(); ++n) {
			pScalarField = static_cast<LBM::ScalarField2D<T2>*>(scalars[n]);
			vtkFile << "<DataArray Name=\"" << scalarFieldName[n] << "\" type=\"Float64\" format=\"ascii\">" << std::endl;
			for (int j = y0; j <= y1; ++j) {
				for (int i = x0; i <= x1; ++i) {
					vtkFile << pScalarField->ps[i][j] << " ";
				}
			}
			vtkFile << std::endl;
			vtkFile << "</DataArray>" << std::endl;
		}
		
		for (int n = 0; n < scalarsSpecial.size(); ++n) {
			pScalarSpecialField = static_cast<LBM::ScalarField2D<T1>*>(scalarsSpecial[n]);
			vtkFile << "<DataArray Name=\"" << scalarSpecialFieldName[n] << "\" type=\"Int32\" format=\"ascii\">" << std::endl;
			for (int j = y0; j <= y1; ++j) {
				for (int i = x0; i <= x1; ++i) {
					vtkFile << pScalarSpecialField->ps[i][j] << " ";
				}
			}
			vtkFile << std::endl;
			vtkFile << "</DataArray>" << std::endl;
		}
		
		for (int n = 0; n < vectors.size(); ++n) {
			pVectorField = static_cast<LBM::VectorField2D<T2>*>(vectors[n]);
			vtkFile << "<DataArray Name=\"" << vectorFieldName[n] << "\" type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
			for (int j = y0; j <= y1; ++j) {
				for (int i = x0; i <= x1; ++i) {
					vtkFile << pVectorField->pv[i][j][0] << " " << pVectorField->pv[i][j][1] << " 0 ";
				}
			}
			vtkFile << std::endl;
			vtkFile << "</DataArray>" << std::endl;
		}
		
		for (int n = 0; n < vectorsSpecial.size(); ++n) {
			pVectorSpecialField = static_cast<LBM::VectorField2D<T1>*>(vectorsSpecial[n]);
			vtkFile << "<DataArray Name=\"" << vectorSpecialFieldName[n] << "\" type=\"Int32\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
			for (int j = y0; j <= y1; ++j) {
				for (int i = x0; i <= x1; ++i) {
					vtkFile << pVectorSpecialField->pv[i][j][0] << " " << pVectorSpecialField->pv[i][j][1] << " 0 ";
				}
			}
			vtkFile << std::endl;
			vtkFile << "</DataArray>" << std::endl;
		}
		
		/**************************************** 最后4行 ****************************************/
		vtkFile << "</PointData>" << std::endl;
		vtkFile << "</Piece>" << std::endl;
		vtkFile << "</ImageData>" << std::endl;
		vtkFile << "</VTKFile>" << std::endl;
		vtkFile.close();
	
	}


};

#endif // !_VTKWRITER_H_
